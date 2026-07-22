import hashlib
import json
import logging
import os
import tempfile
from io import BytesIO
from pathlib import Path
from uuid import uuid4

import dnaio
from fastapi import FastAPI, File, Form, Request, UploadFile
from fastapi.responses import HTMLResponse
from fastapi.templating import Jinja2Templates
from github import Github
from primalbedtools.scheme import Scheme
from primalbedtools.validate import validate_ref_and_bed
from pydantic import BaseModel, ValidationError

from primaschema.schema.info import Checksums, PrimerScheme
from primaschema.util import serialize_fasta_records

logger = logging.getLogger(__name__)

app = FastAPI(title="primaschema submission portal", version="0.1.0")

templates = Jinja2Templates(directory=Path(__file__).parent / "templates")

# ---------------------------------------------------------------------------
# Upload limits — all overridable via environment variables
# ---------------------------------------------------------------------------
MAX_METADATA_BYTES = int(os.environ.get("PRIMASCHEMA_MAX_METADATA_KB", "64")) * 1024
MAX_BED_BYTES = int(os.environ.get("PRIMASCHEMA_MAX_BED_MB", "10")) * 1024 * 1024
MAX_FASTA_BYTES = int(os.environ.get("PRIMASCHEMA_MAX_FASTA_MB", "50")) * 1024 * 1024
MAX_BEDLINES = int(os.environ.get("PRIMASCHEMA_MAX_BEDLINES", "10000"))
MAX_FASTA_RECORDS = int(os.environ.get("PRIMASCHEMA_MAX_FASTA_RECORDS", "100"))
MAX_SEQ_LENGTH = int(os.environ.get("PRIMASCHEMA_MAX_SEQ_LENGTH", "5000000"))


async def _read_upload(upload: UploadFile, max_bytes: int, label: str) -> bytes:
    data = await upload.read(max_bytes + 1)
    if len(data) > max_bytes:
        raise ValueError(
            f"{label} exceeds {max_bytes // (1024 * 1024)} MB upload limit"
        )
    return data


# ---------------------------------------------------------------------------
# Response models
# ---------------------------------------------------------------------------
class ValidationResponse(BaseModel):
    valid: bool
    errors: list[str] = []
    scheme: dict | None = None


class SubmitResponse(BaseModel):
    submitted: bool
    errors: list[str] = []
    branch: str | None = None


# ---------------------------------------------------------------------------
# Core validation logic (shared by /validate/full and /submit)
# ---------------------------------------------------------------------------
async def _run_validation(
    metadata: str,
    primer_bed: UploadFile,
    reference_fasta: UploadFile,
) -> tuple[ValidationResponse, bytes, bytes, PrimerScheme | None]:
    """Run full in-memory validation. Returns (response, bed_bytes, ref_bytes, scheme).

    bed_bytes and ref_bytes are the normalised serialised content.
    scheme is None if validation failed.
    """
    errors: list[str] = []

    if len(metadata.encode()) > MAX_METADATA_BYTES:
        return (
            ValidationResponse(
                valid=False,
                errors=[f"metadata exceeds {MAX_METADATA_BYTES // 1024} KB limit"],
            ),
            b"",
            b"",
            None,
        )

    try:
        scheme = PrimerScheme.model_validate_json(metadata)
    except ValidationError as exc:
        return (
            ValidationResponse(
                valid=False,
                errors=[
                    f"{'.'.join(str(loc) for loc in e['loc'])}: {e['msg']}"
                    for e in exc.errors()
                ],
            ),
            b"",
            b"",
            None,
        )
    except Exception:
        logger.exception("metadata JSON parsing failed")
        return (
            ValidationResponse(valid=False, errors=["metadata: invalid JSON"]),
            b"",
            b"",
            None,
        )

    try:
        bed_bytes = await _read_upload(primer_bed, MAX_BED_BYTES, "primer.bed")
    except ValueError as exc:
        errors.append(str(exc))
        bed_bytes = b""

    try:
        ref_bytes = await _read_upload(
            reference_fasta, MAX_FASTA_BYTES, "reference.fasta"
        )
    except ValueError as exc:
        errors.append(str(exc))
        ref_bytes = b""

    if errors:
        return ValidationResponse(valid=False, errors=errors), b"", b"", None

    scheme_bed: Scheme | None = None
    try:
        scheme_bed = Scheme.from_str(bed_bytes.decode("utf-8"))
        scheme_bed.sort_bedlines()
    except UnicodeDecodeError:
        errors.append("primer.bed: file must be valid UTF-8 text")
    except Exception:
        logger.exception("BED parsing failed")
        errors.append("primer.bed: parsing failed — check file format")

    if scheme_bed is not None and len(scheme_bed.bedlines) > MAX_BEDLINES:
        errors.append(f"primer.bed: exceeds {MAX_BEDLINES} primer limit")
        scheme_bed = None

    ref_records: list = []
    try:
        with dnaio.open(BytesIO(ref_bytes), fileformat="fasta") as reader:
            for record in reader:
                if not record.sequence:
                    errors.append(
                        f"reference.fasta: empty sequence for record {record.id!r}"
                    )
                    break
                if len(record.sequence) > MAX_SEQ_LENGTH:
                    errors.append(
                        f"reference.fasta: sequence {record.id!r} exceeds "
                        f"{MAX_SEQ_LENGTH // 1_000_000} Mbp limit"
                    )
                    break
                ref_records.append(record)
                if len(ref_records) > MAX_FASTA_RECORDS:
                    errors.append(
                        f"reference.fasta: exceeds {MAX_FASTA_RECORDS} record limit"
                    )
                    ref_records = []
                    break
    except Exception:
        logger.exception("FASTA parsing failed")
        errors.append("reference.fasta: parsing failed — check file format")

    if not errors and not ref_records:
        errors.append("reference.fasta: contains no records")

    if errors:
        return ValidationResponse(valid=False, errors=errors), b"", b"", None

    normalised_bed_bytes = scheme_bed.to_str().encode("utf-8")  # type: ignore[union-attr]
    normalised_ref_bytes = serialize_fasta_records(ref_records)

    computed_bed_sha = hashlib.sha256(normalised_bed_bytes).hexdigest()
    computed_ref_sha = hashlib.sha256(normalised_ref_bytes).hexdigest()

    if scheme.checksums:
        if (
            scheme.checksums.primer_sha256
            and computed_bed_sha != scheme.checksums.primer_sha256
        ):
            errors.append("primer.bed: SHA256 does not match checksums in metadata")
        if (
            scheme.checksums.reference_sha256
            and computed_ref_sha != scheme.checksums.reference_sha256
        ):
            errors.append(
                "reference.fasta: SHA256 does not match checksums in metadata"
            )

    if errors:
        return ValidationResponse(valid=False, errors=errors), b"", b"", None

    try:
        with tempfile.TemporaryDirectory() as tmp:
            ref_path = Path(tmp) / "reference.fasta"
            ref_path.write_bytes(normalised_ref_bytes)
            validate_ref_and_bed(scheme_bed.bedlines, str(ref_path.absolute()))  # type: ignore[union-attr]
    except Exception:
        logger.exception("primer/reference cross-validation failed")
        errors.append(
            "primer/reference cross-validation failed — primers may not align to reference"
        )

    if errors:
        return ValidationResponse(valid=False, errors=errors), b"", b"", None

    scheme.checksums = Checksums(
        primer_sha256=computed_bed_sha,
        reference_sha256=computed_ref_sha,
    )
    result = ValidationResponse(
        valid=True,
        scheme=json.loads(scheme.model_dump_json(exclude_none=True)),
    )
    return result, normalised_bed_bytes, normalised_ref_bytes, scheme


# ---------------------------------------------------------------------------
# GitHub helpers
# ---------------------------------------------------------------------------
def _push_scheme_to_branch(
    scheme: PrimerScheme,
    metadata: str,
    bed_bytes: bytes,
    ref_bytes: bytes,
    branch: str,
) -> None:
    g = Github(os.environ["GITHUB_TOKEN"])
    repo = g.get_repo(os.environ["PRIMER_SCHEMES_REPO"])
    base_branch = os.environ.get("PRIMER_SCHEMES_BASE_BRANCH", "main")
    base = repo.get_branch(base_branch)

    repo.create_git_ref(f"refs/heads/{branch}", base.commit.sha)
    prefix = f"{scheme.name}/{scheme.amplicon_size}/{scheme.version}"

    repo.create_file(
        f"{prefix}/info.json", "add info.json", metadata.encode(), branch=branch
    )
    repo.create_file(f"{prefix}/primer.bed", "add primer.bed", bed_bytes, branch=branch)
    repo.create_file(
        f"{prefix}/reference.fasta", "add reference.fasta", ref_bytes, branch=branch
    )


def _trigger_pr_action(scheme: PrimerScheme, branch: str) -> None:
    g = Github(os.environ["GITHUB_TOKEN"])
    submit_repo = g.get_repo(os.environ["SUBMIT_REPO"])
    submit_repo.create_repository_dispatch(
        "create-scheme-pr",
        {
            "branch": branch,
            "name": scheme.name,
            "version": scheme.version,
            "amplicon_size": scheme.amplicon_size,
        },
    )


# ---------------------------------------------------------------------------
# Routes
# ---------------------------------------------------------------------------
@app.get("/health")
def health():
    return {"status": "ok"}


@app.get("/", response_class=HTMLResponse)
def form_page(request: Request):
    return templates.TemplateResponse(request, "index.html")


@app.post("/validate", response_model=ValidationResponse)
def validate(body: dict):
    try:
        scheme = PrimerScheme.model_validate(body)
        return ValidationResponse(
            valid=True,
            scheme=json.loads(scheme.model_dump_json(exclude_none=True)),
        )
    except ValidationError as exc:
        return ValidationResponse(
            valid=False,
            errors=[
                f"{'.'.join(str(loc) for loc in e['loc'])}: {e['msg']}"
                for e in exc.errors()
            ],
        )


@app.post("/validate/full", response_model=ValidationResponse)
async def validate_full(
    metadata: str = Form(...),
    primer_bed: UploadFile = File(...),
    reference_fasta: UploadFile = File(...),
):
    result, _, _, _ = await _run_validation(metadata, primer_bed, reference_fasta)
    return result


@app.post("/submit", response_model=SubmitResponse)
async def submit(
    metadata: str = Form(...),
    primer_bed: UploadFile = File(...),
    reference_fasta: UploadFile = File(...),
):
    result, normalised_bed_bytes, normalised_ref_bytes, scheme = await _run_validation(
        metadata, primer_bed, reference_fasta
    )

    if not result.valid:
        return SubmitResponse(submitted=False, errors=result.errors)

    branch = f"submit/{scheme.name}-{scheme.version}-{uuid4().hex[:8]}"  # type: ignore[union-attr]

    try:
        _push_scheme_to_branch(
            scheme,
            metadata,
            normalised_bed_bytes,
            normalised_ref_bytes,
            branch,  # type: ignore[arg-type]
        )
    except Exception:
        logger.exception("failed to push scheme branch to GitHub")
        return SubmitResponse(
            submitted=False,
            errors=["submission failed — could not push files to GitHub"],
        )

    try:
        _trigger_pr_action(scheme, branch)  # type: ignore[arg-type]
    except Exception:
        logger.exception("failed to trigger PR action")
        return SubmitResponse(
            submitted=False,
            errors=[
                "submission failed — files pushed but could not trigger PR workflow"
            ],
        )

    return SubmitResponse(submitted=True, branch=branch)


@app.get("/schema")
def schema():
    return PrimerScheme.model_json_schema()


def main():
    import uvicorn

    uvicorn.run("primaschema_submit.server:app", host="127.0.0.1", port=8000)
