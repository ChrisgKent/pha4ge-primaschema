from pydantic import model_validator

from primaschema.schema.info import PrimerScheme as _GeneratedPrimerScheme

# info.py is generated straight from info.yml by gen-pydantic and gets
# overwritten by the pre-commit hook, so it can't hand-encode behaviour that
# gen-pydantic doesn't translate from LinkML into the Pydantic model. One
# such gap: primer_scheme_identifier's LinkML `structured_pattern` (its
# composition from name/amplicon_size/version) is dropped entirely during
# generation — the generated field is just a plain optional string with no
# computation or validation behind it. This subclass adds that behaviour
# back on top of the generated model instead.


def compute_primer_scheme_identifier(
    name: str, amplicon_size: int, version: str
) -> str:
    return f"{name}/{amplicon_size}/{version}"


class PrimerScheme(_GeneratedPrimerScheme):
    @model_validator(mode="after")
    def _sync_primer_scheme_identifier(self):
        expected = compute_primer_scheme_identifier(
            self.primer_scheme_name, self.amplicon_size, self.primer_scheme_version
        )
        if self.primer_scheme_identifier != expected:
            self.primer_scheme_identifier = expected
        return self
