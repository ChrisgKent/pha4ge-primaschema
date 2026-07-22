FROM python:3.12-slim

WORKDIR /app

COPY --from=ghcr.io/astral-sh/uv:latest /uv /usr/local/bin/uv

COPY pyproject.toml uv.lock* ./
COPY src/ ./src/

RUN uv sync --no-dev

EXPOSE 8000

CMD [".venv/bin/gunicorn", \
     "--bind", "0.0.0.0:8000", \
     "--workers", "2", \
     "-k", "uvicorn.workers.UvicornWorker", \
     "primaschema_submit.server:app"]
