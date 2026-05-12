FROM ubuntu:24.04

# Based on CP2K's tools/precommit/Dockerfile, adapted for ABACUS.

WORKDIR /opt/abacus-precommit
COPY install_requirements.sh requirements.txt ./
RUN ./install_requirements.sh
ENV PATH="/opt/venv/bin:/opt/abacus-precommit:$PATH"

ARG REVISION
ENV REVISION=${REVISION}

COPY . ./
CMD ["gunicorn", "--bind=:8080", "--workers=1", "--threads=8", "--timeout=0", "precommit_server:app"]

# EOF
