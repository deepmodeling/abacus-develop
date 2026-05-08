# syntax=docker/dockerfile:1.7
# Usage: DOCKER_BUILDKIT=1 docker buildx build -f Dockerfile -t abacus-make-pretty ../.. --target export --progress=plain --output type=local,dest=../..

FROM ubuntu:24.04 AS pretty

ARG DEBIAN_FRONTEND=noninteractive
ARG BUILD_DIR=build
ARG JOBS=16
ARG CMAKE_ARGS=""
ARG CLANG_TIDY_EXTRA_ARGS=""
ARG STRICT_CLANG_TIDY=0

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN apt-get update -qq
RUN apt-get install -qy --no-install-recommends \
    ca-certificates \
    bash \
    git \
    findutils \
    build-essential \
    cmake \
    ninja-build \
    pkg-config \
    python3 \
    python3-pip \
    python3-venv \
    openmpi-bin \
    libopenmpi-dev \
    libfftw3-dev \
    libelpa-dev \
    libopenblas-dev \
    libscalapack-openmpi-dev \
    libxc-dev \
    libcereal-dev \
    libgtest-dev \
    libgmock-dev \
    libomp-dev
RUN apt-get clean

RUN python3 -m venv /opt/precommit-venv \
    && /opt/precommit-venv/bin/pip install --no-cache-dir --upgrade pip \
    && /opt/precommit-venv/bin/pip install --no-cache-dir clang-format clang-tidy
ENV PATH="/opt/precommit-venv/bin:${PATH}"

WORKDIR /repo

COPY . .

ENV BUILD_DIR=${BUILD_DIR}
ENV JOBS=${JOBS}
ENV CMAKE_ARGS=${CMAKE_ARGS}
ENV CLANG_TIDY_EXTRA_ARGS=${CLANG_TIDY_EXTRA_ARGS}
ENV STRICT_CLANG_TIDY=${STRICT_CLANG_TIDY}

RUN python3 /repo/tools/precommit/make_pretty.py

FROM scratch AS export
COPY --from=pretty /out/ /
