#!/usr/bin/env bash
# Provisions ABACUS inside a WSL2 Ubuntu distribution.
# Expected to run as root (installer invokes `wsl -u root`).

set -euo pipefail

MINIFORGE_DIR="/opt/abacus-miniforge"
ENV_NAME="abacus_env"
ENV_BIN="$MINIFORGE_DIR/envs/$ENV_NAME/bin"
CHINA_MIRROR="${ABACUS_CHINA_MIRROR:-0}"

if [ "$CHINA_MIRROR" = "1" ]; then
    MINIFORGE_URL="https://mirrors.tuna.tsinghua.edu.cn/github-release/conda-forge/miniforge/LatestRelease/Miniforge3-Linux-x86_64.sh"
    CONDA_FORGE_CHANNEL="https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge"
else
    MINIFORGE_URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
    CONDA_FORGE_CHANNEL="conda-forge"
fi

log() { printf '[provision] %s\n' "$*"; }

if [ "$(id -u)" -ne 0 ]; then
    echo "[provision] Must run as root." >&2
    exit 1
fi

if [ "$CHINA_MIRROR" = "1" ]; then
    log "China mirror mode: using TUNA (Tsinghua) for apt / Miniforge / conda-forge."
    if [ -f /etc/apt/sources.list ] && grep -qE 'archive\.ubuntu\.com|security\.ubuntu\.com' /etc/apt/sources.list; then
        log "Rewriting /etc/apt/sources.list to TUNA (backup at sources.list.orig)..."
        cp -n /etc/apt/sources.list /etc/apt/sources.list.orig
        sed -i \
            -e 's|http://archive.ubuntu.com/ubuntu|https://mirrors.tuna.tsinghua.edu.cn/ubuntu|g' \
            -e 's|http://security.ubuntu.com/ubuntu|https://mirrors.tuna.tsinghua.edu.cn/ubuntu|g' \
            /etc/apt/sources.list
    fi
fi

log "Installing apt prerequisites (apt-get update + curl/ca-certs/bzip2)..."
export DEBIAN_FRONTEND=noninteractive
apt-get update
apt-get install -y --no-install-recommends curl ca-certificates bzip2

if [ ! -x "$MINIFORGE_DIR/bin/conda" ]; then
    log "Downloading Miniforge installer (~80 MB)..."
    tmp="$(mktemp --suffix=.sh)"
    trap 'rm -f "$tmp"' EXIT
    curl -fL --progress-bar -o "$tmp" "$MINIFORGE_URL"
    log "Installing Miniforge into $MINIFORGE_DIR..."
    # Run the installer via its own shebang (dash on Ubuntu) rather than
    # `bash installer.sh`. Some Miniforge builds have a sourced-check that
    # misfires under bash and aborts with "Please run using bash/dash/...".
    chmod +x "$tmp"
    "$tmp" -b -p "$MINIFORGE_DIR"
    rm -f "$tmp"
    trap - EXIT
else
    log "Miniforge already present at $MINIFORGE_DIR."
fi

# shellcheck disable=SC1091
source "$MINIFORGE_DIR/etc/profile.d/conda.sh"

if conda env list | awk 'NF && $1 !~ /^#/ {print $1}' | grep -qx "$ENV_NAME"; then
    log "Updating existing env '$ENV_NAME' (channel: $CONDA_FORGE_CHANNEL)..."
    conda install -n "$ENV_NAME" -y --override-channels -c "$CONDA_FORGE_CHANNEL" abacus
else
    log "Creating env '$ENV_NAME' (channel: $CONDA_FORGE_CHANNEL)..."
    conda create -n "$ENV_NAME" -y --override-channels -c "$CONDA_FORGE_CHANNEL" abacus
fi

log "Installing system launchers..."

cat > /usr/local/bin/abacus <<EOF
#!/usr/bin/env bash
source "$MINIFORGE_DIR/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"
exec "$ENV_BIN/abacus" "\$@"
EOF
chmod +x /usr/local/bin/abacus

cat > /usr/local/bin/abacus-mpi <<EOF
#!/usr/bin/env bash
source "$MINIFORGE_DIR/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"
exec mpirun "\$@" "$ENV_BIN/abacus"
EOF
chmod +x /usr/local/bin/abacus-mpi

log "Verifying..."
if /usr/local/bin/abacus --version >/dev/null 2>&1; then
    /usr/local/bin/abacus --version || true
elif [ -x "$ENV_BIN/abacus" ]; then
    log "abacus binary present at $ENV_BIN/abacus (no --version flag)."
else
    log "WARNING: abacus binary not found after install."
    exit 1
fi

log "Done."
