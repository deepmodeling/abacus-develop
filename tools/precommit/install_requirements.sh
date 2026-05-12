#!/bin/bash -e

# Based on CP2K's tools/precommit/install_requirements.sh, adapted for ABACUS.

export DEBIAN_FRONTEND=noninteractive DEBCONF_NONINTERACTIVE_SEEN=true
apt-get update -qq
apt-get install -qq --no-install-recommends \
  ca-certificates \
  clang-format \
  git \
  less \
  nano \
  python3 \
  python3-venv \
  python3-pip \
  python3-wheel \
  python3-setuptools \
  shellcheck \
  wget
rm -rf /var/lib/apt/lists/*

python3 -m venv /opt/venv
export PATH="/opt/venv/bin:$PATH"
pip3 install --quiet -r requirements.txt

# Install shfmt.
wget -q https://github.com/mvdan/sh/releases/download/v3.2.2/shfmt_v3.2.2_linux_amd64
echo '3a32a69286a19491a81fcd854154f0d886c379ff28d99e32d5594490b8bbef4b shfmt_v3.2.2_linux_amd64' | sha256sum --check
chmod +x shfmt_v3.2.2_linux_amd64
ln -s /opt/abacus-precommit/shfmt_v3.2.2_linux_amd64 /usr/bin/shfmt

# EOF
