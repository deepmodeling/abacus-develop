#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

./scripts/stage3/install_libxc.sh
./scripts/stage3/install_fftw.sh
./scripts/stage3/install_elpa.sh
./scripts/stage3/install_scalapack.sh

# EOF
