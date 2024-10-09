#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

./scripts/stage4/install_cereal.sh
./scripts/stage4/install_rapidjson.sh
./scripts/stage4/install_libtorch.sh
./scripts/stage4/install_libnpy.sh
./scripts/stage4/install_libri.sh
./scripts/stage4/install_libcomm.sh

# EOF
