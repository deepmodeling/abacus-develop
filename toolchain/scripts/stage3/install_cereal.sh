#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all
# CEREAL is not need any complex setting

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

cereal_ver="1.3.2"
#cereal_sha256="a0f6f1bba7ba5c0c85b2bfe65aca1591025f509a7f11471b4cd651a79491b045"
source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

# [ -f "${BUILDDIR}/setup_cereal" ] && rm "${BUILDDIR}/setup_cereal"

# ! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
# cd "${BUILDDIR}"


echo "==================== Installing CEREAL ===================="
pkg_install_dir="./libs/cereal-${cereal_ver}"
# install_lock_file="$pkg_install_dir/install_successful"
# if verify_checksums "${install_lock_file}"; then
#     echo "cereal-${cereal_ver} is already installed, skipping it."
# else
#     if [ -f cereal-${cereal_ver}.tar.gz ]; then
#     echo "cereal-${cereal_ver}.tar.gz is found"
#     else
#     download_pkg_from_ABACUS_org "${cereal_sha256}" "cereal-${cereal_ver}.tar.gz"
#     fi
#     echo "Installing from scratch into ${pkg_install_dir}"
#     [ -d cereal-${cereal_ver} ] && rm -rf cereal-${cereal_ver}
cd ./libs
tar -xzf cereal-${cereal_ver}.tar.gz

# fi
# ;;

# EOF
# fi

# load "${BUILDDIR}/setup_cereal"
# write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "cereal"
