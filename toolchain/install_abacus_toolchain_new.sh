#!/bin/bash -e

# ============================================================================
# ABACUS Toolchain Installation Script (Modular Version)
# ============================================================================
# This is the new modular version of the ABACUS toolchain installation script.
# It provides the same functionality as the original script but with improved
# modularity, maintainability, and extensibility.
#
# Author: ABACUS Development Team
# Date: 2025-01-12
# ============================================================================

# Set script name for error reporting
[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0

# Set directory variables
export ROOTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export SCRIPTDIR="${ROOTDIR}/scripts"
export BUILDDIR="${ROOTDIR}/build"
export INSTALLDIR="${ROOTDIR}/install"
export SETUPFILE="${INSTALLDIR}/setup"
export SHA256_CHECKSUM="${SCRIPTDIR}/checksums.sha256"

# Make a copy of all options for $SETUPFILE
TOOLCHAIN_OPTIONS="$@"

# Source required modules
source "${SCRIPTDIR}/tool_kit.sh"
source "${SCRIPTDIR}/lib/error_handler.sh"
source "${SCRIPTDIR}/lib/config_manager.sh"
source "${SCRIPTDIR}/lib/version_helper.sh"
source "${SCRIPTDIR}/lib/package_manager.sh"
source "${SCRIPTDIR}/lib/user_interface.sh"
source "${SCRIPTDIR}/lib/config_validator.sh"

# Initialize modules (version helper only - config will be initialized in main)
version_helper_init

# Show help function
show_help() {
    ui_show_help
}

# Main function
main() {
    local args=("$@")
    
    # Initialize configuration with command line arguments
    if ! config_init "${args[@]}"; then
        show_help
        exit 0
    fi
    
    # Handle special version-related requests
    if [[ "$(config_get show_version)" == "__TRUE__" ]]; then
        version_show_available
        exit 0
    fi
    
    local version_info_request="$(config_get show_version_info)"
    if [[ -n "$version_info_request" ]]; then
        if [[ "$version_info_request" == "all" ]]; then
            version_show_available
        else
            version_show_available "$version_info_request"
        fi
        exit 0
    fi
    
    if [[ "$(config_get show_help)" == "__TRUE__" ]]; then
        show_help
        version_show_help
        exit 0
    fi
    
    # Show configuration summary if requested
    if [[ "$(config_get dry_run)" == "__TRUE__" ]]; then
        echo "=== DRY RUN MODE ==="
        echo "Configuration will be written but no packages will be installed."
        echo ""
    fi
    
    if [[ "$(config_get PACK_RUN)" == "__TRUE__" ]]; then
        echo "=== PACK RUN MODE ==="
        echo "Only downloading packages, skipping installation."
        echo ""
    fi
    
    # Show version configuration
    version_show_current
    
    # Show summary
    ui_show_summary
    
    # Validate version configuration
    if ! version_validate_config; then
        echo ""
        echo "Warning: Some version configuration issues were detected."
        echo "Please review the warnings above."
        echo ""
    fi
    
    # Run configuration validation unless skipped
    if ! should_skip_validation; then
        echo ""
        if ! validate_configuration; then
            echo ""
            report_error $LINENO "Configuration validation failed with errors" "CONFIG_ERROR"
            echo "Please fix the configuration issues above and try again."
            echo "You can skip validation with --skip-system-checks if needed."
            exit 1
        fi
    fi
    
    # Skip user confirmation - proceed directly with installation
    
    # Export configuration to environment variables
    config_export_to_env
    
    # Export version configuration for stage scripts
    package_export_version_config
    
    # Preliminaries
    mkdir -p ${INSTALLDIR}

    # Start writing setup file
    cat << EOF > "$SETUPFILE"
#!/bin/bash
source "${SCRIPTDIR}/tool_kit.sh"
export ABACUS_TOOLCHAIN_OPTIONS="${TOOLCHAIN_OPTIONS}"
EOF

    echo "Compiling with $(get_nprocs) processes for target ${TARGET_CPU}."

    write_toolchain_env ${INSTALLDIR}

    # write toolchain config
    echo "tool_list=\"${tool_list}\"" > ${INSTALLDIR}/toolchain.conf
    for ii in ${package_list}; do
      install_mode="$(eval echo \${with_${ii}})"
      echo "with_${ii}=\"${install_mode}\"" >> ${INSTALLDIR}/toolchain.conf
    done
    
    # Exit if dry run
    if [[ "$(config_get dry_run)" == "__TRUE__" ]]; then
        echo "Wrote only configuration files (--dry-run)."
        exit 0
    fi
    
    # Build packages (following original toolchain logic)
    echo "# Leak suppressions" > ${INSTALLDIR}/lsan.supp
    ./scripts/stage0/install_stage0.sh
    ./scripts/stage1/install_stage1.sh
    ./scripts/stage2/install_stage2.sh
    ./scripts/stage3/install_stage3.sh
    ./scripts/stage4/install_stage4.sh

cat << EOF
========================== usage =========================
Done!
To use the installed tools and libraries and ABACUS version
compiled with it you will first need to execute at the prompt:
  source ${SETUPFILE}
To build ABACUS by gnu-toolchain, just use:
    ./build_abacus_gnu.sh
To build ABACUS by intel-toolchain, just use:
    ./build_abacus_intel.sh
To build ABACUS by amd-toolchain in gcc-aocl, just use:
    ./build_abacus_gcc-aocl.sh
To build ABACUS by amd-toolchain in aocc-aocl, just use:
    ./build_abacus_aocc-aocl.sh
or you can modify the builder scripts to suit your needs.
EOF
    
    return 0
}

# Run main function with all arguments
main "$@"