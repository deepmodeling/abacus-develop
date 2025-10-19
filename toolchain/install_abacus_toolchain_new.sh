#!/bin/bash -e

# ============================================================================
# ABACUS Toolchain Installation Script (New Version - Enhanced UI)
# ============================================================================
# This is the new refactored version of the ABACUS toolchain installation script.
# It provides the same functionality as the original script but with improved
# modularity, maintainability, extensibility, and beautiful terminal output.
#
# Author: Quantum Misaka by Trae SOLO
# Date: 2025-10-16
# Enhanced: Beautiful terminal UI with harmonious color scheme
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

# Initialize modules
version_helper_init
ui_init
ui_setup_signals

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
    
    # Show beautiful welcome banner (only when actually running installation)
    ui_welcome_banner
    
    # Show configuration summary with enhanced styling
    if [[ "$(config_get dry_run)" == "__TRUE__" ]]; then
        ui_section "Dry Run Mode"
        ui_warning "Configuration will be written but no packages will be installed"
        echo ""
    fi
    
    if [[ "$(config_get PACK_RUN)" == "__TRUE__" ]]; then
        ui_section "Pack Run Mode"
        ui_info "Only downloading packages, skipping installation"
        echo ""
    fi
    
    # Show version configuration
    version_show_current
    
    # Show enhanced summary
    ui_show_summary
    
    # Validate version configuration
    if ! version_validate_config; then
        echo ""
        ui_warning "Some version configuration issues were detected"
        ui_info "Please review the warnings above"
        echo ""
    fi
    
    # Run configuration validation unless skipped
    if ! should_skip_validation; then
        echo ""
        ui_section "System Validation"
        if ! validate_configuration; then
            echo ""
            report_error $LINENO "Configuration validation failed with errors" "CONFIG_ERROR"
            ui_error "Please fix the configuration issues above and try again"
            ui_info "You can skip validation with --skip-system-checks if needed"
            exit 1
        fi
        ui_success "System validation completed successfully"
    fi
    
    # Skip user confirmation - proceed directly with installation
    
    # Export configuration to environment variables
    config_export_to_env
    
    # Export version configuration for stage scripts
    package_export_version_config
    
    # Preliminaries
    ui_section "Preparing Installation Environment"
    ui_info "Creating installation directories..."
    mkdir -p ${INSTALLDIR}

    # Start writing setup file
    ui_info "Generating setup configuration..."
    cat << EOF > "$SETUPFILE"
#!/bin/bash
source "${SCRIPTDIR}/tool_kit.sh"
export ABACUS_TOOLCHAIN_OPTIONS="${TOOLCHAIN_OPTIONS}"
EOF

    ui_info "Compiling with $(get_nprocs) processes for target ${TARGET_CPU}"

    write_toolchain_env ${INSTALLDIR}

    # write toolchain config
    ui_info "Writing toolchain configuration..."
    echo "tool_list=\"${tool_list}\"" > ${INSTALLDIR}/toolchain.conf
    for ii in ${package_list}; do
      install_mode="$(eval echo \${with_${ii}})"
      echo "with_${ii}=\"${install_mode}\"" >> ${INSTALLDIR}/toolchain.conf
    done
    
    # Exit if dry run
    if [[ "$(config_get dry_run)" == "__TRUE__" ]]; then
        ui_success "Configuration files generated successfully (dry-run mode)"
        ui_info "To proceed with actual installation, run without --dry-run"
        exit 0
    fi
    
    # Build packages (following original toolchain logic)
    ui_section "Starting Package Installation"
    ui_info "Initializing build environment..."
    echo "# Leak suppressions" > ${INSTALLDIR}/lsan.supp
    
    # Install packages stage by stage with enhanced progress display
    ui_stage_progress "0" "Build Tools & Compilers"
    ./scripts/stage0/install_stage0.sh
    
    ui_stage_progress "1" "MPI Implementations"
    ./scripts/stage1/install_stage1.sh
    
    ui_stage_progress "2" "Mathematical Libraries"
    ./scripts/stage2/install_stage2.sh
    
    ui_stage_progress "3" "Scientific Computing Libraries"
    ./scripts/stage3/install_stage3.sh
    
    ui_stage_progress "4" "Advanced Feature Libraries"
    ./scripts/stage4/install_stage4.sh

    # Show beautiful completion message
    ui_section "Installation Complete"
    ui_success "ABACUS Toolchain installation completed successfully!"
    echo ""
    
    # Enhanced usage instructions with beautiful formatting
    ui_print_color "${UI_BLUE}${UI_BOLD}" "🎉 Congratulations! Your ABACUS toolchain installation finished successfully."
    echo ""
    ui_print_color "${UI_CYAN}${UI_BOLD}" "$(printf '═%.0s' $(seq 1 60))"
    ui_print_color "${UI_BLUE}${UI_BOLD}" "                    USAGE INSTRUCTIONS"
    ui_print_color "${UI_CYAN}${UI_BOLD}" "$(printf '═%.0s' $(seq 1 60))"
    echo ""
    
    ui_print_color "${UI_GREEN}${UI_BOLD}" "🔧 To activate the toolchain environment:"
    ui_print_color "${UI_WHITE}" "   source ${SETUPFILE}"
    echo ""
    
    ui_print_color "${UI_GREEN}${UI_BOLD}" "🚀 Build ABACUS with your preferred toolchain:"
    ui_print_color "${UI_WHITE}" "   ./build_abacus_gnu.sh        # GNU toolchain (GCC + OpenBLAS)"
    ui_print_color "${UI_WHITE}" "   ./build_abacus_intel.sh      # Intel toolchain (ICC + MKL)"
    ui_print_color "${UI_WHITE}" "   ./build_abacus_gcc-aocl.sh   # AMD GCC + AOCL"
    ui_print_color "${UI_WHITE}" "   ./build_abacus_aocc-aocl.sh  # AMD AOCC + AOCL"
    echo ""
    
    ui_print_color "${UI_ORANGE}${UI_BOLD}" "💡 Pro Tips:"
    ui_print_color "${UI_GRAY}" "   • Modify the builder scripts to suit your specific needs"
    ui_print_color "${UI_GRAY}" "   • The environment setup is required for each new shell session"
    ui_print_color "${UI_GRAY}" "   • Check the generated setup file for all available environment variables"
    echo ""
    
    ui_print_color "${UI_CYAN}${UI_BOLD}" "$(printf '═%.0s' $(seq 1 60))"
    ui_print_color "${UI_GREEN}${UI_BOLD}" "✨ Happy DFT calculation with ABACUS! ✨"
    ui_print_color "${UI_CYAN}${UI_BOLD}" "$(printf '═%.0s' $(seq 1 60))"
    echo ""
    
    return 0
}

# Run main function with all arguments
main "$@"