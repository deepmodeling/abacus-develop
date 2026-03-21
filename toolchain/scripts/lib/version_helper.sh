#!/bin/bash

# ============================================================================
# ABACUS Toolchain Version Helper
# ============================================================================
# Provides version display, help information, and user interaction for
# dual-version switching functionality
# Author: Quantum Misaka by Trae SOLO
# Date: 2025-10-16
# ============================================================================

# Source required modules
if [[ -f "${SCRIPTDIR}/package_versions.sh" ]]; then
    source "${SCRIPTDIR}/package_versions.sh"
else
    echo "Error: package_versions.sh not found in ${SCRIPTDIR}" >&2
    return 1
fi

# Display available package versions
# Usage: version_show_available [package_name]
version_show_available() {
    local package_name="${1:-all}"
    
    echo "Available Package Versions:"
    echo "=========================="
    echo
    
    if [[ "${package_name}" == "all" ]]; then
        # Show all packages
        for pkg in gcc cmake openmpi mpich openblas elpa fftw libxc scalapack libtorch libnpy; do
            version_show_package_info "${pkg}"
        done
    else
        version_show_package_info "${package_name}"
    fi
}

# Helper function to show package information
version_show_package_info() {
    local pkg="$1"
    local main_var="${pkg}_main_ver"
    local alt_var="${pkg}_alt_ver"
    
    # Get values using indirect variable expansion
    local main_ver="${!main_var}"
    local alt_ver="${!alt_var}"
    
    if [[ -n "$main_ver" || -n "$alt_ver" ]]; then
        printf "%-12s: " "$pkg"
        
        if [[ -n "$main_ver" ]]; then
            printf "%s (main)" "$main_ver"
        fi
        
        if [[ -n "$alt_ver" ]]; then
            if [[ -n "$main_ver" ]]; then
                printf ", %s (alt)" "$alt_ver"
            else
                printf "%s (alt)" "$alt_ver"
            fi
        fi
        
        echo ""
    fi
}

# Display version switching help
# Usage: version_show_help
version_show_help() {
    cat << 'EOF'
ABACUS Toolchain Version Switching Help
=======================================

Command Line Options:
  --use-alt-versions          Use alternative versions for all packages
  --package-version PKG:VER   Set specific package version (PKG:main|alt)
                              Supports two writing styles:
                              1. Multiple independent parameters: --package-version PKG1:VER1 --package-version PKG2:VER2
                              2. Single parameter with multiple key-value pairs: --package-version PKG1:VER1 PKG2:VER2 ...
  --config-file FILE          Load configuration from file
  --version-info [PKG]        Show version info for package (or all)

Configuration File Support:
  VERSION_STRATEGY=main|alt   Global version strategy
  PACKAGE_VERSION_PKG=main|alt  Specific package version

Examples:
  # Use alternative versions for all packages
  ./install_abacus_toolchain_new.sh --use-alt-versions
  
  # Use specific package versions (single)
  ./install_abacus_toolchain_new.sh --package-version openmpi:alt
  
  # Use specific package versions (multiple independent parameters)
  ./install_abacus_toolchain_new.sh --package-version openmpi:alt --package-version gcc:main --package-version cmake:alt
  
  # Use specific package versions (single parameter with multiple key-value pairs)
  ./install_abacus_toolchain_new.sh --package-version openmpi:alt gcc:main cmake:alt
  
  # Mixed with other options
  ./install_abacus_toolchain_new.sh --package-version elpa:alt cmake:alt --dry-run
  
  # Show version information
  ./install_abacus_toolchain_new.sh --version-info openmpi
  ./install_abacus_toolchain_new.sh --version-info

Backward Compatibility:
  --with-4th-openmpi         Legacy flag (use --package-version openmpi=alt)
  OPENMPI_4TH=yes           Legacy variable (use PACKAGE_VERSION_OPENMPI=alt)

Supported Packages:
  gcc, cmake, openmpi, mpich, openblas, elpa, fftw, libxc, scalapack, libtorch, libnpy

EOF
}

# Interactive version selection interface
# Usage: version_interactive_select
version_interactive_select() {
    echo "Interactive Version Selection"
    echo "============================"
    echo ""
    
    local packages=(
        "openmpi" "mpich" "openblas" "elpa" "fftw" 
        "libxc" "scalapack" "libtorch" "libnpy"
    )
    
    echo "Select versions for packages (press Enter for default 'main'):"
    echo ""
    
    for pkg in "${packages[@]}"; do
        local main_var="${pkg^^}_main_ver"
        local alt_var="${pkg^^}_alt_ver"
        
        # Skip if no versions available
        if [[ -z "${!main_var}" && -z "${!alt_var}" ]]; then
            continue
        fi
        
        printf "%-12s: " "$pkg"
        
        # Show available versions
        if [[ -n "${!main_var}" ]]; then
            printf "main(%s)" "${!main_var}"
        fi
        if [[ -n "${!alt_var}" ]]; then
            if [[ -n "${!main_var}" ]]; then
                printf ", alt(%s)" "${!alt_var}"
            else
                printf "alt(%s)" "${!alt_var}"
            fi
        fi
        
        printf " [main]: "
        read -r user_choice
        
        # Set user choice or default to main
        case "$user_choice" in
            alt|alternative)
                if [[ -n "${!alt_var}" ]]; then
                    CONFIG_CACHE["PACKAGE_VERSION_${pkg^^}"]="alt"
                    echo "  → Selected: alt (${!alt_var})"
                else
                    echo "  → Alternative version not available, using main"
                    CONFIG_CACHE["PACKAGE_VERSION_${pkg^^}"]="main"
                fi
                ;;
            main|"")
                CONFIG_CACHE["PACKAGE_VERSION_${pkg^^}"]="main"
                if [[ -n "${!main_var}" ]]; then
                    echo "  → Selected: main (${!main_var})"
                else
                    echo "  → Main version not available"
                fi
                ;;
            *)
                echo "  → Invalid choice, using main"
                CONFIG_CACHE["PACKAGE_VERSION_${pkg^^}"]="main"
                ;;
        esac
        echo ""
    done
    
    echo "Version selection completed!"
    echo ""
}

# Get effective version for a package
# Usage: version_get_effective PACKAGE_NAME
version_get_effective() {
    local package_name="$1"
    local pkg_upper="${package_name^^}"
    
    # Check for specific package version setting
    local specific_version="${CONFIG_CACHE[PACKAGE_VERSION_${pkg_upper}]}"
    if [[ -n "$specific_version" ]]; then
        echo "$specific_version"
        return 0
    fi
    
    # Check for global version strategy
    local global_strategy="${CONFIG_CACHE[VERSION_STRATEGY]}"
    if [[ -n "$global_strategy" ]]; then
        echo "$global_strategy"
        return 0
    fi
    
    # Check for legacy OPENMPI_4TH support
    if [[ "$package_name" == "openmpi" && "${CONFIG_CACHE[OPENMPI_4TH]}" == "yes" ]]; then
        echo "alt"
        return 0
    fi
    
    # Default to main
    echo "main"
    return 0
}

# Load package variables with version selection
# Usage: version_load_package_vars PACKAGE_NAME
version_load_package_vars() {
    local package_name="$1"
    local effective_version
    
    effective_version=$(version_get_effective "$package_name")
    
    # Convert effective version to suffix for load_package_vars
    local version_suffix=""
    if [[ "$effective_version" == "alt" ]]; then
        version_suffix="alt"
    fi
    
    # Call the original load_package_vars function
    load_package_vars "$package_name" "$version_suffix"
    
    # Export version information for debugging
    export "${package_name^^}_EFFECTIVE_VERSION"="$effective_version"
}

# Display current version configuration
# Usage: version_show_current
version_show_current() {
    echo "Current Version Configuration:"
    echo "============================="
    
    # Show global strategy
    local global_strategy="${CONFIG_CACHE[VERSION_STRATEGY]}"
    if [[ -n "$global_strategy" ]]; then
        echo "Global Strategy: $global_strategy"
    else
        echo "Global Strategy: main (default)"
    fi
    echo ""
    
    # Show per-package settings
    echo "Package-specific settings:"
    local packages=(
        "gcc" "cmake" "openmpi" "mpich" "openblas" 
        "elpa" "fftw" "libxc" "scalapack" "libtorch" "libnpy"
    )
    
    local has_specific=false
    for pkg in "${packages[@]}"; do
        local pkg_upper="${pkg^^}"
        local specific_version="${CONFIG_CACHE[PACKAGE_VERSION_${pkg_upper}]}"
        
        if [[ -n "$specific_version" ]]; then
            printf "  %-12s: %s\n" "$pkg" "$specific_version"
            has_specific=true
        fi
    done
    
    if [[ "$has_specific" == false ]]; then
        echo "  (none - using global strategy)"
    fi
    
    echo ""
    
    # Show legacy compatibility
    if [[ "${CONFIG_CACHE[OPENMPI_4TH]}" == "yes" ]]; then
        echo "Legacy compatibility: OPENMPI_4TH=yes (equivalent to openmpi:alt)"
        echo ""
    fi
}

# Validate version configuration
# Usage: version_validate_config
version_validate_config() {
    local validation_errors=0
    
    # Check if requested versions are available
    for key in "${!CONFIG_CACHE[@]}"; do
        if [[ "$key" =~ ^PACKAGE_VERSION_(.+)$ ]]; then
            local pkg_name="${BASH_REMATCH[1]}"
            local requested_version="${CONFIG_CACHE[$key]}"
            local pkg_lower=$(echo "$pkg_name" | tr '[:upper:]' '[:lower:]')
            
            # Check if the requested version exists
            local version_var="${pkg_lower}_${requested_version}_ver"
            if [[ -z "${!version_var}" ]]; then
                echo "Warning: Requested version '$requested_version' for package '$pkg_lower' is not available"
                validation_errors=$((validation_errors + 1))
            fi
        fi
    done
    
    return $validation_errors
}

# Initialize version helper
# Usage: version_helper_init
version_helper_init() {
    # Set default version strategy if not set
    if [[ -z "${CONFIG_CACHE[VERSION_STRATEGY]}" ]]; then
        CONFIG_CACHE["VERSION_STRATEGY"]="main"
    fi
    
    # Handle legacy OPENMPI_4TH environment variable
    if [[ "${OPENMPI_4TH}" == "yes" && -z "${CONFIG_CACHE[PACKAGE_VERSION_OPENMPI]}" ]]; then
        CONFIG_CACHE["PACKAGE_VERSION_OPENMPI"]="alt"
        CONFIG_CACHE["OPENMPI_4TH"]="yes"
        echo "Notice: Legacy OPENMPI_4TH=yes detected, using openmpi:alt"
    fi
    
    return 0
}

# Export functions for use by other modules
export -f version_show_available
export -f version_show_help
export -f version_interactive_select
export -f version_get_effective
export -f version_load_package_vars
export -f version_show_current
export -f version_validate_config
export -f version_helper_init