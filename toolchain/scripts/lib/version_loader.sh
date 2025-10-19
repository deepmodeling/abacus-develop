#!/bin/bash

# ============================================================================
# ABACUS Toolchain Version Loader Helper
# ============================================================================
# Provides a unified version loading mechanism for all stage scripts
# Author: Quantum Misaka by Trae SOLO
# Date: 2025-10-16
# ============================================================================

# Load package variables with version suffix support
# Usage: load_package_with_version "package_name"
load_package_with_version() {
    local package_name="$1"
    
    # Check for version configuration from environment or individual package setting
    local version_suffix=""
    
    if [[ -n "${ABACUS_TOOLCHAIN_PACKAGE_VERSIONS}" ]]; then
        # Check for individual package version override
        if echo "${ABACUS_TOOLCHAIN_PACKAGE_VERSIONS}" | grep -q "${package_name}:alt"; then
            version_suffix="alt"
        elif echo "${ABACUS_TOOLCHAIN_PACKAGE_VERSIONS}" | grep -q "${package_name}:main"; then
            version_suffix="main"
        fi
    fi
    
    # Fall back to global version suffix if no individual setting
    if [[ -z "$version_suffix" && -n "${ABACUS_TOOLCHAIN_VERSION_SUFFIX}" ]]; then
        version_suffix="${ABACUS_TOOLCHAIN_VERSION_SUFFIX}"
    fi
    
    # Load package variables with appropriate version
    load_package_vars "$package_name" "$version_suffix"
    
    # Debug output if verbose mode is enabled
    if [[ "${VERBOSE_MODE}" == "__TRUE__" ]]; then
        echo "DEBUG: Loaded $package_name with version suffix: ${version_suffix:-main}"
    fi
}

# Check if version configuration is available
# Usage: has_version_config
has_version_config() {
    [[ -n "${ABACUS_TOOLCHAIN_VERSION_SUFFIX}" || -n "${ABACUS_TOOLCHAIN_PACKAGE_VERSIONS}" ]]
}

# Get effective version suffix for a package
# Usage: get_package_version_suffix "package_name"
get_package_version_suffix() {
    local package_name="$1"
    local version_suffix=""
    
    if [[ -n "${ABACUS_TOOLCHAIN_PACKAGE_VERSIONS}" ]]; then
        # Check for individual package version override
        if echo "${ABACUS_TOOLCHAIN_PACKAGE_VERSIONS}" | grep -q "${package_name}:alt"; then
            version_suffix="alt"
        elif echo "${ABACUS_TOOLCHAIN_PACKAGE_VERSIONS}" | grep -q "${package_name}:main"; then
            version_suffix="main"
        fi
    fi
    
    # Fall back to global version suffix if no individual setting
    if [[ -z "$version_suffix" && -n "${ABACUS_TOOLCHAIN_VERSION_SUFFIX}" ]]; then
        version_suffix="${ABACUS_TOOLCHAIN_VERSION_SUFFIX}"
    fi
    
    echo "${version_suffix:-main}"
}

# Show version information for debugging
# Usage: show_version_debug "package_name"
show_version_debug() {
    local package_name="$1"
    local version_suffix=$(get_package_version_suffix "$package_name")
    
    echo "Version Debug for $package_name:"
    echo "  Global suffix: ${ABACUS_TOOLCHAIN_VERSION_SUFFIX:-none}"
    echo "  Package versions: ${ABACUS_TOOLCHAIN_PACKAGE_VERSIONS:-none}"
    echo "  Effective suffix: $version_suffix"
}