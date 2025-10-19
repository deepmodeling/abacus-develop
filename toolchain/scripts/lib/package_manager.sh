#!/bin/bash

# ============================================================================
# ABACUS Toolchain Package Manager
# ============================================================================
# Handles package downloading, building, and installation
# Author: Quantum Misaka by Trae SOLO
# Date: 2025-10-16
# ============================================================================

# Global package manager variables
PACKAGE_MANAGER_INITIALIZED=false
declare -A PACKAGE_BUILD_STATUS
declare -A PACKAGE_DEPENDENCIES

# Initialize package manager
# Usage: package_manager_init
package_manager_init() {
    if [[ "$PACKAGE_MANAGER_INITIALIZED" == "true" ]]; then
        return 0
    fi
    
    # Initialize dependencies
    if ! command -v config_manager_init &> /dev/null; then
        source "${SCRIPTDIR}/lib/config_manager.sh"
    fi
    
    config_manager_init
    
    # Define package dependencies
    package_manager_define_dependencies
    
    PACKAGE_MANAGER_INITIALIZED=true
    return 0
}

# Define package dependencies
# Usage: package_manager_define_dependencies
package_manager_define_dependencies() {
    # Stage 0: Compilers and Build Tools
    PACKAGE_DEPENDENCIES["gcc"]=""
    PACKAGE_DEPENDENCIES["intel"]=""
    PACKAGE_DEPENDENCIES["amd"]=""
    PACKAGE_DEPENDENCIES["cmake"]=""
    
    # Stage 1: MPI Implementations (depend on compilers)
    PACKAGE_DEPENDENCIES["openmpi"]="gcc cmake"
    PACKAGE_DEPENDENCIES["mpich"]="gcc cmake"
    PACKAGE_DEPENDENCIES["intelmpi"]="intel"
    
    # Stage 2: Mathematical Libraries (depend on compilers and MPI)
    PACKAGE_DEPENDENCIES["openblas"]="gcc cmake"
    PACKAGE_DEPENDENCIES["mkl"]="intel"
    PACKAGE_DEPENDENCIES["aocl"]="amd"
    
    # Stage 3: Scientific Computing Libraries (depend on math libs and MPI)
    PACKAGE_DEPENDENCIES["fftw"]="gcc cmake"
    PACKAGE_DEPENDENCIES["libxc"]="gcc cmake"
    PACKAGE_DEPENDENCIES["scalapack"]="openblas"
    PACKAGE_DEPENDENCIES["elpa"]="scalapack"
    
    # Stage 4: Advanced Libraries
    PACKAGE_DEPENDENCIES["cereal"]="gcc cmake"
    PACKAGE_DEPENDENCIES["rapidjson"]="gcc cmake"
    PACKAGE_DEPENDENCIES["libtorch"]="gcc cmake"
    PACKAGE_DEPENDENCIES["libnpy"]="gcc cmake"
    PACKAGE_DEPENDENCIES["libri"]="gcc cmake"
    PACKAGE_DEPENDENCIES["libcomm"]="gcc cmake"
    PACKAGE_DEPENDENCIES["nep"]="gcc cmake"
}

# Check if package is enabled for installation
# Usage: package_is_enabled "package_name"
package_is_enabled() {
    local package="$1"
    local status=$(config_get "with_${package}")
    
    case "$status" in
        "__INSTALL__"|"__SYSTEM__")
            return 0
            ;;
        *)
            return 1
            ;;
    esac
}

# Get package installation mode
# Usage: package_get_mode "package_name"
package_get_mode() {
    local package="$1"
    config_get "with_${package}"
}


# Mark package as built
# Usage: package_mark_built "package_name"
package_mark_built() {
    local package="$1"
    PACKAGE_BUILD_STATUS["$package"]="built"
}

# Get package dependencies
# Usage: package_get_dependencies "package_name"
package_get_dependencies() {
    local package="$1"
    echo "${PACKAGE_DEPENDENCIES[$package]}"
}

# Check if all dependencies are satisfied
# Usage: package_check_dependencies "package_name"
package_check_dependencies() {
    local package="$1"
    local deps=$(package_get_dependencies "$package")
    
    for dep in $deps; do
        if package_is_enabled "$dep"; then
            if ! package_is_built "$dep"; then
                echo "Dependency $dep for $package is not built yet"
                return 1
            fi
        fi
    done
    
    return 0
}

# Install package by stage
# Usage: package_install_stage "stage_number"
package_install_stage() {
    local stage="$1"
    local stage_script="${SCRIPTDIR}/stage${stage}/install_stage${stage}.sh"
    
    if [[ ! -f "$stage_script" ]]; then
        report_error ${LINENO} "Stage script not found: $stage_script"
        return 1
    fi
    
    echo "Installing Stage ${stage} packages..."
    
    # Check if dry run mode
    if [[ "$(config_get dry_run)" == "__TRUE__" ]]; then
        echo "Dry run: would execute $stage_script"
        return 0
    fi
    
    # Export version configuration for stage scripts
    package_export_version_config
    
    # Execute stage installation script
    if ! "$stage_script"; then
        report_error ${LINENO} "Failed to install Stage ${stage} packages"
        return 1
    fi
    
    echo "Stage ${stage} installation completed successfully"
    return 0
}

# Install all enabled packages
# Usage: package_install_all
package_install_all() {
    echo "Starting package installation process..."
    
    # Check if pack run mode (only check and install required packages)
    if [[ "$(config_get PACK_RUN)" == "__TRUE__" ]]; then
        echo "Pack run mode: checking and installing required packages only"
        package_check_system_requirements
        return $?
    fi
    
    # Install packages by stages
    for stage in 0 1 2 3 4; do
        if ! package_install_stage "$stage"; then
            report_error ${LINENO} "Failed to install stage $stage"
            return 1
        fi
    done
    
    echo "All package installations completed successfully"
    return 0
}

# Check system requirements
# Usage: package_check_system_requirements
package_check_system_requirements() {
    echo "Checking system requirements..."
    
    # Check for essential system tools
    local required_tools="wget curl tar gzip make"
    local missing_tools=""
    
    for tool in $required_tools; do
        if ! command -v "$tool" &> /dev/null; then
            missing_tools="$missing_tools $tool"
        fi
    done
    
    if [[ -n "$missing_tools" ]]; then
        report_error ${LINENO} "Missing required system tools:$missing_tools"
        echo "Please install these tools using your system package manager"
        return 1
    fi
    
    # Check for compiler availability
    local compiler_mode=$(config_get MATH_MODE)
    case "$compiler_mode" in
        gcc)
            if ! command -v gcc &> /dev/null; then
                report_warning ${LINENO} "GCC compiler not found in system PATH"
            fi
            ;;
        intel)
            if ! command -v icc &> /dev/null && ! command -v icx &> /dev/null; then
                report_warning ${LINENO} "Intel compiler not found in system PATH"
            fi
            ;;
        amd)
            if ! command -v clang &> /dev/null; then
                report_warning ${LINENO} "AMD/Clang compiler not found in system PATH"
            fi
            ;;
    esac
    
    echo "System requirements check completed"
    return 0
}

# Generate package configuration file
# Usage: package_write_config "config_file"
package_write_config() {
    local config_file="$1"
    
    echo "# ABACUS Toolchain Package Configuration" > "$config_file"
    echo "# Generated on $(date)" >> "$config_file"
    echo "" >> "$config_file"
    
    # Write package lists
    echo "tool_list=\"${tool_list}\"" >> "$config_file"
    echo "mpi_list=\"${mpi_list}\"" >> "$config_file"
    echo "math_list=\"${math_list}\"" >> "$config_file"
    echo "lib_list=\"${lib_list}\"" >> "$config_file"
    echo "package_list=\"${package_list}\"" >> "$config_file"
    echo "" >> "$config_file"
    
    # Write package configurations
    for pkg in ${package_list}; do
        local status=$(config_get "with_${pkg}")
        echo "with_${pkg}=\"${status}\"" >> "$config_file"
    done
    
    echo "" >> "$config_file"
    echo "# Configuration modes" >> "$config_file"
    echo "MPI_MODE=\"$(config_get MPI_MODE)\"" >> "$config_file"
    echo "MATH_MODE=\"$(config_get MATH_MODE)\"" >> "$config_file"
    echo "TARGET_CPU=\"$(config_get TARGET_CPU)\"" >> "$config_file"
    echo "GPUVER=\"$(config_get GPUVER)\"" >> "$config_file"
    
    echo "Package configuration written to: $config_file"
}

# Get list of packages to install
# Usage: package_list_to_install
package_list_to_install() {
    local install_packages=""
    
    for pkg in ${package_list}; do
        if [[ "$(config_get with_${pkg})" == "__INSTALL__" ]]; then
            install_packages="$install_packages $pkg"
        fi
    done
    
    echo "$install_packages"
}

# Get list of system packages
# Usage: package_get_system_list
package_get_system_list() {
    local system_packages=""
    
    for pkg in ${package_list}; do
        if [[ "$(config_get with_${pkg})" == "__SYSTEM__" ]]; then
            system_packages="$system_packages $pkg"
        fi
    done
    
    echo "$system_packages"
}

# Print package installation summary
# Usage: package_print_summary
package_print_summary() {
    echo "Package Installation Summary:"
    echo "============================"
    echo ""
    
    local install_list=$(package_get_install_list)
    local system_list=$(package_get_system_list)
    
    if [[ -n "$install_list" ]]; then
        echo "Packages to be installed from source:"
        for pkg in $install_list; do
            echo "  - $pkg"
        done
        echo ""
    fi
    
    if [[ -n "$system_list" ]]; then
        echo "Packages to be used from system:"
        for pkg in $system_list; do
            echo "  - $pkg"
        done
        echo ""
    fi
    
    echo "Installation modes:"
    echo "  MPI Mode: $(config_get MPI_MODE)"
    echo "  Math Mode: $(config_get MATH_MODE)"
    echo "  Target CPU: $(config_get TARGET_CPU)"
    echo "  GPU Version: $(config_get GPUVER)"
    echo ""
}

# Validate package configuration
# Usage: package_validate_config
package_validate_config() {
    local errors=0
    
    # Check for conflicting MPI implementations
    local mpi_count=0
    for mpi in mpich openmpi intelmpi; do
        if [[ "$(config_get with_${mpi})" != "__DONTUSE__" ]]; then
            ((mpi_count++))
        fi
    done
    
    if [[ $mpi_count -gt 1 ]]; then
        report_error ${LINENO} "Multiple MPI implementations selected. Please choose only one."
        ((errors++))
    fi
    
    # Check for conflicting math libraries
    local math_count=0
    for math in mkl openblas aocl; do
        if [[ "$(config_get with_${math})" != "__DONTUSE__" ]]; then
            ((math_count++))
        fi
    done
    
    if [[ $math_count -gt 1 ]]; then
        report_warning ${LINENO} "Multiple math libraries selected. This may cause conflicts."
    fi
    
    # Check for GPU support consistency
    local gpu_ver=$(config_get GPUVER)
    if [[ "$gpu_ver" != "no" ]]; then
        if [[ "$(config_get enable_cuda)" != "__TRUE__" && "$(config_get enable_hip)" != "__TRUE__" ]]; then
            report_warning ${LINENO} "GPU version specified but no GPU support enabled"
        fi
    fi
    
    return $errors
}

# Clean build directories
# Usage: package_clean_build
package_clean_build() {
    if [[ -d "$BUILDDIR" ]]; then
        echo "Cleaning build directory: $BUILDDIR"
        rm -rf "$BUILDDIR"
    fi
    mkdir -p "$BUILDDIR"
}

# Clean install directories
# Usage: package_clean_install
package_clean_install() {
    if [[ -d "$INSTALLDIR" ]]; then
        echo "Cleaning install directory: $INSTALLDIR"
        rm -rf "$INSTALLDIR"
    fi
    mkdir -p "$INSTALLDIR"
}

# Install all packages in pack run mode
# Usage: package_install_all_pack_run
package_install_all_pack_run() {
    echo "Pack run mode - installing all packages in sequence"
    
    # Execute all stages in sequence
    for stage in 0 1 2 3 4; do
        local stage_script="${SCRIPTDIR}/stage${stage}/install_stage${stage}.sh"
        if [[ -f "$stage_script" ]]; then
            echo "Executing stage $stage in pack run mode..."
            if ! "$stage_script"; then
                echo "Stage $stage failed in pack run mode"
                return 1
            fi
        fi
    done
    
    return 0
}

# Export version configuration for stage scripts
# Usage: package_export_version_config
package_export_version_config() {
    # Export version suffix configuration
    local use_alt_versions=$(config_get use_alt_versions)
    
    if [[ "$use_alt_versions" == "__TRUE__" ]]; then
        export ABACUS_TOOLCHAIN_VERSION_SUFFIX="alt"
    else
        export ABACUS_TOOLCHAIN_VERSION_SUFFIX="main"
    fi
    
    # Build package versions string from CONFIG_CACHE PACKAGE_VERSION_* entries
    local package_versions=""
    for key in "${!CONFIG_CACHE[@]}"; do
        if [[ "$key" =~ ^PACKAGE_VERSION_(.+)$ ]]; then
            local pkg_name="${BASH_REMATCH[1]}"
            local pkg_version="${CONFIG_CACHE[$key]}"
            local pkg_lower=$(echo "$pkg_name" | tr '[:upper:]' '[:lower:]')
            
            if [[ -n "$package_versions" ]]; then
                package_versions="${package_versions} ${pkg_lower}:${pkg_version}"
            else
                package_versions="${pkg_lower}:${pkg_version}"
            fi
        fi
    done
    
    # Export individual package version overrides
    if [[ -n "$package_versions" ]]; then
        export ABACUS_TOOLCHAIN_PACKAGE_VERSIONS="$package_versions"
    fi
    
    # Export for backward compatibility
    export VERSION_SUFFIX="$ABACUS_TOOLCHAIN_VERSION_SUFFIX"
}

# Execute stage-based installation (stage0 -> stage1 -> stage2 -> stage3 -> stage4)
# This follows the original toolchain logic exactly
package_install_all_stages() {
    echo "Starting stage-based installation..."
    
    # Export version configuration before starting stages
    package_export_version_config
    
    # Execute stages in order, just like the original fallback version
    local stages=("stage0" "stage1" "stage2" "stage3" "stage4")
    
    for stage in "${stages[@]}"; do
        local stage_script="${SCRIPTDIR}/${stage}/install_${stage}.sh"
        
        if [[ -f "$stage_script" ]]; then
            echo "Executing ${stage}..."
            if ! "$stage_script"; then
                report_error ${LINENO} "Failed to execute ${stage}"
                return 1
            fi
        else
            echo "Warning: Stage script not found: $stage_script"
        fi
    done
    
    echo "All stages completed successfully."
    return 0
}