#!/bin/bash

# ============================================================================
# ABACUS Toolchain Configuration Validator
# ============================================================================
# Validates configuration settings and detects potential conflicts
# Author: Quantum Misaka by Trae SOLO
# Date: 2025-10-16
# ============================================================================

# Global validation state
VALIDATION_ERRORS=()
VALIDATION_WARNINGS=()
VALIDATION_ERROR_GROUPS=()
VALIDATION_WARNING_GROUPS=()
VALIDATION_INITIALIZED=false

# Initialize configuration validator
# Usage: config_validator_init
config_validator_init() {
    if [[ "$VALIDATION_INITIALIZED" == "true" ]]; then
        return 0
    fi
    
    VALIDATION_ERRORS=()
    VALIDATION_WARNINGS=()
    VALIDATION_ERROR_GROUPS=()
    VALIDATION_WARNING_GROUPS=()
    VALIDATION_INITIALIZED=true
    
    return 0
}

# Add validation error
# Usage: add_validation_error "error message"
add_validation_error() {
    local message="$1"
    VALIDATION_ERRORS+=("ERROR: $message")
}

# Add validation warning
# Usage: add_validation_warning "warning message"
add_validation_warning() {
    local message="$1"
    VALIDATION_WARNINGS+=("WARNING: $message")
}

# Start a new error group (represents one logical error)
# Usage: start_error_group "group_name"
start_error_group() {
    local group_name="$1"
    VALIDATION_ERROR_GROUPS+=("$group_name")
}

# Start a new warning group (represents one logical warning)
# Usage: start_warning_group "group_name"
start_warning_group() {
    local group_name="$1"
    VALIDATION_WARNING_GROUPS+=("$group_name")
}

# Add validation error group (complete error with all details)
# Usage: add_validation_error_group "group_name" "line1" "line2" ...
add_validation_error_group() {
    local group_name="$1"
    shift
    
    # Add the group to track unique errors
    VALIDATION_ERROR_GROUPS+=("$group_name")
    
    # Add all error lines
    for line in "$@"; do
        VALIDATION_ERRORS+=("ERROR: $line")
    done
}

# Add validation warning group (complete warning with all details)
# Usage: add_validation_warning_group "group_name" "line1" "line2" ...
add_validation_warning_group() {
    local group_name="$1"
    shift
    
    # Add the group to track unique warnings
    VALIDATION_WARNING_GROUPS+=("$group_name")
    
    # Add all warning lines
    for line in "$@"; do
        VALIDATION_WARNINGS+=("WARNING: $line")
    done
}

# Add validation info (for successful validations)
# Usage: add_validation_info "info message"
add_validation_info() {
    local message="$1"
    # Info messages are displayed immediately during validation
    echo "INFO: $message"
}

# Check for conflicting math libraries
# Usage: validate_math_libraries
validate_math_libraries() {
    local math_libs_enabled=0
    local enabled_libs=()
    
    # Check which math libraries are enabled
    for lib in mkl aocl openblas; do
        if [[ "${CONFIG_CACHE[with_${lib}]}" == "__INSTALL__" || "${CONFIG_CACHE[with_${lib}]}" == "__SYSTEM__" ]]; then
            math_libs_enabled=$((math_libs_enabled + 1))
            enabled_libs+=("$lib")
        fi
    done
    
    # Validate math library configuration
    if [[ $math_libs_enabled -gt 1 ]]; then
        add_validation_error "Multiple math libraries enabled: ${enabled_libs[*]}. Only one should be active."
    elif [[ $math_libs_enabled -eq 0 ]]; then
        add_validation_warning "No math library enabled. This may cause compilation issues."
    fi
    
    # Check MKL-specific conflicts
    if [[ "${CONFIG_CACHE[with_mkl]}" == "__SYSTEM__" || "${CONFIG_CACHE[with_mkl]}" == "__INSTALL__" ]]; then
        if [[ "${CONFIG_CACHE[with_fftw]}" == "__INSTALL__" ]]; then
            add_validation_warning "MKL includes FFTW interface. Consider setting FFTW to __DONTUSE__ or __SYSTEM__."
        fi
        if [[ "${CONFIG_CACHE[with_scalapack]}" == "__INSTALL__" ]]; then
            add_validation_warning "MKL includes ScaLAPACK. Consider setting ScaLAPACK to __DONTUSE__ or __SYSTEM__."
        fi
    fi
}

# Check for conflicting MPI implementations
# Usage: validate_mpi_implementations
validate_mpi_implementations() {
    local mpi_libs_enabled=0
    local enabled_mpis=()
    
    # Check which MPI implementations are enabled
    for mpi in mpich openmpi intelmpi; do
        if [[ "${CONFIG_CACHE[with_${mpi}]}" == "__INSTALL__" || "${CONFIG_CACHE[with_${mpi}]}" == "__SYSTEM__" ]]; then
            mpi_libs_enabled=$((mpi_libs_enabled + 1))
            enabled_mpis+=("$mpi")
        fi
    done
    
    # Validate MPI configuration
    if [[ $mpi_libs_enabled -gt 1 ]]; then
        add_validation_error "Multiple MPI implementations enabled: ${enabled_mpis[*]}. Only one should be active."
    elif [[ $mpi_libs_enabled -eq 0 ]]; then
        add_validation_warning "No MPI implementation enabled. This may limit parallel functionality."
    fi
}

# Check for compiler consistency
# Usage: validate_compiler_consistency
validate_compiler_consistency() {
    local compilers_enabled=0
    local enabled_compilers=()
    
    # Check which compilers are enabled
    for compiler in gcc intel amd; do
        if [[ "${CONFIG_CACHE[with_${compiler}]}" == "__INSTALL__" || "${CONFIG_CACHE[with_${compiler}]}" == "__SYSTEM__" ]]; then
            compilers_enabled=$((compilers_enabled + 1))
            enabled_compilers+=("$compiler")
        fi
    done
    
    # Validate compiler configuration
    if [[ $compilers_enabled -gt 1 ]]; then
        add_validation_warning "Multiple compilers enabled: ${enabled_compilers[*]}. This may cause compatibility issues."
    elif [[ $compilers_enabled -eq 0 ]]; then
        add_validation_error "No compiler enabled. At least one compiler must be available."
    fi
    
    # Check Intel-specific dependencies
    if [[ "${CONFIG_CACHE[with_intel]}" == "__SYSTEM__" || "${CONFIG_CACHE[with_intel]}" == "__INSTALL__" ]]; then
        if [[ "${CONFIG_CACHE[with_mkl]}" == "__DONTUSE__" ]]; then
            add_validation_warning "Intel compiler is enabled but MKL is disabled. Consider enabling MKL for optimal performance."
        fi
        if [[ "${CONFIG_CACHE[with_intelmpi]}" == "__DONTUSE__" ]]; then
            add_validation_warning "Intel compiler is enabled but Intel MPI is disabled. Consider using Intel MPI for better integration."
        fi
    fi
}

# Check system requirements
# Usage: validate_system_requirements
validate_system_requirements() {
    # ABACUS itself and some dependencies require cmake.
    # Check cmake requirement based on configuration (mirrors original L768-772)
    if [[ "${CONFIG_CACHE[with_cmake]}" == "__DONTUSE__" ]]; then
        add_validation_error_group "cmake_required" \
            "CMake is required for ABACUS and some dependencies. Please enable it." \
            "" \
            "SOLUTION: Use '--with-cmake=install' to automatically install CMake:" \
            "  ./install_abacus_toolchain_new.sh --with-cmake=install [other options]"
        return
    fi
    
    # Check for required system tools (excluding cmake which is handled above)
    local required_tools=("make" "git" "wget" "tar" "gzip")
    local missing_tools=()
    
    # Only check for system cmake if using system cmake
    if [[ "${CONFIG_CACHE[with_cmake]}" == "__SYSTEM__" ]]; then
        if ! command -v "cmake" &> /dev/null; then
            missing_tools+=("cmake")
        else
            # Check cmake version if it exists
            local cmake_min_version="3.16"
            local cmake_version=$(cmake --version 2>/dev/null | head -n 1 | awk '{print $3}')
            
            if [[ -z "$cmake_version" ]]; then
                add_validation_error "Failed to determine CMake version"
                add_validation_error "This usually indicates a corrupted or non-standard CMake installation."
                add_validation_error ""
                add_validation_error "SOLUTION: Use '--with-cmake=install' to install a known-good CMake version:"
                add_validation_error "  ./install_abacus_toolchain_new.sh --with-cmake=install [other options]"
                add_validation_error ""
                add_validation_error "This will install CMake ${cmake_ver:-3.29.6} with proper version information."
                return
            fi
            
            # Extract major and minor version numbers for comparison
            local cmake_major=$(echo "$cmake_version" | awk -F. '{print $1}')
            local cmake_minor=$(echo "$cmake_version" | awk -F. '{print $2}')
            local min_major=$(echo "$cmake_min_version" | awk -F. '{print $1}')
            local min_minor=$(echo "$cmake_min_version" | awk -F. '{print $2}')
            
            # Validate version format
            if ! [[ "$cmake_major" =~ ^[0-9]+$ ]] || ! [[ "$cmake_minor" =~ ^[0-9]+$ ]]; then
                add_validation_error "Unable to parse CMake version from: $cmake_version"
                add_validation_error "Expected format: X.Y.Z (e.g., 3.29.6), but got: $cmake_version"
                add_validation_error ""
                add_validation_error "SOLUTION: Use '--with-cmake=install' to install a standard CMake version:"
                add_validation_error "  ./install_abacus_toolchain_new.sh --with-cmake=install [other options]"
                add_validation_error ""
                add_validation_error "This will install CMake ${cmake_ver:-3.29.6} with standard version format."
                return
            fi
            
            # Compare versions (major.minor comparison)
            local version_too_old=false
            if [[ "$cmake_major" -lt "$min_major" ]]; then
                version_too_old=true
            elif [[ "$cmake_major" -eq "$min_major" ]] && [[ "$cmake_minor" -lt "$min_minor" ]]; then
                version_too_old=true
            fi
            
            if [[ "$version_too_old" == "true" ]]; then
                add_validation_error "CMake version $cmake_version is too old"
                add_validation_error "Minimum required: CMake $cmake_min_version or newer (ABACUS requirement)"
                add_validation_error ""
                add_validation_error_group "cmake_outdated" \
                    "Your system CMake is outdated and may cause build failures." \
                    "" \
                    "SOLUTION: Use '--with-cmake=install' to install a modern CMake version:" \
                    "  ./install_abacus_toolchain_new.sh --with-cmake=install [other options]" \
                    "" \
                    "This will install CMake ${cmake_ver:-3.29.6} (>= $cmake_min_version) with full feature support."
                return
            fi
            
            # Success - add informational message
            add_validation_info "System CMake validated: version $cmake_version (>= $cmake_min_version required)"
        fi
    fi
    
    for tool in "${required_tools[@]}"; do
        if ! command -v "$tool" &> /dev/null; then
            missing_tools+=("$tool")
        fi
    done
    
    if [[ ${#missing_tools[@]} -gt 0 ]]; then
        # Check if cmake is in the missing tools and provide specific guidance
        if [[ " ${missing_tools[*]} " =~ " cmake " ]]; then
            add_validation_error_group "cmake_missing" \
                "CMake is not installed on your system" \
                "CMake is required for building ABACUS and its dependencies." \
                "" \
                "SOLUTION: Use '--with-cmake=install' to automatically install CMake:" \
                "  ./install_abacus_toolchain_new.sh --with-cmake=install [other options]" \
                "" \
                "This will install CMake ${cmake_ver:-3.29.6} with full feature support."
            
            # Remove cmake from missing_tools to avoid duplicate error messages
            local filtered_tools=()
            for tool in "${missing_tools[@]}"; do
                if [[ "$tool" != "cmake" ]]; then
                    filtered_tools+=("$tool")
                fi
            done
            missing_tools=("${filtered_tools[@]}")
        fi
        
        # Report other missing tools if any
        if [[ ${#missing_tools[@]} -gt 0 ]]; then
            add_validation_error "Missing required system tools: ${missing_tools[*]}"
        fi
    fi
    
    # Check for development packages (common names)
    local dev_packages=("build-essential" "gcc" "g++" "gfortran")
    local missing_dev=()
    
    for pkg in "${dev_packages[@]}"; do
        if ! command -v "$pkg" &> /dev/null && ! dpkg -l | grep -q "$pkg" 2>/dev/null; then
            missing_dev+=("$pkg")
        fi
    done
    
    if [[ ${#missing_dev[@]} -gt 0 ]]; then
        add_validation_warning "Potentially missing development packages: ${missing_dev[*]}"
    fi
    
    # Check GCC version consistency and minimum requirements
    # This mirrors the logic from original install_abacus_toolchain.sh L636-667
    if [[ "${CONFIG_CACHE[with_gcc]}" != "__INSTALL__" ]]; then
        local gcc_min_version=5
        
        # Check if GCC tools are available
        if ! command -v gcc &> /dev/null || ! command -v g++ &> /dev/null || ! command -v gfortran &> /dev/null; then
            local missing_tools=()
            ! command -v gcc &> /dev/null && missing_tools+=("gcc")
            ! command -v g++ &> /dev/null && missing_tools+=("g++")
            ! command -v gfortran &> /dev/null && missing_tools+=("gfortran")
            
            add_validation_error_group "gcc_incomplete" \
                "System GCC toolchain incomplete. Missing: ${missing_tools[*]}" \
                "" \
                "SOLUTION: Use '--with-gcc=install' to automatically download and install a complete GCC toolchain:" \
                "  ./install_abacus_toolchain_new.sh --with-gcc=install [other options]" \
                "" \
                "This will install GCC ${gcc_ver:-13.2.0} with all required components (gcc, g++, gfortran)."
            return
        fi
        
        # Get versions of GCC components
        local gcc_version=$(gcc --version 2>/dev/null | head -n 1 | awk '{print $NF}')
        local gxx_version=$(g++ --version 2>/dev/null | head -n 1 | awk '{print $NF}')
        local gfc_version=$(gfortran --version 2>/dev/null | head -n 1 | awk '{print $NF}')
        
        # Check if version extraction was successful
        if [[ -z "$gcc_version" || -z "$gxx_version" || -z "$gfc_version" ]]; then
            local failed_tools=()
            [[ -z "$gcc_version" ]] && failed_tools+=("gcc")
            [[ -z "$gxx_version" ]] && failed_tools+=("g++")
            [[ -z "$gfc_version" ]] && failed_tools+=("gfortran")
            
            add_validation_error_group "gcc_version_failed" \
                "Failed to determine GCC toolchain versions for: ${failed_tools[*]}" \
                "This usually indicates corrupted or non-standard GCC installation." \
                "" \
                "SOLUTION: Use '--with-gcc=install' to install a known-good GCC version:" \
                "  ./install_abacus_toolchain_new.sh --with-gcc=install [other options]" \
                "" \
                "This will install GCC ${gcc_ver:-13.2.0} with proper version information."
            return
        fi
        
        # Check version consistency
        if [[ "$gcc_version" != "$gxx_version" ]] || [[ "$gcc_version" != "$gfc_version" ]]; then
            add_validation_error_group "gcc_version_inconsistent" \
                "GCC toolchain versions are inconsistent:" \
                "  gcc:      $gcc_version" \
                "  g++:      $gxx_version" \
                "  gfortran: $gfc_version" \
                "" \
                "All GCC components must have the same version for proper compilation." \
                "SOLUTION: Use '--with-gcc=install' to install a consistent GCC toolchain:" \
                "  ./install_abacus_toolchain_new.sh --with-gcc=install [other options]" \
                "" \
                "This will install GCC ${gcc_ver:-13.2.0} with all components at the same version."
            return
        fi
        
        # Extract major version number
        local gcc_major=$(echo "$gcc_version" | awk -F. '{print $1}')
        
        # Validate major version is numeric and meets minimum requirement
        if ! [[ "$gcc_major" =~ ^[0-9]+$ ]]; then
            add_validation_error_group "gcc_version_parse_failed" \
                "Unable to parse GCC major version from: $gcc_version" \
                "Expected format: X.Y.Z (e.g., 13.2.0), but got: $gcc_version" \
                "" \
                "SOLUTION: Use '--with-gcc=install' to install a standard GCC version:" \
                "  ./install_abacus_toolchain_new.sh --with-gcc=install [other options]" \
                "" \
                "This will install GCC ${gcc_ver:-13.2.0} with standard version format."
            return
        fi
        
        if [[ "$gcc_major" -lt "$gcc_min_version" ]]; then
            add_validation_error_group "gcc_version_too_old" \
                "GCC version $gcc_version is too old (major version: $gcc_major)" \
                "Minimum required: GCC $gcc_min_version.x or newer (ABACUS requirement)" \
                "" \
                "Your system GCC is outdated and may cause compilation failures." \
                "" \
                "SOLUTION: Use '--with-gcc=install' to install a modern GCC version:" \
                "  ./install_abacus_toolchain_new.sh --with-gcc=install [other options]" \
                "" \
                "This will install GCC ${gcc_ver:-13.2.0} (>= $gcc_min_version.x) with full C++17/C++20 support."
            return
        fi
        
        # Success - add informational message
        add_validation_info "System GCC toolchain validated: version $gcc_version (>= $gcc_min_version.x required)"
    fi
}

# Check for logical inconsistencies
# Usage: validate_logical_consistency
validate_logical_consistency() {
    # Check if ELPA is enabled without ScaLAPACK
    if [[ "${CONFIG_CACHE[with_elpa]}" == "__INSTALL__" || "${CONFIG_CACHE[with_elpa]}" == "__SYSTEM__" ]]; then
        if [[ "${CONFIG_CACHE[with_scalapack]}" == "__DONTUSE__" ]]; then
            # Check if MKL is enabled (MKL provides ScaLAPACK functionality)
            if [[ "${CONFIG_CACHE[MATH_MODE]}" == "mkl" || "${CONFIG_CACHE[with_mkl]}" == "__SYSTEM__" || "${CONFIG_CACHE[with_mkl]}" == "__INSTALL__" ]]; then
                # MKL provides ScaLAPACK, so this is acceptable
                :
            else
                add_validation_error "ELPA requires ScaLAPACK but ScaLAPACK is disabled."
            fi
        fi
    fi
    
    # Check if ScaLAPACK is enabled without MPI
    if [[ "${CONFIG_CACHE[with_scalapack]}" == "__INSTALL__" || "${CONFIG_CACHE[with_scalapack]}" == "__SYSTEM__" ]]; then
        local mpi_enabled=false
        for mpi in mpich openmpi intelmpi; do
            if [[ "${CONFIG_CACHE[with_${mpi}]}" == "__INSTALL__" || "${CONFIG_CACHE[with_${mpi}]}" == "__SYSTEM__" ]]; then
                mpi_enabled=true
                break
            fi
        done
        
        if [[ "$mpi_enabled" == "false" ]]; then
            add_validation_warning "ScaLAPACK is enabled but no MPI implementation is active."
        fi
    fi
    
    # Check for GPU-related inconsistencies
    if [[ "${CONFIG_CACHE[GPUVER]}" != "no" ]]; then
        if [[ "${CONFIG_CACHE[with_elpa]}" == "__INSTALL__" ]]; then
            add_validation_warning "GPU support is enabled. Ensure ELPA is compiled with GPU support."
        fi
    fi
}

# Validate package versions compatibility
# Usage: validate_package_versions
validate_package_versions() {
    # This is a placeholder for version compatibility checks
    # In a real implementation, you would check for known incompatible version combinations
    
    # Example: Check for known problematic combinations
    if [[ "${CONFIG_CACHE[with_gcc]}" == "__INSTALL__" ]]; then
        local gcc_version="${gcc_ver:-unknown}"
        if [[ "$gcc_version" == "unknown" ]]; then
            add_validation_warning "GCC version not specified. Using default version."
        fi
    fi
}

# Run all validation checks
# Usage: validate_configuration
validate_configuration() {
    config_validator_init
    
    echo "Running configuration validation..."
    
    # Run all validation checks
    validate_math_libraries
    validate_mpi_implementations
    validate_compiler_consistency
    validate_system_requirements
    validate_logical_consistency
    validate_package_versions
    
    # Report results using error groups for accurate counting
    local error_groups=${#VALIDATION_ERROR_GROUPS[@]}
    local warning_groups=${#VALIDATION_WARNING_GROUPS[@]}
    local total_issues=$((error_groups + warning_groups))
    
    if [[ ${#VALIDATION_ERRORS[@]} -gt 0 ]]; then
        echo ""
        echo "Configuration Errors Found:"
        echo "=========================="
        for error in "${VALIDATION_ERRORS[@]}"; do
            echo "  $error"
        done
    fi
    
    if [[ ${#VALIDATION_WARNINGS[@]} -gt 0 ]]; then
        echo ""
        echo "Configuration Warnings:"
        echo "======================"
        for warning in "${VALIDATION_WARNINGS[@]}"; do
            echo "  $warning"
        done
    fi
    
    if [[ $total_issues -eq 0 ]]; then
        echo "✓ Configuration validation passed with no issues."
        return 0
    else
        echo ""
        echo "Configuration validation completed with $total_issues issue(s)."
        echo "  Errors: $error_groups"
        echo "  Warnings: $warning_groups"
        
        # Return error code if there are validation errors
        if [[ $error_groups -gt 0 ]]; then
            return 1
        else
            return 0
        fi
    fi
}

# Check if validation should be skipped
# Usage: should_skip_validation
should_skip_validation() {
    if [[ "${CONFIG_CACHE[SKIP_SYSTEM_CHECKS]}" == "true" ]]; then
        echo "Skipping configuration validation (SKIP_SYSTEM_CHECKS=true)"
        return 0
    fi
    return 1
}
