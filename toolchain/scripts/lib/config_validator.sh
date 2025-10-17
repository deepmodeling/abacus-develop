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
VALIDATION_INITIALIZED=false

# Initialize configuration validator
# Usage: config_validator_init
config_validator_init() {
    if [[ "$VALIDATION_INITIALIZED" == "true" ]]; then
        return 0
    fi
    
    VALIDATION_ERRORS=()
    VALIDATION_WARNINGS=()
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
    # Check for required system tools
    local required_tools=("make" "cmake" "git" "wget" "tar" "gzip")
    local missing_tools=()
    
    for tool in "${required_tools[@]}"; do
        if ! command -v "$tool" &> /dev/null; then
            missing_tools+=("$tool")
        fi
    done
    
    if [[ ${#missing_tools[@]} -gt 0 ]]; then
        add_validation_error "Missing required system tools: ${missing_tools[*]}"
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
            
            add_validation_error "System GCC toolchain incomplete. Missing: ${missing_tools[*]}"
            add_validation_error ""
            add_validation_error "SOLUTION: Use '--with-gcc=install' to automatically download and install a complete GCC toolchain:"
            add_validation_error "  ./install_abacus_toolchain_new.sh --with-gcc=install [other options]"
            add_validation_error ""
            add_validation_error "This will install GCC ${gcc_ver:-13.2.0} with all required components (gcc, g++, gfortran)."
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
            
            add_validation_error "Failed to determine GCC toolchain versions for: ${failed_tools[*]}"
            add_validation_error "This usually indicates corrupted or non-standard GCC installation."
            add_validation_error ""
            add_validation_error "SOLUTION: Use '--with-gcc=install' to install a known-good GCC version:"
            add_validation_error "  ./install_abacus_toolchain_new.sh --with-gcc=install [other options]"
            add_validation_error ""
            add_validation_error "This will install GCC ${gcc_ver:-13.2.0} with proper version information."
            return
        fi
        
        # Check version consistency
        if [[ "$gcc_version" != "$gxx_version" ]] || [[ "$gcc_version" != "$gfc_version" ]]; then
            add_validation_error "GCC toolchain versions are inconsistent:"
            add_validation_error "  gcc:      $gcc_version"
            add_validation_error "  g++:      $gxx_version"
            add_validation_error "  gfortran: $gfc_version"
            add_validation_error ""
            add_validation_error "All GCC components must have the same version for proper compilation."
            add_validation_error ""
            add_validation_error "SOLUTION: Use '--with-gcc=install' to install a consistent GCC toolchain:"
            add_validation_error "  ./install_abacus_toolchain_new.sh --with-gcc=install [other options]"
            add_validation_error ""
            add_validation_error "This will install GCC ${gcc_ver:-13.2.0} with all components at the same version."
            return
        fi
        
        # Extract major version number
        local gcc_major=$(echo "$gcc_version" | awk -F. '{print $1}')
        
        # Validate major version is numeric and meets minimum requirement
        if ! [[ "$gcc_major" =~ ^[0-9]+$ ]]; then
            add_validation_error "Unable to parse GCC major version from: $gcc_version"
            add_validation_error "Expected format: X.Y.Z (e.g., 13.2.0), but got: $gcc_version"
            add_validation_error ""
            add_validation_error "SOLUTION: Use '--with-gcc=install' to install a standard GCC version:"
            add_validation_error "  ./install_abacus_toolchain_new.sh --with-gcc=install [other options]"
            add_validation_error ""
            add_validation_error "This will install GCC ${gcc_ver:-13.2.0} with standard version format."
            return
        fi
        
        if [[ "$gcc_major" -lt "$gcc_min_version" ]]; then
            add_validation_error "GCC version $gcc_version is too old (major version: $gcc_major)"
            add_validation_error "Minimum required: GCC $gcc_min_version.x or newer"
            add_validation_error ""
            add_validation_error "Your system GCC is outdated and may cause compilation failures."
            add_validation_error ""
            add_validation_error "SOLUTION: Use '--with-gcc=install' to install a modern GCC version:"
            add_validation_error "  ./install_abacus_toolchain_new.sh --with-gcc=install [other options]"
            add_validation_error ""
            add_validation_error "This will install GCC ${gcc_ver:-13.2.0} (>= $gcc_min_version.x) with full C++17/C++20 support."
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
    
    # Report results
    local total_issues=$((${#VALIDATION_ERRORS[@]} + ${#VALIDATION_WARNINGS[@]}))
    
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
        echo "  Errors: ${#VALIDATION_ERRORS[@]}"
        echo "  Warnings: ${#VALIDATION_WARNINGS[@]}"
        
        # Return error code if there are validation errors
        if [[ ${#VALIDATION_ERRORS[@]} -gt 0 ]]; then
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
