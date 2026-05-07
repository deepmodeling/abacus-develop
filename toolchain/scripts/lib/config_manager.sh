#!/bin/bash

# ============================================================================
# ABACUS Toolchain Configuration Manager
# ============================================================================
# Handles configuration parsing, validation, and management
# Author: Quantum Misaka by Trae SOLO
# Date: 2025-10-16
# ============================================================================

# Global configuration variables
declare -A CONFIG_CACHE
declare -A USER_EXPLICIT_MPI
declare -A USER_EXPLICIT_MATH
CONFIG_INITIALIZED=false
CONFIG_FILE_LOADED=false

# Package lists (from original script)
tool_list="gcc intel amd cmake"
mpi_list="mpich openmpi intelmpi"
math_list="mkl aocl openblas"
lib_list="fftw libxc scalapack elpa cereal rapidjson libtorch libnpy libri libcomm nep"
package_list="${tool_list} ${mpi_list} ${math_list} ${lib_list}"

# Configuration file paths for loading in advance
# Not in use now
CONFIG_FILE_PATHS=(
    "${ROOTDIR}/toolchain_config.conf"
)

# Load configuration from file
# Usage: config_load_from_file [config_file_path]
config_load_from_file() {
    local config_file="$1"
    
    # If no file specified, search for default locations
    if [[ -z "$config_file" ]]; then
        for path in "${CONFIG_FILE_PATHS[@]}"; do
            if [[ -f "$path" && -r "$path" ]]; then
                config_file="$path"
                break
            fi
        done
    fi
    
    # If still no config file found, return silently (not an error)
    if [[ -z "$config_file" || ! -f "$config_file" ]]; then
        return 0
    fi
    
    echo "Loading configuration from: $config_file"
    
    # Read configuration file line by line
    local line_num=0
    while IFS= read -r line || [[ -n "$line" ]]; do
        line_num=$((line_num + 1))
        
        # Skip empty lines and comments
        [[ -z "$line" || "$line" =~ ^[[:space:]]*# ]] && continue
        
        # Remove leading/trailing whitespace
        line=$(echo "$line" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
        
        # Parse variable assignments
        if [[ "$line" =~ ^([A-Z_][A-Z0-9_]*)=(.*)$ ]]; then
            local var_name="${BASH_REMATCH[1]}"
            local var_value="${BASH_REMATCH[2]}"
            
            # Remove quotes if present
            var_value=$(echo "$var_value" | sed 's/^"//;s/"$//')
            
            # Validate and store configuration
            case "$var_name" in
                # Basic configuration
                INSTALL_PREFIX|BUILD_PREFIX)
                    CONFIG_CACHE["$var_name"]="$var_value"
                    ;;
                NPROCS)
                    if [[ "$var_value" =~ ^[0-9]+$ ]]; then
                        CONFIG_CACHE["NPROCS_OVERWRITE"]="$var_value"
                    else
                        echo "Warning: Invalid NPROCS value '$var_value' in config file line $line_num"
                    fi
                    ;;
                    
                # Mode selections
                MPI_MODE)
                    case "$var_value" in
                        mpich|openmpi|intelmpi|none)
                            CONFIG_CACHE["MPI_MODE"]="$var_value"
                            ;;
                        *)
                            echo "Warning: Invalid MPI_MODE value '$var_value' in config file line $line_num"
                            ;;
                    esac
                    ;;
                MATH_MODE)
                    case "$var_value" in
                        mkl|aocl|openblas|none)
                            CONFIG_CACHE["MATH_MODE"]="$var_value"
                            ;;
                        *)
                            echo "Warning: Invalid MATH_MODE value '$var_value' in config file line $line_num"
                            ;;
                    esac
                    ;;
                    
                # Package configurations
                WITH_*)
                    local package_name=$(echo "${var_name#WITH_}" | tr '[:upper:]' '[:lower:]')
                    # Handle special case for 4th-openmpi (convert to OPENMPI_4TH)
                    if [[ "$package_name" == "4th-openmpi" ]]; then
                        case "$var_value" in
                            yes|no|__DONTUSE__)
                                # Map __DONTUSE__ to "no" for OPENMPI_4TH
                                if [[ "$var_value" == "__DONTUSE__" ]]; then
                                    CONFIG_CACHE["OPENMPI_4TH"]="no"
                                else
                                    CONFIG_CACHE["OPENMPI_4TH"]="$var_value"
                                fi
                                ;;
                            *)
                                echo "Warning: Invalid value '$var_value' for WITH_4TH_OPENMPI in config file line $line_num"
                                ;;
                        esac
                    else
                        case "$var_value" in
                            __SYSTEM__|__INSTALL__|__DONTUSE__)
                                CONFIG_CACHE["with_${package_name}"]="$var_value"
                                ;;
                            *)
                                # Allow custom paths for some packages
                                CONFIG_CACHE["with_${package_name}"]="$var_value"
                                ;;
                        esac
                    fi
                    ;;
                    
                # Version strategy options (NEW)
                VERSION_STRATEGY)
                    case "$var_value" in
                        main|alt)
                            CONFIG_CACHE["VERSION_STRATEGY"]="$var_value"
                            ;;
                        *)
                            echo "Warning: Invalid VERSION_STRATEGY value '$var_value' in config file line $line_num. Valid options: main, alt"
                            ;;
                    esac
                    ;;
                PACKAGE_VERSION_*)
                    # Extract package name from PACKAGE_VERSION_PACKAGENAME
                    local pkg_name="${var_name#PACKAGE_VERSION_}"
                    case "$var_value" in
                        main|alt)
                            CONFIG_CACHE["PACKAGE_VERSION_${pkg_name}"]="$var_value"
                            ;;
                        *)
                            echo "Warning: Invalid package version value '$var_value' for $var_name in config file line $line_num. Valid options: main, alt"
                            ;;
                    esac
                    ;;
                    
                # Advanced options
                # Not in use, leave it for future
                DEBUG_MODE|VERBOSE_MODE|SKIP_SYSTEM_CHECKS|FORCE_REINSTALL)
                    case "$var_value" in
                        true|false)
                            CONFIG_CACHE["$var_name"]="$var_value"
                            ;;
                        *)
                            echo "Warning: Invalid boolean value '$var_value' for $var_name in config file line $line_num"
                            ;;
                    esac
                    ;;
                    
                *)
                    echo "Warning: Unknown configuration option '$var_name' in config file line $line_num"
                    ;;
            esac
        elif [[ -n "$line" ]]; then
            echo "Warning: Invalid configuration syntax in config file line $line_num: $line"
        fi
    done < "$config_file"
    
    CONFIG_FILE_LOADED=true
    echo "Configuration loaded successfully from: $config_file"
    return 0
}

# Apply mode-based configuration after loading from file
# Usage: config_apply_modes_from_file
config_apply_modes_from_file() {
    # Apply MPI mode if set in config file (silent mode to avoid duplicate output)
    if [[ -n "${CONFIG_CACHE[MPI_MODE]}" ]]; then
        case "${CONFIG_CACHE[MPI_MODE]}" in
            mpich)
                CONFIG_CACHE["with_mpich"]="__INSTALL__"
                CONFIG_CACHE["with_openmpi"]="__DONTUSE__"
                CONFIG_CACHE["with_intelmpi"]="__DONTUSE__"
                ;;
            openmpi)
                CONFIG_CACHE["with_mpich"]="__DONTUSE__"
                CONFIG_CACHE["with_openmpi"]="__INSTALL__"
                CONFIG_CACHE["with_intelmpi"]="__DONTUSE__"
                ;;
            intelmpi)
                CONFIG_CACHE["with_mpich"]="__DONTUSE__"
                CONFIG_CACHE["with_openmpi"]="__DONTUSE__"
                CONFIG_CACHE["with_intelmpi"]="__INSTALL__"
                ;;
            none)
                CONFIG_CACHE["with_mpich"]="__DONTUSE__"
                CONFIG_CACHE["with_openmpi"]="__DONTUSE__"
                CONFIG_CACHE["with_intelmpi"]="__DONTUSE__"
                ;;
        esac
    fi
    
    # Apply Math mode if set in config file (silent mode to avoid duplicate output)
    if [[ -n "${CONFIG_CACHE[MATH_MODE]}" ]]; then
        case "${CONFIG_CACHE[MATH_MODE]}" in
            mkl)
                CONFIG_CACHE["with_mkl"]="__SYSTEM__"
                CONFIG_CACHE["with_aocl"]="__DONTUSE__"
                CONFIG_CACHE["with_openblas"]="__DONTUSE__"
                # MKL provides FFTW and ScaLAPACK, so disable them (original logic)
                CONFIG_CACHE["with_fftw"]="__DONTUSE__"
                CONFIG_CACHE["with_scalapack"]="__DONTUSE__"
                ;;
            aocl)
                CONFIG_CACHE["with_mkl"]="__DONTUSE__"
                CONFIG_CACHE["with_aocl"]="__INSTALL__"
                CONFIG_CACHE["with_openblas"]="__DONTUSE__"
                ;;
            openblas)
                CONFIG_CACHE["with_mkl"]="__DONTUSE__"
                CONFIG_CACHE["with_aocl"]="__DONTUSE__"
                CONFIG_CACHE["with_openblas"]="__INSTALL__"
                ;;
            none)
                CONFIG_CACHE["with_mkl"]="__DONTUSE__"
                CONFIG_CACHE["with_aocl"]="__DONTUSE__"
                CONFIG_CACHE["with_openblas"]="__DONTUSE__"
                ;;
        esac
    fi
}

# Initialize configuration manager
# Usage: config_manager_init
config_manager_init() {
    if [[ "$CONFIG_INITIALIZED" == "true" ]]; then
        return 0
    fi
    
    # Version management is now handled directly in individual scripts via package_versions.sh
    
    # Set default configuration first
    config_set_defaults
    
    # Load configuration from file (if available) - this will override defaults
    config_load_from_file
    
    # Apply mode-based configurations from file - this will override defaults
    config_apply_modes_from_file
    
    CONFIG_INITIALIZED=true
    
    return 0
}

# Set default configuration values
# Usage: config_set_defaults
config_set_defaults() {
    # First set everything to __DONTUSE__
    for ii in ${package_list}; do
        CONFIG_CACHE["with_${ii}"]="__DONTUSE__"
    done

    # Tools to turn on by default
    CONFIG_CACHE["with_gcc"]="__SYSTEM__"
    CONFIG_CACHE["with_cmake"]="__INSTALL__"
    
    # Set default parallel jobs (NEW: Fix for parallel jobs parameter)
    if [[ -z "${NPROCS_OVERWRITE}" ]]; then
        if command -v nproc > /dev/null 2>&1; then
            CONFIG_CACHE["NPROCS_OVERWRITE"]="$(nproc --all)"
        elif command -v sysctl > /dev/null 2>&1; then
            CONFIG_CACHE["NPROCS_OVERWRITE"]="$(sysctl -n hw.ncpu)"
        else
            CONFIG_CACHE["NPROCS_OVERWRITE"]="1"
        fi
    else
        CONFIG_CACHE["NPROCS_OVERWRITE"]="${NPROCS_OVERWRITE}"
    fi
    
    # Default MPI and Math modes (following original script logic)
    # Default math library settings to openblas
    CONFIG_CACHE["MATH_MODE"]="openblas"
    
    # Set math library defaults based on MATH_MODE (before user input processing)
    case "${CONFIG_CACHE[MATH_MODE]}" in
        mkl)
            CONFIG_CACHE["with_mkl"]="__SYSTEM__"
            ;;
        aocl)
            CONFIG_CACHE["with_aocl"]="__SYSTEM__"
            ;;
        openblas)
            CONFIG_CACHE["with_openblas"]="__INSTALL__"
            ;;
    esac
    
    # For MPI, we try to detect system MPI variant (following original script logic)
    # Only detect if MPI_MODE is not already set (to avoid duplicate detection)
    if [[ -z "${CONFIG_CACHE[MPI_MODE]}" ]] && command -v mpiexec > /dev/null 2>&1; then
        # check if we are dealing with openmpi, mpich or intelmpi
        if mpiexec --version 2>&1 | grep -s -q "HYDRA"; then
            if command -v ui_info &> /dev/null; then
                ui_info "🔍 Detected system MPI: MPICH"
            else
                echo "MPI is detected and it appears to be MPICH"
            fi
            CONFIG_CACHE["MPI_MODE"]="mpich"
            CONFIG_CACHE["with_mpich"]="__SYSTEM__"
        elif mpiexec --version 2>&1 | grep -s -q "OpenRTE"; then
            if command -v ui_info &> /dev/null; then
                ui_info "🔍 Detected system MPI: OpenMPI"
            else
                echo "MPI is detected and it appears to be OpenMPI"
            fi
            CONFIG_CACHE["MPI_MODE"]="openmpi"
            CONFIG_CACHE["with_openmpi"]="__SYSTEM__"
        elif mpiexec --version 2>&1 | grep -s -q "Intel"; then
            if command -v ui_info &> /dev/null; then
                ui_info "🔍 Detected system MPI: Intel MPI"
            else
                echo "MPI is detected and it appears to be Intel MPI"
            fi
            CONFIG_CACHE["with_gcc"]="__DONTUSE__"
            CONFIG_CACHE["with_amd"]="__DONTUSE__"
            CONFIG_CACHE["with_aocl"]="__DONTUSE__"
            CONFIG_CACHE["with_intel"]="__SYSTEM__"
            CONFIG_CACHE["with_intelmpi"]="__SYSTEM__"
            CONFIG_CACHE["MPI_MODE"]="intelmpi"
        else # default to mpich
            if command -v ui_info &> /dev/null; then
                ui_info "🔍 MPI detected, defaulting to MPICH configuration"
            else
                echo "MPI is detected and defaults to MPICH"
            fi
            CONFIG_CACHE["MPI_MODE"]="mpich"
            CONFIG_CACHE["with_mpich"]="__SYSTEM__"
        fi
    else
        if command -v report_warning &> /dev/null; then
            report_warning ${LINENO} "No MPI installation detected (ignore this message in Cray Linux Environment or when MPI installation was requested)."
        else
            echo "Warning: No MPI installation detected (ignore this message in Cray Linux Environment or when MPI installation was requested)."
        fi
        CONFIG_CACHE["MPI_MODE"]="no"
    fi
    
    # Default libraries to install
    CONFIG_CACHE["with_fftw"]="__INSTALL__"
    CONFIG_CACHE["with_libxc"]="__INSTALL__"
    CONFIG_CACHE["with_scalapack"]="__INSTALL__"
    CONFIG_CACHE["with_elpa"]="__INSTALL__"
    CONFIG_CACHE["with_cereal"]="__INSTALL__"
    CONFIG_CACHE["with_rapidjson"]="__INSTALL__"
    CONFIG_CACHE["with_libri"]="__INSTALL__"
    CONFIG_CACHE["with_libcomm"]="__INSTALL__"
    CONFIG_CACHE["with_nep"]="__DONTUSE__"
    
    # Default enable options (following original script logic)
    CONFIG_CACHE["dry_run"]="__FALSE__"
    CONFIG_CACHE["enable_tsan"]="__FALSE__"
    CONFIG_CACHE["enable_opencl"]="__FALSE__"
    CONFIG_CACHE["enable_cuda"]="__FALSE__"
    CONFIG_CACHE["enable_hip"]="__FALSE__"
    CONFIG_CACHE["intel_classic"]="no"
    CONFIG_CACHE["PACK_RUN"]="__FALSE__"
    CONFIG_CACHE["INTELMPI_CLASSIC"]="no"
    CONFIG_CACHE["WITH_IFX"]="yes"
    CONFIG_CACHE["WITH_FLANG"]="no"
    CONFIG_CACHE["OPENMPI_4TH"]="no"
    CONFIG_CACHE["GPUVER"]="no"
    CONFIG_CACHE["MPICH_DEVICE"]="ch4"
    CONFIG_CACHE["TARGET_CPU"]="native"
    CONFIG_CACHE["LOG_LINES"]="200"
    CONFIG_CACHE["show_help"]="false"
    CONFIG_CACHE["DOWNLOADER_FLAGS"]=""
    
    # Version strategy defaults (NEW)
    CONFIG_CACHE["VERSION_STRATEGY"]="main"
    
    # Defaults for CRAY Linux Environment (following original script logic)
    if [[ -n "${CRAY_LD_LIBRARY_PATH}" ]]; then
        CONFIG_CACHE["enable_cray"]="__TRUE__"
        CONFIG_CACHE["MATH_MODE"]="cray"
        # Default MPI used by CLE is assumed to be MPICH, in any case
        # do not use the installers for the MPI libraries
        CONFIG_CACHE["with_mpich"]="__DONTUSE__"
        CONFIG_CACHE["with_openmpi"]="__DONTUSE__"
        CONFIG_CACHE["with_intelmpi"]="__DONTUSE__"
        CONFIG_CACHE["MPI_MODE"]="mpich"
        # set default value for some installers appropriate for CLE
        CONFIG_CACHE["with_gcc"]="__DONTUSE__"
        # Reset math library defaults for CRAY environment
        CONFIG_CACHE["with_mkl"]="__DONTUSE__"
        CONFIG_CACHE["with_aocl"]="__DONTUSE__"
        CONFIG_CACHE["with_openblas"]="__DONTUSE__"
    fi
    
    return 0
}

# Helper function to read --enable-* options
# Usage: read_enable_option "--enable-feature=value"
read_enable_option() {
    local arg="$1"
    if [[ "$arg" == *"="* ]]; then
        local value="${arg#*=}"
        case "$value" in
            yes|true|1) echo "__TRUE__" ;;
            no|false|0) echo "__FALSE__" ;;
            *) echo "__INVALID__" ;;
        esac
    else
        echo "__TRUE__"
    fi
}

# Helper function to read --with-* options
# Usage: read_with_option "--with-package=value"
read_with_option() {
    local arg="$1"
    if [[ "$arg" == *"="* ]]; then
        local value="${arg#*=}"
        case "$value" in
            install) echo "__INSTALL__" ;;
            system) echo "__SYSTEM__" ;;
            no) echo "__DONTUSE__" ;;
            *) echo "$value" ;;  # Allow custom paths
        esac
    else
        echo "__INSTALL__"
    fi
}

# Validate configuration
# Usage: config_validate
config_validate() {
    # Validate numeric values
    if [[ -n "${CONFIG_CACHE[NPROCS_OVERWRITE]}" ]]; then
        if ! [[ "${CONFIG_CACHE[NPROCS_OVERWRITE]}" =~ ^[0-9]+$ ]]; then
            report_error ${LINENO} "Invalid number of processes: ${CONFIG_CACHE[NPROCS_OVERWRITE]}"
            exit 1
        fi
    fi
    
    if ! [[ "${CONFIG_CACHE[LOG_LINES]}" =~ ^[0-9]+$ ]]; then
        report_error ${LINENO} "Invalid log lines value: ${CONFIG_CACHE[LOG_LINES]}"
        exit 1
    fi
    
    # Validate GPU version - support only numeric formats
    local gpu_ver="${CONFIG_CACHE[GPUVER]}"
    if [[ "$gpu_ver" != "no" ]]; then
        # Set ARCH_NUM based on GPUVER (remove decimal point for numeric versions)
        local arch_num="${gpu_ver//.}"
        
        # Check if it's a valid numeric format (like 8.0, 70, 80, etc.)
        if [[ "$arch_num" =~ ^[1-9][0-9]*$ ]]; then
            CONFIG_CACHE["ARCH_NUM"]="$arch_num"
        else
            report_error ${LINENO} "Invalid GPU version: $gpu_ver. Supported formats: numeric with decimal (6.0, 7.0, 8.0, 8.9, etc.) or numeric without decimal (60, 70, 80, 89, etc.)"
            exit 1
        fi
    else
        CONFIG_CACHE["ARCH_NUM"]="no"
    fi
    
    # Backward compatibility: also export ARCH_NUM to environment when set
    if [[ -n "${CONFIG_CACHE[ARCH_NUM]}" ]]; then
        export ARCH_NUM="${CONFIG_CACHE[ARCH_NUM]}"
    fi
    
    return 0
}

# Apply environment variable processing logic (from original script lines 668-744)
# Usage: config_apply_env_logic
config_apply_env_logic() {
    # Compiler conflicts (from original script L668-677)
    if [ "${CONFIG_CACHE[with_intel]}" != "__DONTUSE__" ] && [ "${CONFIG_CACHE[with_gcc]}" = "__INSTALL__" ]; then
        echo "You have chosen to use the Intel compiler, therefore the installation of the GNU compiler will be skipped."
        CONFIG_CACHE["with_gcc"]="__SYSTEM__"
    fi
    if [ "${CONFIG_CACHE[with_amd]}" != "__DONTUSE__" ] && [ "${CONFIG_CACHE[with_gcc]}" = "__INSTALL__" ]; then
        echo "You have chosen to use the AMD compiler, therefore the installation of the GNU compiler will be skipped."
        CONFIG_CACHE["with_gcc"]="__SYSTEM__"
    fi
    if [ "${CONFIG_CACHE[with_amd]}" != "__DONTUSE__" ] && [ "${CONFIG_CACHE[with_intel]}" != "__DONTUSE__" ]; then
        report_error "You have chosen to use the AMD and the Intel compiler to compile dependent packages. Select only one compiler."
        exit 1
    fi

    # MPI library conflicts (from original script L679-710)
    if [ "${CONFIG_CACHE[MPI_MODE]}" = "no" ]; then
        if [ "${CONFIG_CACHE[with_scalapack]}" != "__DONTUSE__" ]; then
            echo "Not using MPI, so scalapack is disabled."
            CONFIG_CACHE["with_scalapack"]="__DONTUSE__"
        fi
        if [ "${CONFIG_CACHE[with_elpa]}" != "__DONTUSE__" ]; then
            echo "Not using MPI, so ELPA is disabled."
            CONFIG_CACHE["with_elpa"]="__DONTUSE__"
        fi
    else
        # if gcc is installed, then mpi needs to be installed too
        if [ "${CONFIG_CACHE[with_gcc]}" = "__INSTALL__" ]; then
            echo "You have chosen to install the GNU compiler, therefore MPI libraries have to be installed too"
            case ${CONFIG_CACHE[MPI_MODE]} in
                mpich)
                    CONFIG_CACHE["with_mpich"]="__INSTALL__"
                    CONFIG_CACHE["with_openmpi"]="__DONTUSE__"
                    ;;
                openmpi)
                    CONFIG_CACHE["with_mpich"]="__DONTUSE__"
                    CONFIG_CACHE["with_openmpi"]="__INSTALL__"
                    ;;
            esac
            echo "and the use of the Intel compiler and Intel MPI will be disabled."
            CONFIG_CACHE["with_intel"]="__DONTUSE__"
            CONFIG_CACHE["with_intelmpi"]="__DONTUSE__"
        fi
        # Enable only one MPI implementation
        case ${CONFIG_CACHE[MPI_MODE]} in
            mpich)
                CONFIG_CACHE["with_openmpi"]="__DONTUSE__"
                CONFIG_CACHE["with_intelmpi"]="__DONTUSE__"
                ;;
            openmpi)
                CONFIG_CACHE["with_mpich"]="__DONTUSE__"
                CONFIG_CACHE["with_intelmpi"]="__DONTUSE__"
                ;;
            intelmpi)
                CONFIG_CACHE["with_mpich"]="__DONTUSE__"
                CONFIG_CACHE["with_openmpi"]="__DONTUSE__"
                ;;
        esac
    fi

    # If MATH_MODE is mkl, then openblas, scalapack and fftw is not needed
    # QuantumMisaka in 2023-09-17
    if [ "${CONFIG_CACHE[MATH_MODE]}" = "mkl" ]; then
        if [ "${CONFIG_CACHE[with_openblas]}" != "__DONTUSE__" ]; then
            echo "Using MKL, so openblas is disabled."
            CONFIG_CACHE["with_openblas"]="__DONTUSE__"
        fi
        if [ "${CONFIG_CACHE[with_scalapack]}" != "__DONTUSE__" ]; then
            echo "Using MKL, so scalapack is disabled."
            CONFIG_CACHE["with_scalapack"]="__DONTUSE__"
        fi
        if [ "${CONFIG_CACHE[with_fftw]}" != "__DONTUSE__" ]; then
            echo "Using MKL, so fftw is disabled."
            CONFIG_CACHE["with_fftw"]="__DONTUSE__"
        fi
    fi

    # Select the correct compute number based on the GPU architecture
    # QuantumMisaka in 2025-03-19
    local gpuver="${CONFIG_CACHE[GPUVER]}"
    if [ "$gpuver" != "no" ]; then
        export ARCH_NUM="${gpuver//.}"
        CONFIG_CACHE["ARCH_NUM"]="$ARCH_NUM"
    else
        export ARCH_NUM="no"
        CONFIG_CACHE["ARCH_NUM"]="no"
    fi

    # If CUDA or HIP are enabled, make sure the GPU version has been defined.
    if [ "${CONFIG_CACHE[enable_cuda]}" = "__TRUE__" ] || [ "${CONFIG_CACHE[enable_hip]}" = "__TRUE__" ]; then
        if [ "${CONFIG_CACHE[GPUVER]}" = "no" ]; then
            report_error "Please choose GPU architecture to compile for with --gpu-ver"
            exit 1
        fi
        if [[ "$ARCH_NUM" =~ ^[1-9][0-9]*$ ]] || [ "$ARCH_NUM" = "no" ]; then
            echo "Notice: GPU compilation is enabled, and GPU compatibility is set via --gpu-ver to sm_${ARCH_NUM}."
        else
            report_error ${LINENO} \
                "When GPU compilation is enabled, the --gpu-ver variable should be properly set regarding to GPU compatibility. For check your GPU compatibility, visit https://developer.nvidia.com/cuda-gpus. For example: A100 -> 8.0 (or 80), V100 -> 7.0 (or 70), 4090 -> 8.9 (or 89)"
            exit 1
        fi
    fi

    # ABACUS itself and some dependencies require cmake.
    if [ "${CONFIG_CACHE[with_cmake]}" = "__DONTUSE__" ]; then
        report_error "CMake is required for ABACUS and some dependencies. Please enable it."
        exit 1
    fi
}

# Export configuration to environment variables
# Usage: config_export_to_env
config_export_to_env() {
    # Apply environment variable processing logic first
    config_apply_env_logic
    
    # Export all configuration values as environment variables
    for key in "${!CONFIG_CACHE[@]}"; do
        export "$key"="${CONFIG_CACHE[$key]}"
    done

    # Backward compatibility for stage scripts expecting uppercase GPU flags
    # Installers (e.g., stage3/install_elpa.sh) read ENABLE_CUDA, not enable_cuda
    export ENABLE_CUDA="${CONFIG_CACHE[enable_cuda]}"
    
    # Export package list variables
    export tool_list
    export mpi_list
    export math_list
    export lib_list
    export package_list
}

# Get configuration value
# Usage: config_get "key"
config_get() {
    local key="$1"
    echo "${CONFIG_CACHE[$key]}"
}

# Set configuration value
# Usage: config_set "key" "value"
config_set() {
    local key="$1"
    local value="$2"
    CONFIG_CACHE["$key"]="$value"
}

# Check if configuration key exists
# Usage: config_has "key"
config_has() {
    local key="$1"
    [[ -n "${CONFIG_CACHE[$key]}" ]]
}

# Print configuration summary
# Usage: config_print_summary
config_print_summary() {
    echo "Configuration Summary:"
    echo "====================="
    echo "MPI Mode: ${CONFIG_CACHE[MPI_MODE]}"
    echo "Math Mode: ${CONFIG_CACHE[MATH_MODE]}"
    echo "Target CPU: ${CONFIG_CACHE[TARGET_CPU]}"
    echo "GPU Version: ${CONFIG_CACHE[GPUVER]}"
    echo "Parallel Jobs: ${CONFIG_CACHE[NPROCS_OVERWRITE]:-$(get_nprocs)}"
    echo "Dry Run: ${CONFIG_CACHE[dry_run]}"
    echo ""
    
    echo "Package Configuration:"
    echo "====================="
    for pkg in ${package_list}; do
        local status="${CONFIG_CACHE[with_${pkg}]}"
        if [[ "$status" != "__DONTUSE__" ]]; then
            printf "%-15s: %s\n" "$pkg" "$status"
        fi
    done
    echo ""
}

# Parse command line arguments (New unified version)
# Usage: config_parse_arguments "$@"
config_parse_arguments() {
    local show_help=false
    local show_version=false
    
    while [[ $# -gt 0 ]]; do
        case $1 in
            # Help and version
            -h|--help)
                show_help=true
                shift
                ;;
            --version)
                show_version=true
                shift
                ;;
            --version-info)
                if [[ -n "$2" && "$2" != -* ]]; then
                    CONFIG_CACHE["show_version_info"]="$2"
                    shift 2
                else
                    CONFIG_CACHE["show_version_info"]="all"
                    shift
                fi
                ;;
                
            # Build options
            -j)
                if [[ -n "$2" && "$2" =~ ^[0-9]+$ ]]; then
                    CONFIG_CACHE["NPROCS_OVERWRITE"]="$2"
                    shift 2
                else
                    report_error $LINENO "-j requires a number argument"
                    return 1
                fi
                ;;
            --dry-run)
                CONFIG_CACHE["dry_run"]="__TRUE__"
                shift
                ;;
            --pack-run)
                CONFIG_CACHE["PACK_RUN"]="__TRUE__"
                shift
                ;;
                
            # Package version selection - Support multiple package:version pairs
            --package-version)
                local package_version_args=""
                local processed_count=0
                
                # Handle --package-version=value format
                if [[ "$1" =~ ^--package-version=(.+)$ ]]; then
                    package_version_args="${BASH_REMATCH[1]}"
                    shift
                    processed_count=1
                else
                    # Collect all consecutive non-option arguments
                    shift  # Skip --package-version
                    while [[ $# -gt 0 && "$1" != -* ]]; do
                        if [[ -n "$package_version_args" ]]; then
                            package_version_args="$package_version_args $1"
                        else
                            package_version_args="$1"
                        fi
                        shift
                        ((processed_count++))
                    done
                fi
                
                # Validate we have at least one argument
                if [[ -z "$package_version_args" ]]; then
                    report_error $LINENO "--package-version requires at least one package:version argument"
                    return 1
                fi
                
                # Process each package:version pair
                local pair_count=0
                for pair in $package_version_args; do
                    if [[ "$pair" =~ ^([a-zA-Z0-9_]+):(main|alt)$ ]]; then
                        local pkg="${BASH_REMATCH[1]}"
                        local ver="${BASH_REMATCH[2]}"
                        CONFIG_CACHE["PACKAGE_VERSION_${pkg^^}"]="$ver"
                        ((pair_count++))
                    else
                        report_error $LINENO "Invalid package version format: '$pair'. Use format 'package:version' (e.g., openmpi:alt, openblas:main)"
                        return 1
                    fi
                done
                
                # Report successful processing
                if [[ $pair_count -gt 1 ]]; then
                    echo "INFO: Processed $pair_count package version overrides"
                fi
                ;;
                
            # Configuration file
            --config-file)
                if [[ -n "$2" && "$2" != -* ]]; then
                    config_load_from_file "$2"
                    shift 2
                else
                    report_error $LINENO "--config-file requires an argument"
                    return 1
                fi
                ;;
                
            # Mode selections - Support both --mode=value and --mode value formats
            --mpi-mode|--mpi-mode=*)
                local mpi_value=""
                # Handle --mpi-mode=value format
                if [[ "$1" =~ ^--mpi-mode=(.+)$ ]]; then
                    mpi_value="${BASH_REMATCH[1]}"
                    shift
                elif [[ -n "$2" && "$2" != -* ]]; then
                    mpi_value="$2"
                    shift 2
                else
                    report_error $LINENO "--mpi-mode requires an argument"
                    return 1
                fi
                
                case "$mpi_value" in
                    mpich|openmpi|intelmpi|no)
                        CONFIG_CACHE["MPI_MODE"]="$mpi_value"
                        ;;
                    *)
                        report_error $LINENO "Invalid MPI mode: $mpi_value. Valid options: mpich, openmpi, intelmpi, no"
                        return 1
                        ;;
                esac
                ;;
            --math-mode|--math-mode=*)
                local math_value=""
                # Handle --math-mode=value format
                if [[ "$1" =~ ^--math-mode=(.+)$ ]]; then
                    math_value="${BASH_REMATCH[1]}"
                    shift
                elif [[ -n "$2" && "$2" != -* ]]; then
                    math_value="$2"
                    shift 2
                else
                    report_error $LINENO "--math-mode requires an argument"
                    return 1
                fi
                
                case "$math_value" in
                    mkl|aocl|openblas|cray|no)
                        CONFIG_CACHE["MATH_MODE"]="$math_value"
                        # Apply automatic parameter settings for specific math modes
                        case "$math_value" in
                            aocl)
                                CONFIG_CACHE["with_aocl"]="__SYSTEM__"
                                CONFIG_CACHE["with_fftw"]="__SYSTEM__"
                                CONFIG_CACHE["with_scalapack"]="__SYSTEM__"
                                ;;
                        esac
                        ;;
                    *)
                        report_error $LINENO "Invalid math mode: $math_value. Valid options: mkl, aocl, openblas, cray, no"
                        return 1
                        ;;
                esac
                ;;
                
            # Package options
            --with-*)
                local option="${1#--with-}"
                local value="__INSTALL__"
                
                # Handle --with-package=value format
                if [[ "$option" =~ ^(.+)=(.+)$ ]]; then
                    option="${BASH_REMATCH[1]}"
                    value="${BASH_REMATCH[2]}"
                fi
                
                # Special handling for --with-mpich-device=* (must be handled before general mpich)
                if [[ "$option" == "mpich-device" ]]; then
                    CONFIG_CACHE["MPICH_DEVICE"]="$value"
                    CONFIG_CACHE["MPI_MODE"]="mpich"
                    shift
                    continue
                fi
                
                # Convert value to internal format (matching original read_with function)
                case "$value" in
                    system)
                        value="__SYSTEM__"
                        ;;
                    install)
                        value="__INSTALL__"
                        ;;
                    no)
                        value="__DONTUSE__"
                        ;;
                    __SYSTEM__|__INSTALL__|__DONTUSE__)
                        # Already in internal format, keep as is
                        ;;
                    *)
                        # For custom paths, expand ~ to $HOME (matching original behavior)
                        value="${value//\~/$HOME}"
                        ;;
                esac
                
                # Handle special cases with specific processing logic
                case "$option" in
                    4th-openmpi)
                        # Handle --with-4th-openmpi parameter (only yes/no options)
                        case "$value" in
                            "yes"|"install"|"__INSTALL__")
                                CONFIG_CACHE["OPENMPI_4TH"]="yes"
                                ;;
                            "no"|"__DONTUSE__")
                                CONFIG_CACHE["OPENMPI_4TH"]="no"
                                ;;
                            *)
                                # Default to "no" for any other values
                                CONFIG_CACHE["OPENMPI_4TH"]="no"
                                ;;
                        esac
                        ;;
                    mpich)
                        CONFIG_CACHE["with_mpich"]="$value"
                        USER_EXPLICIT_MPI["with_mpich"]="true"  # Mark as explicitly set by user
                        # Set MPI_MODE if not disabled
                        if [ "$value" != "__DONTUSE__" ]; then
                            CONFIG_CACHE["MPI_MODE"]="mpich"
                        fi
                        ;;
                    openmpi)
                        CONFIG_CACHE["with_openmpi"]="$value"
                        USER_EXPLICIT_MPI["with_openmpi"]="true"  # Mark as explicitly set by user
                        # Set MPI_MODE if not disabled
                        if [ "$value" != "__DONTUSE__" ]; then
                            CONFIG_CACHE["MPI_MODE"]="openmpi"
                        fi
                        ;;
                    intelmpi)
                        CONFIG_CACHE["with_intelmpi"]="$value"
                        USER_EXPLICIT_MPI["with_intelmpi"]="true"  # Mark as explicitly set by user
                        # Set MPI_MODE if not disabled
                        if [ "$value" != "__DONTUSE__" ]; then
                            CONFIG_CACHE["MPI_MODE"]="intelmpi"
                        fi
                        ;;
                    intel-classic)
                        # Special handling for intel-classic: only accepts yes/no values
                        case "$value" in
                            "__INSTALL__"|""|"__DONTUSE__")
                                value="no"  # Default to "no"
                                ;;
                            "yes")
                                value="yes"
                                ;;
                            "no")
                                value="no"
                                ;;
                            *)
                                report_error $LINENO "Invalid value '$value' for --with-intel-classic. Only 'yes' or 'no' are allowed." "CONFIG_ERROR"
                                return 1
                                ;;
                        esac
                        CONFIG_CACHE["intel_classic"]="$value"
                        ;;
                    intel-mpi-clas*)
                        # Special handling for intel-mpi-classic: only accepts yes/no values
                        case "$value" in
                            "__INSTALL__"|""|"__DONTUSE__")
                                value="no"  # Default to "no"
                                ;;
                            "yes")
                                value="yes"
                                ;;
                            "no")
                                value="no"
                                ;;
                            *)
                                report_error $LINENO "Invalid value '$value' for --with-intel-mpi-classic. Only 'yes' or 'no' are allowed." "CONFIG_ERROR"
                                return 1
                                ;;
                        esac
                        CONFIG_CACHE["INTELMPI_CLASSIC"]="$value"
                        ;;
                    intel)
                        CONFIG_CACHE["with_intel"]="$value"
                        ;;
                    ifx)
                        # Special handling for ifx: only accepts yes/no values
                        case "$value" in
                            "__INSTALL__"|""|"__DONTUSE__")
                                value="no"  # Default to "no"
                                ;;
                            "yes")
                                value="yes"
                                ;;
                            "no")
                                value="no"
                                ;;
                            *)
                                report_error $LINENO "Invalid value '$value' for --with-ifx. Only 'yes' or 'no' are allowed." "CONFIG_ERROR"
                                return 1
                                ;;
                        esac
                        CONFIG_CACHE["WITH_IFX"]="$value"
                        ;;
                    amd)
                        CONFIG_CACHE["with_amd"]="$value"
                        ;;
                    flang)
                        # Special handling for flang: only accepts yes/no values
                        case "$value" in
                            "__INSTALL__"|""|"__DONTUSE__")
                                value="no"  # Default to "no"
                                ;;
                            "yes")
                                value="yes"
                                ;;
                            "no")
                                value="no"
                                ;;
                            *)
                                report_error $LINENO "Invalid value '$value' for --with-flang. Only 'yes' or 'no' are allowed." "CONFIG_ERROR"
                                return 1
                                ;;
                        esac
                        CONFIG_CACHE["WITH_FLANG"]="$value"
                        ;;
                    aocl)
                        CONFIG_CACHE["with_aocl"]="$value"
                        USER_EXPLICIT_MATH["with_aocl"]="true"  # Mark as explicitly set by user
                        ;;
                    mkl)
                        CONFIG_CACHE["with_mkl"]="$value"
                        USER_EXPLICIT_MATH["with_mkl"]="true"  # Mark as explicitly set by user
                        # Set MATH_MODE if not disabled
                        if [ "$value" != "__DONTUSE__" ]; then
                            CONFIG_CACHE["MATH_MODE"]="mkl"
                        fi
                        ;;
                    openblas)
                        CONFIG_CACHE["with_openblas"]="$value"
                        USER_EXPLICIT_MATH["with_openblas"]="true"  # Mark as explicitly set by user
                        # Set MATH_MODE if not disabled
                        if [ "$value" != "__DONTUSE__" ]; then
                            CONFIG_CACHE["MATH_MODE"]="openblas"
                        fi
                        ;;
                    fftw)
                        CONFIG_CACHE["with_fftw"]="$value"
                        USER_EXPLICIT_MATH["with_fftw"]="true"  # Mark as explicitly set by user
                        ;;
                    scalapack)
                        CONFIG_CACHE["with_scalapack"]="$value"
                        USER_EXPLICIT_MATH["with_scalapack"]="true"  # Mark as explicitly set by user
                        ;;
                    *)
                        # Convert to standard format for all other options
                        option="${option,,}"  # Convert to lowercase
                        CONFIG_CACHE["with_${option}"]="$value"
                        ;;
                esac
                shift
                ;;
                
            # Enable options
            --enable-*)
                local option="${1#--enable-}"
                local value="__TRUE__"
                
                # Handle --enable-feature=value format
                if [[ "$option" =~ ^(.+)=(.+)$ ]]; then
                    option="${BASH_REMATCH[1]}"
                    value="${BASH_REMATCH[2]}"
                    case "$value" in
                        yes|true|1) value="__TRUE__" ;;
                        no|false|0) value="__FALSE__" ;;
                        *) 
                            report_error $LINENO "Invalid value for --enable-${option}: $value. Use yes/no"
                            return 1
                            ;;
                    esac
                fi
                
                CONFIG_CACHE["enable_${option}"]="$value"
                shift
                ;;
                
            # GPU version
            --gpu-ver|--gpu-ver=*)
                local gpu_value=""
                # Handle --gpu-ver=value format
                if [[ "$1" =~ ^--gpu-ver=(.+)$ ]]; then
                    gpu_value="${BASH_REMATCH[1]}"
                    shift
                elif [[ -n "$2" && "$2" != -* ]]; then
                    gpu_value="$2"
                    shift 2
                else
                    report_error $LINENO "--gpu-ver requires an argument"
                    return 1
                fi
                CONFIG_CACHE["GPUVER"]="$gpu_value"
                ;;
                
            # Target CPU
            --target-cpu|--target-cpu=*)
                local cpu_value=""
                # Handle --target-cpu=value format
                if [[ "$1" =~ ^--target-cpu=(.+)$ ]]; then
                    cpu_value="${BASH_REMATCH[1]}"
                    shift
                elif [[ -n "$2" && "$2" != -* ]]; then
                    cpu_value="$2"
                    shift 2
                else
                    report_error $LINENO "--target-cpu requires an argument"
                    return 1
                fi
                CONFIG_CACHE["TARGET_CPU"]="$cpu_value"
                ;;
                
            # Log lines
            --log-lines|--log-lines=*)
                local log_value=""
                # Handle --log-lines=value format
                if [[ "$1" =~ ^--log-lines=(.+)$ ]]; then
                    log_value="${BASH_REMATCH[1]}"
                    shift
                elif [[ -n "$2" && "$2" =~ ^[0-9]+$ ]]; then
                    log_value="$2"
                    shift 2
                else
                    report_error $LINENO "--log-lines requires a number argument"
                    return 1
                fi
                
                if [[ "$log_value" =~ ^[0-9]+$ ]]; then
                    CONFIG_CACHE["LOG_LINES"]="$log_value"
                else
                    report_error $LINENO "--log-lines requires a number argument"
                    return 1
                fi
                ;;
                
            # Skip system checks
            --skip-system-checks)
                CONFIG_CACHE["skip_system_checks"]="__TRUE__"
                shift
                ;;
                
            # Positional arguments (ignored for now)
            *)
                shift
                ;;
        esac
    done
    
    # Handle help and version requests
    if [[ "$show_help" == true ]]; then
        CONFIG_CACHE["show_help"]="__TRUE__"
        return 0
    fi
    
    if [[ "$show_version" == true ]]; then
        CONFIG_CACHE["show_version"]="__TRUE__"
        return 0
    fi
    
    return 0
}

# Initialize configuration with command line arguments
# Usage: config_init "$@"
config_init() {
    # Set defaults first
    config_set_defaults
    
    # Initialize version helper to ensure VERSION_STRATEGY defaults are set
    if command -v version_helper_init > /dev/null 2>&1; then
        version_helper_init
    fi
    
    # Load configuration from file (if available) - this will override defaults
    config_load_from_file
    
    # Apply mode-based configurations from file - this will override defaults
    config_apply_modes_from_file
    
    # Parse command line arguments - this will override file settings
    config_parse_arguments "$@"
    
    # Apply mode-based configurations from command line
    config_apply_modes
    
    # Validate configuration
    config_validate
    
    return 0
}

# Apply configuration based on modes
# Usage: config_apply_modes
config_apply_modes() {
    # Apply MPI mode settings (with output for user feedback)
    local mpi_mode="${CONFIG_CACHE[MPI_MODE]}"
    if [[ -n "$mpi_mode" ]]; then
        if command -v ui_info &> /dev/null; then
            ui_info "⚙️ Configuring MPI mode: $mpi_mode"
        else
            echo "Applying MPI mode: $mpi_mode"
        fi
        case "$mpi_mode" in
            mpich)
                # Only override if user hasn't explicitly set these values via command line
                [[ -z "${USER_EXPLICIT_MPI[with_mpich]}" ]] && CONFIG_CACHE["with_mpich"]="__INSTALL__"
                [[ -z "${USER_EXPLICIT_MPI[with_openmpi]}" ]] && CONFIG_CACHE["with_openmpi"]="__DONTUSE__"
                [[ -z "${USER_EXPLICIT_MPI[with_intelmpi]}" ]] && CONFIG_CACHE["with_intelmpi"]="__DONTUSE__"
                ;;
            openmpi)
                [[ -z "${USER_EXPLICIT_MPI[with_mpich]}" ]] && CONFIG_CACHE["with_mpich"]="__DONTUSE__"
                [[ -z "${USER_EXPLICIT_MPI[with_openmpi]}" ]] && CONFIG_CACHE["with_openmpi"]="__INSTALL__"
                [[ -z "${USER_EXPLICIT_MPI[with_intelmpi]}" ]] && CONFIG_CACHE["with_intelmpi"]="__DONTUSE__"
                ;;
            intelmpi)
                [[ -z "${USER_EXPLICIT_MPI[with_mpich]}" ]] && CONFIG_CACHE["with_mpich"]="__DONTUSE__"
                [[ -z "${USER_EXPLICIT_MPI[with_openmpi]}" ]] && CONFIG_CACHE["with_openmpi"]="__DONTUSE__"
                [[ -z "${USER_EXPLICIT_MPI[with_intelmpi]}" ]] && CONFIG_CACHE["with_intelmpi"]="__INSTALL__"
                ;;
            no)
                [[ -z "${USER_EXPLICIT_MPI[with_mpich]}" ]] && CONFIG_CACHE["with_mpich"]="__DONTUSE__"
                [[ -z "${USER_EXPLICIT_MPI[with_openmpi]}" ]] && CONFIG_CACHE["with_openmpi"]="__DONTUSE__"
                [[ -z "${USER_EXPLICIT_MPI[with_intelmpi]}" ]] && CONFIG_CACHE["with_intelmpi"]="__DONTUSE__"
                ;;
        esac
    fi
    
    # Apply math mode settings (with output for user feedback)
    local math_mode="${CONFIG_CACHE[MATH_MODE]}"
    if [[ -n "$math_mode" ]]; then
        if command -v ui_info &> /dev/null; then
            ui_info "🧮 Configuring Math mode: $math_mode"
        else
            echo "Applying Math mode: $math_mode"
        fi
        case "$math_mode" in
            mkl)
                [[ -z "${USER_EXPLICIT_MATH[with_mkl]}" ]] && CONFIG_CACHE["with_mkl"]="__SYSTEM__"
                [[ -z "${USER_EXPLICIT_MATH[with_aocl]}" ]] && CONFIG_CACHE["with_aocl"]="__DONTUSE__"
                [[ -z "${USER_EXPLICIT_MATH[with_openblas]}" ]] && CONFIG_CACHE["with_openblas"]="__DONTUSE__"
                # MKL provides FFTW and ScaLAPACK, so disable them (original logic)
                [[ -z "${USER_EXPLICIT_MATH[with_fftw]}" ]] && CONFIG_CACHE["with_fftw"]="__DONTUSE__"
                [[ -z "${USER_EXPLICIT_MATH[with_scalapack]}" ]] && CONFIG_CACHE["with_scalapack"]="__DONTUSE__"
                ;;
            aocl)
                [[ -z "${USER_EXPLICIT_MATH[with_mkl]}" ]] && CONFIG_CACHE["with_mkl"]="__DONTUSE__"
                [[ -z "${USER_EXPLICIT_MATH[with_aocl]}" ]] && CONFIG_CACHE["with_aocl"]="__SYSTEM__"
                [[ -z "${USER_EXPLICIT_MATH[with_openblas]}" ]] && CONFIG_CACHE["with_openblas"]="__DONTUSE__"
                [[ -z "${USER_EXPLICIT_MATH[with_fftw]}" ]] && CONFIG_CACHE["with_fftw"]="__SYSTEM__"
                [[ -z "${USER_EXPLICIT_MATH[with_scalapack]}" ]] && CONFIG_CACHE["with_scalapack"]="__SYSTEM__"
                ;;
            openblas)
                [[ -z "${USER_EXPLICIT_MATH[with_mkl]}" ]] && CONFIG_CACHE["with_mkl"]="__DONTUSE__"
                [[ -z "${USER_EXPLICIT_MATH[with_aocl]}" ]] && CONFIG_CACHE["with_aocl"]="__DONTUSE__"
                [[ -z "${USER_EXPLICIT_MATH[with_openblas]}" ]] && CONFIG_CACHE["with_openblas"]="__INSTALL__"
                ;;
            cray)
                [[ -z "${USER_EXPLICIT_MATH[with_mkl]}" ]] && CONFIG_CACHE["with_mkl"]="__DONTUSE__"
                [[ -z "${USER_EXPLICIT_MATH[with_aocl]}" ]] && CONFIG_CACHE["with_aocl"]="__DONTUSE__"
                [[ -z "${USER_EXPLICIT_MATH[with_openblas]}" ]] && CONFIG_CACHE["with_openblas"]="__DONTUSE__"
                ;;
        esac
    fi
}