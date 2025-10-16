#!/bin/bash

# ============================================================================
# ABACUS Toolchain User Interface Module
# ============================================================================
# Provides consistent user interaction, help messages, and progress display
# Author: Quantum Misaka by Trae SOLO
# Date: 2025-10-16
# ============================================================================

# Global UI variables
UI_INITIALIZED=false
UI_VERBOSE=false
UI_QUIET=false
UI_LOG_FILE=""

# Color definitions for better user experience
if [[ -t 1 ]]; then
    # Terminal supports colors
    readonly UI_RED='\033[0;31m'
    readonly UI_GREEN='\033[0;32m'
    readonly UI_YELLOW='\033[0;33m'
    readonly UI_BLUE='\033[0;34m'
    readonly UI_PURPLE='\033[0;35m'
    readonly UI_CYAN='\033[0;36m'
    readonly UI_WHITE='\033[0;37m'
    readonly UI_BOLD='\033[1m'
    readonly UI_NC='\033[0m' # No Color
else
    # No color support
    readonly UI_RED=''
    readonly UI_GREEN=''
    readonly UI_YELLOW=''
    readonly UI_BLUE=''
    readonly UI_PURPLE=''
    readonly UI_CYAN=''
    readonly UI_WHITE=''
    readonly UI_BOLD=''
    readonly UI_NC=''
fi

# Initialize user interface
# Usage: ui_init [verbose] [quiet] [log_file]
ui_init() {
    if [[ "$UI_INITIALIZED" == "true" ]]; then
        return 0
    fi
    
    # Parse options
    while [[ $# -gt 0 ]]; do
        case "$1" in
            verbose)
                UI_VERBOSE=true
                ;;
            quiet)
                UI_QUIET=true
                ;;
            --log-file=*)
                UI_LOG_FILE="${1#*=}"
                ;;
            *)
                # Assume it's a log file path
                UI_LOG_FILE="$1"
                ;;
        esac
        shift
    done
    
    UI_INITIALIZED=true
    return 0
}

# Print colored message
# Usage: ui_print_color "color" "message"
ui_print_color() {
    local color="$1"
    local message="$2"
    
    if [[ "$UI_QUIET" == "true" ]]; then
        return 0
    fi
    
    echo -e "${color}${message}${UI_NC}"
    
    # Log to file if specified
    if [[ -n "$UI_LOG_FILE" ]]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] $message" >> "$UI_LOG_FILE"
    fi
}

# Print info message
# Usage: ui_info "message"
ui_info() {
    ui_print_color "$UI_BLUE" "INFO: $1"
}

# Print success message
# Usage: ui_success "message"
ui_success() {
    ui_print_color "$UI_GREEN" "SUCCESS: $1"
}

# Print warning message
# Usage: ui_warning "message"
ui_warning() {
    ui_print_color "$UI_YELLOW" "WARNING: $1"
}

# Print error message
# Usage: ui_error "message"
ui_error() {
    ui_print_color "$UI_RED" "ERROR: $1" >&2
}

# Print debug message (only in verbose mode)
# Usage: ui_debug "message"
ui_debug() {
    if [[ "$UI_VERBOSE" == "true" ]]; then
        ui_print_color "$UI_PURPLE" "DEBUG: $1"
    fi
}

# Print section header
# Usage: ui_section "title"
ui_section() {
    local title="$1"
    local line=$(printf '=%.0s' {1..60})
    
    if [[ "$UI_QUIET" != "true" ]]; then
        echo ""
        ui_print_color "$UI_BOLD$UI_CYAN" "$line"
        ui_print_color "$UI_BOLD$UI_CYAN" "$title"
        ui_print_color "$UI_BOLD$UI_CYAN" "$line"
        echo ""
    fi
}

# Print progress bar
# Usage: ui_progress "current" "total" "description"
ui_progress() {
    local current="$1"
    local total="$2"
    local description="$3"
    
    if [[ "$UI_QUIET" == "true" ]]; then
        return 0
    fi
    
    local percent=$((current * 100 / total))
    local filled=$((percent / 2))
    local empty=$((50 - filled))
    
    printf "\r${UI_BLUE}Progress: [%s%s] %d%% - %s${UI_NC}" \
        "$(printf '#%.0s' $(seq 1 $filled))" \
        "$(printf ' %.0s' $(seq 1 $empty))" \
        "$percent" \
        "$description"
    
    if [[ $current -eq $total ]]; then
        echo ""
    fi
}

# Show help message
# Usage: ui_show_help
ui_show_help() {
    cat << 'EOF'
ABACUS Toolchain Installation Script

USAGE:
    install_abacus_toolchain.sh [OPTIONS]

DESCRIPTION:
    This script installs the ABACUS toolchain and its dependencies.
    It supports various compilers, MPI implementations, and mathematical libraries.

RECOMMENDED WORKFLOW:
    1. Run with --help to see all options
    2. Use --dry-run to preview what will be installed
    3. Run the actual installation
    4. Source the setup file before building ABACUS

OPTIONS:
    -j <N>                    Number of parallel compilation processes (default: auto-detect)
    --help                    Show this help message
    --dry-run                 Show what would be done without actually doing it
    --pack-run                Only check and install required packages
    --install-all             Install all packages from source (except Intel/AMD)
    
    --mpi-mode <MODE>         MPI implementation to use
                              Options: mpich, openmpi, intelmpi, no
                              Default: auto-detect or mpich
    
    --math-mode <MODE>        Mathematical library to use
                              Options: cray, mkl, openblas, aocl
                              Default: openblas
    
    --gpu-ver <VERSION>       GPU version for optimizations
                              Options: Numeric (6.0, 7.0, 8.0, etc.) or Models (K20X, K40, K80, P100, V100, Mi50, Mi100, Mi250)
                              Default: no
    
    --target-cpu <CPU>        Target CPU architecture
                              Default: native
    
    --log-lines <N>           Number of log lines to show during compilation
                              Default: 20
    
    --enable-<FEATURE>        Enable specific features
                              Features: tsan, cuda, hip, opencl, cray
    
    --with-<PACKAGE>=<MODE>   Package installation mode
                              Modes: install, system, no, <path>
                              Packages: gcc, intel, intel-classic, ifx, amd, flang,
                                       cmake, openmpi, mpich, intelmpi, libxc, fftw,
                                       mkl, aocl, openblas, scalapack, elpa, cereal,
                                       rapidjson, libtorch, libnpy, libri, libcomm

EXAMPLES:
    # Basic installation with OpenMPI and OpenBLAS
    ./install_abacus_toolchain.sh --mpi-mode openmpi --math-mode openblas
    
    # Dry run to see what would be installed
    ./install_abacus_toolchain.sh --dry-run --install-all
    
    # Install with Intel compiler and MKL
    ./install_abacus_toolchain.sh --with-intel=install --math-mode mkl
    
    # GPU-enabled installation
    ./install_abacus_toolchain.sh --enable-cuda --gpu-ver 75

NOTES:
    - The build and install directories can be safely deleted after installation
    - Source the setup file (setup or abacus_env.sh) before running ABACUS
    - Use --dry-run to preview changes before actual installation
    - Check the log files in case of compilation errors

For more information, see the documentation in Details.md
EOF
}

# Show package installation summary
# Usage: ui_show_summary
ui_show_summary() {
    ui_section "Installation Summary"
    
    # Check system requirements first
    if command -v check_system_requirements &> /dev/null; then
        echo "System Requirements Check:"
        echo "========================="
        if check_system_requirements; then
            echo "✅ All system requirements met"
        else
            echo "⚠️  Some system requirements not met (see above)"
        fi
        echo ""
    fi
    
    echo "Configuration:"
    echo "============="
    echo "MPI Mode: $(config_get MPI_MODE)"
    echo "Math Mode: $(config_get MATH_MODE)"
    echo "Target CPU: $(config_get TARGET_CPU)"
    echo "GPU Version: $(config_get GPUVER)"
    echo "Parallel Jobs: $(config_get NPROCS_OVERWRITE)"
    # Only show dry run status if it's enabled (avoid showing "false" or empty)
    if [[ "$(config_get dry_run)" == "__TRUE__" ]]; then
        echo "Dry Run: enabled"
    fi
    # Only show pack run status if it's enabled (avoid showing "false" or empty)
    if [[ "$(config_get PACK_RUN)" == "__TRUE__" ]]; then
        echo "Pack Run: enabled"
    fi
    echo ""
    
    echo "Package Configuration:"
    echo "====================="
    for pkg in ${package_list}; do
        local status=$(config_get "with_${pkg}")
        if [[ "$status" != "__DONTUSE__" ]]; then
            # Convert __INSTALL__ to install for display consistency
            local display_status="$status"
            if [[ "$status" == "__INSTALL__" ]]; then
                display_status="install"
            elif [[ "$status" == "__SYSTEM__" ]]; then
                display_status="system"
            fi
            printf "%-15s: %s\n" "$pkg" "$display_status"
        fi
    done
    echo ""
    
    # Show packages to install
    local install_list=$(package_list_to_install)
    if [[ -n "$install_list" ]]; then
        echo "Packages to be installed from source:"
        for pkg in $install_list; do
            echo "  - $pkg"
        done
        echo ""
    fi
    
    # Show system packages
    local system_packages=""
    for pkg in ${package_list}; do
        if [[ "$(config_get with_${pkg})" == "__SYSTEM__" ]]; then
            system_packages="$system_packages $pkg"
        fi
    done
    
    if [[ -n "$system_packages" ]]; then
        echo "Packages to be used from system:"
        for pkg in $system_packages; do
            echo "  - $pkg"
        done
        echo ""
    fi
}

# Confirm installation with user
# Usage: ui_confirm_installation
ui_confirm_installation() {
    # Skip confirmation in dry-run mode
    if [[ "$(config_get dry_run)" == "__TRUE__" ]]; then
        ui_info "Dry run mode - no actual installation will be performed"
        return 0
    fi
    
    # Skip confirmation in quiet mode
    if [[ "$UI_QUIET" == "true" ]]; then
        return 0
    fi
    
    # Simple confirmation prompt
    echo ""
    echo -n "Proceed with installation? [y/N]: "
    read -r response
    
    case "$response" in
        [yY]|[yY][eE][sS])
            return 0
            ;;
        *)
            ui_info "Installation cancelled"
            return 1
            ;;
    esac
}

# Show installation progress for a stage
# Usage: ui_stage_progress "stage_number" "stage_name"
ui_stage_progress() {
    local stage="$1"
    local name="$2"
    
    ui_section "Stage $stage: $name"
    ui_info "Installing packages for stage $stage..."
}

# Show package build progress
# Usage: ui_package_progress "package_name" "action"
ui_package_progress() {
    local package="$1"
    local action="$2"
    
    case "$action" in
        start)
            ui_info "Building $package..."
            ;;
        download)
            ui_debug "Downloading $package..."
            ;;
        extract)
            ui_debug "Extracting $package..."
            ;;
        configure)
            ui_debug "Configuring $package..."
            ;;
        compile)
            ui_debug "Compiling $package..."
            ;;
        install)
            ui_debug "Installing $package..."
            ;;
        success)
            ui_success "Successfully built $package"
            ;;
        skip)
            ui_info "Skipping $package (already built or disabled)"
            ;;
        error)
            ui_error "Failed to build $package"
            ;;
    esac
}

# Show final installation results
# Usage: ui_show_results "success_count" "total_count" "failed_packages"
ui_show_results() {
    local success_count="$1"
    local total_count="$2"
    local failed_packages="$3"
    
    ui_section "Installation Results"
    
    if [[ $success_count -eq $total_count ]]; then
        ui_success "All packages installed successfully! ($success_count/$total_count)"
        echo ""
        ui_info "To use the installed toolchain:"
        echo "  source ${SETUPFILE:-setup}"
        echo "  # or"
        echo "  source ${INSTALLDIR:-install}/abacus_env.sh"
        echo ""
        ui_info "You can now build ABACUS with the installed toolchain."
    else
        local failed_count=$((total_count - success_count))
        ui_error "Installation completed with errors ($success_count/$total_count successful)"
        echo ""
        if [[ -n "$failed_packages" ]]; then
            ui_error "Failed packages: $failed_packages"
        fi
        echo ""
        ui_info "Check the log files for detailed error information."
        ui_info "You may need to install missing dependencies or fix configuration issues."
    fi
}

# Show environment setup instructions
# Usage: ui_show_env_setup
ui_show_env_setup() {
    ui_section "Environment Setup"
    
    echo "To use the installed ABACUS toolchain, run one of the following:"
    echo ""
    ui_print_color "$UI_GREEN" "  source ${SETUPFILE:-setup}"
    echo "    or"
    ui_print_color "$UI_GREEN" "  source ${INSTALLDIR:-install}/abacus_env.sh"
    echo ""
    echo "Then you can build ABACUS with:"
    ui_print_color "$UI_GREEN" "  make -j\$(nproc)"
    echo ""
    ui_info "The environment setup needs to be done in each new shell session."
}

# Handle user interruption (Ctrl+C)
# Usage: ui_handle_interrupt
ui_handle_interrupt() {
    echo ""
    ui_warning "Installation interrupted by user"
    ui_info "Cleaning up temporary files..."
    
    # Clean up any temporary files or processes
    if [[ -n "$BUILDDIR" && -d "$BUILDDIR" ]]; then
        ui_debug "Cleaning build directory: $BUILDDIR"
        # Don't remove the entire build directory, just mark as interrupted
        touch "$BUILDDIR/.interrupted"
    fi
    
    ui_info "Installation cancelled"
    exit 130
}

# Set up signal handlers
# Usage: ui_setup_signals
ui_setup_signals() {
    trap ui_handle_interrupt SIGINT SIGTERM
}

# Validate user input
# Usage: ui_validate_input "input" "type"
ui_validate_input() {
    local input="$1"
    local type="$2"
    
    case "$type" in
        number)
            if [[ ! "$input" =~ ^[0-9]+$ ]]; then
                ui_error "Invalid number: $input"
                return 1
            fi
            ;;
        path)
            if [[ ! -d "$input" ]]; then
                ui_error "Directory does not exist: $input"
                return 1
            fi
            ;;
        file)
            if [[ ! -f "$input" ]]; then
                ui_error "File does not exist: $input"
                return 1
            fi
            ;;
        mpi_mode)
            case "$input" in
                mpich|openmpi|intelmpi|no)
                    return 0
                    ;;
                *)
                    ui_error "Invalid MPI mode: $input"
                    ui_info "Valid options: mpich, openmpi, intelmpi, no"
                    return 1
                    ;;
            esac
            ;;
        math_mode)
            case "$input" in
                cray|mkl|openblas|aocl)
                    return 0
                    ;;
                *)
                    ui_error "Invalid math mode: $input"
                    ui_info "Valid options: cray, mkl, openblas, aocl"
                    return 1
                    ;;
            esac
            ;;
        gpu_version)
            # Support only numeric formats for GPU versions
            if [[ "$input" == "no" ]]; then
                return 0
            fi
            
            # Check if it's a valid numeric format (like 8.0, 70, 80, etc.)
            local arch_num="${input//.}"
            if [[ "$arch_num" =~ ^[1-9][0-9]*$ ]]; then
                return 0
            fi
            
            # Invalid format - show error message
            ui_error "Invalid GPU version: $input"
            ui_info "Valid formats:"
            ui_info "  - Numeric with decimal: 6.0, 7.0, 8.0, 8.9, etc. (CUDA compute capability)"
            ui_info "  - Numeric without decimal: 60, 70, 80, 89, etc."
            ui_info "  - Use 'no' to disable GPU support"
            ui_info "Examples: 8.0 or 80 for compute capability 8.0"
            return 1
            ;;
        *)
            ui_error "Unknown validation type: $type"
            return 1
            ;;
    esac
    
    return 0
}

# Print system information
# Usage: ui_show_system_info
ui_show_system_info() {
    ui_section "System Information"
    
    echo "Operating System: $(uname -s)"
    echo "Architecture: $(uname -m)"
    echo "Kernel: $(uname -r)"
    echo "CPU Cores: $(nproc 2>/dev/null || echo "unknown")"
    
    if command -v free &> /dev/null; then
        local mem_gb=$(free -g | awk '/^Mem:/ {print $2}')
        echo "Memory: ${mem_gb}GB"
    fi
    
    echo "Shell: $SHELL"
    echo "User: $(whoami)"
    echo "Working Directory: $(pwd)"
    echo ""
}

# Check system requirements
# Usage: ui_check_system_requirements
ui_check_system_requirements() {
    ui_section "System Requirements Check"
    
    local missing_tools=""
    local required_tools="wget curl tar gzip make"
    
    for tool in $required_tools; do
        if ! command -v "$tool" &> /dev/null; then
            missing_tools="$missing_tools $tool"
        else
            ui_debug "Found: $tool"
        fi
    done
    
    if [[ -n "$missing_tools" ]]; then
        ui_error "Missing required system tools:$missing_tools"
        ui_info "Please install these tools using your system package manager:"
        ui_info "  Ubuntu/Debian: sudo apt-get install$missing_tools"
        ui_info "  CentOS/RHEL: sudo yum install$missing_tools"
        ui_info "  Fedora: sudo dnf install$missing_tools"
        return 1
    else
        ui_success "All required system tools are available"
        return 0
    fi
}