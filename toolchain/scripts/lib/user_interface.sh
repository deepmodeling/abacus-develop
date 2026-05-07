#!/bin/bash

# ============================================================================
# ABACUS Toolchain User Interface Module (Enhanced)
# ============================================================================
# Provides beautiful and consistent user interaction, help messages, and progress display
# Author: Quantum Misaka by Trae SOLO
# Date: 2025-10-16
# Enhanced: Beautiful terminal output with harmonious color scheme
# ============================================================================

# Global UI variables
UI_INITIALIZED=false
UI_VERBOSE=false
UI_QUIET=false
UI_LOG_FILE=""

# Enhanced color definitions with harmonious blue-based theme
if [[ -t 1 ]]; then
    # Terminal supports colors - Enhanced color palette
    readonly UI_RED='\033[38;5;196m'        # Bright red for errors
    readonly UI_GREEN='\033[38;5;46m'       # Bright green for success
    readonly UI_YELLOW='\033[38;5;226m'     # Bright yellow for warnings
    readonly UI_BLUE='\033[38;5;39m'        # Bright blue for info (main theme)
    readonly UI_PURPLE='\033[38;5;141m'     # Purple for debug
    readonly UI_CYAN='\033[38;5;51m'        # Cyan for progress
    readonly UI_WHITE='\033[38;5;255m'      # Pure white
    readonly UI_GRAY='\033[38;5;244m'       # Gray for secondary text
    readonly UI_ORANGE='\033[38;5;208m'     # Orange for highlights
    
    # Style definitions
    readonly UI_BOLD='\033[1m'
    readonly UI_DIM='\033[2m'
    readonly UI_ITALIC='\033[3m'
    readonly UI_UNDERLINE='\033[4m'
    readonly UI_BLINK='\033[5m'
    readonly UI_REVERSE='\033[7m'
    readonly UI_NC='\033[0m'                # No Color
    
    # Background colors for special effects
    readonly UI_BG_RED='\033[48;5;196m'
    readonly UI_BG_GREEN='\033[48;5;46m'
    readonly UI_BG_BLUE='\033[48;5;39m'
    readonly UI_BG_YELLOW='\033[48;5;226m'
else
    # No color support - fallback to empty strings
    readonly UI_RED=''
    readonly UI_GREEN=''
    readonly UI_YELLOW=''
    readonly UI_BLUE=''
    readonly UI_PURPLE=''
    readonly UI_CYAN=''
    readonly UI_WHITE=''
    readonly UI_GRAY=''
    readonly UI_ORANGE=''
    readonly UI_BOLD=''
    readonly UI_DIM=''
    readonly UI_ITALIC=''
    readonly UI_UNDERLINE=''
    readonly UI_BLINK=''
    readonly UI_REVERSE=''
    readonly UI_NC=''
    readonly UI_BG_RED=''
    readonly UI_BG_GREEN=''
    readonly UI_BG_BLUE=''
    readonly UI_BG_YELLOW=''
fi

# Terminal Unicode support detection
UI_UNICODE_SUPPORT=""

# Detect terminal Unicode support capability
# Returns: "full", "basic", or "none"
ui_detect_unicode_support() {
    # Check if already detected
    if [[ -n "$UI_UNICODE_SUPPORT" ]]; then
        echo "$UI_UNICODE_SUPPORT"
        return 0
    fi
    
    # Manual override via environment variable
    if [[ "${ABACUS_UI_SIMPLE:-}" == "1" ]] || [[ "${ABACUS_UI_ASCII:-}" == "1" ]]; then
        UI_UNICODE_SUPPORT="none"
        echo "none"
        return 0
    fi
    
    # Force Unicode mode if explicitly requested
    if [[ "${ABACUS_UI_UNICODE:-}" == "1" ]]; then
        UI_UNICODE_SUPPORT="full"
        echo "full"
        return 0
    fi
    
    # Check if we're in a non-interactive environment
    if [[ ! -t 1 ]]; then
        UI_UNICODE_SUPPORT="none"
        echo "none"
        return 0
    fi
    
    # Check locale settings for UTF-8 support
    local locale_utf8=false
    
    # First check environment variables
    if [[ "${LC_ALL:-}" =~ [Uu][Tt][Ff]-?8 ]] || \
       [[ "${LC_CTYPE:-}" =~ [Uu][Tt][Ff]-?8 ]] || \
       [[ "${LANG:-}" =~ [Uu][Tt][Ff]-?8 ]]; then
        locale_utf8=true
    fi
    
    # Additional check using locale command if available
    if [[ "$locale_utf8" == "false" ]] && command -v locale &>/dev/null; then
        local charset=$(locale charmap 2>/dev/null || echo "")
        if [[ "$charset" =~ UTF-?8 ]]; then
            locale_utf8=true
        fi
    fi
    
    # Check terminal type and capabilities
    local term_support="none"
    case "${TERM:-}" in
        xterm*|screen*|tmux*|rxvt*|gnome*|konsole*|alacritty*|kitty*|iterm*)
            if [[ "$locale_utf8" == "true" ]]; then
                term_support="full"
            else
                term_support="basic"
            fi
            ;;
        linux|vt*)
            # Linux console - basic Unicode support
            if [[ "$locale_utf8" == "true" ]]; then
                term_support="basic"
            else
                term_support="none"
            fi
            ;;
        dumb|unknown)
            term_support="none"
            ;;
        *)
            # Unknown terminal - be conservative
            if [[ "$locale_utf8" == "true" ]]; then
                term_support="basic"
            else
                term_support="none"
            fi
            ;;
    esac
    
    # Final validation - ensure we have both locale and terminal support
    if [[ "$term_support" != "none" && "$locale_utf8" == "false" ]]; then
        # Terminal claims support but locale doesn't - be conservative
        term_support="none"
    fi
    
    UI_UNICODE_SUPPORT="$term_support"
    echo "$term_support"
}

# Unicode symbols with ASCII fallbacks
# These will be set based on terminal capability
UI_ICON_SUCCESS=""
UI_ICON_ERROR=""
UI_ICON_WARNING=""
UI_ICON_INFO=""
UI_ICON_PROGRESS=""
UI_ICON_ROCKET=""
UI_ICON_GEAR=""
UI_ICON_PACKAGE=""
UI_ICON_DOWNLOAD=""
UI_ICON_BUILD=""
UI_ICON_INSTALL=""
UI_ICON_CHECK=""
UI_ICON_CROSS=""
UI_ICON_ARROW=""
UI_ICON_STAR=""

# Progress bar characters
UI_PROGRESS_FULL=""
UI_PROGRESS_PARTIAL1=""
UI_PROGRESS_PARTIAL2=""
UI_PROGRESS_EMPTY=""

# Initialize Unicode icons based on terminal capability
ui_init_icons() {
    local support=$(ui_detect_unicode_support)
    
    case "$support" in
        "full")
            # Full Unicode support - use emoji and special symbols
            readonly UI_ICON_SUCCESS="✅"
            readonly UI_ICON_ERROR="❌"
            readonly UI_ICON_WARNING="⚠️"
            readonly UI_ICON_INFO="ℹ️"
            readonly UI_ICON_PROGRESS="🔄"
            readonly UI_ICON_ROCKET="🚀"
            readonly UI_ICON_GEAR="⚙️"
            readonly UI_ICON_PACKAGE="📦"
            readonly UI_ICON_DOWNLOAD="⬇️"
            readonly UI_ICON_BUILD="🔨"
            readonly UI_ICON_INSTALL="📥"
            readonly UI_ICON_CHECK="✓"
            readonly UI_ICON_CROSS="✗"
            readonly UI_ICON_ARROW="→"
            readonly UI_ICON_STAR="⭐"
            # Full Unicode progress bar
            readonly UI_PROGRESS_FULL="█"
            readonly UI_PROGRESS_PARTIAL1="▓"
            readonly UI_PROGRESS_PARTIAL2="▒"
            readonly UI_PROGRESS_EMPTY="░"
            ;;
        "basic")
            # Basic Unicode support - use simple Unicode symbols
            readonly UI_ICON_SUCCESS="✓"
            readonly UI_ICON_ERROR="✗"
            readonly UI_ICON_WARNING="!"
            readonly UI_ICON_INFO="i"
            readonly UI_ICON_PROGRESS="*"
            readonly UI_ICON_ROCKET="^"
            readonly UI_ICON_GEAR="*"
            readonly UI_ICON_PACKAGE="+"
            readonly UI_ICON_DOWNLOAD="v"
            readonly UI_ICON_BUILD="#"
            readonly UI_ICON_INSTALL="+"
            readonly UI_ICON_CHECK="✓"
            readonly UI_ICON_CROSS="✗"
            readonly UI_ICON_ARROW="→"
            readonly UI_ICON_STAR="*"
            # Basic Unicode progress bar
            readonly UI_PROGRESS_FULL="█"
            readonly UI_PROGRESS_PARTIAL1="▓"
            readonly UI_PROGRESS_PARTIAL2="▒"
            readonly UI_PROGRESS_EMPTY="░"
            ;;
        *)
            # ASCII fallback - maximum compatibility
            readonly UI_ICON_SUCCESS="[OK]"
            readonly UI_ICON_ERROR="[ERR]"
            readonly UI_ICON_WARNING="[WARN]"
            readonly UI_ICON_INFO="[INFO]"
            readonly UI_ICON_PROGRESS="[*]"
            readonly UI_ICON_ROCKET="[^]"
            readonly UI_ICON_GEAR="[*]"
            readonly UI_ICON_PACKAGE="[+]"
            readonly UI_ICON_DOWNLOAD="[v]"
            readonly UI_ICON_BUILD="[#]"
            readonly UI_ICON_INSTALL="[+]"
            readonly UI_ICON_CHECK="[+]"
            readonly UI_ICON_CROSS="[-]"
            readonly UI_ICON_ARROW="->"
            readonly UI_ICON_STAR="[*]"
            # ASCII progress bar
            readonly UI_PROGRESS_FULL="#"
            readonly UI_PROGRESS_PARTIAL1="#"
            readonly UI_PROGRESS_PARTIAL2="-"
            readonly UI_PROGRESS_EMPTY="."
            ;;
    esac
}

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
    
    # Initialize icons based on terminal capability
    ui_init_icons
    
    UI_INITIALIZED=true
    return 0
}

# Print colored message with enhanced formatting
# Usage: ui_print_color "color" "message" [icon]
ui_print_color() {
    local color="$1"
    local message="$2"
    local icon="${3:-}"
    
    if [[ "$UI_QUIET" == "true" ]]; then
        return 0
    fi
    
    if [[ -n "$icon" ]]; then
        echo -e "${color}${icon} ${message}${UI_NC}"
    else
        echo -e "${color}${message}${UI_NC}"
    fi
    
    # Log to file if specified
    if [[ -n "$UI_LOG_FILE" ]]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] $message" >> "$UI_LOG_FILE"
    fi
}

# Print info message with enhanced styling
# Usage: ui_info "message"
ui_info() {
    ui_print_color "${UI_BLUE}${UI_BOLD}" "$1" "$UI_ICON_INFO"
}

# Print success message with enhanced styling
# Usage: ui_success "message"
ui_success() {
    ui_print_color "${UI_GREEN}${UI_BOLD}" "$1" "$UI_ICON_SUCCESS"
}

# Print warning message with enhanced styling
# Usage: ui_warning "message"
ui_warning() {
    ui_print_color "${UI_YELLOW}${UI_BOLD}" "$1" "$UI_ICON_WARNING"
}

# Print error message with enhanced styling
# Usage: ui_error "message"
ui_error() {
    ui_print_color "${UI_RED}${UI_BOLD}" "$1" "$UI_ICON_ERROR" >&2
}

# Print debug message (only in verbose mode) with enhanced styling
# Usage: ui_debug "message"
ui_debug() {
    if [[ "$UI_VERBOSE" == "true" ]]; then
        ui_print_color "${UI_PURPLE}${UI_DIM}" "DEBUG: $1" "🔍"
    fi
}

# Detect glibc version using multiple methods for maximum compatibility
# Usage: ui_get_glibc_version
ui_get_glibc_version() {
    local glibc_version="unknown"
    
    # Skip on non-Linux systems
    if [[ "$(uname -s)" != "Linux" ]]; then
        echo "unknown"
        return 0
    fi
    
    # Method 1: ldd --version (most reliable)
    if command -v ldd &> /dev/null; then
        local ldd_output=$(ldd --version 2>/dev/null | head -n1)
        if [[ $? -eq 0 && -n "$ldd_output" ]]; then
            # Extract version from ldd output (e.g., "ldd (Ubuntu GLIBC 2.35-0ubuntu3.1) 2.35")
            glibc_version=$(echo "$ldd_output" | grep -oE '[0-9]+\.[0-9]+(\.[0-9]+)?' | head -n1)
            if [[ -n "$glibc_version" ]]; then
                echo "$glibc_version"
                return 0
            fi
        fi
    fi
    
    # Method 2: Direct library call (fallback)
    local libc_paths=(
        "/lib/x86_64-linux-gnu/libc.so.6"
        "/lib64/libc.so.6"
        "/lib/libc.so.6"
        "/usr/lib/x86_64-linux-gnu/libc.so.6"
        "/usr/lib64/libc.so.6"
    )
    
    for libc_path in "${libc_paths[@]}"; do
        if [[ -x "$libc_path" ]]; then
            local libc_output=$("$libc_path" 2>/dev/null | head -n1)
            if [[ $? -eq 0 && -n "$libc_output" ]]; then
                # Extract version from libc output
                glibc_version=$(echo "$libc_output" | grep -oE '[0-9]+\.[0-9]+(\.[0-9]+)?' | head -n1)
                if [[ -n "$glibc_version" ]]; then
                    echo "$glibc_version"
                    return 0
                fi
            fi
        fi
    done
    
    # Method 3: getconf GNU_LIBC_VERSION (alternative)
    if command -v getconf &> /dev/null; then
        local getconf_output=$(getconf GNU_LIBC_VERSION 2>/dev/null)
        if [[ $? -eq 0 && -n "$getconf_output" ]]; then
            # Extract version from getconf output (e.g., "glibc 2.35")
            glibc_version=$(echo "$getconf_output" | grep -oE '[0-9]+\.[0-9]+(\.[0-9]+)?' | head -n1)
            if [[ -n "$glibc_version" ]]; then
                echo "$glibc_version"
                return 0
            fi
        fi
    fi
    
    # Method 4: Check /proc/version_signature or similar (last resort)
    if [[ -r "/proc/version" ]]; then
        local proc_version=$(cat /proc/version 2>/dev/null)
        if [[ -n "$proc_version" ]]; then
            # Look for glibc version in kernel version string (rare but possible)
            glibc_version=$(echo "$proc_version" | grep -oE 'glibc[[:space:]]*[0-9]+\.[0-9]+(\.[0-9]+)?' | grep -oE '[0-9]+\.[0-9]+(\.[0-9]+)?' | head -n1)
            if [[ -n "$glibc_version" ]]; then
                echo "$glibc_version"
                return 0
            fi
        fi
    fi
    
    # All methods failed - return unknown
    echo "unknown"
    return 0
}

# Read version from VERSION file
# Usage: ui_get_version
ui_get_version() {
    local version_file
    local version="Unknown"
    
    # Try to find VERSION file in multiple locations
    if [[ -n "${SCRIPT_DIR:-}" && -f "${SCRIPT_DIR}/VERSION" ]]; then
        version_file="${SCRIPT_DIR}/VERSION"
    elif [[ -f "$(dirname "${BASH_SOURCE[0]}")/../VERSION" ]]; then
        version_file="$(dirname "${BASH_SOURCE[0]}")/../VERSION"
    elif [[ -f "./scripts/VERSION" ]]; then
        version_file="./scripts/VERSION"
    elif [[ -f "./VERSION" ]]; then
        version_file="./VERSION"
    fi
    
    if [[ -n "$version_file" && -f "$version_file" ]]; then
        # Extract version from VERSION file, handling different formats
        version=$(grep -E '^VERSION=' "$version_file" 2>/dev/null | cut -d'=' -f2 | tr -d '"' | tr -d "'" | xargs)
        if [[ -z "$version" ]]; then
            # Fallback: try to read the entire file content if no VERSION= line
            version=$(head -n 1 "$version_file" 2>/dev/null | grep -v '^#' | xargs)
        fi
        if [[ -z "$version" ]]; then
            version="Unknown"
        fi
    fi
    
    echo "$version"
}

# Print beautiful welcome banner
# Usage: ui_welcome_banner
ui_welcome_banner() {
    if [[ "$UI_QUIET" == "true" ]]; then
        return 0
    fi
    
    local banner_width=80
    local title="ABACUS Toolchain Installation"
    local subtitle="DFT Calculation Made Simple"
    local version=$(ui_get_version)
    local version_text="Version $version"
    
    echo ""
    ui_print_color "${UI_CYAN}${UI_BOLD}" "$(printf '═%.0s' $(seq 1 $banner_width))"
    echo ""
    ui_print_color "${UI_BLUE}${UI_BOLD}" "$(printf '%*s' $(((banner_width + ${#title}) / 2)) "$title")"
    ui_print_color "${UI_GRAY}${UI_ITALIC}" "$(printf '%*s' $(((banner_width + ${#subtitle}) / 2)) "$subtitle")"
    ui_print_color "${UI_GRAY}${UI_DIM}" "$(printf '%*s' $(((banner_width + ${#version_text}) / 2)) "$version_text")"
    echo ""
    ui_print_color "${UI_CYAN}${UI_BOLD}" "$(printf '═%.0s' $(seq 1 $banner_width))"
    echo ""
}

# Print beautiful section header with enhanced styling
# Usage: ui_section "title"
ui_section() {
    local title="$1"
    local line_width=60
    
    if [[ "$UI_QUIET" != "true" ]]; then
        echo ""
        ui_print_color "${UI_CYAN}${UI_BOLD}" "$(printf '─%.0s' $(seq 1 $line_width))"
        ui_print_color "${UI_BLUE}${UI_BOLD}" "${UI_ICON_GEAR} $title"
        ui_print_color "${UI_CYAN}${UI_BOLD}" "$(printf '─%.0s' $(seq 1 $line_width))"
        echo ""
    fi
}

# Print enhanced progress bar with percentage and ETA
# Usage: ui_progress "current" "total" "description" [eta_seconds]
ui_progress() {
    local current="$1"
    local total="$2"
    local description="$3"
    local eta_seconds="${4:-}"
    
    if [[ "$UI_QUIET" == "true" ]]; then
        return 0
    fi
    
    local percent=$((current * 100 / total))
    local filled=$((percent / 2))
    local empty=$((50 - filled))
    
    # Create progress bar with gradient effect using detected characters
    local progress_bar=""
    for ((i=1; i<=filled; i++)); do
        if [[ $i -le 10 ]]; then
            progress_bar+="$UI_PROGRESS_FULL"
        elif [[ $i -le 30 ]]; then
            progress_bar+="$UI_PROGRESS_PARTIAL1"
        else
            progress_bar+="$UI_PROGRESS_PARTIAL2"
        fi
    done
    
    for ((i=1; i<=empty; i++)); do
        progress_bar+="$UI_PROGRESS_EMPTY"
    done
    
    # Format ETA if provided
    local eta_text=""
    if [[ -n "$eta_seconds" && "$eta_seconds" -gt 0 ]]; then
        local eta_min=$((eta_seconds / 60))
        local eta_sec=$((eta_seconds % 60))
        eta_text=" (ETA: ${eta_min}m${eta_sec}s)"
    fi
    
    printf "\r${UI_CYAN}${UI_ICON_PROGRESS} Progress: ${UI_BOLD}[${progress_bar}]${UI_NC} ${UI_GREEN}${percent}%%${UI_NC} - ${UI_BLUE}%s${UI_NC}${eta_text}" "$description"
    
    if [[ $current -eq $total ]]; then
        echo ""
        ui_success "Task completed successfully!"
    fi
}

# Show enhanced help message with beautiful formatting
# Usage: ui_show_help
ui_show_help() {
    ui_welcome_banner
    
    cat << 'EOF'
USAGE:
    install_abacus_toolchain_new.sh [OPTIONS]

DESCRIPTION:
    This script installs the ABACUS toolchain and its dependencies with a beautiful
    and user-friendly interface. It supports various compilers, MPI implementations,
    and mathematical libraries for density functional theory calculations.

RECOMMENDED WORKFLOW:
    1. 🔍 Run with --help to see all available options
    2. 🧪 Use --dry-run to preview what will be installed
    3. 🚀 Run the actual installation
    4. ⚙️  Source the setup file before building ABACUS

BASIC OPTIONS:
    -h, --help                Show this comprehensive help message
    --version                 Show version information
    --version-info [PACKAGE]  Show version information for specific package or all
    --dry-run                 Preview installation without executing (recommended first)
    --pack-run                Only download packages without building
    
CONFIGURATION OPTIONS:
    --config-file <FILE>      Load configuration from specified file
                              🎯 Default: ./install_abacus_toolchain.conf
    
    --mpi-mode <MODE>         MPI implementation to use
                              📋 Options: mpich, openmpi, intelmpi, no
                              🎯 Default: mpich (auto-detect available)
    
    --math-mode <MODE>        Mathematical library to use
                              📋 Options: mkl, aocl, openblas, cray, no
                              🎯 Default: openblas
    
    --gpu-ver <VERSION>       GPU version for ELPA (CUDA compute capability)
                              📋 Options: Numeric (7.0, 7.5, 8.0, etc.) or (70, 75, 80, etc.)
                              🎯 Default: no (CPU-only)
    
    --target-cpu <CPU>        Target CPU architecture for optimizations
                              🎯 Default: native (auto-detect)
    
    --log-lines <N>           Number of log lines to show during compilation
                              🎯 Default: 200

PACKAGE CONTROL OPTIONS:
    --package-version <PKG:VER>  Set package version strategy
                              📋 Format: package:version (e.g., openmpi:alt, openblas:main)
                              📋 Versions: main (latest stable), alt (alternative/legacy)
    
    --with-<PACKAGE>=<MODE>   Fine-tune package installation modes
                              📋 Modes: install, system, no, <custom_path>
                              📦 Build Tools: gcc, cmake
                              📦 Compilers: intel, intel-classic, ifx, amd, flang
                              📦 MPI Libraries: openmpi, mpich, intelmpi
                              📦 Math Libraries: mkl, aocl, openblas, scalapack
                              📦 Scientific: libxc, fftw, elpa
                              📦 Advanced: cereal, rapidjson, libri, libcomm, libtorch, libnpy, nep

ADVANCED OPTIONS:
    --enable-<FEATURE>[=yes/no]  Enable specific advanced features
                              📋 Features: tsan, cuda, hip, opencl, cray
                              🎯 Default: no for all features
    
    --with-intel-classic[=yes/no]     Use Intel Classic compiler (icc/ifort)
                              🎯 Default: no (uses OneAPI icx/ifx)
    
    --with-intel-mpi-classic[=yes/no] Use Intel MPI Classic
                              🎯 Default: no
    
    --with-ifx[=yes/no]       Use Intel Fortran compiler (ifx)
                              🎯 Default: yes (when Intel is enabled)
    
    --with-flang[=yes/no]     Use AMD Flang Fortran compiler
                              🎯 Default: no
    
    --with-4th-openmpi[=yes/no]      Use OpenMPI 4th generation (v4.x)
                              🎯 Default: no (uses v5.x)
    
    --with-mpich-device=<DEV> MPICH device type
                              📋 Options: ch3, ch4
                              🎯 Default: ch4
    
    --skip-system-checks      Skip system validation checks

EXAMPLES:
    # 🎯 Basic installation with OpenMPI and OpenBLAS
    ./install_abacus_toolchain_new.sh --mpi-mode openmpi --math-mode openblas
    
    # 🧪 Preview installation with all packages
    ./install_abacus_toolchain_new.sh --dry-run --mpi-mode mpich
    
    # 🏭 Intel compiler with MKL (high performance)
    ./install_abacus_toolchain_new.sh --with-intel=install --math-mode mkl
    
    # 🎮 GPU-enabled installation for CUDA compute capability 8.0
    ./install_abacus_toolchain_new.sh --enable-cuda --gpu-ver 8.0
    
    # 🔧 Custom configuration with specific package versions
    ./install_abacus_toolchain_new.sh --package-version openmpi:alt --with-fftw=system
    
    # 📁 Load configuration from file
    ./install_abacus_toolchain_new.sh --config-file my_config.conf --dry-run
    
    # 🚀 Use pre-configured toolchain scripts (recommended)
    ./toolchain_gnu.sh         # GNU toolchain (GCC + OpenMPI + OpenBLAS)
    ./toolchain_intel.sh       # Intel toolchain (Intel + MPI + MKL)
    ./toolchain_gcc-aocl.sh    # GCC + AMD AOCL
    ./toolchain_aocc-aocl.sh   # AMD AOCC + AOCL

ENVIRONMENT VARIABLES:
    NPROCS_OVERWRITE=N        Override parallel compilation processes
    DOWNLOAD_CERT_POLICY      Certificate verification policy (strict/smart/skip)
                              🎯 Default: smart (try secure, fallback if needed)

NOTES:
    📁 Build and install directories can be safely deleted after installation
    🔧 Source the setup file (install/setup) before building ABACUS
    🧪 Always use --dry-run first to preview changes
    📋 Check log files in build/PKG_NAME/make.log for compilation errors
    💡 For detailed information, see the documentation in README.md
    🎛️  Configuration files allow saving and reusing complex setups
    🚀 Use toolchain_*.sh scripts for easier pre-configured installations

EOF
    
    ui_print_color "${UI_ORANGE}${UI_BOLD}" "${UI_ICON_STAR} For the best experience, start with: ${UI_WHITE}./install_abacus_toolchain_new.sh --dry-run${UI_NC}"
    echo ""
}

# Show enhanced package installation summary with beautiful table formatting
# Usage: ui_show_summary
ui_show_summary() {
    ui_section "Installation Configuration Summary"
    
    # System information box
    ui_print_color "${UI_BLUE}${UI_BOLD}" "🖥️  System Information:"
    echo "   ├─ OS: $(uname -s) $(uname -m)"
    echo "   ├─ Kernel: $(uname -r)"
    echo "   ├─ glibc: $(ui_get_glibc_version)"
    echo "   ├─ CPU Cores: $(nproc 2>/dev/null || echo "unknown")"
    # CPU model detection
    local cpu_model="unknown"
    if command -v lscpu &> /dev/null; then
        cpu_model=$(LC_ALL=C lscpu | awk -F: '/^Model name/{print $2}' | sed 's/^[[:space:]]*//')
    fi
    if [[ -z "$cpu_model" || "$cpu_model" == "unknown" ]] && [[ -r /proc/cpuinfo ]]; then
        cpu_model=$(awk -F: '/model name|Hardware|Processor/{print $2; exit}' /proc/cpuinfo | sed 's/^[[:space:]]*//')
    fi
    echo "   ├─ CPU Model: ${cpu_model}"
    if command -v free &> /dev/null; then
        local mem_gb=$(free -g | awk '/^Mem:/ {print $2}')
        echo "   ├─ Memory: ${mem_gb}GB"
    else
        echo "   ├─ Memory: unknown"
    fi
    
    # GPU detection with comprehensive vendor support
    local gpu_info="no GPU detected"
    local gpu_count=0
    local gpu_models=""
    
    # Method 1: Try nvidia-smi for NVIDIA GPUs
    if command -v nvidia-smi &> /dev/null; then
        local nvidia_output=$(nvidia-smi --query-gpu=name --format=csv,noheader,nounits 2>/dev/null)
        if [[ $? -eq 0 && -n "$nvidia_output" ]]; then
            gpu_count=$(echo "$nvidia_output" | wc -l)
            if [[ $gpu_count -gt 0 ]]; then
                local first_gpu=$(echo "$nvidia_output" | head -n1 | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
                if [[ $gpu_count -eq 1 ]]; then
                    gpu_info="$first_gpu"
                else
                    gpu_info="$first_gpu (${gpu_count} devices)"
                fi
            fi
        fi
    fi
    
    # Method 2: Try rocm-smi for AMD GPUs (if no NVIDIA found)
    if [[ "$gpu_info" == "no GPU detected" ]] && command -v rocm-smi &> /dev/null; then
        local amd_output=$(rocm-smi --showproductname 2>/dev/null | grep -E "GPU\[" | head -n1)
        if [[ $? -eq 0 && -n "$amd_output" ]]; then
            local amd_name=$(echo "$amd_output" | sed -n 's/.*: \(.*\)/\1/p' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
            if [[ -n "$amd_name" ]]; then
                local amd_count=$(rocm-smi --showproductname 2>/dev/null | grep -c "GPU\[" || echo "1")
                if [[ $amd_count -eq 1 ]]; then
                    gpu_info="$amd_name"
                else
                    gpu_info="$amd_name (${amd_count} devices)"
                fi
            fi
        fi
    fi
    
    # Method 3: Try lspci as fallback for any GPU vendor
    if [[ "$gpu_info" == "no GPU detected" ]] && command -v lspci &> /dev/null; then
        local pci_gpus=$(lspci 2>/dev/null | grep -i "vga\|3d\|display" | grep -v "audio")
        if [[ -n "$pci_gpus" ]]; then
            gpu_count=$(echo "$pci_gpus" | wc -l)
            local first_gpu_line=$(echo "$pci_gpus" | head -n1)
            
            # Extract GPU name from lspci output
            local gpu_name=""
            if echo "$first_gpu_line" | grep -qi "nvidia"; then
                gpu_name=$(echo "$first_gpu_line" | sed -n 's/.*NVIDIA Corporation \(.*\) (rev.*/\1/p' | sed 's/\[.*\]//g' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
                [[ -z "$gpu_name" ]] && gpu_name=$(echo "$first_gpu_line" | sed -n 's/.*NVIDIA \(.*\)/\1/p' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
                [[ -z "$gpu_name" ]] && gpu_name="NVIDIA GPU"
            elif echo "$first_gpu_line" | grep -qi "amd\|ati\|radeon"; then
                gpu_name=$(echo "$first_gpu_line" | sed -n 's/.*Advanced Micro Devices.*\[\(.*\)\].*/\1/p' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
                [[ -z "$gpu_name" ]] && gpu_name=$(echo "$first_gpu_line" | sed -n 's/.*AMD\/ATI \(.*\)/\1/p' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
                [[ -z "$gpu_name" ]] && gpu_name="AMD GPU"
            elif echo "$first_gpu_line" | grep -qi "intel"; then
                gpu_name=$(echo "$first_gpu_line" | sed -n 's/.*Intel Corporation \(.*\) (rev.*/\1/p' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
                [[ -z "$gpu_name" ]] && gpu_name="Intel GPU"
            else
                gpu_name="Unknown GPU"
            fi
            
            if [[ -n "$gpu_name" && "$gpu_name" != "Unknown GPU" ]]; then
                if [[ $gpu_count -eq 1 ]]; then
                    gpu_info="$gpu_name"
                else
                    gpu_info="$gpu_name (${gpu_count} devices)"
                fi
            fi
        fi
    fi
    
    echo "   └─ GPU: $gpu_info"
    echo ""
    
    # Configuration box with aligned formatting
    ui_print_color "${UI_BLUE}${UI_BOLD}" "⚙️  Build Configuration:"
    printf "   ├─ %-15s %s\n" "MPI Mode:" "$(config_get MPI_MODE)"
    printf "   ├─ %-15s %s\n" "Math Mode:" "$(config_get MATH_MODE)"
    printf "   ├─ %-15s %s\n" "Target CPU:" "$(config_get TARGET_CPU)"
    printf "   ├─ %-15s %s\n" "GPU Version:" "$(config_get GPUVER)"
    printf "   └─ %-15s %s\n" "Parallel Jobs:" "$(config_get NPROCS_OVERWRITE)"
    
    # Special modes
    if [[ "$(config_get dry_run)" == "__TRUE__" ]]; then
        ui_print_color "${UI_YELLOW}${UI_BOLD}" "   🧪 Dry Run Mode: Preview only (no actual installation)"
    fi
    if [[ "$(config_get PACK_RUN)" == "__TRUE__" ]]; then
        ui_print_color "${UI_YELLOW}${UI_BOLD}" "   📦 Pack Run Mode: System packages only"
    fi
    echo ""
    
    # Package configuration table
    ui_print_color "${UI_BLUE}${UI_BOLD}" "📦 Package Configuration:"
    local has_packages=false
    for pkg in ${package_list}; do
        local status=$(config_get "with_${pkg}")
        if [[ "$status" != "__DONTUSE__" ]]; then
            has_packages=true
            # Convert status for display
            local display_status="$status"
            local status_icon=""
            case "$status" in
                "__INSTALL__")
                    display_status="${UI_GREEN}install${UI_NC}"
                    status_icon="📥"
                    ;;
                "__SYSTEM__")
                    display_status="${UI_BLUE}system${UI_NC}"
                    status_icon="🔗"
                    ;;
                "__DONTUSE__")
                    display_status="${UI_GRAY}disabled${UI_NC}"
                    status_icon="❌"
                    ;;
                *)
                    display_status="${UI_ORANGE}${status}${UI_NC}"
                    status_icon="📁"
                    ;;
            esac
            printf "   ├─ %s %-12s %b %b\n" "$status_icon" "$pkg:" "$UI_ARROW" "$display_status"
        fi
    done
    
    if [[ "$has_packages" == "false" ]]; then
        ui_print_color "${UI_GRAY}" "   └─ No packages configured"
    else
        echo "   └─ Configuration complete"
    fi
    echo ""
    
    # Installation summary
    local install_list=$(package_list_to_install 2>/dev/null || echo "")
    if [[ -n "$install_list" ]]; then
        ui_print_color "${UI_GREEN}${UI_BOLD}" "🚀 Packages to be built from source:"
        for pkg in $install_list; do
            ui_print_color "${UI_GREEN}" "   ${UI_ICON_BUILD} $pkg"
        done
        echo ""
    fi
    
    # System packages
    local system_packages=""
    for pkg in ${package_list}; do
        if [[ "$(config_get with_${pkg})" == "__SYSTEM__" ]]; then
            system_packages="$system_packages $pkg"
        fi
    done
    
    if [[ -n "$system_packages" ]]; then
        ui_print_color "${UI_BLUE}${UI_BOLD}" "🔗 System packages to be used:"
        for pkg in $system_packages; do
            ui_print_color "${UI_BLUE}" "   ${UI_ICON_CHECK} $pkg"
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

# Show installation progress for a stage with enhanced styling
# Usage: ui_stage_progress "stage_number" "stage_name"
ui_stage_progress() {
    local stage="$1"
    local name="$2"
    
    ui_section "Stage $stage: $name"
    ui_print_color "${UI_CYAN}${UI_BOLD}" "${UI_ICON_ROCKET} Installing packages for stage $stage..."
}

# Show package build progress with enhanced visual feedback
# Usage: ui_package_progress "package_name" "action"
ui_package_progress() {
    local package="$1"
    local action="$2"
    
    case "$action" in
        start)
            ui_print_color "${UI_BLUE}${UI_BOLD}" "${UI_ICON_BUILD} Building $package..." 
            ;;
        download)
            ui_debug "${UI_ICON_DOWNLOAD} Downloading $package..."
            ;;
        extract)
            ui_debug "${UI_ICON_PACKAGE} Extracting $package..."
            ;;
        configure)
            ui_debug "${UI_ICON_GEAR} Configuring $package..."
            ;;
        compile)
            ui_debug "${UI_ICON_BUILD} Compiling $package..."
            ;;
        install)
            ui_debug "${UI_ICON_INSTALL} Installing $package..."
            ;;
        success)
            ui_success "Successfully built $package"
            ;;
        skip)
            ui_print_color "${UI_GRAY}${UI_BOLD}" "${UI_ICON_INFO} Skipping $package (already built or disabled)"
            ;;
        error)
            ui_error "Failed to build $package"
            ;;
    esac
}

# Show enhanced final installation results
# Usage: ui_show_results "success_count" "total_count" "failed_packages"
ui_show_results() {
    local success_count="$1"
    local total_count="$2"
    local failed_packages="$3"
    
    ui_section "Installation Results"
    
    if [[ $success_count -eq $total_count ]]; then
        ui_success "All packages installed successfully! ($success_count/$total_count)"
        echo ""
        ui_print_color "${UI_GREEN}${UI_BOLD}" "${UI_ICON_ROCKET} Ready to use ABACUS toolchain!"
        echo ""
        ui_print_color "${UI_BLUE}${UI_BOLD}" "🔧 To activate the toolchain environment:"
        ui_print_color "${UI_WHITE}" "   source ${SETUPFILE:-setup}"
        ui_print_color "${UI_GRAY}" "   # or alternatively:"
        ui_print_color "${UI_WHITE}" "   source ${INSTALLDIR:-install}/abacus_env.sh"
        echo ""
        ui_print_color "${UI_BLUE}${UI_BOLD}" "🚀 Build ABACUS with:"
        ui_print_color "${UI_WHITE}" "   ./build_abacus_gnu.sh      # GNU toolchain"
        ui_print_color "${UI_WHITE}" "   ./build_abacus_intel.sh    # Intel toolchain"
        ui_print_color "${UI_WHITE}" "   ./build_abacus_gcc-aocl.sh # AMD GCC+AOCL"
        ui_print_color "${UI_WHITE}" "   ./build_abacus_aocc-aocl.sh # AMD AOCC+AOCL"
    else
        local failed_count=$((total_count - success_count))
        ui_error "Installation completed with errors ($success_count/$total_count successful)"
        echo ""
        if [[ -n "$failed_packages" ]]; then
            ui_print_color "${UI_RED}${UI_BOLD}" "❌ Failed packages:"
            for pkg in $failed_packages; do
                ui_print_color "${UI_RED}" "   ${UI_ICON_CROSS} $pkg"
            done
        fi
        echo ""
        ui_print_color "${UI_YELLOW}${UI_BOLD}" "🔍 Troubleshooting tips:"
        ui_print_color "${UI_YELLOW}" "   • Check log files for detailed error information"
        ui_print_color "${UI_YELLOW}" "   • Install missing system dependencies"
        ui_print_color "${UI_YELLOW}" "   • Verify network connectivity for downloads"
        ui_print_color "${UI_YELLOW}" "   • Try with --install-all for problematic packages"
    fi
}

# Show enhanced environment setup instructions
# Usage: ui_show_env_setup
ui_show_env_setup() {
    ui_section "Environment Setup Instructions"
    
    ui_print_color "${UI_BLUE}${UI_BOLD}" "🔧 To use the installed ABACUS toolchain:"
    echo ""
    ui_print_color "${UI_GREEN}${UI_BOLD}" "   source ${SETUPFILE:-setup}"
    ui_print_color "${UI_GRAY}" "   # or"
    ui_print_color "${UI_GREEN}${UI_BOLD}" "   source ${INSTALLDIR:-install}/abacus_env.sh"
    echo ""
    ui_print_color "${UI_BLUE}${UI_BOLD}" "🚀 Then build ABACUS with:"
    ui_print_color "${UI_GREEN}${UI_BOLD}" "   make -j\$(nproc)"
    echo ""
    ui_warning "The environment setup needs to be done in each new shell session"
}

# Handle user interruption with graceful cleanup
# Usage: ui_handle_interrupt
ui_handle_interrupt() {
    echo ""
    ui_warning "Installation interrupted by user (Ctrl+C)"
    ui_print_color "${UI_BLUE}${UI_BOLD}" "${UI_ICON_GEAR} Cleaning up temporary files..."
    
    # Clean up any temporary files or processes
    if [[ -n "$BUILDDIR" && -d "$BUILDDIR" ]]; then
        ui_debug "Cleaning build directory: $BUILDDIR"
        # Don't remove the entire build directory, just mark as interrupted
        touch "$BUILDDIR/.interrupted"
    fi
    
    ui_print_color "${UI_YELLOW}${UI_BOLD}" "${UI_ICON_INFO} Installation cancelled - you can resume later"
    exit 130
}

# Set up signal handlers
# Usage: ui_setup_signals
ui_setup_signals() {
    trap ui_handle_interrupt SIGINT SIGTERM
}

# Enhanced input validation with helpful error messages
# Usage: ui_validate_input "input" "type"
ui_validate_input() {
    local input="$1"
    local type="$2"
    
    case "$type" in
        number)
            if [[ ! "$input" =~ ^[0-9]+$ ]]; then
                ui_error "Invalid number: $input"
                ui_info "Please provide a positive integer"
                return 1
            fi
            ;;
        path)
            if [[ ! -d "$input" ]]; then
                ui_error "Directory does not exist: $input"
                ui_info "Please provide a valid directory path"
                return 1
            fi
            ;;
        file)
            if [[ ! -f "$input" ]]; then
                ui_error "File does not exist: $input"
                ui_info "Please provide a valid file path"
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
                    ui_print_color "${UI_BLUE}${UI_BOLD}" "📋 Valid options:"
                    ui_print_color "${UI_GREEN}" "   • mpich    - MPICH implementation (recommended)"
                    ui_print_color "${UI_GREEN}" "   • openmpi  - Open MPI implementation"
                    ui_print_color "${UI_GREEN}" "   • intelmpi - Intel MPI (requires Intel compiler)"
                    ui_print_color "${UI_GREEN}" "   • no       - Disable MPI support"
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
                    ui_print_color "${UI_BLUE}${UI_BOLD}" "📋 Valid options:"
                    ui_print_color "${UI_GREEN}" "   • openblas - OpenBLAS (open source, recommended)"
                    ui_print_color "${UI_GREEN}" "   • mkl      - Intel Math Kernel Library (high performance)"
                    ui_print_color "${UI_GREEN}" "   • aocl     - AMD Optimizing CPU Libraries"
                    ui_print_color "${UI_GREEN}" "   • cray     - Cray LibSci (for Cray systems)"
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
            
            # Invalid format - show enhanced error message
            ui_error "Invalid GPU version: $input"
            ui_print_color "${UI_BLUE}${UI_BOLD}" "📋 Valid formats:"
            ui_print_color "${UI_GREEN}" "   • Decimal format: 6.0, 7.0, 8.0, 8.9, etc. (CUDA compute capability)"
            ui_print_color "${UI_GREEN}" "   • Integer format: 60, 70, 80, 89, etc."
            ui_print_color "${UI_GREEN}" "   • Disable GPU: no"
            ui_print_color "${UI_ORANGE}${UI_BOLD}" "💡 Examples:"
            ui_print_color "${UI_WHITE}" "   • 8.0 or 80 for compute capability 8.0 (RTX 30xx series)"
            ui_print_color "${UI_WHITE}" "   • 7.5 or 75 for compute capability 7.5 (RTX 20xx series)"
            return 1
            ;;
        *)
            ui_error "Unknown validation type: $type"
            return 1
            ;;
    esac
    
    return 0
}

# Enhanced system information display
# Usage: ui_show_system_info
ui_show_system_info() {
    ui_section "System Information"
    
    ui_print_color "${UI_BLUE}${UI_BOLD}" "🖥️  Hardware & OS:"
    printf "   ├─ %-20s %s\n" "Operating System:" "$(uname -s)"
    printf "   ├─ %-20s %s\n" "Architecture:" "$(uname -m)"
    printf "   ├─ %-20s %s\n" "Kernel Version:" "$(uname -r)"
    printf "   ├─ %-20s %s\n" "CPU Cores:" "$(nproc 2>/dev/null || echo "unknown")"
    
    if command -v free &> /dev/null; then
        local mem_gb=$(free -g | awk '/^Mem:/ {print $2}')
        printf "   ├─ %-20s %sGB\n" "Total Memory:" "$mem_gb"
    fi
    
    printf "   ├─ %-20s %s\n" "Shell:" "$SHELL"
    printf "   ├─ %-20s %s\n" "User:" "$(whoami)"
    printf "   └─ %-20s %s\n" "Working Directory:" "$(pwd)"
    echo ""
}

# Enhanced system requirements check
# Usage: ui_check_system_requirements
ui_check_system_requirements() {
    ui_section "System Requirements Check"
    
    local missing_tools=""
    local required_tools="wget curl tar gzip make"
    local found_tools=""
    
    ui_print_color "${UI_BLUE}${UI_BOLD}" "🔍 Checking required system tools..."
    
    for tool in $required_tools; do
        if ! command -v "$tool" &> /dev/null; then
            missing_tools="$missing_tools $tool"
            ui_print_color "${UI_RED}" "   ${UI_ICON_CROSS} $tool - not found"
        else
            found_tools="$found_tools $tool"
            ui_print_color "${UI_GREEN}" "   ${UI_ICON_CHECK} $tool - found"
        fi
    done
    
    echo ""
    
    if [[ -n "$missing_tools" ]]; then
        ui_error "Missing required system tools:$missing_tools"
        echo ""
        ui_print_color "${UI_YELLOW}${UI_BOLD}" "📦 Install missing tools using your package manager:"
        ui_print_color "${UI_WHITE}" "   # Ubuntu/Debian:"
        ui_print_color "${UI_GREEN}" "   sudo apt-get install$missing_tools"
        ui_print_color "${UI_WHITE}" "   # CentOS/RHEL:"
        ui_print_color "${UI_GREEN}" "   sudo yum install$missing_tools"
        ui_print_color "${UI_WHITE}" "   # Fedora:"
        ui_print_color "${UI_GREEN}" "   sudo dnf install$missing_tools"
        return 1
    else
        ui_success "All required system tools are available"
        return 0
    fi
}