#!/bin/bash

# ============================================================================
# ABACUS Toolchain Error Handler
# ============================================================================
# Provides error handling and reporting functions for the toolchain
# Author: Quantum Misaka by Trae SOLO
# Date: 2025-10-16
# ============================================================================

# Global error handling state
ERROR_HANDLER_INITIALIZED=false

# Initialize error handler
# Usage: error_handler_init
error_handler_init() {
    if [[ "$ERROR_HANDLER_INITIALIZED" == "true" ]]; then
        return 0
    fi
    
    # Set up error trapping
    set -e
    trap 'error_handler ${LINENO}' ERR
    
    ERROR_HANDLER_INITIALIZED=true
    return 0
}

# Main error handler function
# Usage: error_handler line_number [error_code]
error_handler() {
    local line_number="$1"
    local error_code="${2:-$?}"
    local script_name="${SCRIPT_NAME:-${BASH_SOURCE[1]}}"
    
    echo "ERROR: Script failed at line $line_number in $script_name (exit code: $error_code)" >&2
    echo "Command that failed: ${BASH_COMMAND}" >&2
    
    # Print stack trace
    echo "Call stack:" >&2
    local frame=0
    while caller $frame; do
        ((frame++))
    done | while read line func file; do
        echo "  at $func ($file:$line)" >&2
    done
    
    exit $error_code
}

# Report error with context
# Usage: report_error line_number "error message" [context]
report_error() {
    local line_number=""
    local message=""
    local context=""
    
    if [[ $# -gt 1 ]]; then
        line_number="$1"
        message="$2"
        context="${3:-}"
    else
        message="$1"
    fi
    
    local location=""
    if [[ -n "$line_number" ]]; then
        location=", line $line_number"
    fi
    
    echo "ERROR: (${SCRIPT_NAME:-$0}${location}) $message" >&2
    if [[ -n "$context" ]]; then
        echo "Context: $context" >&2
    fi
}

# Report warning with context
# Usage: report_warning line_number "warning message" [context]
report_warning() {
    local line_number=""
    local message=""
    local context=""
    
    if [[ $# -gt 1 ]]; then
        line_number="$1"
        message="$2"
        context="${3:-}"
    else
        message="$1"
    fi
    
    local location=""
    if [[ -n "$line_number" ]]; then
        location=", line $line_number"
    fi
    
    echo "WARNING: (${SCRIPT_NAME:-$0}${location}) $message" >&2
    if [[ -n "$context" ]]; then
        echo "Context: $context" >&2
    fi
}

# Enhanced error reporting with error codes
# Usage: report_error_enhanced error_code "error message" [context] [line_number]
report_error_enhanced() {
    local error_code="$1"
    local error_message="$2"
    local context="${3:-}"
    local line_number="${4:-}"
    
    local location=""
    if [[ -n "$line_number" ]]; then
        location=", line $line_number"
    fi
    
    echo "ERROR [$error_code]: (${SCRIPT_NAME:-$0}${location}) $error_message" >&2
    if [[ -n "$context" ]]; then
        echo "Context: $context" >&2
    fi
    
    return $error_code
}

# Enhanced warning reporting with severity levels
# Usage: report_warning_enhanced severity "warning message" [context]
report_warning_enhanced() {
    local severity="$1"
    local warning_message="$2"
    local context="${3:-}"
    
    echo "$severity: (${SCRIPT_NAME:-$0}) $warning_message" >&2
    if [[ -n "$context" ]]; then
        echo "Context: $context" >&2
    fi
}