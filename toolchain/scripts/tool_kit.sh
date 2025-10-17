# ABACUS Toolchain Utility Functions
# A set of tools used in the toolchain installer, intended to be used
# by sourcing this file inside other scripts.
# Enhanced with modular architecture support
# Author: ABACUS Development Team
# Date: 2025-01-12

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all
# shellcheck shell=bash

SYS_INCLUDE_PATH=${SYS_INCLUDE_PATH:-"/usr/local/include:/usr/include"}
SYS_LIB_PATH=${SYS_LIB_PATH:-"/usr/local/lib64:/usr/local/lib:/usr/lib64:/usr/lib:/lib64:/lib"}
INCLUDE_PATHS=${INCLUDE_PATHS:-"CPATH SYS_INCLUDE_PATH"}
LIB_PATHS=${LIB_PATHS:-"LIBRARY_PATH LD_LIBRARY_PATH LD_RUN_PATH SYS_LIB_PATH"}
time_start=$(date +%s)

# ============================================================================
# Timing and Reporting Functions
# ============================================================================

# report timing
report_timing() {
  time_stop=$(date +%s)
  printf "Step %s took %0.2f seconds.\n" $1 $((time_stop - time_start))
}

# Enhanced timing function with better formatting
report_timing_enhanced() {
  local step_name="$1"
  local start_time="$2"
  local end_time="${3:-$(date +%s)}"
  local duration=$((end_time - start_time))
  
  if [[ $duration -gt 3600 ]]; then
    printf "Step '%s' completed in %d hours, %d minutes, %d seconds.\n" \
      "$step_name" $((duration / 3600)) $(((duration % 3600) / 60)) $((duration % 60))
  elif [[ $duration -gt 60 ]]; then
    printf "Step '%s' completed in %d minutes, %d seconds.\n" \
      "$step_name" $((duration / 60)) $((duration % 60))
  else
    printf "Step '%s' completed in %d seconds.\n" "$step_name" $duration
  fi
}

# ============================================================================
# Error and Warning Reporting Functions
# ============================================================================

# report a warning message with script name and line number
report_warning() {
  if [ $# -gt 1 ]; then
    local __lineno=", line $1"
    local __message="$2"
  else
    local __lineno=''
    local __message="$1"
  fi
  echo "WARNING: (${SCRIPT_NAME}${__lineno}) $__message" >&2
}

# report an error message with script name and line number
report_error() {
  if [ $# -gt 1 ]; then
    local __lineno=", line $1"
    local __message="$2"
  else
    local __lineno=''
    local __message="$1"
  fi
  echo "ERROR: (${SCRIPT_NAME}${__lineno}) $__message" >&2
}

# Enhanced error reporting with context
report_error_enhanced() {
  local error_code="$1"
  local error_message="$2"
  local context="${3:-}"
  local lineno="${4:-}"
  
  local location=""
  if [[ -n "$lineno" ]]; then
    location=", line $lineno"
  fi
  
  echo "ERROR [$error_code]: (${SCRIPT_NAME}${location}) $error_message" >&2
  if [[ -n "$context" ]]; then
    echo "Context: $context" >&2
  fi
  
  return $error_code
}

# Enhanced warning reporting with severity levels
report_warning_enhanced() {
  local severity="${1:-INFO}"  # INFO, WARNING, CRITICAL
  local message="$2"
  local context="${3:-}"
  local lineno="${4:-}"
  
  local location=""
  if [[ -n "$lineno" ]]; then
    location=", line $lineno"
  fi
  
  echo "$severity: (${SCRIPT_NAME}${location}) $message" >&2
  if [[ -n "$context" ]]; then
    echo "Context: $context" >&2
  fi
}

# recommend users to use offline installation when download failed
# zhaoqing in 2025.10.15
recommend_offline_installation(){
  __filename=$1
  __url=$2
  cat << EOF
========================== NOTICE =========================
You can use OFFLINE installation method manually
By download $__filename from $__url, 
Rename it as $__filename and put it into ${BUILDDIR},
And re-run toolchain installation script.

Also. the --pack-run option can help in OFFLINE installation.

You can manually install requirements packages via:
0. Download by 'wget $__url -O $__filename' manually
1. Download from www.cp2k.org/static/downloads (for OpenBLAS, OpenMPI and Others)
2. Download from github.com (especially for CEREAL, RapidJSON, libnpy, LibRI and other stage4 packages)
3. for Intel-oneAPI and AMD AOCC/AOCL, please contact your server manager or visit their official website
EOF
}

# error handler for line trap from set -e
error_handler() {
  local __lineno="$1"
  report_error $1 "Non-zero exit code detected."
  exit 1
}

# source a file if it exists, otherwise do nothing
load() {
  if [ -f "$1" ]; then
    source "$1"
  fi
}

# A more portable command that will give the full path, removing
# symlinks, of a given path. This is more portable than readlink -f
# which does not work on Mac OS X
realpath() {
  local __path="$1"
  if [ "x$__path" = x ]; then
    return 0
  fi
  local __basename=$(basename "$__path")
  if [ -e "$__path" ]; then
    echo $(
      cd "$(dirname "$__path")"
      pwd -P
    )/"$__basename"
    return 0
  else
    return 1
  fi
}

# given a list, outputs a list with duplicated items filtered out
unique() (
  # given a list, outputs a list with duplicated items filtered out.
  # If -d <delimiter> option exists, then output the list delimited
  # by <delimiter>; note that this option does not effect the input.
  local __result=''
  local __delimiter=' '
  local __item=''
  if [ "$1" = "-d" ]; then
    shift
    __delimiter="$1"
    shift
  fi
  # It is essential that we quote $@, which makes it equivalent to
  # "$1" "$2" ...  So this works if any of the arguments contains
  # space.  And we use \n to separate the fields in the
  # __result for now, so that fields that contain spaces are
  # correctly grepped.
  for __item in "$@"; do
    if [ x"$__result" = x ]; then
      __result="${__item}"
    # Note that quoting $__result after echo is essential to
    # retain the \n in the variable from the output of echo.  Also
    # remember grep only works on a line by line basis, so if
    # items are delimited by newlines, then for grep search it
    # should be delimited by ^ and $ (beginning and end of line)
    elif ! (echo "$__result" |
      grep -s -q -e "^$__item\$"); then
      __result="${__result}
${__item}"
    fi
  done
  __result="$(echo "$__result" | paste -s -d "$__delimiter" -)"
  # quoting $__result below is again essential for correct
  # behaviour if IFS is set to be the same $__delimiter in the
  # parent shell calling this macro
  echo "$__result"
)

# reverse a list
reverse() (
  # given a list, output a list with reversed order. If -d
  # <delimiter> option exists, then output the list delimited by
  # <delimiter>; note that this option does not effect the input.
  local __result=''
  local __delimiter=' '
  local __item=''
  if [ "$1" = "-d" ]; then
    shift
    __delimiter="$1"
    shift
  fi
  for __item in "$@"; do
    if [ x"$__result" = x ]; then
      __result="$__item"
    else
      __result="${__item}${__delimiter}${__result}"
    fi
  done
  echo "$__result"
)

# get the number of processes available for compilation
get_nprocs() {
  if [ -n "${NPROCS_OVERWRITE}" ]; then
    echo ${NPROCS_OVERWRITE} | sed 's/^0*//'
  elif $(command -v nproc > /dev/null 2>&1); then
    echo $(nproc --all)
  elif $(command -v sysctl > /dev/null 2>&1); then
    echo $(sysctl -n hw.ncpu)
  else
    echo 1
  fi
}

# convert a list of paths to -L<dir> ... used by ld
paths_to_ld() {
  # need to define the POSIX default IFS values here, cannot just do
  # __ifs=$IFS first, because IFS can be unset, and if so __ifs will
  # becomes an empty string (null) and NOT unset, so later when IFS
  # is set to __ifs it becomes null rather than unset, and thus
  # causing wrong behaviour.  So if IFS is unset, __ifs should be
  # the POSIX default value.  Further more, due to shell
  # automatically remove the tailing "\n" in a string during
  # variable assignment, we need to add x after \n and then remove
  # it.
  local __paths=$@
  local __name=''
  local __raw_path=''
  local __dir=''
  local __lib_dirs=''
  # set default IFS first
  local __ifs=$(printf " \t\nx")
  __ifs="${__ifs%x}"
  [ "$IFS" ] && __ifs="$IFS"
  for __name in $__paths; do
    eval __raw_path=\$"$__name"
    # change internal field separator to :
    IFS=':'
    # loop over all dirs in path, and filter out duplications
    for __dir in $__raw_path; do
      if ! [ x"$__dir" = x ]; then
        if ! [[ "$__lib_dirs" =~ (^|[[:space:]])"-L'$__dir'"($|[[:space:]]) ]]; then
          __lib_dirs="$__lib_dirs -L'$__dir'"
        fi
      fi
    done
    IFS="$__ifs"
  done
  echo $__lib_dirs
}

# Find a file from directories given in a list of paths, each has the
# same format as env variable PATH. If the file is found, then echoes
# the full path of the file. If the file is not found, then echoes
# __FALSE__. The file name can also contain wildcards that are
# acceptable for bash, and in that case the full path of the first
# matching file will be echoed.
find_in_paths() {
  local __target=$1
  shift
  local __paths=$@
  local __name=''
  local __raw_path=''
  local __dir=''
  local __file=''
  local __files=''
  # use the IFS variable to take care of possible spaces in file/dir names
  local __ifs="$(printf " \t\nx")"
  __ifs="${__ifs%x}"
  [ "$IFS" ] && __ifs="$IFS"
  for __name in $__paths; do
    eval __raw_path=\$"$__name"
    # fields in paths are separated by :
    IFS=':'
    for __dir in $__raw_path; do
      # files in possible glob expansion are to be delimited by "\n\b"
      IFS="$(printf "\nx")"
      IFS="${IFS%x}"
      for __file in $__dir/$__target; do
        if [ -e "$__file" ]; then
          echo $(realpath "$__file")
          # must remember to change IFS back when exiting
          IFS="$__ifs"
          return 0
        fi
      done
      IFS=':'
    done
    IFS=$__ifs
  done
  echo "__FALSE__"
}

# search through a list of given paths, try to find the required file
# or directory, and if found then add full path of dirname file, or
# directory, to the -I include list for CFLAGS and append to a user
# specified variable (__cflags_name). If not found, then nothing is
# done. If the option -p is present, then if the search target is a
# directory, then the parent directory of the directory is used for -I
# instead.  The search target accepts bash wildcards, and in this case
# the first match will be used.
add_include_from_paths() {
  local __parent_dir_only=false
  if [ $1 = "-p" ]; then
    __parent_dir_only=true
    shift
  fi
  local __cflags_name=$1
  shift
  local __search_target=$1
  shift
  local __paths=$@
  local __found_target=""
  local __cflags=""
  __found_target="$(find_in_paths "$__search_target" \
    $__paths)"
  if [ "$__found_target" != "__FALSE__" ]; then
    if [ -f "$__found_target" ] || $__parent_dir_only; then
      __found_target="$(dirname "$__found_target")"
    fi
    echo "Found include directory $__found_target"
    eval __cflags=\$"${__cflags_name}"
    __cflags="${__cflags} -I'${__found_target}'"
    # remove possible duplicates
    __cflags="$(unique $__cflags)"
    # must escape all quotes again before the last eval, as
    # otherwise all quotes gets interpreted by the shell when
    # assigning to variable because eval will reduce one escape
    # level
    __cflags="${__cflags//'/\\'/}"
    eval $__cflags_name=\"$__cflags\"
  fi
}

# search through a list of given paths, try to find the required file
# or directory, and if found then add full path of dirname file, or
# directory, to the -L library list (including -Wl,-rpath) for LDFLAGS
# and append to a user specified variable (__ldflags_name). If not
# found, then nothing is done. If the option -p is present, then if
# the search target is a directory, then the parent directory of the
# directory is used for -L instead.  The search target accepts bash
# wildcards, and in this case the first match will be used.
add_lib_from_paths() {
  local __parent_dir_only=false
  if [ $1 = "-p" ]; then
    __parent_dir_only=true
    shift
  fi
  local __ldflags_name=$1
  shift
  local __search_target=$1
  shift
  local __paths=$@
  local __found_target=""
  local __ldflags=""
  __found_target="$(find_in_paths "$__search_target" \
    $__paths)"
  if [ "$__found_target" != "__FALSE__" ]; then
    if [ -f "$__found_target" ] || $__parent_dir_only; then
      __found_target="$(dirname "$__found_target")"
    fi
    echo "Found lib directory $__found_target"
    eval __ldflags=\$"${__ldflags_name}"
    __ldflags="${__ldflags} -L'${__found_target}' -Wl,-rpath,'${__found_target}'"
    # remove possible duplicates
    __ldflags="$(unique $__ldflags)"
    # must escape all quotes again before the last eval, as
    # otherwise all quotes gets interpreted by the shell when
    # assigning to variable because eval will reduce one escape
    # level
    __ldflags="${__ldflags//'/\\'/}"
    eval $__ldflags_name=\"$__ldflags\"
  fi
}

# ============================================================================
# Environment and System Validation Functions
# ============================================================================

# check if environment variable is assigned and non-empty
# https://serverfault.com/questions/7503/how-to-determine-if-a-bash-variable-is-empty
require_env() {
  local __env_var_name=$1
  local __env_var="$(eval echo \"\$$__env_var_name\")"
  if [ -z "${__env_var+set}" ]; then
    report_error "requires environment variable $__env_var_name to work"
    return 1
  fi
}

# Enhanced environment variable validation
validate_env_vars() {
  local required_vars=("$@")
  local missing_vars=()
  
  for var in "${required_vars[@]}"; do
    if [[ -z "${!var:-}" ]]; then
      missing_vars+=("$var")
    fi
  done
  
  if [[ ${#missing_vars[@]} -gt 0 ]]; then
    report_error_enhanced 1 "Missing required environment variables" "$(printf '%s ' "${missing_vars[@]}")"
    return 1
  fi
  
  return 0
}

# Check system requirements
check_system_requirements() {
  local min_disk_space_gb="${1:-10}"  # Default 10GB
  local required_commands=("${@:2}")
  
  # Check disk space
  local available_space=$(df . | awk 'NR==2 {print int($4/1024/1024)}')
  if [[ $available_space -lt $min_disk_space_gb ]]; then
    report_warning_enhanced "WARNING" "Low disk space: ${available_space}GB available, ${min_disk_space_gb}GB recommended"
  fi
  
  # Check required commands
  for cmd in "${required_commands[@]}"; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
      report_error_enhanced 1 "Required command not found: $cmd"
      return 1
    fi
  done
  
  return 0
}

# ============================================================================
# String and Configuration Processing Functions
# ============================================================================

resolve_string() {
  local __to_resolve=$1
  shift
  local __flags=$@

  echo $("${SCRIPTDIR}/parse_if.py" $__flags <<< "${__to_resolve}")
}

# Enhanced string processing with validation
resolve_string_safe() {
  local __to_resolve="$1"
  shift
  local __flags=("$@")
  
  if [[ -z "$__to_resolve" ]]; then
    report_error_enhanced 1 "Empty string provided for resolution"
    return 1
  fi
  
  if [[ ! -f "${SCRIPTDIR}/parse_if.py" ]]; then
    report_error_enhanced 1 "Parser script not found" "${SCRIPTDIR}/parse_if.py"
    return 1
  fi
  
  local result
  result=$("${SCRIPTDIR}/parse_if.py" "${__flags[@]}" <<< "${__to_resolve}" 2>/dev/null)
  local exit_code=$?
  
  if [[ $exit_code -ne 0 ]]; then
    report_error_enhanced $exit_code "String resolution failed" "$__to_resolve"
    return $exit_code
  fi
  
  echo "$result"
  return 0
}

# ============================================================================
# Command and Package Validation Functions
# ============================================================================

# check if a command is available
check_command() {
  local __command=${1}
  if [ $# -eq 1 ]; then
    local __package=${1}
  elif [ $# -gt 1 ]; then
    local __package=${2}
  fi
  if $(command -v ${__command} > /dev/null 2>&1); then
    echo "path to ${__command} is $(realpath $(command -v ${__command}))"
  else
    report_error "Cannot find ${__command}, please check if the package ${__package} is installed or in system search path"
    return 1
  fi
}

# Enhanced command checking with version validation
check_command_version() {
  local command_name="$1"
  local min_version="$2"
  local package_name="${3:-$command_name}"
  
  if ! command -v "$command_name" >/dev/null 2>&1; then
    report_error_enhanced 1 "Command not found: $command_name" "Package: $package_name"
    return 1
  fi
  
  local command_path
  command_path=$(command -v "$command_name")
  echo "Found $command_name at: $(realpath "$command_path")"
  
  if [[ -n "$min_version" ]]; then
    # This is a placeholder for version checking logic
    # Different commands have different version output formats
    echo "Version check for $command_name (minimum: $min_version) - implement as needed"
  fi
  
  return 0
}

# Batch command checking
check_commands_batch() {
  local commands=("$@")
  local failed_commands=()
  
  for cmd in "${commands[@]}"; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
      failed_commands+=("$cmd")
    else
      echo "✓ Found: $cmd at $(command -v "$cmd")"
    fi
  done
  
  if [[ ${#failed_commands[@]} -gt 0 ]]; then
    report_error_enhanced 1 "Missing commands" "$(printf '%s ' "${failed_commands[@]}")"
    return 1
  fi
  
  return 0
}

# check if directory exists
# add more error msg by QuantumMisaka in 2025.03.19
check_dir() {
  local __dir=$1
  if [ -d "$__dir" ]; then
    echo "Found directory $__dir"
  else
    report_error "Cannot find $__dir, please check your --with-PKG input to march options: [system|install|no|(path/to/pkg)]"
    return 1
  fi
}

# check if a command has been installed correctly
check_install() {
  local __command=${1}
  if [ $# -eq 1 ]; then
    local __package=${1}
  elif [ $# -gt 1 ]; then
    local __package=${2}
  fi
  if $(command -v ${__command} > /dev/null 2>&1); then
    echo "$(basename ${__command}) is installed as $(command -v ${__command})"
  else
    report_error "cannot find ${__command}, please check if the package ${__package} has been installed correctly"
    return 1
  fi
}

# check if a library can be found by ld, library names should in the
# format -lname, which would then referred to libname.a or libname.so
# by ld
check_lib() {
  local __libname="${1#-l}"
  if [ $# -eq 1 ]; then
    local __package=lib"$__libname"
  elif [ $# -gt 1 ]; then
    local __package=$2
  fi
  # Note that LD_LIBRARY_PATH is NOT used by ld during linking
  # stage, and is only used for searching to the shared libraries
  # required by the executable AFTER it has already been compiled, to
  # override its internal search paths built into the binary when it
  # was compiled. Here, we explicitly include the commonly defined
  # library search paths---including LD_LIBRARY_PATH---in the -L
  # search paths of ld.  This is the only way ld can include
  # non-standard directories in its search path.  If we use gcc
  # instead of ld for linker then we can use LIBRARY_PATH, which IS
  # used during link stage. However, I think using ld is more
  # general, as in most systems LIBRARY_PATH is rarely defined, and
  # we would have to rely on gcc.
  local __search_engine="ld -o /dev/null"
  local __search_paths="$LIB_PATHS"
  # convert a list of paths to -L<dir> list used by ld
  __search_engine="$__search_engine $(paths_to_ld $__search_paths)"
  # needed the eval to interpret the quoted directories correctly (somehow)
  if (eval $__search_engine -l$__libname 2>&1 | grep -q -s "\-l$__libname"); then
    # if library not found, ld will return error message
    # containing the library name
    report_error \
      "ld cannot find -l$__libname, please check if $__package is installed or in system search path"
    return 1
  else
    # if library is found, then ld will return error message about
    # not able to find _start or _main symbol
    echo "lib$__libname is found in ld search path"
  fi
}

# check if a module is available for the current version of gfortran,
# returns 0 if available and 1 if not
check_gfortran_module() {
  local __module_name=$1
  local __FC=${FC:-gfortran}
  cat << EOF | $__FC -c -o /dev/null -xf95 -ffree-form - > /dev/null 2>&1
PROGRAM check_gfortran_module
USE ${__module_name}
IMPLICIT NONE
PRINT *, "PASS"
END PROGRAM check_gfortran_module
EOF
}

# check if a flag is allowed for the current version of
# gfortran. returns 0 if allowed and 1 if not
check_gfortran_flag() {
  local __flag=$1
  local __FC=${FC:-gfortran}
  # no need to do a full compilation, just -E -cpp would do for
  # checking flags
  cat << EOF | $__FC -E -cpp $__flag -xf95 -ffree-form - > /dev/null 2>&1
PROGRAM test_code
  IMPLICIT NONE
  PRINT *, "PASS"
END PROGRAM test_code
EOF
}

# check if a flag is allowed for the current version of
# gcc. returns 0 if allowed and 1 if not
check_gcc_flag() {
  local __flag=$1
  local __CC=${CC:-gcc}
  # no need to do a full compilation, just -E -cpp would do for
  # checking flags
  cat << EOF | $__CC -E -cpp $__flag -xc - > /dev/null 2>&1
#include <stdio.h>
int main() {
  printf("PASS\n");
}
EOF
}

# check if a flag is allowed for the current version of
# g++. returns 0 if allowed and 1 if not
check_gxx_flag() {
  local __flag=$1
  local __CXX=${CXX:-g++}
  # no need to do a full compilation, just -E -cpp would do for
  # checking flags
  cat << EOF | $__CXX -E -cpp $__flag -xc - > /dev/null 2>&1
#include <stdio.h>
int main() {
  printf("PASS\n");
}
EOF
}

# given a list of flags, only print out what is allowed by the current
# version of gfortran
allowed_gfortran_flags() {
  local __flags=$@
  local __flag=''
  local __result=''
  for __flag in $__flags; do
    if (check_gfortran_flag $__flag); then
      [ -z "$__result" ] && __result="$__flag" || __result="$__result $__flag"
    fi
  done
  echo $__result
}

# given a list of flags, only print out what is allowed by the current
# version of gcc
allowed_gcc_flags() {
  local __flags=$@
  local __flag=''
  local __result=''
  for __flag in $__flags; do
    if (check_gcc_flag $__flag); then
      [ -z "$__result" ] && __result="$__flag" || __result="$__result $__flag"
    fi
  done
  echo $__result
}

# given a list of flags, only print out what is allowed by the current
# version of g++
allowed_gxx_flags() {
  local __flags=$@
  local __flag=''
  local __result=''
  for __flag in $__flags; do
    if (check_gxx_flag $__flag); then
      [ -z "$__result" ] && __result="$__flag" || __result="$__result $__flag"
    fi
  done
  echo $__result
}

# remove a directory to a given path
remove_path() {
  local __path_name=$1
  local __directory=$2
  local __path="$(eval echo \$$__path_name)"
  # must remove all the middle ones first before treating two ends,
  # otherwise there can be cases where not all __directory are
  # removed.
  __path=${__path//:$__directory:/:}
  __path=${__path#$__directory:}
  __path=${__path%:$__directory}
  __path=$(echo "$__path" | sed "s:^$__directory\$::g")
  eval $__path_name=\"$__path\"
  export $__path_name
}

# prepend a directory to a given path
prepend_path() {
  # prepend directory to $path_name and then export path_name. If
  # the directory already exists in path, bring the directory to the
  # front of the list.
  # $1 is path name
  # $2 is directory
  remove_path "$1" "$2"
  eval $1=\"$2\${$1:+\":\$$1\"}\"
  eval export $1
}

# append a directory to a given path
append_path() {
  # append directory to $path_name and then export path_name. If
  # the directory already exists in path, bring the directory to the
  # back of the list.
  # $1 is path name
  # $2 is directory
  remove_path "$1" "$2"
  eval $1=\"\${$1:+\"\$$1:\"}$2\"
  eval export $1
}

# ============================================================================
# Configuration Parsing Functions
# ============================================================================

# helper routine for reading --enable=* input options
read_enable() {
  local __input_var="${1#*=}"
  case $__input_var in
    "$1")
      # if there is no "=" then treat as "yes"
      echo "__TRUE__"
      ;;
    yes)
      echo "__TRUE__"
      ;;
    no)
      echo "__FALSE__"
      ;;
    *)
      echo "__INVALID__"
      ;;
  esac
}

# helper routine for reading --with=* input options
read_with() {
  local __input_var="${1#--with*=}"
  case $__input_var in
    "${1}")
      # if there is no "=" then treat as "install"
      if [ ${#} -gt 1 ]; then
        echo "${2}"
      else
        echo "__INSTALL__"
      fi
      ;;
    install)
      echo "__INSTALL__"
      ;;
    system)
      echo "__SYSTEM__"
      ;;
    no)
      echo "__DONTUSE__"
      ;;
    *)
      echo "${__input_var//\~/$HOME}"
      ;;
  esac
}

# Enhanced configuration parsing with validation
parse_config_option() {
  local option="$1"
  local default_value="${2:-}"
  local valid_values="${3:-}"  # Space-separated list of valid values
  
  local parsed_value
  if [[ "$option" =~ ^--enable ]]; then
    parsed_value=$(read_enable "$option")
  elif [[ "$option" =~ ^--with ]]; then
    parsed_value=$(read_with "$option" "$default_value")
  else
    parsed_value="$option"
  fi
  
  # Validate against allowed values if provided
  if [[ -n "$valid_values" && "$parsed_value" != "__INVALID__" ]]; then
    local is_valid=false
    for valid in $valid_values; do
      if [[ "$parsed_value" == "$valid" ]]; then
        is_valid=true
        break
      fi
    done
    
    if [[ "$is_valid" == false ]]; then
      report_warning_enhanced "WARNING" "Invalid configuration value: $parsed_value" "Valid options: $valid_values"
      parsed_value="__INVALID__"
    fi
  fi
  
  echo "$parsed_value"
}

# ============================================================================
# File Download and Integrity Functions
# ============================================================================

# helper routine to check integrity of downloaded files
checksum() {
  local __filename=$1
  local __sha256=$2
  local __shasum_command='sha256sum'
  # check if we have sha256sum command, Mac OS X does not have
  # sha256sum, but has an equivalent with shasum -a 256
  command -v "$__shasum_command" > /dev/null 2>&1 ||
    __shasum_command="shasum -a 256"
  if echo "$__sha256  $__filename" | ${__shasum_command} --check; then
    echo "Checksum of $__filename Ok"
  else
    rm -v ${__filename}
    report_error "Checksum of $__filename could not be verified, abort."
    return 1
  fi
}

# Enhanced checksum verification with multiple hash algorithms
verify_file_integrity() {
  local filename="$1"
  local expected_hash="$2"
  local hash_type="${3:-sha256}"  # Default to SHA256
  
  if [[ ! -f "$filename" ]]; then
    report_error_enhanced 1 "File not found for integrity check" "$filename"
    return 1
  fi
  
  local hash_command
  case "$hash_type" in
    sha256)
      hash_command="sha256sum"
      command -v "$hash_command" >/dev/null 2>&1 || hash_command="shasum -a 256"
      ;;
    sha1)
      hash_command="sha1sum"
      command -v "$hash_command" >/dev/null 2>&1 || hash_command="shasum -a 1"
      ;;
    md5)
      hash_command="md5sum"
      command -v "$hash_command" >/dev/null 2>&1 || hash_command="md5"
      ;;
    *)
      report_error_enhanced 1 "Unsupported hash type" "$hash_type"
      return 1
      ;;
  esac
  
  local computed_hash
  computed_hash=$($hash_command "$filename" | cut -d' ' -f1)
  
  if [[ "$computed_hash" == "$expected_hash" ]]; then
    echo "✓ Integrity check passed for $filename ($hash_type)"
    return 0
  else
    report_error_enhanced 1 "Integrity check failed for $filename" "Expected: $expected_hash, Got: $computed_hash"
    return 1
  fi
}

# downloader for the package tars, includes checksum
# backup and deprecated
download_pkg_from_org() {
  # usage: download_pkg_from_org sha256 filename
  echo "use cp2k mirror to download $__filename"
  local __sha256="$1"
  local __filename="$2"
  local __url="https://www.cp2k.org/static/downloads/$__filename"
  # download
  #echo "wget ${DOWNLOADER_FLAGS} --quiet $__url"
  #if ! wget ${DOWNLOADER_FLAGS} --quiet $__url; then
  echo "wget ${DOWNLOADER_FLAGS} $__url"
  if ! wget ${DOWNLOADER_FLAGS} $__url; then
    report_error "failed to download $__url"
    recommend_offline_installation $__filename $__url
    if [ "${PACK_RUN}" != "__TRUE__" ]; then
        return 1
    fi
  fi
  # checksum
  checksum "$__filename" "$__sha256"
}

download_pkg_from_url() {
  # usage: download_pkg_from_url sha256 filename url
  local __sha256="$1" # if set to "--no-checksum", do not check checksum
  local __filename="$2"
  local __url="$3"
  
  # Smart certificate validation strategy
  case "${DOWNLOAD_CERT_POLICY:-smart}" in
    "strict")
      echo "Downloading with strict certificate validation: $__url"
      if ! wget ${DOWNLOADER_FLAGS} "$__url" -O "$__filename"; then
        report_error "failed to download $__url (strict certificate validation)"
        recommend_offline_installation "$__filename" "$__url"
        if [ "${PACK_RUN}" != "__TRUE__" ]; then
          return 1
        fi
      fi
      ;;
    "skip")
      echo "Downloading with certificate validation disabled: $__url"
      if ! wget ${DOWNLOADER_FLAGS} "$__url" -O "$__filename" --no-check-certificate; then
        report_error "failed to download $__url"
        recommend_offline_installation "$__filename" "$__url"
        if [ "${PACK_RUN}" != "__TRUE__" ]; then
          return 1
        fi
      fi
      ;;
    "smart"|*)
      # Smart fallback: try with certificate validation first, then without
      echo "Attempting secure download: $__url"
      if wget ${DOWNLOADER_FLAGS} "$__url" -O "$__filename" 2>/dev/null; then
        echo "Download successful with certificate validation"
      else
        echo "Certificate validation failed, retrying without certificate check..."
        if ! wget ${DOWNLOADER_FLAGS} "$__url" -O "$__filename" --no-check-certificate; then
          report_error "failed to download $__url (both secure and insecure attempts failed)"
          recommend_offline_installation "$__filename" "$__url"
          if [ "${PACK_RUN}" != "__TRUE__" ]; then
            return 1
          fi
        else
          echo "Download successful without certificate validation"
        fi
      fi
      ;;
  esac
  
  # checksum validation (unchanged)
  if [ "$__sha256" != "--no-checksum" ]; then
    checksum "$__filename" "$__sha256"
  fi
}

# verify the checksums inside the given checksum file
verify_checksums() {
  local __checksum_file=$1
  local __shasum_command='sha256sum'

  # check if we have sha256sum command, Mac OS X does not have
  # sha256sum, but has an equivalent with shasum -a 256
  command -v "$__shasum_command" > /dev/null 2>&1 ||
    __shasum_command="shasum -a 256"

  ${__shasum_command} --check "${__checksum_file}" > /dev/null 2>&1
}

# write a checksum file $1 containing checksums for each given file $2, $3, ... (plus the $VERSION_FILE)
write_checksums() {
  local __checksum_file=$1
  shift # remove output file from arguments to be able to pass them along properly quoted
  local __shasum_command='sha256sum'

  # check if we have sha256sum command, Mac OS X does not have
  # sha256sum, but has an equivalent with shasum -a 256
  command -v "$__shasum_command" > /dev/null 2>&1 ||
    __shasum_command="shasum -a 256"

  ${__shasum_command} "${VERSION_FILE}" "$@" > "${__checksum_file}"
}

# generate a filtered toolchain.env
write_toolchain_env() {
  local __installdir=$1

  # run the following in a subshell to not affect the currently running shell
  # we do not need to achieve complete filtering, it is sufficient to
  # remove problematic variables (TERM/TERMCAP/COLORTERM) which may trigger
  # 'too many arguments' (since the environment vars are stored in the same memory block as command line arguments)
  # or which may not be valid anymore the next time the user runs the toolchain scripts,
  # like the proxy vars which may affect fetching tarballs
  (
    unset COLORTERM DISPLAY EDITOR LESS LESSOPEN LOGNAME LS_COLORS PAGER
    unset TERM TERMCAP USER
    unset ftp_proxy http_proxy no_proxy
    unset GPG_AGENT_INFO SSH_AGENT_PID SSH_AUTH_SOCK SSH_CLIENT SSH_CONNECTION SSH_TTY
    unset LS_COLORS LS_OPTIONS
    unset STY WINDOW XAUTHORITY
    unset XDG_CURRENT_DESKTOP XDG_RUNTIME_DIR XDG_SEAT XDG_SESSION_CLASS XDG_SESSION_DESKTOP XDG_SESSION_ID XDG_SESSION_TYPE XDG_VTNR XDG_CONFIG_DIRS XDG_DATA_DIRS
    unset DBUS_SESSION_BUS_ADDRESS

    export -p
  ) > "${__installdir}/toolchain.env"
}
