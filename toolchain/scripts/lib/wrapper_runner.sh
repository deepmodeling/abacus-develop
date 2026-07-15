#!/usr/bin/env bash

run_toolchain_with_log() {
  local log_file="$1"
  shift

  "$@" | tee "$log_file"
  local installer_status=${PIPESTATUS[0]}

  return "$installer_status"
}
