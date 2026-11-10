#!/bin/bash

# Property collectors for the PW DFPT esolver (esolver_type dfpt).
#
# The frequency blocks are printed to stdout by DFPT_PW (captured in the
# test's log file): one "DFPT phonon frequencies at q #<i>" block per q
# point and, at Gamma, a "DFPT LO-TO corrected frequencies" block along
# the q->0 direction. Each block lists "<mode j> : <freq> cm^-1" lines.
# Both blocks are collected so result.ref pins the analytic (TO) and the
# direction-dependent (LO) spectrum separately.
#
# Keys emitted:
#   dfpt_freq_q<i>_m<j>      plain (analytic) frequency of q point i, mode j
#   dfpt_loto_q<i>_m<j>      LO-TO corrected frequency (Gamma only)

run_dfpt_props() {
  if [ -z "$esolver_type" ] || [ "$esolver_type" != "dfpt" ]; then
    return
  fi
  local dfpt_log="${1:-log.txt}"
  if [ ! -f "$dfpt_log" ]; then
    echo "WARNING: run_dfpt_props: $dfpt_log not found, DFPT frequencies not collected" >&2
    return
  fi
  # emit one "key value" line per frequency; awk rewrites the block header
  # into the key prefix and drops everything else. The mode index keeps its
  # zero-based printing so keys stay stable across cell sizes.
  awk '
    /^ *DFPT phonon frequencies at q #/ {
      match($0, /q #[0-9]+/); qi = substr($0, RSTART + 3, RLENGTH - 3);
      block = "freq"; next
    }
    /^ *DFPT LO-TO corrected frequencies at q #/ {
      match($0, /q #[0-9]+/); qi = substr($0, RSTART + 3, RLENGTH - 3);
      block = "loto"; next
    }
    /^ *mode +[0-9]+ :/ {
      if (block == "") next
      mi = $2
      val = $4
      printf "dfpt_%s_q%s_m%s %s\n", block, qi, mi, val
      next
    }
    { if ($0 !~ /^ *mode/ && block != "" && $0 ~ /[^ ]/) block = "" }
  ' "$dfpt_log" >> "$props_result_file"
}
