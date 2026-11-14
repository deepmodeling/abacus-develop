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
#
# Completeness guard: the collection fails (returns 1) when the log is
# missing, contains no analytic frequency block, any block's mode set
# does not start at mode 0 and form a dense 0..N sequence, any freq
# block's total does not equal 3*natom (read from STRU), or a freq
# block header is followed by no mode rows at all. These checks apply
# on first-time reference generation as well as comparison runs.

run_dfpt_props() {
  if [ -z "$esolver_type" ] || [ "$esolver_type" != "dfpt" ]; then
    return
  fi
  local dfpt_log="${1:-log.txt}"
  if [ ! -f "$dfpt_log" ]; then
    echo "ERROR: run_dfpt_props: $dfpt_log not found, DFPT frequencies not collected" >&2
    return 1
  fi
  # count atoms from STRU: walk each species block in ATOMIC_POSITIONS
  # (element name → magnetic moment → atom count → coordinate lines) and
  # sum the counts; handles multi-species cells and single-letter elements
  local natom=0
  if [ -f "STRU" ]; then
    natom=$(awk '
      /ATOMIC_POSITIONS/{found=1; next}
      !found{next}
      /^[ \t]*$/{next}
      found && /^(Direct|Cartesian|direct|cartesian)/{next}
      found && /^[A-Za-z][A-Za-z0-9]*[ \t]*$/{
        getline
        getline
        natom += $1 + 0
        for(i=0; i<($1+0); i++) getline
      }
      END{print natom+0}
    ' STRU 2>/dev/null)
  fi
  local expected_total=$((natom * 3))
  # emit one "key value" line per frequency; awk rewrites the block header
  # into the key prefix and drops everything else. Mode 0 must open every
  # block; indices are checked for dense 0..N; every freq block must end
  # at exactly expected_total (3*natom from STRU); empty freq blocks fail.
  if ! awk -v expected_total="$expected_total" '
    function fail(msg) { printf "ERROR: run_dfpt_props: %s\n", msg > "/dev/stderr"; bad = 1; exit 1 }
    function close_block() {
      if (seq == "") return  # no previous block to close
      if (nextmode[seq] == 0) {
        fail(block_type[seq] " block " seq " has no mode rows")
      }
      if (expected_total > 0 && nextmode[seq] != expected_total) {
        fail(block_type[seq] " block " seq " has " nextmode[seq] " modes; expected " expected_total " (3*natom from STRU)")
      }
      if (expected_total == 0 && expected_modes > 0 && nextmode[seq] != expected_modes) {
        fail(block_type[seq] " block " seq " has " nextmode[seq] " modes; expected " expected_modes " (from first block)")
      }
      if (expected_modes == 0) expected_modes = nextmode[seq]
    }
    /^ *DFPT phonon frequencies at q #/ {
      close_block()
      match($0, /q #[0-9]+/); qi = substr($0, RSTART + 3, RLENGTH - 3);
      block = "freq"; block_qi = qi; seq = "freq:" qi; nextmode[seq] = 0; block_type[seq] = "freq"; next
    }
    /^ *DFPT LO-TO corrected frequencies at q #/ {
      close_block()
      match($0, /q #[0-9]+/); qi = substr($0, RSTART + 3, RLENGTH - 3);
      block = "loto"; block_qi = qi; seq = "loto:" qi; nextmode[seq] = 0; block_type[seq] = "loto"; next
    }
    /^ *mode +[0-9]+ :/ {
      if (block == "") next
      mi = $2
      if (mi != nextmode[seq]) {
        fail("mode " mi " out of sequence in " seq " (expected " nextmode[seq] "): missing or dropped rows")
      }
      nextmode[seq] = mi + 1
      val = $4
      # reject empty or nonnumeric frequency tokens (format drift)
      if (val !~ /^-?[0-9]+\.?[0-9]*([eE][-+]?[0-9]+)?$/) {
        fail("nonnumeric frequency token [" val "] in " seq " mode " mi)
      }
      printf "dfpt_%s_q%s_m%s %s\n", block, block_qi, mi, val
      next
    }
    { if ($0 !~ /^ *mode/ && block != "" && $0 ~ /[^ ]/) block = "" }
    END {
      if (!bad) {
        close_block()
        if (NR == 0) fail("empty log")
        if (expected_modes == 0) fail("no complete frequency block found")
      }
    }
  ' "$dfpt_log" >> "$props_result_file"; then
    echo "ERROR: run_dfpt_props: incomplete phonon frequency spectrum in $dfpt_log (truncated output or format drift); failing the collection" >&2
    return 1
  fi
  # zero-key guard
  if ! grep -q "^dfpt_freq_" "$props_result_file"; then
    echo "ERROR: run_dfpt_props: no complete phonon frequency block found in $dfpt_log (truncated output or format drift); failing the collection" >&2
    return 1
  fi
  # cross-reference guard (comparison mode only; first-time generation
  # relies on the STRU-derived 3*natom check above)
  local ref_file="result.ref"
  if [ -f "$ref_file" ] && grep -q "^dfpt_freq_" "$ref_file"; then
    local missing
    missing=$(comm -23 <(grep -o "^dfpt_[a-z]*_q[0-9]*_m[0-9]*" "$ref_file" | sort -u) \
                      <(grep -o "^dfpt_[a-z]*_q[0-9]*_m[0-9]*" "$props_result_file" | sort -u))
    if [ -n "$missing" ]; then
      echo "ERROR: run_dfpt_props: collected spectrum is incomplete; missing keys (present in result.ref): $missing" >&2
      return 1
    fi
  fi
}
