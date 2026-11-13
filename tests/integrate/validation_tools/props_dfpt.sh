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
# does not form a dense 0..N sequence starting at 0, or any freq block
# has a different total mode count than the first complete block (every
# q point has exactly 3*natom modes — a truncated tail at EOF that is
# internally dense would otherwise become an incomplete reference on
# first-time generation, where no old result.ref exists to cross-check).

run_dfpt_props() {
  if [ -z "$esolver_type" ] || [ "$esolver_type" != "dfpt" ]; then
    return
  fi
  local dfpt_log="${1:-log.txt}"
  if [ ! -f "$dfpt_log" ]; then
    echo "ERROR: run_dfpt_props: $dfpt_log not found, DFPT frequencies not collected" >&2
    return 1
  fi
  # count atoms from STRU to derive the expected 3*natom mode count;
  # this catches a single-block truncation that the cross-block consistency
  # check cannot detect (1-atom complete vs 2-atom truncated look identical)
  local natom=0
  if [ -f "STRU" ]; then
    natom=$(awk '/^ *Direct|^ *Cartesian/{found=1;next} found && /^[A-Z][a-z]/{getline; getline; natom+=$1} END{print natom}' STRU 2>/dev/null)
  fi
  local expected_total=$((natom * 3))
  # emit one "key value" line per frequency; awk rewrites the block header
  # into the key prefix and drops everything else. The mode index keeps its
  # zero-based printing so keys stay stable across cell sizes. Per-block
  # mode indices are checked for density (0..N with no gaps); every freq
  # block must have exactly expected_total modes (3*natom from STRU).
  if ! awk -v expected_total="$expected_total" '
    function fail(msg) { printf "ERROR: run_dfpt_props: %s\n", msg > "/dev/stderr"; bad = 1; exit 1 }
    function close_block() {
      if (seq != "" && nextmode[seq] > 0) {
        if (expected_total > 0 && block_type[seq] == "freq" && nextmode[seq] != expected_total) {
          fail("freq block " seq " has " nextmode[seq] " modes; expected " expected_total " (3*natom from STRU): truncated at EOF or mid-block")
        } else if (expected_modes == 0) {
          expected_modes = nextmode[seq]
          if (expected_modes % 3 != 0) {
            fail("first freq block has " expected_modes " modes; not a multiple of 3 (expected 3*natom)")
          }
        } else if (block_type[seq] == "freq" && nextmode[seq] != expected_modes) {
          fail("freq block " seq " has " nextmode[seq] " modes; expected " expected_modes " (truncated at EOF or mid-block)")
        }
      }
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
      if (mi != nextmode[seq] && nextmode[seq] > 0) {
        fail("mode " mi " out of sequence in " seq " (expected " nextmode[seq] "): truncated or dropped rows")
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
  # zero-key guard: a DFPT log with no analytic block at all
  if ! grep -q "^dfpt_freq_" "$props_result_file"; then
    echo "ERROR: run_dfpt_props: no complete phonon frequency block found in $dfpt_log (truncated output or format drift); failing the collection" >&2
    return 1
  fi
  # cross-reference guard: when generating result.out against an existing
  # result.ref (comparison mode), every dfpt_* key in the reference must
  # appear in the collected output — a log truncated mid-block would
  # otherwise pass on its surviving prefix
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
