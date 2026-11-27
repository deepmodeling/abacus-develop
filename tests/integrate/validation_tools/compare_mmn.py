#!/usr/bin/env python3
"""Gauge-invariant comparison of Wannier90 .mmn files.

M_mn(k,b) = <u_mk|u_{n,k+b}> is covariant under the band gauge:
M(k,b) -> U_k^dagger M(k,b) U_{k+b}.  The basis of the selected bands is
solver-dependent: the CPU (LAPACK) and GPU (cuSOLVER) diagonalizations of
the same subspace differ by a per-band phase at each k.  For this diamond
case the gauge recovered from the .amn files is diagonal to ~1e-7 (median
~1e-10), and the elementwise magnitudes of A agree to ~1e-7.  Under such a
diagonal gauge |M_mn| is invariant, so the magnitudes are compared here.
The raw complex elements cannot be compared because their arbitrary phases
are solver-dependent.  The diagonal-phase picture assumes the selected
bands are non-degenerate; a mixed degenerate subspace would require the
full gauge-invariant content (loop products of M around closed paths).

File layout (as written by toW90_PW::cal_Mmn): an optional preamble such
as " Created on ...", the header 'nbnd nk nntot', then nk*nntot blocks.
Each block starts with 'ik ikb G1 G2 G3' (ikb is the 1-based neighbour
k-point index from the .nnkp file), followed by nbnd*nbnd lines 'Re Im'
with m varying fastest.

Usage:
    compare_mmn.py ref.mmn cal.mmn [accuracy]

accuracy is the number of decimal digits: the tolerance is 10^-accuracy
(default 6).  Exit code is 0 when the files agree up to the band gauge and
1 otherwise.
"""

import math
import sys

MAX_REPORTED = 10


def fail(message):
    print("CompareMMN Error: " + message)
    sys.exit(1)


def is_int(token):
    try:
        int(token)
        return True
    except ValueError:
        return False


def is_float(token):
    try:
        float(token)
        return True
    except ValueError:
        return False


def read_mmn(path):
    """Return ((nbnd, nk, nntot), blocks) from a Wannier90 .mmn file.

    blocks maps (ik, ikb) -> (G, matrix) with matrix[m][n] = M_mn.
    Lines before the 'nbnd nk nntot' header (e.g. ABACUS's
    'Created on ...' preamble) are ignored.
    """
    try:
        with open(path) as handle:
            lines = handle.readlines()
    except OSError:
        fail("can not find file %s" % path)

    header = None
    lineno = 0
    for index, line in enumerate(lines):
        tokens = line.split()
        if not tokens:
            continue
        if len(tokens) == 3 and all(is_int(t) for t in tokens):
            nbnd, nk, nntot = (int(t) for t in tokens)
            if nbnd <= 0 or nk <= 0 or nntot <= 0:
                fail("%s:%d: invalid header '%s'" % (path, index + 1, line.strip()))
            header = (nbnd, nk, nntot)
            lineno = index + 1
            break
        if len(tokens) == 5 and all(is_int(t) for t in tokens):
            fail("%s:%d: block header before the nbnd nk nntot header"
                 % (path, index + 1))
        if len(tokens) == 2 and all(is_float(t) for t in tokens):
            fail("%s:%d: data line before the header" % (path, index + 1))
        # Otherwise this is a preamble line such as "Created on ...".
    if header is None:
        fail("%s: no valid header 'nbnd nk nntot' found" % path)

    nbnd, nk, nntot = header
    blocks = {}
    counts = [0] * nk
    total = nk * nntot

    def next_content_line(describe):
        nonlocal lineno
        while lineno < len(lines) and not lines[lineno].strip():
            lineno += 1
        if lineno >= len(lines):
            fail("%s: unexpected end of file while %s" % (path, describe))
        line = lines[lineno]
        lineno += 1
        return line

    for _ in range(total):
        line = next_content_line("reading block header %d of %d"
                                 % (len(blocks) + 1, total))
        tokens = line.split()
        if len(tokens) != 5 or not all(is_int(t) for t in tokens):
            fail("%s:%d: expected block header 'ik ikb G1 G2 G3', got '%s'"
                 % (path, lineno, line.strip()))
        ik, ikb, gx, gy, gz = (int(t) for t in tokens)
        if not (1 <= ik <= nk and 1 <= ikb <= nk):
            fail("%s:%d: ik/ikb out of range (ik=%d ikb=%d, nk=%d)"
                 % (path, lineno, ik, ikb, nk))
        if (ik, ikb) in blocks:
            fail("%s:%d: duplicate block ik=%d ikb=%d" % (path, lineno, ik, ikb))

        matrix = [[0j] * nbnd for _ in range(nbnd)]
        for n in range(nbnd):
            for m in range(nbnd):
                line = next_content_line(
                    "reading M(%d,%d) of block ik=%d ikb=%d" % (m + 1, n + 1, ik, ikb))
                tokens = line.split()
                if len(tokens) != 2 or not all(is_float(t) for t in tokens):
                    fail("%s:%d: expected 'Re Im' in block ik=%d ikb=%d, got '%s'"
                         % (path, lineno, ik, ikb, line.strip()))
                real = float(tokens[0])
                imag = float(tokens[1])
                if not (math.isfinite(real) and math.isfinite(imag)):
                    fail("%s:%d: non-finite matrix element" % (path, lineno))
                matrix[m][n] = complex(real, imag)

        blocks[(ik, ikb)] = ((gx, gy, gz), matrix)
        counts[ik - 1] += 1

    for index in range(lineno, len(lines)):
        if lines[index].strip():
            fail("%s:%d: unexpected trailing content '%s'"
                 % (path, index + 1, lines[index].strip()))
    for ik in range(nk):
        if counts[ik] != nntot:
            fail("%s: k-point %d has %d blocks, expected %d"
                 % (path, ik + 1, counts[ik], nntot))
    return header, blocks


def compare(ref_path, cal_path, tolerance):
    ref_header, ref_blocks = read_mmn(ref_path)
    cal_header, cal_blocks = read_mmn(cal_path)
    if ref_header != cal_header:
        fail("header mismatch: %s has %s, %s has %s"
             % (ref_path, ref_header, cal_path, cal_header))

    missing = sorted(set(ref_blocks) - set(cal_blocks))
    extra = sorted(set(cal_blocks) - set(ref_blocks))
    if missing or extra:
        fail("block set mismatch: %d missing (e.g. %s), %d extra (e.g. %s)"
             % (len(missing), missing[:1], len(extra), extra[:1]))

    nbnd, nk, nntot = ref_header
    failures = 0
    reported = 0
    max_diff = 0.0
    for ik, ikb in sorted(ref_blocks):
        ref_g, ref_matrix = ref_blocks[(ik, ikb)]
        cal_g, cal_matrix = cal_blocks[(ik, ikb)]
        if ref_g != cal_g:
            fail("block ik=%d ikb=%d: G vector differs, ref=%s cal=%s"
                 % (ik, ikb, ref_g, cal_g))
        for m in range(nbnd):
            for n in range(nbnd):
                diff = abs(abs(ref_matrix[m][n]) - abs(cal_matrix[m][n]))
                max_diff = max(max_diff, diff)
                if diff > tolerance:
                    failures += 1
                    if reported < MAX_REPORTED:
                        print("Mismatch at ik=%d ikb=%d, |M[%d,%d]|: "
                              "ref=%.12g cal=%.12g diff=%.3g"
                              % (ik, ikb, m + 1, n + 1,
                                 abs(ref_matrix[m][n]), abs(cal_matrix[m][n]), diff))
                        reported += 1

    if failures:
        print("MMN check FAILED: %d magnitude components exceed %g "
              "(max diff %g)" % (failures, tolerance, max_diff))
        sys.exit(1)

    print("MMN check PASSED: %d k-points, %d bands, %d b-vectors per k-point; "
          "max |d|M|| = %.3g (tolerance %g)"
          % (nk, nbnd, nntot, max_diff, tolerance))


def main():
    if len(sys.argv) not in (3, 4):
        print(__doc__.strip())
        sys.exit(1)

    accuracy = 6
    if len(sys.argv) == 4:
        try:
            accuracy = int(sys.argv[3])
        except ValueError:
            fail("accuracy must be an integer, got '%s'" % sys.argv[3])
    if accuracy < 0:
        fail("accuracy must be non-negative")
    compare(sys.argv[1], sys.argv[2], 10.0 ** (-accuracy))


if __name__ == "__main__":
    main()
