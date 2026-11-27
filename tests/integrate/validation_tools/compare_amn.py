#!/usr/bin/env python3
"""Gauge-invariant comparison of Wannier90 .amn files.

The overlap matrix elements A_mn(k) = <u_mk|g_n> are covariant under the
band gauge: an arbitrary unitary mixing U of the selected bands maps
A -> U^dagger A.  The Gram matrix A^dagger A is invariant under this
transformation, so it is compared here instead of the raw elements.  This
lets runs whose eigensolvers pick different gauges (e.g. LAPACK on CPU and
cuSOLVER on GPU) be validated against the same reference.

Usage:
    compare_amn.py ref.amn cal.amn [accuracy]

accuracy is the number of decimal digits: the tolerance is 10^-accuracy
(default 6).  Exit code is 0 when the files agree up to the band gauge and
1 otherwise.
"""

import math
import sys

MAX_REPORTED = 10


def fail(message):
    print("CompareAMN Error: " + message)
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


def read_amn(path):
    """Return ((nbnd, nk, nwann), A[k][m][n]) from a Wannier90 .amn file.

    Lines before the 'nbnd nk nwann' header (e.g. ABACUS's
    'Created on ...' line) are ignored.  Every entry is checked for range,
    duplicates and completeness so that a truncated or corrupted output is
    reported instead of silently producing a small Gram matrix.
    """
    try:
        with open(path) as handle:
            lines = handle.readlines()
    except OSError:
        fail("can not find file %s" % path)

    header = None
    data = None
    filled = None
    nentries = 0

    for lineno, line in enumerate(lines, 1):
        tokens = line.split()
        if not tokens:
            continue

        if header is None:
            if len(tokens) == 3 and all(is_int(t) for t in tokens):
                nbnd, nk, nwann = (int(t) for t in tokens)
                if nbnd <= 0 or nk <= 0 or nwann <= 0:
                    fail("%s:%d: invalid header '%s'" % (path, lineno, line.strip()))
                header = (nbnd, nk, nwann)
                data = [[[0j] * nwann for _ in range(nbnd)] for _ in range(nk)]
                filled = [[[False] * nwann for _ in range(nbnd)] for _ in range(nk)]
            elif len(tokens) == 5 and is_int(tokens[0]) \
                    and is_int(tokens[1]) and is_int(tokens[2]) \
                    and is_float(tokens[3]) and is_float(tokens[4]):
                fail("%s:%d: data line before the header" % (path, lineno))
            # Otherwise this is a preamble line such as "Created on ...".
            continue

        if len(tokens) != 5:
            fail("%s:%d: expected 5 columns 'm n k Re Im', got %d"
                 % (path, lineno, len(tokens)))
        if not (is_int(tokens[0]) and is_int(tokens[1]) and is_int(tokens[2])
                and is_float(tokens[3]) and is_float(tokens[4])):
            fail("%s:%d: malformed data line '%s'" % (path, lineno, line.strip()))

        m = int(tokens[0])
        n = int(tokens[1])
        k = int(tokens[2])
        real = float(tokens[3])
        imag = float(tokens[4])
        if not (math.isfinite(real) and math.isfinite(imag)):
            fail("%s:%d: non-finite matrix element" % (path, lineno))

        nbnd, nk, nwann = header
        if not (1 <= m <= nbnd and 1 <= n <= nwann and 1 <= k <= nk):
            fail("%s:%d: index out of range (m=%d n=%d k=%d, header %d %d %d)"
                 % (path, lineno, m, n, k, nbnd, nk, nwann))
        if filled[k - 1][m - 1][n - 1]:
            fail("%s:%d: duplicate entry m=%d n=%d k=%d" % (path, lineno, m, n, k))
        data[k - 1][m - 1][n - 1] = complex(real, imag)
        filled[k - 1][m - 1][n - 1] = True
        nentries += 1

    if header is None:
        fail("%s: no valid header 'nbnd nk nwann' found" % path)

    expected = header[0] * header[1] * header[2]
    if nentries != expected:
        fail("%s: found %d entries, expected %d" % (path, nentries, expected))
    for k in range(header[1]):
        for m in range(header[0]):
            for n in range(header[2]):
                if not filled[k][m][n]:
                    fail("%s: missing entry m=%d n=%d k=%d" % (path, m + 1, n + 1, k + 1))
    return header, data


def gram_matrix(matrix):
    """Return A^dagger A for a (nbnd x nwann) matrix."""
    nbnd = len(matrix)
    nwann = len(matrix[0])
    gram = [[0j] * nwann for _ in range(nwann)]
    for m in range(nbnd):
        for i in range(nwann):
            ai_conj = matrix[m][i].conjugate()
            for j in range(nwann):
                gram[i][j] += ai_conj * matrix[m][j]
    return gram


def main():
    if len(sys.argv) not in (3, 4):
        print(__doc__.strip())
        sys.exit(1)

    ref_path = sys.argv[1]
    cal_path = sys.argv[2]
    accuracy = 6
    if len(sys.argv) == 4:
        try:
            accuracy = int(sys.argv[3])
        except ValueError:
            fail("accuracy must be an integer, got '%s'" % sys.argv[3])
    if accuracy < 0:
        fail("accuracy must be non-negative")
    tolerance = 10.0 ** (-accuracy)

    ref_header, ref_data = read_amn(ref_path)
    cal_header, cal_data = read_amn(cal_path)
    if ref_header != cal_header:
        fail("header mismatch: %s has %s, %s has %s"
             % (ref_path, ref_header, cal_path, cal_header))

    nbnd, nk, nwann = ref_header
    failures = 0
    reported = 0
    max_diff = 0.0
    for k in range(nk):
        ref_gram = gram_matrix(ref_data[k])
        cal_gram = gram_matrix(cal_data[k])
        for i in range(nwann):
            for j in range(nwann):
                diff_real = abs(ref_gram[i][j].real - cal_gram[i][j].real)
                diff_imag = abs(ref_gram[i][j].imag - cal_gram[i][j].imag)
                max_diff = max(max_diff, diff_real, diff_imag)
                if diff_real > tolerance or diff_imag > tolerance:
                    failures += 1
                    if reported < MAX_REPORTED:
                        print("Mismatch at k=%d, (A^dagger A)[%d,%d]: "
                              "ref=(%.12g,%.12g) cal=(%.12g,%.12g) diff=(%.3g,%.3g)"
                              % (k + 1, i + 1, j + 1,
                                 ref_gram[i][j].real, ref_gram[i][j].imag,
                                 cal_gram[i][j].real, cal_gram[i][j].imag,
                                 diff_real, diff_imag))
                        reported += 1

    if failures:
        print("AMN check FAILED: %d Gram-matrix components exceed %g "
              "(max diff %g)" % (failures, tolerance, max_diff))
        sys.exit(1)

    print("AMN check PASSED: %d k-points, %d bands, %d Wannier functions; "
          "max |d(A^dagger A)| = %.3g (tolerance %g)"
          % (nk, nbnd, nwann, max_diff, tolerance))


if __name__ == "__main__":
    main()
