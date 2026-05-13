# CI Test Failure Summary

**Date:** 2026-05-12
**Branch:** feat/dftu-pw-port-v2
**Status:** Failed

## Critical Errors

The following test cases failed with energy deviations exceeding the threshold. These errors indicate significant discrepancies in the DeltaSpin implementation for specific configurations (Collinear, DFT+U, ReadLam).

| # | Case | Mode | Deviation (eV) | Calculated | Reference |
|---|------|------|----------------|------------|-----------|
| 12 | `12_PW_DS_S2_Z` | DeltaSpin (S2) | ~904.56 | -36051.41 | -36955.97 |
| 18 | `18_PW_DFTU_DS_S2_Z` | DFT+U+DS (S2) | ~4301.41 | -61644.94 | -57343.53 |
| 36 | `36_PW_DS_S2_ReadLam_Z` | ReadLam (S2) | ~1035.27 | -6368.96 | -5333.69 |
| 37 | `37_PW_DS_S4_ReadLam_XY` | ReadLam (S4) | ~1035.55 | -6370.63 | -5335.07 |

## Observations

- **ReadLam Cases (36/37):** Both show ~1035 eV deviation. This suggests a consistent issue with the lambda reading/restart logic or the energy contribution calculation in this mode.
- **Collinear DeltaSpin (12):** Significant energy mismatch (~900 eV).
- **DFT+U + DeltaSpin (18):** The largest deviation (~4300 eV), suggesting the interaction between DFT+U correction and DeltaSpin potential is problematic in collinear mode.

## Warnings

The following cases passed but showed deviations in the range of 1e-4 to 1e-2 eV. These are currently within tolerance but should be monitored.

| # | Case | Deviation |
|---|------|-----------|
| 14 | `14_PW_DS_S4_XYZ` | ~1.1e-4 |
| 17 | `17_PW_DS_S4_XYZ` | ~1.9e-5 |
| 22 | `22_PW_DFTU_DS_S4_XY` | ~1.2e-7 |
| 25 | `25_LCAO_DS_S4_XY` | ~1.3e-5 |
| 26 | `26_LCAO_DS_S4_XYZ` | ~7.6e-6 |
| 27 | `27_LCAO_DS_S4_Z` | ~4.6e-5 |
| 28 | `28_LCAO_DS_S4_XY` | ~3.7e-5 |
| 29 | `29_LCAO_DS_S4_XYZ` | ~5.6e-5 |
| 31 | `31_LCAO_DFTU_DS_S4_XY` | ~3.8e-4 |
| 32 | `32_LCAO_DFTU_DS_S4_XYZ` | ~4.4e-3 |
| 33 | `33_LCAO_DFTU_DS_S4_Z` | ~3.8e-4 |
| 34 | `34_LCAO_DFTU_DS_S4_XY` | ~3.8e-4 |
| 35 | `35_LCAO_DFTU_DS_S4_XYZ` | ~4.4e-3 |

## New Test Cases

The following test cases were added in commit `795454a1b`:

- **53-54:** Magnetic Moment & Lambda Verification
- **55:** NSCF Mode
- **56-59:** Direction Only Constraint
- **60:** NSCF Band Structure

*Note: CI log for new cases (53-60) was truncated in the provided report.*
