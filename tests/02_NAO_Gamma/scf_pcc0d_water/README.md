Vacuum PCC 0D regression with `assume_isolated pcc_0d`.

The integration harness compares total energy, ionic force magnitudes and `E_pcc`
(in eV). References were generated with two MPI ranks, one OpenMP thread and
the repository pseudopotentials/orbitals. Numerical kernels have independent
analytic and finite-difference unit tests; PW neutral/charged results were also
compared against the original full PCC implementation.
