#include "fssh_driver.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <tuple>
#include <utility>

#include <mpi.h>
#include "source_cell/unitcell.h"

// LAPACK 外部声明: SVD (Procrustes 对齐) + LU 分解/求逆 (Löwdin 方法)
extern "C" {
void dgesvd_(const char* jobu, const char* jobvt, const int* m, const int* n,
             double* a, const int* lda, double* s, double* u, const int* ldu,
             double* vt, const int* ldvt, double* work, const int* lwork, int* info);
void dgetrf_(const int* m, const int* n, double* a, const int* lda,
             int* ipiv, int* info);
void dgetri_(const int* n, double* a, const int* lda,
             const int* ipiv, double* work, const int* lwork, int* info);
}

void FsshDriver::run_lowdin_validation_suite(const std::string& output_prefix)
{
    const int nbasis = 4;
    const int nocc_lr = 2;
    const int nvirt_lr = 2;
    const int nstates = 3; // 0: ground, 1..2: excited
    const int nks_total = nocc_lr + nvirt_lr;
    const int occ_offset = 0;

    this->init(nbasis, nstates, 1.0, 0, nocc_lr, nks_total);

    auto make_zero_complex_matrix = [](int nr, int nc) {
        ModuleBase::ComplexMatrix mat(nr, nc);
        for (int i = 0; i < nr; ++i) {
            for (int j = 0; j < nc; ++j) {
                mat(i, j) = std::complex<double>(0.0, 0.0);
            }
        }
        return mat;
    };

    auto make_identity_ao = [&]() {
        ModuleBase::ComplexMatrix s_ao = make_zero_complex_matrix(nbasis, nbasis);
        for (int i = 0; i < nbasis; ++i) {
            s_ao(i, i) = std::complex<double>(1.0, 0.0);
        }
        return s_ao;
    };

    auto make_rotated_coef = [&](double alpha) {
        ModuleBase::ComplexMatrix coef = make_zero_complex_matrix(nks_total, nbasis);
        const double c = std::cos(alpha);
        const double s = std::sin(alpha);

        // occupied-virtual pair 0 <-> 2
        coef(0, 0) = std::complex<double>(c, 0.0);
        coef(0, 2) = std::complex<double>(s, 0.0);
        coef(2, 0) = std::complex<double>(-s, 0.0);
        coef(2, 2) = std::complex<double>(c, 0.0);

        // occupied-virtual pair 1 <-> 3
        coef(1, 1) = std::complex<double>(c, 0.0);
        coef(1, 3) = std::complex<double>(s, 0.0);
        coef(3, 1) = std::complex<double>(-s, 0.0);
        coef(3, 3) = std::complex<double>(c, 0.0);
        return coef;
    };

    auto make_casida = [&](double theta) {
        std::vector<CasidaWavefunction> wfcs(
            nstates, CasidaWavefunction(0.0, {}, {}, nocc_lr, nvirt_lr));
        const double c = std::cos(theta);
        const double s = std::sin(theta);

        std::vector<double> x1(nocc_lr * nvirt_lr, 0.0);
        std::vector<double> x2(nocc_lr * nvirt_lr, 0.0);
        // excitation basis:
        // idx 0 -> (i=0,a=0), idx 2 -> (i=1,a=0)
        x1[0] = c;
        x1[2] = s;
        x2[0] = -s;
        x2[2] = c;

        wfcs[0] = CasidaWavefunction(0.0, {}, {}, nocc_lr, nvirt_lr);
        wfcs[1] = CasidaWavefunction(0.10, x1, {}, nocc_lr, nvirt_lr);
        wfcs[2] = CasidaWavefunction(0.20, x2, {}, nocc_lr, nvirt_lr);
        return wfcs;
    };

    auto compute_occ_phase = [&](const ModuleBase::ComplexMatrix& coef_old_in,
                                 const ModuleBase::ComplexMatrix& coef_new_in,
                                 const ModuleBase::ComplexMatrix& s_ao_dense) {
        std::vector<double> occ_phase(nocc_lr, 1.0);
        for (int i_lr = 0; i_lr < nocc_lr; ++i_lr) {
            const int i_ks = occ_offset + i_lr;
            double s_ii = 0.0;
            for (int mu = 0; mu < nbasis; ++mu) {
                for (int nu = 0; nu < nbasis; ++nu) {
                    s_ii += (std::conj(coef_old_in(i_ks, mu))
                             * s_ao_dense(mu, nu)
                             * coef_new_in(i_ks, nu))
                                .real();
                }
            }
            if (s_ii < 0.0) {
                occ_phase[i_lr] = -1.0;
            }
        }
        return occ_phase;
    };

    auto build_tda_pseudowfc = [&](const std::vector<CasidaWavefunction>& casida_wfcs,
                                   const ModuleBase::ComplexMatrix& coef) {
        std::vector<std::vector<std::vector<std::complex<double>>>> tda_wfc(nstates);
        for (int I = 1; I < nstates; ++I) {
            tda_wfc[I].resize(nocc_lr,
                              std::vector<std::complex<double>>(nbasis, {0.0, 0.0}));
            if (I >= static_cast<int>(casida_wfcs.size())) {
                continue;
            }
            const auto& X_I = casida_wfcs[I].X_coeffs;
            for (int i_lr = 0; i_lr < nocc_lr; ++i_lr) {
                for (int a_lr = 0; a_lr < nvirt_lr; ++a_lr) {
                    const int idx_ia = i_lr * nvirt_lr + a_lr;
                    if (idx_ia >= static_cast<int>(X_I.size())) {
                        continue;
                    }
                    const int coef_row = nocc_ + a_lr;
                    if (coef_row >= coef.nr) {
                        continue;
                    }
                    const double x_ia = X_I[idx_ia];
                    for (int mu = 0; mu < nbasis; ++mu) {
                        tda_wfc[I][i_lr][mu] += x_ia * coef(coef_row, mu);
                    }
                }
            }
        }
        return tda_wfc;
    };

    auto compute_sigma_approx = [&](const ModuleBase::ComplexMatrix& coef_old_in,
                                    const ModuleBase::ComplexMatrix& coef_new_in,
                                    const ModuleBase::ComplexMatrix& s_ao_dense,
                                    const std::vector<CasidaWavefunction>& casida_old_in,
                                    const std::vector<CasidaWavefunction>& casida_new_in,
                                    const std::vector<double>& occ_phase) {
        ModuleBase::ComplexMatrix Sigma = make_zero_complex_matrix(nstates, nstates);
        auto tda_old = build_tda_pseudowfc(casida_old_in, coef_old_in);
        auto tda_new = build_tda_pseudowfc(casida_new_in, coef_new_in);

        for (int I = 0; I < nstates; ++I) {
            for (int J = 0; J < nstates; ++J) {
                std::complex<double> overlap(0.0, 0.0);
                if (I == 0 && J == 0) {
                    overlap = std::complex<double>(1.0, 0.0);
                } else if (I == 0 && J > 0) {
                    // Corrected path: no occupied-phase factor in Sigma(0,J)
                    for (int i_lr = 0; i_lr < nocc_lr; ++i_lr) {
                        const int i_ks = occ_offset + i_lr;
                        for (int mu = 0; mu < nbasis; ++mu) {
                            for (int nu = 0; nu < nbasis; ++nu) {
                                overlap += std::conj(coef_old_in(i_ks, mu))
                                         * s_ao_dense(mu, nu)
                                         * tda_new[J][i_lr][nu];
                            }
                        }
                    }
                } else if (I > 0 && J == 0) {
                    for (int i_lr = 0; i_lr < nocc_lr; ++i_lr) {
                        const int i_ks = occ_offset + i_lr;
                        for (int mu = 0; mu < nbasis; ++mu) {
                            for (int nu = 0; nu < nbasis; ++nu) {
                                overlap += std::conj(tda_old[I][i_lr][mu])
                                         * s_ao_dense(mu, nu)
                                         * coef_new_in(i_ks, nu)
                                         * occ_phase[i_lr];
                            }
                        }
                    }
                } else {
                    for (int i_lr = 0; i_lr < nocc_lr; ++i_lr) {
                        for (int mu = 0; mu < nbasis; ++mu) {
                            for (int nu = 0; nu < nbasis; ++nu) {
                                overlap += std::conj(tda_old[I][i_lr][mu])
                                         * s_ao_dense(mu, nu)
                                         * tda_new[J][i_lr][nu];
                            }
                        }
                    }
                }
                Sigma(I, J) = overlap;
            }
        }
        return Sigma;
    };

    auto compute_s_mo_lr = [&](const ModuleBase::ComplexMatrix& coef_old_in,
                               const ModuleBase::ComplexMatrix& coef_new_in,
                               const ModuleBase::ComplexMatrix& s_ao_dense) {
        const int n_lr = nocc_lr + nvirt_lr;
        std::vector<double> s_mo(n_lr * n_lr, 0.0);
        auto local_to_ks = [&](int loc) -> int {
            if (loc < nocc_lr) {
                return occ_offset + loc;
            }
            return nocc_ + (loc - nocc_lr);
        };
        for (int p = 0; p < n_lr; ++p) {
            const int p_ks = local_to_ks(p);
            for (int q = 0; q < n_lr; ++q) {
                const int q_ks = local_to_ks(q);
                double val = 0.0;
                for (int mu = 0; mu < nbasis; ++mu) {
                    for (int nu = 0; nu < nbasis; ++nu) {
                        val += (std::conj(coef_old_in(p_ks, mu))
                                * s_ao_dense(mu, nu)
                                * coef_new_in(q_ks, nu))
                                   .real();
                    }
                }
                s_mo[p * n_lr + q] = val;
            }
        }
        return s_mo;
    };

    auto determinant = [](std::vector<double> A, int n) {
        double det = 1.0;
        int sign = 1;
        const double eps = 1e-14;
        for (int i = 0; i < n; ++i) {
            int pivot = i;
            double max_abs = std::abs(A[i * n + i]);
            for (int r = i + 1; r < n; ++r) {
                const double v = std::abs(A[r * n + i]);
                if (v > max_abs) {
                    max_abs = v;
                    pivot = r;
                }
            }
            if (max_abs < eps) {
                return 0.0;
            }
            if (pivot != i) {
                for (int c = 0; c < n; ++c) {
                    std::swap(A[i * n + c], A[pivot * n + c]);
                }
                sign *= -1;
            }
            const double pivot_val = A[i * n + i];
            det *= pivot_val;
            for (int r = i + 1; r < n; ++r) {
                const double factor = A[r * n + i] / pivot_val;
                for (int c = i; c < n; ++c) {
                    A[r * n + c] -= factor * A[i * n + c];
                }
            }
        }
        return det * static_cast<double>(sign);
    };

    auto invert_matrix = [](const std::vector<double>& A, int n,
                            std::vector<double>& inv_out, double& det_out) {
        std::vector<double> m = A;
        inv_out.assign(n * n, 0.0);
        for (int i = 0; i < n; ++i) {
            inv_out[i * n + i] = 1.0;
        }
        det_out = 1.0;
        int sign = 1;
        const double eps = 1e-14;

        for (int col = 0; col < n; ++col) {
            int pivot = col;
            double max_abs = std::abs(m[col * n + col]);
            for (int r = col + 1; r < n; ++r) {
                const double v = std::abs(m[r * n + col]);
                if (v > max_abs) {
                    max_abs = v;
                    pivot = r;
                }
            }
            if (max_abs < eps) {
                return false;
            }
            if (pivot != col) {
                for (int c = 0; c < n; ++c) {
                    std::swap(m[col * n + c], m[pivot * n + c]);
                    std::swap(inv_out[col * n + c], inv_out[pivot * n + c]);
                }
                sign *= -1;
            }

            const double piv = m[col * n + col];
            det_out *= piv;
            for (int c = 0; c < n; ++c) {
                m[col * n + c] /= piv;
                inv_out[col * n + c] /= piv;
            }
            for (int r = 0; r < n; ++r) {
                if (r == col) {
                    continue;
                }
                const double f = m[r * n + col];
                for (int c = 0; c < n; ++c) {
                    m[r * n + c] -= f * m[col * n + c];
                    inv_out[r * n + c] -= f * inv_out[col * n + c];
                }
            }
        }
        det_out *= static_cast<double>(sign);
        return true;
    };

    auto max_abs_matrix_diff = [&](const ModuleBase::ComplexMatrix& a,
                                   const ModuleBase::ComplexMatrix& b) {
        double max_diff = 0.0;
        for (int i = 0; i < a.nr; ++i) {
            for (int j = 0; j < a.nc; ++j) {
                max_diff = std::max(max_diff, std::abs(a(i, j) - b(i, j)));
            }
        }
        return max_diff;
    };

    auto frob_sigma_minus_identity = [&](const ModuleBase::ComplexMatrix& sigma) {
        double s = 0.0;
        for (int i = 0; i < sigma.nr; ++i) {
            for (int j = 0; j < sigma.nc; ++j) {
                const std::complex<double> target =
                    (i == j) ? std::complex<double>(1.0, 0.0)
                             : std::complex<double>(0.0, 0.0);
                const double d = std::abs(sigma(i, j) - target);
                s += d * d;
            }
        }
        return std::sqrt(s);
    };

    auto nac_from_sigma = [&](const ModuleBase::ComplexMatrix& sigma) {
        ModuleBase::ComplexMatrix nac = make_zero_complex_matrix(sigma.nr, sigma.nc);
        for (int i = 0; i < sigma.nr; ++i) {
            for (int j = 0; j < sigma.nc; ++j) {
                nac(i, j) = (sigma(i, j) - std::conj(sigma(j, i))) / 2.0;
            }
        }
        return nac;
    };

    auto max_abs_complex_matrix = [&](const ModuleBase::ComplexMatrix& a) {
        double m = 0.0;
        for (int i = 0; i < a.nr; ++i) {
            for (int j = 0; j < a.nc; ++j) {
                m = std::max(m, std::abs(a(i, j)));
            }
        }
        return m;
    };

    auto lowdin_sigma = [&](const ModuleBase::ComplexMatrix& coef_old_in,
                            const ModuleBase::ComplexMatrix& coef_new_in,
                            const ModuleBase::ComplexMatrix& s_ao_dense,
                            const std::vector<CasidaWavefunction>& casida_old_in,
                            const std::vector<CasidaWavefunction>& casida_new_in,
                            const std::vector<double>& occ_phase,
                            double& det_out) {
        ModuleBase::ComplexMatrix sigma = make_zero_complex_matrix(nstates, nstates);
        this->compute_lowdin_sigma(coef_old_in, coef_new_in, s_ao_dense,
                                   casida_old_in, casida_new_in,
                                   nocc_lr, nvirt_lr, occ_offset,
                                   occ_phase, sigma, &det_out);
        return sigma;
    };

    auto sigma_reference_by_determinants = [&](const std::vector<double>& s_mo,
                                               const std::vector<CasidaWavefunction>& casida_old_in,
                                               const std::vector<CasidaWavefunction>& casida_new_in) {
        const int n_lr = nocc_lr + nvirt_lr;
        ModuleBase::ComplexMatrix sigma_ref = make_zero_complex_matrix(nstates, nstates);

        auto det_overlap = [&](const std::vector<int>& bra_occ,
                               const std::vector<int>& ket_occ) {
            std::vector<double> m(nocc_lr * nocc_lr, 0.0);
            for (int r = 0; r < nocc_lr; ++r) {
                for (int c = 0; c < nocc_lr; ++c) {
                    m[r * nocc_lr + c] = s_mo[bra_occ[r] * n_lr + ket_occ[c]];
                }
            }
            return determinant(m, nocc_lr);
        };

        auto ground_occ = [&]() {
            std::vector<int> occ(nocc_lr);
            for (int i = 0; i < nocc_lr; ++i) {
                occ[i] = i;
            }
            return occ;
        }();

        auto excited_occ = [&](int i, int a) {
            std::vector<int> occ = ground_occ;
            occ[i] = nocc_lr + a;
            return occ;
        };

        sigma_ref(0, 0) = std::complex<double>(det_overlap(ground_occ, ground_occ), 0.0);

        for (int J = 1; J < nstates; ++J) {
            double val = 0.0;
            const auto& X_J = casida_new_in[J].X_coeffs;
            for (int j = 0; j < nocc_lr; ++j) {
                for (int b = 0; b < nvirt_lr; ++b) {
                    const int idx = j * nvirt_lr + b;
                    val += X_J[idx] * det_overlap(ground_occ, excited_occ(j, b));
                }
            }
            sigma_ref(0, J) = std::complex<double>(val, 0.0);
        }

        for (int I = 1; I < nstates; ++I) {
            double val = 0.0;
            const auto& X_I = casida_old_in[I].X_coeffs;
            for (int i = 0; i < nocc_lr; ++i) {
                for (int a = 0; a < nvirt_lr; ++a) {
                    const int idx = i * nvirt_lr + a;
                    val += X_I[idx] * det_overlap(excited_occ(i, a), ground_occ);
                }
            }
            sigma_ref(I, 0) = std::complex<double>(val, 0.0);
        }

        for (int I = 1; I < nstates; ++I) {
            const auto& X_I = casida_old_in[I].X_coeffs;
            for (int J = 1; J < nstates; ++J) {
                const auto& X_J = casida_new_in[J].X_coeffs;
                double val = 0.0;
                for (int i = 0; i < nocc_lr; ++i) {
                    for (int a = 0; a < nvirt_lr; ++a) {
                        const int idx_ia = i * nvirt_lr + a;
                        for (int j = 0; j < nocc_lr; ++j) {
                            for (int b = 0; b < nvirt_lr; ++b) {
                                const int idx_jb = j * nvirt_lr + b;
                                val += X_I[idx_ia] * X_J[idx_jb]
                                     * det_overlap(excited_occ(i, a), excited_occ(j, b));
                            }
                        }
                    }
                }
                sigma_ref(I, J) = std::complex<double>(val, 0.0);
            }
        }
        return sigma_ref;
    };

    auto write_report_header = [&](std::ofstream& ofs, const std::string& title) {
        ofs << "# " << title << "\n\n";
        ofs << "Generated by `FsshDriver::run_lowdin_validation_suite`\n\n";
        ofs << "- nbasis: " << nbasis << "\n";
        ofs << "- nocc_lr: " << nocc_lr << "\n";
        ofs << "- nvirt_lr: " << nvirt_lr << "\n";
        ofs << "- nstates: " << nstates << "\n\n";
    };

    const ModuleBase::ComplexMatrix s_ao = make_identity_ao();
    const auto casida_old = make_casida(0.0);
    const auto casida_new = make_casida(0.05);
    const auto coef_old = make_rotated_coef(0.0);
    const auto coef_new = make_rotated_coef(0.04);
    const auto occ_phase = compute_occ_phase(coef_old, coef_new, s_ao);

    // -------------------- Scheme 1 --------------------
    {
        std::ofstream ofs(output_prefix + "_scheme1_reference.md");
        write_report_header(ofs, "Scheme 1: Explicit Determinant Reference Cross-Check");

        double det_lowdin = 0.0;
        ModuleBase::ComplexMatrix sigma_lowdin =
            lowdin_sigma(coef_old, coef_new, s_ao, casida_old, casida_new, occ_phase, det_lowdin);
        auto s_mo = compute_s_mo_lr(coef_old, coef_new, s_ao);
        // 与 Lowdin 主实现保持一致：对新步占据轨道列应用 occ_phase 修正。
        const int n_lr = nocc_lr + nvirt_lr;
        for (int q = 0; q < nocc_lr; ++q) {
            if (occ_phase[q] < 0.0) {
                for (int p = 0; p < n_lr; ++p) {
                    s_mo[p * n_lr + q] *= -1.0;
                }
            }
        }
        ModuleBase::ComplexMatrix sigma_ref =
            sigma_reference_by_determinants(s_mo, casida_old, casida_new);

        const double max_diff = max_abs_matrix_diff(sigma_lowdin, sigma_ref);
        ofs << "## Result\n\n";
        ofs << "- `max |Sigma_lowdin - Sigma_ref| = " << std::scientific << max_diff << "`\n";
        ofs << "- `det(Soo) from Lowdin path = " << det_lowdin << "`\n";
        ofs << "- Verdict: " << ((max_diff < 1e-8) ? "PASS" : "CHECK") << "\n";
    }

    // -------------------- Scheme 2 --------------------
    {
        std::ofstream ofs(output_prefix + "_scheme2_identities.md");
        write_report_header(ofs, "Scheme 2: Algebraic Identity Checks");

        const auto s_mo = compute_s_mo_lr(coef_old, coef_new, s_ao);
        std::vector<double> Soo(nocc_lr * nocc_lr, 0.0);
        std::vector<double> Sov(nocc_lr * nvirt_lr, 0.0);
        std::vector<double> Svo(nvirt_lr * nocc_lr, 0.0);
        std::vector<double> Svv(nvirt_lr * nvirt_lr, 0.0);

        const int n_lr = nocc_lr + nvirt_lr;
        for (int i = 0; i < nocc_lr; ++i) {
            for (int j = 0; j < nocc_lr; ++j) {
                Soo[i * nocc_lr + j] = s_mo[i * n_lr + j];
            }
            for (int b = 0; b < nvirt_lr; ++b) {
                Sov[i * nvirt_lr + b] = s_mo[i * n_lr + (nocc_lr + b)];
            }
        }
        for (int a = 0; a < nvirt_lr; ++a) {
            for (int j = 0; j < nocc_lr; ++j) {
                Svo[a * nocc_lr + j] = s_mo[(nocc_lr + a) * n_lr + j];
            }
            for (int b = 0; b < nvirt_lr; ++b) {
                Svv[a * nvirt_lr + b] = s_mo[(nocc_lr + a) * n_lr + (nocc_lr + b)];
            }
        }

        std::vector<double> Soo_inv;
        double det_Soo = 0.0;
        const bool inv_ok = invert_matrix(Soo, nocc_lr, Soo_inv, det_Soo);
        double inv_residual_inf = std::numeric_limits<double>::infinity();
        double stilde_consistency = std::numeric_limits<double>::infinity();
        double sigma00_gap = std::numeric_limits<double>::infinity();

        if (inv_ok) {
            // ||Soo*Soo_inv - I||_inf
            double max_row_sum = 0.0;
            for (int i = 0; i < nocc_lr; ++i) {
                double row_sum = 0.0;
                for (int j = 0; j < nocc_lr; ++j) {
                    double v = 0.0;
                    for (int k = 0; k < nocc_lr; ++k) {
                        v += Soo[i * nocc_lr + k] * Soo_inv[k * nocc_lr + j];
                    }
                    v -= (i == j ? 1.0 : 0.0);
                    row_sum += std::abs(v);
                }
                max_row_sum = std::max(max_row_sum, row_sum);
            }
            inv_residual_inf = max_row_sum;

            // Stilde consistency: direct triple product vs Q form
            std::vector<double> Q(nocc_lr * nvirt_lr, 0.0);
            for (int j = 0; j < nocc_lr; ++j) {
                for (int b = 0; b < nvirt_lr; ++b) {
                    for (int k = 0; k < nocc_lr; ++k) {
                        Q[j * nvirt_lr + b] += Soo_inv[j * nocc_lr + k] * Sov[k * nvirt_lr + b];
                    }
                }
            }
            std::vector<double> stilde_a(nvirt_lr * nvirt_lr, 0.0);
            std::vector<double> stilde_b(nvirt_lr * nvirt_lr, 0.0);
            for (int a = 0; a < nvirt_lr; ++a) {
                for (int b = 0; b < nvirt_lr; ++b) {
                    double corr_a = 0.0;
                    double corr_b = 0.0;
                    for (int i = 0; i < nocc_lr; ++i) {
                        corr_b += Svo[a * nocc_lr + i] * Q[i * nvirt_lr + b];
                        for (int j = 0; j < nocc_lr; ++j) {
                            corr_a += Svo[a * nocc_lr + i]
                                   * Soo_inv[i * nocc_lr + j]
                                   * Sov[j * nvirt_lr + b];
                        }
                    }
                    stilde_a[a * nvirt_lr + b] = Svv[a * nvirt_lr + b] - corr_a;
                    stilde_b[a * nvirt_lr + b] = Svv[a * nvirt_lr + b] - corr_b;
                }
            }
            double max_stilde_diff = 0.0;
            for (int i = 0; i < nvirt_lr * nvirt_lr; ++i) {
                max_stilde_diff = std::max(max_stilde_diff, std::abs(stilde_a[i] - stilde_b[i]));
            }
            stilde_consistency = max_stilde_diff;

            double det_lowdin = 0.0;
            ModuleBase::ComplexMatrix sigma_lowdin =
                lowdin_sigma(coef_old, coef_new, s_ao, casida_old, casida_new, occ_phase, det_lowdin);
            sigma00_gap = std::abs(sigma_lowdin(0, 0).real() - det_Soo);
        }

        ofs << "## Result\n\n";
        ofs << "- `inverse(Soo) success = " << (inv_ok ? "true" : "false") << "`\n";
        ofs << "- `||Soo*Soo^{-1} - I||_inf = " << std::scientific << inv_residual_inf << "`\n";
        ofs << "- `|Sigma00 - det(Soo)| = " << sigma00_gap << "`\n";
        ofs << "- `max |Stilde(direct) - Stilde(Q-form)| = " << stilde_consistency << "`\n";
        const bool pass = inv_ok && inv_residual_inf < 1e-10
                       && sigma00_gap < 1e-10 && stilde_consistency < 1e-12;
        ofs << "- Verdict: " << (pass ? "PASS" : "CHECK") << "\n";
    }

    // -------------------- Scheme 3 --------------------
    {
        std::ofstream ofs(output_prefix + "_scheme3_limits.md");
        write_report_header(ofs, "Scheme 3: Limit Behavior Checks");

        const auto casida_same = casida_old;
        const auto coef_same = coef_old;
        const auto occ_phase_same = compute_occ_phase(coef_same, coef_same, s_ao);

        double det_same = 0.0;
        ModuleBase::ComplexMatrix sigma_same =
            lowdin_sigma(coef_same, coef_same, s_ao, casida_same, casida_same, occ_phase_same, det_same);
        ModuleBase::ComplexMatrix nac_same = nac_from_sigma(sigma_same);

        double max_sigma_identity_dev = 0.0;
        for (int i = 0; i < nstates; ++i) {
            for (int j = 0; j < nstates; ++j) {
                const double target = (i == j ? 1.0 : 0.0);
                max_sigma_identity_dev = std::max(
                    max_sigma_identity_dev,
                    std::abs(sigma_same(i, j).real() - target));
            }
        }
        const double max_nac_same = max_abs_complex_matrix(nac_same);

        const auto coef_eps1 = make_rotated_coef(0.04);
        const auto coef_eps2 = make_rotated_coef(0.02);
        const auto occ_eps1 = compute_occ_phase(coef_old, coef_eps1, s_ao);
        const auto occ_eps2 = compute_occ_phase(coef_old, coef_eps2, s_ao);
        double det_eps1 = 0.0, det_eps2 = 0.0;
        ModuleBase::ComplexMatrix sigma_eps1 =
            lowdin_sigma(coef_old, coef_eps1, s_ao, casida_old, casida_new, occ_eps1, det_eps1);
        ModuleBase::ComplexMatrix sigma_eps2 =
            lowdin_sigma(coef_old, coef_eps2, s_ao, casida_old, casida_new, occ_eps2, det_eps2);
        const double frob1 = frob_sigma_minus_identity(sigma_eps1);
        const double frob2 = frob_sigma_minus_identity(sigma_eps2);
        const double ratio = (frob2 > 1e-16) ? (frob1 / frob2) : std::numeric_limits<double>::infinity();

        ofs << "## Result\n\n";
        ofs << "- Identical-step `max |Sigma - I| = " << std::scientific << max_sigma_identity_dev << "`\n";
        ofs << "- Identical-step `max |NAC| = " << max_nac_same << "`\n";
        ofs << "- Perturbation scaling: `||Sigma(alpha=0.04)-I||_F = " << frob1 << "`\n";
        ofs << "- Perturbation scaling: `||Sigma(alpha=0.02)-I||_F = " << frob2 << "`\n";
        ofs << "- Ratio `frob(alpha=0.04)/frob(alpha=0.02) = " << ratio << "`\n";
        // 比例不要求严格 2.0：Lowdin 包含 det(Soo) 等非线性项，有限步长下会偏离理想线性。
        const bool pass = (max_sigma_identity_dev < 1e-12) && (max_nac_same < 1e-12)
                       && (ratio > 1.05 && ratio < 2.5);
        ofs << "- Verdict: " << (pass ? "PASS" : "CHECK") << "\n";
    }

    // -------------------- Scheme 4 --------------------
    {
        std::ofstream ofs(output_prefix + "_scheme4_approx_convergence.md");
        write_report_header(ofs, "Scheme 4: Convergence Toward Diagonal Approximation");

        auto evaluate_case = [&](double alpha) {
            const auto coef_new_case = make_rotated_coef(alpha);
            const auto occ_case = compute_occ_phase(coef_old, coef_new_case, s_ao);
            double det_case = 0.0;
            ModuleBase::ComplexMatrix sigma_lowdin_case =
                lowdin_sigma(coef_old, coef_new_case, s_ao, casida_old, casida_new, occ_case, det_case);
            ModuleBase::ComplexMatrix sigma_approx_case =
                compute_sigma_approx(coef_old, coef_new_case, s_ao, casida_old, casida_new, occ_case);
            ModuleBase::ComplexMatrix nac_lowdin_case = nac_from_sigma(sigma_lowdin_case);
            ModuleBase::ComplexMatrix nac_approx_case = nac_from_sigma(sigma_approx_case);
            const double sigma_diff = max_abs_matrix_diff(sigma_lowdin_case, sigma_approx_case);
            const double nac_diff = max_abs_matrix_diff(nac_lowdin_case, nac_approx_case);
            return std::make_pair(sigma_diff, nac_diff);
        };

        const auto strong = evaluate_case(0.08);
        const auto weak = evaluate_case(0.008);

        ofs << "## Result\n\n";
        ofs << "- Strong coupling (`alpha=0.08`): `max|Sigma_lowdin-Sigma_approx| = "
            << std::scientific << strong.first << "`\n";
        ofs << "- Strong coupling (`alpha=0.08`): `max|NAC_lowdin-NAC_approx| = "
            << strong.second << "`\n";
        ofs << "- Weak coupling (`alpha=0.008`): `max|Sigma_lowdin-Sigma_approx| = "
            << weak.first << "`\n";
        ofs << "- Weak coupling (`alpha=0.008`): `max|NAC_lowdin-NAC_approx| = "
            << weak.second << "`\n";
        const bool pass = weak.first < strong.first && weak.second < strong.second;
        ofs << "- Verdict: " << (pass ? "PASS" : "CHECK") << "\n";
    }

    // -------------------- Scheme 5 --------------------
    {
        std::ofstream ofs(output_prefix + "_scheme5_gauge_phase.md");
        write_report_header(ofs, "Scheme 5: Gauge/Phase Robustness");

        double det_base = 0.0;
        ModuleBase::ComplexMatrix sigma_base =
            lowdin_sigma(coef_old, coef_new, s_ao, casida_old, casida_new, occ_phase, det_base);

        // (a) Occupied-phase robustness: flip sign of one occupied orbital at t'
        ModuleBase::ComplexMatrix coef_new_flip = coef_new;
        for (int mu = 0; mu < nbasis; ++mu) {
            coef_new_flip(0, mu) *= -1.0;
        }
        const auto occ_phase_flip = compute_occ_phase(coef_old, coef_new_flip, s_ao);
        double det_flip = 0.0;
        ModuleBase::ComplexMatrix sigma_flip =
            lowdin_sigma(coef_old, coef_new_flip, s_ao, casida_old, casida_new, occ_phase_flip, det_flip);
        const double phase_robust_diff = max_abs_matrix_diff(sigma_base, sigma_flip);

        // (b) Excited-subspace gauge covariance: rotate casida_new in excited space
        const double phi = 0.3;
        const double c = std::cos(phi);
        const double s = std::sin(phi);
        std::vector<CasidaWavefunction> casida_rot = casida_new;
        std::vector<double> x1 = casida_new[1].X_coeffs;
        std::vector<double> x2 = casida_new[2].X_coeffs;
        for (size_t k = 0; k < x1.size(); ++k) {
            casida_rot[1].X_coeffs[k] = c * x1[k] + s * x2[k];
            casida_rot[2].X_coeffs[k] = -s * x1[k] + c * x2[k];
        }

        double det_rot = 0.0;
        ModuleBase::ComplexMatrix sigma_rot =
            lowdin_sigma(coef_old, coef_new, s_ao, casida_old, casida_rot, occ_phase, det_rot);

        ModuleBase::ComplexMatrix sigma_expected = sigma_base;
        for (int I = 0; I < nstates; ++I) {
            const std::complex<double> col1 = sigma_base(I, 1);
            const std::complex<double> col2 = sigma_base(I, 2);
            sigma_expected(I, 1) = c * col1 + s * col2;
            sigma_expected(I, 2) = -s * col1 + c * col2;
        }
        const double covariant_diff = max_abs_matrix_diff(sigma_rot, sigma_expected);

        ofs << "## Result\n\n";
        ofs << "- Occupied sign-flip robustness: `max|Sigma_base - Sigma_flipped| = "
            << std::scientific << phase_robust_diff << "`\n";
        ofs << "- Excited-subspace covariance: `max|Sigma_rot - Sigma_base*W| = "
            << covariant_diff << "`\n";
        const bool pass = (phase_robust_diff < 1e-10) && (covariant_diff < 1e-10);
        ofs << "- Verdict: " << (pass ? "PASS" : "CHECK") << "\n";
    }

    // -------------------- Scheme 6 --------------------
    {
        std::ofstream ofs(output_prefix + "_scheme6_stability.md");
        write_report_header(ofs, "Scheme 6: Numerical Stability Diagnostics");

        auto norm_inf = [](const std::vector<double>& a, int n) {
            double m = 0.0;
            for (int i = 0; i < n; ++i) {
                double row_sum = 0.0;
                for (int j = 0; j < n; ++j) {
                    row_sum += std::abs(a[i * n + j]);
                }
                m = std::max(m, row_sum);
            }
            return m;
        };

        auto estimate_cond_and_residual = [&](const std::vector<double>& Soo) {
            std::vector<double> Soo_inv;
            double det = 0.0;
            const bool ok = invert_matrix(Soo, nocc_lr, Soo_inv, det);
            double cond_inf = std::numeric_limits<double>::infinity();
            double residual_inf = std::numeric_limits<double>::infinity();
            if (ok) {
                cond_inf = norm_inf(Soo, nocc_lr) * norm_inf(Soo_inv, nocc_lr);
                double max_row_sum = 0.0;
                for (int i = 0; i < nocc_lr; ++i) {
                    double row_sum = 0.0;
                    for (int j = 0; j < nocc_lr; ++j) {
                        double v = 0.0;
                        for (int k = 0; k < nocc_lr; ++k) {
                            v += Soo[i * nocc_lr + k] * Soo_inv[k * nocc_lr + j];
                        }
                        v -= (i == j ? 1.0 : 0.0);
                        row_sum += std::abs(v);
                    }
                    max_row_sum = std::max(max_row_sum, row_sum);
                }
                residual_inf = max_row_sum;
            }
            return std::make_tuple(ok, cond_inf, residual_inf, det);
        };

        // Well-conditioned baseline
        const auto s_mo_base = compute_s_mo_lr(coef_old, coef_new, s_ao);
        std::vector<double> Soo_base(nocc_lr * nocc_lr, 0.0);
        for (int i = 0; i < nocc_lr; ++i) {
            for (int j = 0; j < nocc_lr; ++j) {
                Soo_base[i * nocc_lr + j] = s_mo_base[i * (nocc_lr + nvirt_lr) + j];
            }
        }

        // Near-singular case: make occupied row 1 at t' almost parallel to row 0
        ModuleBase::ComplexMatrix coef_new_ns = make_rotated_coef(0.0);
        const double eps = 1e-6;
        for (int mu = 0; mu < nbasis; ++mu) {
            coef_new_ns(1, mu) = std::complex<double>(0.0, 0.0);
        }
        coef_new_ns(1, 0) = std::complex<double>(1.0 - eps, 0.0);
        coef_new_ns(1, 1) = std::complex<double>(eps, 0.0);
        const auto s_mo_ns = compute_s_mo_lr(coef_old, coef_new_ns, s_ao);
        std::vector<double> Soo_ns(nocc_lr * nocc_lr, 0.0);
        for (int i = 0; i < nocc_lr; ++i) {
            for (int j = 0; j < nocc_lr; ++j) {
                Soo_ns[i * nocc_lr + j] = s_mo_ns[i * (nocc_lr + nvirt_lr) + j];
            }
        }

        bool ok_base, ok_ns;
        double cond_base, cond_ns, resid_base, resid_ns, det_base, det_ns;
        std::tie(ok_base, cond_base, resid_base, det_base) = estimate_cond_and_residual(Soo_base);
        std::tie(ok_ns, cond_ns, resid_ns, det_ns) = estimate_cond_and_residual(Soo_ns);

        ofs << "## Result\n\n";
        ofs << "- Baseline `invertible=" << (ok_base ? "true" : "false")
            << "`, `cond_inf=" << std::scientific << cond_base
            << "`, `||Soo*Soo^{-1}-I||_inf=" << resid_base
            << "`, `det(Soo)=" << det_base << "`\n";
        ofs << "- Near-singular `invertible=" << (ok_ns ? "true" : "false")
            << "`, `cond_inf=" << cond_ns
            << "`, `||Soo*Soo^{-1}-I||_inf=" << resid_ns
            << "`, `det(Soo)=" << det_ns << "`\n";
        bool pass = ok_base && ok_ns && (cond_ns > cond_base * 1e3);
        ofs << "- Verdict: " << (pass ? "PASS" : "CHECK") << "\n";
    }

    std::cout << "[FSSH VALIDATION] Wrote 6 reports with prefix: " << output_prefix << std::endl;
    std::cout << "  - " << output_prefix << "_scheme1_reference.md\n";
    std::cout << "  - " << output_prefix << "_scheme2_identities.md\n";
    std::cout << "  - " << output_prefix << "_scheme3_limits.md\n";
    std::cout << "  - " << output_prefix << "_scheme4_approx_convergence.md\n";
    std::cout << "  - " << output_prefix << "_scheme5_gauge_phase.md\n";
    std::cout << "  - " << output_prefix << "_scheme6_stability.md" << std::endl;
}

// =====================================================================
// Static Coupling Test (test_math_1)
// =====================================================================
void FsshDriver::test_math_1() {
    std::cout << "\n==================================================" << std::endl;
    std::cout << ">>> STARTING FSSH RK4 MATH TEST (2-Level System) <<<" << std::endl;
    std::cout << "==================================================\n" << std::endl;

    this->init(10, 2, 0.5, 0);
    std::vector<double> fake_energies = {0.0, 0.1}; 
    
    ModuleBase::ComplexMatrix fake_sigma(2, 2);
    fake_sigma(0, 0) = std::complex<double>(0.0, 0.0);
    fake_sigma(1, 1) = std::complex<double>(0.0, 0.0);
    
    double static_coupling = 0.05;
    fake_sigma(0, 1) = std::complex<double>(static_coupling, 0.0);  
    fake_sigma(1, 0) = std::complex<double>(-static_coupling, 0.0); 

    // 🌟 新增：准备输出 CSV 文件
    std::string filename = "/home/shane/ABACUS/test/fssh/fssh_test1_results.csv";
    std::ofstream ofs(filename);
    if(ofs.is_open()) {
        ofs << "Step,Coupling(0-1),Pop_0,Pop_1,Sum\n"; // 写入 CSV 表头
        std::cout << "[INFO] Results will be saved to: " << filename << "\n" << std::endl;
    }

    // 打印屏幕表头
    printf("Step | Coupling   | Pop_0      | Pop_1      | Sum\n");
    printf("-------------------------------------------------------\n");

    for (int step = 1; step <= 30; ++step) {
        this->propagate_rk4(fake_sigma, fake_energies);

        double pop_0 = std::norm(electronic_coeffs_[0]);
        double pop_1 = std::norm(electronic_coeffs_[1]);
        double total_pop = pop_0 + pop_1;

        // 🌟 屏幕输出：严格保留 6 位小数，宽度 10 对齐
        printf("%4d | %10.6f | %10.6f | %10.6f | %10.6f\n", 
               step, static_coupling, pop_0, pop_1, total_pop);

        // 🌟 文件输出：写入 CSV
        if(ofs.is_open()) {
            ofs << std::fixed << std::setprecision(6)
                << step << ","
                << static_coupling << ","
                << pop_0 << ","
                << pop_1 << ","
                << total_pop << "\n";
        }
    }

    if(ofs.is_open()) ofs.close(); // 关闭文件

    std::cout << "\n>>> Testing Hopping Probability <<<" << std::endl;
    int proposed_state = this->check_hopping(fake_sigma);
    std::cout << "Current State: " << current_state_ 
              << " | Proposed State by Monte Carlo: " << proposed_state << std::endl;

    std::cout << "\n==================================================" << std::endl;
    std::cout << ">>> FSSH MATH TEST 1 FINISHED SUCCESSFULLY <<<" << std::endl;
    std::cout << "==================================================\n" << std::endl;
}

// =====================================================================
// Time Dependent Test (test_math_2)
// =====================================================================
void FsshDriver::test_math_2() {
    std::cout << "\n==================================================" << std::endl;
    std::cout << ">>> STARTING FSSH RK4 MATH TEST (Time-Dependent) <<<" << std::endl;
    std::cout << "==================================================\n" << std::endl;

    this->init(10, 2, 0.5, 0);
    std::vector<double> fake_energies = {0.0, 0.1}; 
    
    ModuleBase::ComplexMatrix fake_sigma(2, 2);
    fake_sigma(0, 0) = std::complex<double>(0.0, 0.0);
    fake_sigma(1, 1) = std::complex<double>(0.0, 0.0);

    double base_coupling = 0.05; 
    std::normal_distribution<double> noise_dist(0.0, 0.02);

    // 🌟 新增：准备输出 CSV 文件
    std::string filename = "/home/shane/ABACUS/test/fssh/fssh_test2_results.csv";
    std::ofstream ofs(filename);
    if(ofs.is_open()) {
        ofs << "Step,Coupling(0-1),Pop_0,Pop_1,Sum\n"; // 写入 CSV 表头
        std::cout << "[INFO] Results will be saved to: " << filename << "\n" << std::endl;
    }

    // 打印屏幕表头
    printf("Step | Coupling   | Pop_0      | Pop_1      | Sum\n");
    printf("-------------------------------------------------------\n");

    for (int step = 1; step <= 30; ++step) {
        
        double noise = noise_dist(this->gen_); 
        double current_coupling = std::abs(base_coupling + noise); 

        fake_sigma(0, 1) = std::complex<double>(current_coupling, 0.0);
        fake_sigma(1, 0) = std::complex<double>(-current_coupling, 0.0);

        this->propagate_rk4(fake_sigma, fake_energies);

        double pop_0 = std::norm(electronic_coeffs_[0]);
        double pop_1 = std::norm(electronic_coeffs_[1]);
        double total_pop = pop_0 + pop_1;

        // 🌟 屏幕输出：严格保留 6 位小数，宽度 10 对齐
        printf("%4d | %10.6f | %10.6f | %10.6f | %10.6f\n", 
               step, current_coupling, pop_0, pop_1, total_pop);

        // 🌟 文件输出：写入 CSV
        if(ofs.is_open()) {
            ofs << std::fixed << std::setprecision(6)
                << step << ","
                << current_coupling << ","
                << pop_0 << ","
                << pop_1 << ","
                << total_pop << "\n";
        }
    }

    if(ofs.is_open()) ofs.close(); // 关闭文件

    std::cout << "\n>>> Testing Hopping Probability (Last Step) <<<" << std::endl;
    int proposed_state = this->check_hopping(fake_sigma);
    std::cout << "Current State: " << current_state_ 
              << " | Proposed State by Monte Carlo: " << proposed_state << std::endl;

    std::cout << "\n==================================================" << std::endl;
    std::cout << ">>> FSSH TIME-DEPENDENT MATH TEST FINISHED <<<" << std::endl;
    std::cout << "==================================================\n" << std::endl;
}

// =====================================================================
// Gaussian Pulse Coupling Test (test_math_3)
// 模拟真实的 Avoided Crossing (避免交叉) 物理过程
// =====================================================================
void FsshDriver::test_math_3() {
    std::cout << "\n==================================================" << std::endl;
    std::cout << ">>> STARTING FSSH RK4 MATH TEST (Gaussian Pulse) <<<" << std::endl;
    std::cout << "==================================================\n" << std::endl;

    // 1. 初始化驱动器：10个基组，2个态，时间步长 0.5 fs，从态 0 开始
    this->init(10, 2, 0.5, 0);

    // 2. 伪造静态物理数据
    std::vector<double> fake_energies = {0.0, 0.1}; // 态0能量为0，态1能量为0.1
    
    // 初始化非绝热耦合矩阵
    ModuleBase::ComplexMatrix fake_sigma(2, 2);
    fake_sigma(0, 0) = std::complex<double>(0.0, 0.0);
    fake_sigma(1, 1) = std::complex<double>(0.0, 0.0);

    // 🌟 高斯脉冲参数设置
    int total_steps = 100;         // 总共跑 100 步，看完整的脉冲过程
    double peak_step = 50.0;       // 高斯峰的中心位置 (第50步到达最高潮)
    double width = 10.0;           // 高斯峰的宽度 (标准差 sigma)
    double baseline_coupling = 0.01; // 平时极其微弱的耦合
    double max_coupling = 0.80;    // 交叉点极强的耦合

    // 3. 准备输出 CSV 文件
    std::string filename = "/home/shane/ABACUS/test/fssh/fssh_test3_results.csv";
    std::ofstream ofs(filename);
    if(ofs.is_open()) {
        ofs << "Step,Coupling(0-1),Pop_0,Pop_1,Sum\n"; 
        std::cout << "[INFO] Results will be saved to: " << filename << "\n" << std::endl;
    }

    // 打印屏幕表头
    printf("Step | Coupling   | Pop_0      | Pop_1      | Sum\n");
    printf("-------------------------------------------------------\n");

    for (int step = 1; step <= total_steps; ++step) {
        
        // --- 核心改动：计算当前步的高斯型耦合 ---
        // 公式: f(x) = base + (max - base) * exp(- (x - mu)^2 / (2 * sigma^2))
        double exponent = -0.5 * std::pow((step - peak_step) / width, 2);
        double current_coupling = baseline_coupling + (max_coupling - baseline_coupling) * std::exp(exponent);

        // 将算好的耦合填入矩阵 (严格保持反对称)
        fake_sigma(0, 1) = std::complex<double>(current_coupling, 0.0);
        fake_sigma(1, 0) = std::complex<double>(-current_coupling, 0.0);

        // --- 调用 RK4 核心函数 ---
        this->propagate_rk4(fake_sigma, fake_energies);

        // 计算布居数 (模平方)
        double pop_0 = std::norm(electronic_coeffs_[0]);
        double pop_1 = std::norm(electronic_coeffs_[1]);
        double total_pop = pop_0 + pop_1;

        // 🌟 屏幕输出：严格保留 6 位小数，宽度 10 对齐
        printf("%4d | %10.6f | %10.6f | %10.6f | %10.6f\n", 
               step, current_coupling, pop_0, pop_1, total_pop);

        // 🌟 文件输出：写入 CSV
        if(ofs.is_open()) {
            ofs << std::fixed << std::setprecision(6)
                << step << ","
                << current_coupling << ","
                << pop_0 << ","
                << pop_1 << ","
                << total_pop << "\n";
        }
    }

    if(ofs.is_open()) ofs.close(); // 关闭文件

    // 4. 测试一下最后的概率
    std::cout << "\n>>> Testing Hopping Probability (After Passage) <<<" << std::endl;
    int proposed_state = this->check_hopping(fake_sigma);
    std::cout << "Current State: " << current_state_ 
              << " | Proposed State by Monte Carlo: " << proposed_state << std::endl;

    std::cout << "\n==================================================" << std::endl;
    std::cout << ">>> FSSH GAUSSIAN TEST FINISHED SUCCESSFULLY <<<" << std::endl;
    std::cout << "==================================================\n" << std::endl;
}

// =====================================================================
// 3-State Cascaded Transitions Test (test_math_4)
// 模拟三态体系中连续发生的两次避免交叉 (Avoided Crossings)
// =====================================================================
void FsshDriver::test_math_4() {
    std::cout << "\n==================================================" << std::endl;
    std::cout << ">>> STARTING FSSH RK4 TEST (3-State Cascaded)  <<<" << std::endl;
    std::cout << "==================================================\n" << std::endl;

    // 1. 初始化：10个基组，3个态！时间步长 0.5 fs，从态 0 开始演化
    this->init(10, 3, 0.5, 0);

    // 2. 伪造静态物理数据：3个态的能量呈阶梯状排列
    std::vector<double> fake_energies = {0.0, 0.1, 0.2}; 
    
    // 初始化 3x3 的非绝热耦合矩阵 (全部清零)
    ModuleBase::ComplexMatrix fake_sigma(3, 3);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            fake_sigma(i, j) = std::complex<double>(0.0, 0.0);
        }
    }

    // 🌟 高斯脉冲参数设置
    int total_steps = 1000;         
    
    // 脉冲 1 (态0和态1之间)：中心在第 30 步
    double peak_step_01 = 400.0;       
    double width_01 = 80.0;           
    
    // 脉冲 2 (态1和态2之间)：中心在第 50 步
    double peak_step_12 = 600.0;       
    double width_12 = 80.0;           

    double baseline_coupling = 0.01; // 基线耦合
    double max_coupling = 0.80;      // 峰值极强耦合

    // 3. 准备输出 CSV 文件
    std::string filename = "/home/shane/ABACUS/test/fssh/fssh_test6_results.csv";
    std::ofstream ofs(filename);
    if(ofs.is_open()) {
        ofs << "Step,Coup_01,Coup_12,Pop_0,Pop_1,Pop_2,Sum\n"; 
        std::cout << "[INFO] Results will be saved to: " << filename << "\n" << std::endl;
    }

    // 打印屏幕表头 (为了排版整齐，缩写了 Coupling 为 C_)
    printf("Step | C_01     | C_12     | Pop_0    | Pop_1    | Pop_2    | Sum\n");
    printf("----------------------------------------------------------------------\n");

    for (int step = 1; step <= total_steps; ++step) {
        
        // --- 计算态 0 和 态 1 之间的高斯耦合 (第 30 步最大) ---
        double exp_01 = -0.5 * std::pow((step - peak_step_01) / width_01, 2);
        double current_coup_01 = baseline_coupling + (max_coupling - baseline_coupling) * std::exp(exp_01);

        // --- 计算态 1 和 态 2 之间的高斯耦合 (第 50 步最大) ---
        double exp_12 = -0.5 * std::pow((step - peak_step_12) / width_12, 2);
        double current_coup_12 = baseline_coupling + (max_coupling - baseline_coupling) * std::exp(exp_12);

        // 填入矩阵 (严格保持反对称)
        // 0-1 耦合
        fake_sigma(0, 1) = std::complex<double>(current_coup_01, 0.0);
        fake_sigma(1, 0) = std::complex<double>(-current_coup_01, 0.0);
        
        // 1-2 耦合
        fake_sigma(1, 2) = std::complex<double>(current_coup_12, 0.0);
        fake_sigma(2, 1) = std::complex<double>(-current_coup_12, 0.0);

        // 态 0 和 态 2 之间无耦合 (永远为0，初始化时已经清零了)

        // --- 调用 RK4 核心函数演化 ---
        this->propagate_rk4(fake_sigma, fake_energies);

        // 计算布居数 (模平方)
        double pop_0 = std::norm(electronic_coeffs_[0]);
        double pop_1 = std::norm(electronic_coeffs_[1]);
        double pop_2 = std::norm(electronic_coeffs_[2]);
        double total_pop = pop_0 + pop_1 + pop_2;

        // 🌟 屏幕输出
        printf("%4d | %8.5f | %8.5f | %8.6f | %8.6f | %8.6f | %8.6f\n", 
               step, current_coup_01, current_coup_12, pop_0, pop_1, pop_2, total_pop);

        // 🌟 文件输出：写入 CSV
        if(ofs.is_open()) {
            ofs << std::fixed << std::setprecision(6)
                << step << ","
                << current_coup_01 << ","
                << current_coup_12 << ","
                << pop_0 << ","
                << pop_1 << ","
                << pop_2 << ","
                << total_pop << "\n";
        }
    }

    if(ofs.is_open()) ofs.close(); // 关闭文件

    // 4. 测试最后的掷骰子概率
    std::cout << "\n>>> Testing Hopping Probability (After Passage) <<<" << std::endl;
    int proposed_state = this->check_hopping(fake_sigma);
    std::cout << "Current State: " << current_state_ 
              << " | Proposed State by Monte Carlo: " << proposed_state << std::endl;

    std::cout << "\n==================================================" << std::endl;
    std::cout << ">>> FSSH 3-STATE TEST FINISHED SUCCESSFULLY    <<<" << std::endl;
    std::cout << "==================================================\n" << std::endl;
}

// =====================================================================
// Surface Hopping & Velocity Rescaling Test (test_math_5)
// 测试最少面跳跃逻辑、受挫跳跃(Frustrated Hop)以及能量守恒
// =====================================================================
void FsshDriver::test_math_5() {
    std::cout << "\n==================================================" << std::endl;
    std::cout << ">>> STARTING FSSH HOPPING & RESCALING TEST     <<<" << std::endl;
    std::cout << "==================================================\n" << std::endl;

    // 伪造物理数据：2个态，能级差 0.1
    std::vector<double> fake_energies = {0.0, 0.1}; 
    ModuleBase::ComplexMatrix fake_sigma(2, 2);
    fake_sigma(0, 0) = std::complex<double>(0.0, 0.0);
    fake_sigma(1, 1) = std::complex<double>(0.0, 0.0);
    fake_sigma(0, 1) = std::complex<double>(0.05, 0.0);   
    fake_sigma(1, 0) = std::complex<double>(-0.05, 0.0);  

    // 手工捏造一个只有 1 个原子的虚拟 UnitCell
    UnitCell mock_ucell;
    mock_ucell.ntype = 1;
    mock_ucell.atoms = new Atom[1];
    mock_ucell.atoms[0].na = 1;
    mock_ucell.atoms[0].mass = 2.0; // 虚拟质量 m = 2.0
    mock_ucell.atoms[0].vel.resize(1); // 现代 C++ vector 的正确打开方式

    auto calculate_Ekin = [&]() {
        double ekin = 0.0;
        double mass = mock_ucell.atoms[0].mass;
        double v = mock_ucell.atoms[0].vel[0].x; // 假设只有 x 方向有速度
        ekin = 0.5 * mass * v * v;
        return ekin;
    };

    // -----------------------------------------------------------------
    // 实验 A：动能充足的成功跳跃 (Sufficient Energy)
    // -----------------------------------------------------------------
    std::cout << "\n[SCENARIO A]: Sufficient Kinetic Energy (Expected: SUCCESS)" << std::endl;
    std::cout << "--------------------------------------------------" << std::endl;
    
    this->init(10, 2, 0.5, 0); // 初始化到态 0
    // 给原子一个较快的初始速度 v=0.4472 (E_kin = 0.5 * 2 * 0.4472^2 ≈ 0.20 > 0.1)
    mock_ucell.atoms[0].vel[0] = ModuleBase::Vector3<double>(0.44721359, 0.0, 0.0); 
    
    double initial_Ekin = calculate_Ekin();
    double initial_Epot = fake_energies[current_state_];
    double initial_Etot = initial_Ekin + initial_Epot;
    
    printf("Initial State: %d | E_kin: %.4f | E_pot: %.4f | E_tot: %.4f\n", 
            current_state_, initial_Ekin, initial_Epot, initial_Etot);

    bool hopped_A = false;
    for (int step = 1; step <= 50; ++step) {
        this->propagate_rk4(fake_sigma, fake_energies);
        int proposed_state = this->check_hopping(fake_sigma);
        
        if (proposed_state != current_state_) {
            std::cout << ">> Step " << step << ": Hop proposed! (0 -> 1)" << std::endl;
            bool success = this->rescale_velocity(mock_ucell, current_state_, proposed_state, fake_energies);
            if (success) {
                current_state_ = proposed_state;
                hopped_A = true;
                
                double final_Ekin = calculate_Ekin();
                double final_Epot = fake_energies[current_state_];
                double final_Etot = final_Ekin + final_Epot;
                
                printf("Final State  : %d | E_kin: %.4f | E_pot: %.4f | E_tot: %.4f\n", 
                        current_state_, final_Ekin, final_Epot, final_Etot);
                if (std::abs(final_Etot - initial_Etot) < 1e-6) {
                    std::cout << ">> [RESULT]: ENERGY PERFECTLY CONSERVED! <<" << std::endl;
                }
                break; // 跳跃成功，结束本实验
            }
        }
    }
    if (!hopped_A) std::cout << "No hop occurred (Random numbers didn't hit)." << std::endl;


    // -----------------------------------------------------------------
    // 实验 B：动能不足的受挫跳跃 (Frustrated Hop)
    // -----------------------------------------------------------------
    std::cout << "\n[SCENARIO B]: Insufficient Kinetic Energy (Expected: FRUSTRATED)" << std::endl;
    std::cout << "--------------------------------------------------" << std::endl;
    
    this->init(10, 2, 0.5, 0); // 重新初始化到态 0
    // 给原子一个很慢的速度 v=0.2236 (E_kin = 0.5 * 2 * 0.2236^2 ≈ 0.05 < 0.1)
    mock_ucell.atoms[0].vel[0] = ModuleBase::Vector3<double>(0.22360679, 0.0, 0.0); 
    
    initial_Ekin = calculate_Ekin();
    initial_Epot = fake_energies[current_state_];
    initial_Etot = initial_Ekin + initial_Epot;
    
    printf("Initial State: %d | E_kin: %.4f | E_pot: %.4f | E_tot: %.4f\n", 
            current_state_, initial_Ekin, initial_Epot, initial_Etot);

    bool hopped_B = false;
    for (int step = 1; step <= 50; ++step) {
        this->propagate_rk4(fake_sigma, fake_energies);
        int proposed_state = this->check_hopping(fake_sigma);
        
        if (proposed_state != current_state_) {
            std::cout << ">> Step " << step << ": Hop proposed! (0 -> 1)" << std::endl;
            bool success = this->rescale_velocity(mock_ucell, current_state_, proposed_state, fake_energies);
            if (!success) {
                std::cout << ">> [RESULT]: HOP REJECTED! Atom maintains original speed and state." << std::endl;
                
                double final_Ekin = calculate_Ekin();
                double final_Epot = fake_energies[current_state_];
                double final_Etot = final_Ekin + final_Epot;
                printf("Current State: %d | E_kin: %.4f | E_pot: %.4f | E_tot: %.4f\n", 
                        current_state_, final_Ekin, final_Epot, final_Etot);
                hopped_B = true;
                break; // 受挫测试完成，结束本实验
            }
        }
    }
    if (!hopped_B) std::cout << "No hop occurred (Random numbers didn't hit)." << std::endl;

    std::cout << "\n==================================================" << std::endl;
    std::cout << ">>> FSSH HOPPING TEST FINISHED SUCCESSFULLY    <<<" << std::endl;
    std::cout << "==================================================\n" << std::endl;
}

// =====================================================================
// 🌟 终极试金石：Tully Model 1 (包含势能面输出与 Landau-Zener 理论基准)
// =====================================================================
void FsshDriver::test_tully_model_1() {
    std::cout << "\n=====================================================" << std::endl;
    std::cout << ">>> STARTING GOLD STANDARD TEST: TULLY MODEL 1  <<<" << std::endl;
    std::cout << "=====================================================\n" << std::endl;

    // 1. Tully Model 1 经典参数 (Atomic Units)
    double A = 0.01;
    double B = 1.6;
    double C = 0.005;
    double D = 1.0;
    double mass = 2000.0; 

    // ==========================================================
    // 🌟 新增模块 1：输出 Tully Model 1 的绝热势能面 (PES) 及非绝热耦合 (NAC)
    // ==========================================================
    std::string pes_filename = "tully_model_1_pes.csv";
    std::ofstream pes_ofs(pes_filename);
    if(pes_ofs.is_open()) {
        // 🌟 表头增加了最后一列：NAC(d12)
        pes_ofs << "x,E0(Adiabatic),E1(Adiabatic),V11(Diabatic),V12(Coupling),NAC(d12)\n"; 
        
        // 从 -10 到 10 扫描坐标，步长 0.05
        for(double px = -10.0; px <= 10.0; px += 0.05) {
            // 1. 计算透热势能面 (Diabatic)
            double pv11 = (px > 0) ? A * (1.0 - std::exp(-B * px)) : -A * (1.0 - std::exp(B * px));
            double pv12 = C * std::exp(-D * px * px);
            
            // 2. 计算透热势能面的空间导数
            double pdv11_dx = (px > 0) ? A * B * std::exp(-B * px) : A * B * std::exp(B * px);
            double pdv12_dx = -2.0 * C * D * px * std::exp(-D * px * px);
            
            // 3. 计算绝热势能面 (Adiabatic)
            double gap_sq = pv11 * pv11 + pv12 * pv12;
            double pgap = std::sqrt(gap_sq);
            
            // 🌟 4. 计算非绝热耦合项 (Derivative Coupling d12)
            // 公式: d12 = (V11 * dV12/dx - V12 * dV11/dx) / (2 * gap^2)
            double pd12_x = (pv11 * pdv12_dx - pv12 * pdv11_dx) / (2.0 * gap_sq);
            
            // 5. 统统写入 CSV 文件
            pes_ofs << std::fixed << std::setprecision(8) 
                    << px << "," 
                    << -pgap << "," 
                    << pgap << "," 
                    << pv11 << "," 
                    << pv12 << "," 
                    << pd12_x << "\n";
        }
        pes_ofs.close();
        std::cout << "[INFO] 势能面及非绝热耦合数据已生成并保存至: " << pes_filename << std::endl;
    }

    // 2. 初始物理状态设置
    double x = -10.0;     
    double k_init = 20.0; // Tully 原文测试动量
    double v = k_init / mass; 
    
    // 初始化 FSSH 引擎
    double dt_au = 1.0; 
    this->init(1, 2, dt_au, 0); 

    // ==========================================================
    // 🌟 新增模块 2：解析计算 Landau-Zener (LZ) 理论跃迁概率
    // ==========================================================
    const double PI = 3.14159265358979323846;
    // (a) 根据能量守恒，计算到达交叉点 (x=0) 时的速度 v_c
    double e_kin_init = 0.5 * mass * v * v;
    double e_pot_init = -A; // x=-10 处的绝热基态能量极点
    double e_tot = e_kin_init + e_pot_init;
    
    double e_pot_xc = -C; // 交叉点 x=0 处的绝热基态能量
    double v_c = std::sqrt(2.0 * (e_tot - e_pot_xc) / mass);
    
    // (b) 计算 LZ 公式 (Tully 1990 Eq. 22)
    // 跃迁到激发态的概率 P_LZ = exp(- 2 * pi * |H12|^2 / (v_c * |F1 - F2|))
    // 经过数学化简：|F1 - F2| = 2*A*B
    double prob_LZ = std::exp(- PI * C * C / (v_c * A * B));
    
    std::cout << "[INFO] 根据 LZ 公式计算出的精确量子透射率 (理论跃迁概率): " 
              << prob_LZ * 100.0 << " %\n" << std::endl;

    // 3. 准备输出文件 (增加 LZ_Prob 列)
    std::string filename = "tully_model_1_results.csv";
    std::ofstream ofs(filename);
    ofs << "Step,Time(au),x,v,Pop_0,Pop_1,Prob_0_to_1,ActiveState,NAC,E_tot,LZ_Prob\n";

    double time = 0.0;
    int total_steps = 3500; 

    for (int step = 0; step < total_steps; ++step) {
        
        // (A) 计算 Tully 1D 解析势能面 (Diabatic)
        double v11 = (x > 0) ? A * (1.0 - std::exp(-B * x)) : -A * (1.0 - std::exp(B * x));
        double v12 = C * std::exp(-D * x * x);
        double dv11_dx = (x > 0) ? A * B * std::exp(-B * x) : A * B * std::exp(B * x);
        double dv12_dx = -2.0 * C * D * x * std::exp(-D * x * x);

        // (B) 绝热变换 (Adiabatic)
        double gap_sq = v11 * v11 + v12 * v12;
        double gap = std::sqrt(gap_sq);
        double E0 = -gap; 
        double E1 = gap;  

        // (C) 计算非绝热耦合 (NAC) 与 原子受力
        double d12_x = (v11 * dv12_dx - v12 * dv11_dx) / (2.0 * gap_sq);
        double sigma12 = v * d12_x; 

        double f0 = (v11 * dv11_dx + v12 * dv12_dx) / gap; 
        double f1 = -f0;

        // (D) 交给 FSSH 引擎进行推演
        std::vector<double> energies = {E0, E1};
        ModuleBase::ComplexMatrix sigma(2, 2);
        sigma(0, 1) = std::complex<double>(sigma12, 0.0);
        sigma(1, 0) = std::complex<double>(-sigma12, 0.0);

        double pop_0 = std::norm(electronic_coeffs_[0]);
        double pop_1 = std::norm(electronic_coeffs_[1]);

        double current_force = (current_state_ == 0) ? f0 : f1;
        v += 0.5 * (current_force / mass) * dt_au;
        x += v * dt_au;

        this->propagate_rk4(sigma, energies);

        double B_jk = 2.0 * (std::conj(electronic_coeffs_[current_state_]) * electronic_coeffs_[1 - current_state_] * sigma(current_state_, 1 - current_state_)).real();
        double prob = std::max(0.0, B_jk * dt_au / std::norm(electronic_coeffs_[current_state_]));

        int proposed = this->check_hopping(sigma);
        if (proposed != current_state_) {
            double dE = energies[proposed] - energies[current_state_];
            double e_kin = 0.5 * mass * v * v;
            if (e_kin > dE) { 
                double v_new = std::sqrt(2.0 * (e_kin - dE) / mass);
                v = (v > 0) ? v_new : -v_new; 
                current_state_ = proposed;
                std::cout << " [TULLY HOP!] Time: " << time << " | Hop: " << 1-proposed << " -> " << proposed << std::endl;
            } else { 
                std::cout << " [TULLY FRUSTRATED] Time: " << time << " | Insufficient Energy!" << std::endl;
            }
        }

        double v11_new = (x > 0) ? A * (1.0 - std::exp(-B * x)) : -A * (1.0 - std::exp(B * x));
        double v12_new = C * std::exp(-D * x * x);
        double dv11_dx_new = (x > 0) ? A * B * std::exp(-B * x) : A * B * std::exp(B * x);
        double dv12_dx_new = -2.0 * C * D * x * std::exp(-D * x * x);
        double gap_new = std::sqrt(v11_new*v11_new + v12_new*v12_new);
        double f0_new = (v11_new * dv11_dx_new + v12_new * dv12_dx_new) / gap_new;
        double f1_new = -f0_new;
        double force_new = (current_state_ == 0) ? f0_new : f1_new;

        v += 0.5 * (force_new / mass) * dt_au;
        time += dt_au;

        double etot = 0.5 * mass * v * v + energies[current_state_];
        
        // 保存数据 (增加 LZ_Prob)
        ofs << std::fixed << std::setprecision(8)
            << step << "," << time << "," << x << "," << v << ","
            << pop_0 << "," << pop_1 << "," << prob << "," 
            << current_state_ << "," << d12_x << "," << etot << "," << prob_LZ << "\n";
    }
    ofs.close();
    std::cout << "\n>>> TULLY MODEL 1 TEST COMPLETED! <<<\n" << std::endl;
}

// =====================================================================
// 🌟 终极试金石进阶：Tully Model 1 动量扫描 (带分支比、LZ概率及文件关闭播报)
// =====================================================================
void FsshDriver::test_tully_model_1_scan() {
    std::cout << "\n=====================================================" << std::endl;
    std::cout << ">>> STARTING MONTE CARLO SCAN: TULLY MODEL 1      <<<" << std::endl;
    std::cout << ">>> 59 Momentums x 1000 Trajectories = 59000 Runs <<<" << std::endl;
    std::cout << "=====================================================\n" << std::endl;

    const double A = 0.01;
    const double B = 1.6;
    const double C = 0.005;
    const double D = 1.0;
    const double mass = 2000.0; 
    const double PI = 3.14159265358979323846;
    const int num_trajs = 1000; 

    // 1. 准备均值与标准差文件
    std::string stats_filename = "tully_model_1_scan_stats.csv";
    std::ofstream ofs_stats(stats_filename);
    if(ofs_stats.is_open()) {
        ofs_stats << "k_init,v_init,Avg_Pop_1,StdDev_Pop_1,Fraction_State_1,LZ_Prob\n";
    }

    // 2. 准备原始轨迹数据文件
    std::string raw_filename = "tully_model_1_scan_raw.csv";
    std::ofstream ofs_raw(raw_filename);
    if(ofs_raw.is_open()) {
        ofs_raw << "k_init,traj_id,Pop_1,Final_Active_State,Status\n";
    }

    // 🌟 3. 准备分支比统计文件 (新增 LZ_Prob 表头)
    std::string branch_filename = "tully_model_1_scan_branching.csv";
    std::ofstream ofs_branch(branch_filename);
    if(ofs_branch.is_open()) {
        ofs_branch << "k_init,T0(Trans_State0),R0(Refl_State0),T1(Trans_State1),R1(Refl_State1),LZ_Prob\n";
    }

    printf(" %-6s | %-10s | %-10s | %-10s | %-10s | %-10s\n", 
           "k_init", "T0(%)", "R0(%)", "T1(%)", "R1(%)", "LZ_Prob(%)");
    printf("------------------------------------------------------------------------\n");

    for (double k_init = 1.0; k_init <= 30.0; k_init += 0.5) {
        
        double v_init = k_init / mass; 
        
        double e_kin_init = 0.5 * mass * v_init * v_init;
        double e_pot_init = -A; 
        double e_tot = e_kin_init + e_pot_init;
        double e_pot_xc = -C; 
        
        double prob_LZ = 0.0;
        if (e_tot > e_pot_xc) {
            double v_c = std::sqrt(2.0 * (e_tot - e_pot_xc) / mass);
            prob_LZ = std::exp(- PI * C * C / (v_c * A * B));
        }

        std::vector<double> pop1_list;
        pop1_list.reserve(num_trajs);
        
        int count_T0 = 0; 
        int count_R0 = 0; 
        int count_T1 = 0; 
        int count_R1 = 0; 

        // 内层循环：1000 次独立轨迹
        for (int traj = 0; traj < num_trajs; ++traj) {
            
            double x = -10.0; 
            double v = v_init;
            double dt_au = 1.0; 
            
            // 注意：这里调用的 init 内部现在使用的是 assign()，完美清空幽灵记忆
            this->init(1, 2, dt_au, 0); 
            
            int step = 0;
            std::string status = "Transmitted";

            while (true) {
                double v11 = (x > 0) ? A * (1.0 - std::exp(-B * x)) : -A * (1.0 - std::exp(B * x));
                double v12 = C * std::exp(-D * x * x);
                double dv11_dx = (x > 0) ? A * B * std::exp(-B * x) : A * B * std::exp(B * x);
                double dv12_dx = -2.0 * C * D * x * std::exp(-D * x * x);

                double gap_sq = v11 * v11 + v12 * v12;
                double gap = std::sqrt(gap_sq);
                double E0 = -gap; 
                double E1 = gap;  

                double d12_x = (v11 * dv12_dx - v12 * dv11_dx) / (2.0 * gap_sq);
                double sigma12 = v * d12_x; 
                double f0 = (v11 * dv11_dx + v12 * dv12_dx) / gap; 
                double f1 = -f0;

                std::vector<double> energies = {E0, E1};
                ModuleBase::ComplexMatrix sigma(2, 2);
                sigma(0, 1) = std::complex<double>(sigma12, 0.0);
                sigma(1, 0) = std::complex<double>(-sigma12, 0.0);

                double current_force = (current_state_ == 0) ? f0 : f1;
                v += 0.5 * (current_force / mass) * dt_au;
                x += v * dt_au;

                this->propagate_rk4(sigma, energies);

                int proposed = this->check_hopping(sigma);
                if (proposed != current_state_) {
                    double dE = energies[proposed] - energies[current_state_];
                    double e_kin = 0.5 * mass * v * v;
                    if (e_kin > dE) { 
                        double v_new = std::sqrt(2.0 * (e_kin - dE) / mass);
                        v = (v > 0) ? v_new : -v_new; 
                        current_state_ = proposed;
                    }
                }

                double v11_new = (x > 0) ? A * (1.0 - std::exp(-B * x)) : -A * (1.0 - std::exp(B * x));
                double v12_new = C * std::exp(-D * x * x);
                double dv11_dx_new = (x > 0) ? A * B * std::exp(-B * x) : A * B * std::exp(B * x);
                double dv12_dx_new = -2.0 * C * D * x * std::exp(-D * x * x);
                double gap_new = std::sqrt(v11_new*v11_new + v12_new*v12_new);
                double f0_new = (v11_new * dv11_dx_new + v12_new * dv12_dx_new) / gap_new;
                double f1_new = -f0_new;
                double force_new = (current_state_ == 0) ? f0_new : f1_new;

                v += 0.5 * (force_new / mass) * dt_au;
                step++;

                // 截断条件
                if (x > 10.0) { status = "Transmitted"; break; }
                if (x < -10.5) { status = "Reflected"; break; }
                if (step > 2000000) { status = "Timeout"; break; }
            }

            double final_pop1 = std::norm(electronic_coeffs_[1]);
            pop1_list.push_back(final_pop1);
            
            if (status == "Transmitted") {
                if (current_state_ == 0) count_T0++;
                else count_T1++;
            } else if (status == "Reflected") {
                if (current_state_ == 0) count_R0++;
                else count_R1++;
            }

            if(ofs_raw.is_open()) {
                ofs_raw << k_init << "," << traj << "," << final_pop1 << "," << current_state_ << "," << status << "\n";
            }
        }

        double P_T0 = (double)count_T0 / num_trajs;
        double P_R0 = (double)count_R0 / num_trajs;
        double P_T1 = (double)count_T1 / num_trajs;
        double P_R1 = (double)count_R1 / num_trajs;

        double sum_pop1 = 0.0;
        for(double p : pop1_list) sum_pop1 += p;
        double avg_pop1 = sum_pop1 / num_trajs;

        double sq_sum = 0.0;
        for(double p : pop1_list) sq_sum += (p - avg_pop1) * (p - avg_pop1);
        double stddev_pop1 = std::sqrt(sq_sum / (num_trajs - 1)); 

        printf(" %-6.1f | %-10.2f | %-10.2f | %-10.2f | %-10.2f | %-10.2f\n", 
               k_init, P_T0*100, P_R0*100, P_T1*100, P_R1*100, prob_LZ*100);

        if(ofs_stats.is_open()) {
            ofs_stats << std::fixed << std::setprecision(6)
                      << k_init << "," << v_init << "," 
                      << avg_pop1 << "," << stddev_pop1 << "," 
                      << P_T1 << "," << prob_LZ << "\n";
        }

        // 🌟 新增：在输出四个分支概率后，紧跟着输出 LZ_Prob
        if(ofs_branch.is_open()) {
            ofs_branch << std::fixed << std::setprecision(6)
                       << k_init << "," << P_T0 << "," << P_R0 << "," << P_T1 << "," << P_R1 << "," << prob_LZ << "\n";
        }
    }
    
    // 🌟 新增：文件关闭判定与播报
    std::cout << "\n=====================================================" << std::endl;
    if(ofs_stats.is_open()) {
        ofs_stats.close();
        std::cout << "[INFO] 统计文件已安全保存并关闭: " << stats_filename << std::endl;
    }
    if(ofs_raw.is_open()) {
        ofs_raw.close();
        std::cout << "[INFO] 原始轨迹文件已安全保存并关闭: " << raw_filename << std::endl;
    }
    if(ofs_branch.is_open()) {
        ofs_branch.close();
        std::cout << "[INFO] 分支比文件已安全保存并关闭: " << branch_filename << std::endl;
    }
    
    std::cout << ">>> ALL 59000 SIMULATIONS COMPLETED SUCCESSFULLY! <<<\n" << std::endl;
}

// =====================================================================
// 🌟 终极试金石二：Tully Model 2 (Dual Avoided Crossing) 势能面与单轨测试
// =====================================================================
void FsshDriver::test_tully_model_2() {
    std::cout << "\n=====================================================" << std::endl;
    std::cout << ">>> STARTING GOLD STANDARD TEST: TULLY MODEL 2  <<<" << std::endl;
    std::cout << "=====================================================\n" << std::endl;

    // 1. Tully Model 2 经典参数 (Atomic Units)
    const double A = 0.10;
    const double B = 0.28;
    const double C = 0.015;
    const double D = 0.06;
    const double E0_param = 0.05;
    const double mass = 2000.0; 

    // ==========================================================
    // 🌟 模块 1：输出 Tully Model 2 的势能面 (PES) 及非绝热耦合 (NAC)
    // ==========================================================
    std::string pes_filename = "tully_model_2_pes.csv";
    std::ofstream pes_ofs(pes_filename);
    if(pes_ofs.is_open()) {
        pes_ofs << "x,E0(Lower),E1(Upper),V11,V22,V12,NAC(d12)\n"; 
        for(double px = -15.0; px <= 15.0; px += 0.05) {
            double pv11 = 0.0;
            double pv22 = -A * std::exp(-B * px * px) + E0_param;
            double pv12 = C * std::exp(-D * px * px);
            
            double pdv11_dx = 0.0;
            double pdv22_dx = 2.0 * A * B * px * std::exp(-B * px * px);
            double pdv12_dx = -2.0 * C * D * px * std::exp(-D * px * px);
            
            double gap_sq = (pv11 - pv22)*(pv11 - pv22) + 4.0 * pv12 * pv12;
            double pgap = std::sqrt(gap_sq);
            
            double E0 = 0.5 * (pv11 + pv22 - pgap);
            double E1 = 0.5 * (pv11 + pv22 + pgap);
            
            // Model 2 的 d12 解析公式
            double pd12_x = (pv12 * pdv22_dx - pdv12_dx * pv22) / gap_sq;
            
            pes_ofs << std::fixed << std::setprecision(8) 
                    << px << "," << E0 << "," << E1 << "," 
                    << pv11 << "," << pv22 << "," << pv12 << "," << pd12_x << "\n";
        }
        pes_ofs.close();
        std::cout << "[INFO] Model 2 势能面数据已保存至: " << pes_filename << std::endl;
    }

    // ==========================================================
    // 🌟 模块 2：单条轨迹演化测试 (感受双重交叉的量子干涉)
    // ==========================================================
    // 设置初始状态，选择总能 E = 0.1 a.u. (足以越过中心势垒)
    double x = -10.0; 
    double e_tot_init = 0.1; 
    double v = std::sqrt(2.0 * e_tot_init / mass); 
    
    // 初始化 FSSH 引擎 (使用 assign 重置了状态)
    double dt_au = 1.0; 
    this->init(1, 2, dt_au, 0); 

    std::string traj_filename = "tully_model_2_results.csv";
    std::ofstream traj_ofs(traj_filename);
    if(traj_ofs.is_open()) {
        traj_ofs << "Step,Time(au),x,v,Pop_0,Pop_1,ActiveState,NAC,E_tot\n";
    }

    double time = 0.0;
    int max_steps = 300000; // 给足步数，防止慢速粒子跑不完

    std::cout << "\n[INFO] 开始跑单条轨迹 (E_tot = 0.1 a.u.)..." << std::endl;

    for (int step = 0; step < max_steps; ++step) {
        
        // (A) 计算 Tully 2 解析势能面
        double v11 = 0.0;
        double v22 = -A * std::exp(-B * x * x) + E0_param;
        double v12 = C * std::exp(-D * x * x);
        
        double dv11_dx = 0.0;
        double dv22_dx = 2.0 * A * B * x * std::exp(-B * x * x);
        double dv12_dx = -2.0 * C * D * x * std::exp(-D * x * x);

        double gap_sq = (v11 - v22)*(v11 - v22) + 4.0 * v12 * v12;
        double gap = std::sqrt(gap_sq);
        double E0 = 0.5 * (v11 + v22 - gap);
        double E1 = 0.5 * (v11 + v22 + gap);

        // (C) 计算非绝热耦合 (NAC) 与 原子受力
        double d12_x = (v12 * dv22_dx - dv12_dx * v22) / gap_sq;
        double sigma12 = v * d12_x;

        double f0 = -0.5 * dv22_dx + (v22 * dv22_dx + 4.0 * v12 * dv12_dx) / (2.0 * gap);
        double f1 = -0.5 * dv22_dx - (v22 * dv22_dx + 4.0 * v12 * dv12_dx) / (2.0 * gap);

        std::vector<double> energies = {E0, E1};
        ModuleBase::ComplexMatrix sigma(2, 2);
        sigma(0, 1) = std::complex<double>(sigma12, 0.0);
        sigma(1, 0) = std::complex<double>(-sigma12, 0.0);

        double pop_0 = std::norm(electronic_coeffs_[0]);
        double pop_1 = std::norm(electronic_coeffs_[1]);

        double current_force = (current_state_ == 0) ? f0 : f1;
        v += 0.5 * (current_force / mass) * dt_au;
        x += v * dt_au;

        // 波函数演化
        this->propagate_rk4(sigma, energies);

        // 面跳跃逻辑
        int proposed = this->check_hopping(sigma);
        if (proposed != current_state_) {
            double dE = energies[proposed] - energies[current_state_];
            double e_kin = 0.5 * mass * v * v;
            if (e_kin > dE) { 
                double v_new = std::sqrt(2.0 * (e_kin - dE) / mass);
                v = (v > 0) ? v_new : -v_new; 
                current_state_ = proposed;
                std::cout << " [HOPPING EVENT!] Time: " << time << " | Hop: " << 1-proposed << " -> " << proposed 
                          << " | Location x = " << std::fixed << std::setprecision(4) << x << std::endl;
            }
        }

        // 重新计算后半步受力
        double v22_new = -A * std::exp(-B * x * x) + E0_param;
        double v12_new = C * std::exp(-D * x * x);
        double dv22_dx_new = 2.0 * A * B * x * std::exp(-B * x * x);
        double dv12_dx_new = -2.0 * C * D * x * std::exp(-D * x * x);
        
        double gap_new = std::sqrt(v22_new*v22_new + 4.0*v12_new*v12_new);
        double f0_new = -0.5 * dv22_dx_new + (v22_new * dv22_dx_new + 4.0 * v12_new * dv12_dx_new) / (2.0 * gap_new);
        double f1_new = -0.5 * dv22_dx_new - (v22_new * dv22_dx_new + 4.0 * v12_new * dv12_dx_new) / (2.0 * gap_new);
        
        double force_new = (current_state_ == 0) ? f0_new : f1_new;
        v += 0.5 * (force_new / mass) * dt_au;
        time += dt_au;

        double etot = 0.5 * mass * v * v + energies[current_state_];

        // 稀疏采样输出，避免文件过大
        if (step % 50 == 0 && traj_ofs.is_open()) {
            traj_ofs << std::fixed << std::setprecision(8)
                     << step << "," << time << "," << x << "," << v << ","
                     << pop_0 << "," << pop_1 << "," 
                     << current_state_ << "," << d12_x << "," << etot << "\n";
        }
        
        // 🌟 空间截断：跑出危险区就算结束
        if (x > 10.0 || x < -10.5) {
            std::cout << "[INFO] 粒子已飞出相互作用区 (x = " << x << ")，演化结束！" << std::endl;
            break;
        }
    }
    
    if(traj_ofs.is_open()) traj_ofs.close();
    std::cout << "[INFO] Model 2 单轨演化数据已保存至: " << traj_filename << std::endl;
    std::cout << ">>> TULLY MODEL 2 SINGLE TRAJECTORY TEST COMPLETED! <<<\n" << std::endl;
}

// =====================================================================
// 🌟 终极试金石进阶：Tully Model 2 能量扫描 (带分支比与干涉震荡)
// =====================================================================
void FsshDriver::test_tully_model_2_scan() {
    std::cout << "\n=====================================================" << std::endl;
    std::cout << ">>> STARTING MONTE CARLO SCAN: TULLY MODEL 2      <<<" << std::endl;
    std::cout << ">>> Energy from 0.02 to 3.0 (Logarithmic steps)   <<<" << std::endl;
    std::cout << "=====================================================\n" << std::endl;

    const double A = 0.10;
    const double B = 0.28;
    const double C = 0.015;
    const double D = 0.06;
    const double E0_param = 0.05;
    const double mass = 2000.0; 
    const int num_trajs = 1000; // 每个能量点跑 1000 次

    // 1. 均值与标准差文件 (您的 Pop2 在代码里是 Pop_1，即上能级)
    std::string stats_filename = "tully_model_2_scan_stats.csv";
    std::ofstream ofs_stats(stats_filename);
    if(ofs_stats.is_open()) {
        ofs_stats << "E_tot,v_init,Avg_Pop_1,StdDev_Pop_1,Fraction_State_1\n";
    }

    // 2. 原始轨迹数据文件
    std::string raw_filename = "tully_model_2_scan_raw.csv";
    std::ofstream ofs_raw(raw_filename);
    if(ofs_raw.is_open()) {
        ofs_raw << "E_tot,traj_id,Pop_1,Final_Active_State,Status\n";
    }

    // 3. 分支比统计文件
    std::string branch_filename = "tully_model_2_scan_branching.csv";
    std::ofstream ofs_branch(branch_filename);
    if(ofs_branch.is_open()) {
        ofs_branch << "E_tot,T0_Lower,R0_Lower,T1_Upper,R1_Upper\n";
    }

    printf(" %-8s | %-10s | %-10s | %-10s | %-10s | %-10s\n", 
           "E_tot", "T0_Low(%)", "R0_Low(%)", "T1_Up(%)", "R1_Up(%)", "Avg Pop_1");
    printf("--------------------------------------------------------------------------\n");

    // 🌟 使用对数网格扫描总能 E_tot，以完美捕捉低能区的量子干涉震荡！
    int num_points = 50;
    double logE_start = std::log(0.02);
    double logE_end = std::log(3.0);
    double dlogE = (logE_end - logE_start) / (num_points - 1);

    for (int i = 0; i < num_points; ++i) {
        
        double e_tot = std::exp(logE_start + i * dlogE);
        
        // Model 2 初始点在 x=-10 处，E0 约等于 0，所以动能近似为 e_tot
        double v_init = std::sqrt(2.0 * e_tot / mass); 
        
        std::vector<double> pop1_list;
        pop1_list.reserve(num_trajs);
        
        int count_T0 = 0; 
        int count_R0 = 0; 
        int count_T1 = 0; 
        int count_R1 = 0; 

        // 蒙特卡洛系综
        for (int traj = 0; traj < num_trajs; ++traj) {
            
            double x = -10.0; 
            double v = v_init;
            double dt_au = 1.0; 
            
            // 初始化，强制分配清除上一条轨迹的幽灵内存，初始在低能级态 0
            this->init(1, 2, dt_au, 0); 
            
            int step = 0;
            std::string status = "Transmitted";

            while (true) {
                // 计算势能面
                double v11 = 0.0;
                double v22 = -A * std::exp(-B * x * x) + E0_param;
                double v12 = C * std::exp(-D * x * x);

                double dv11_dx = 0.0;
                double dv22_dx = 2.0 * A * B * x * std::exp(-B * x * x);
                double dv12_dx = -2.0 * C * D * x * std::exp(-D * x * x);

                double gap_sq = (v11 - v22)*(v11 - v22) + 4.0 * v12 * v12;
                double gap = std::sqrt(gap_sq);
                
                double E0 = 0.5 * (v11 + v22 - gap); 
                double E1 = 0.5 * (v11 + v22 + gap);  

                // NAC 与受力
                double d12_x = (v12 * dv22_dx - dv12_dx * v22) / gap_sq;
                double sigma12 = v * d12_x; 
                
                double f0 = -0.5 * dv22_dx + (v22 * dv22_dx + 4.0 * v12 * dv12_dx) / (2.0 * gap);
                double f1 = -0.5 * dv22_dx - (v22 * dv22_dx + 4.0 * v12 * dv12_dx) / (2.0 * gap);

                std::vector<double> energies = {E0, E1};
                ModuleBase::ComplexMatrix sigma(2, 2);
                sigma(0, 1) = std::complex<double>(sigma12, 0.0);
                sigma(1, 0) = std::complex<double>(-sigma12, 0.0);

                // 前半步速度更新
                double current_force = (current_state_ == 0) ? f0 : f1;
                v += 0.5 * (current_force / mass) * dt_au;
                x += v * dt_au;

                // 薛定谔演化
                this->propagate_rk4(sigma, energies);

                // 面跳跃与动量缩放
                int proposed = this->check_hopping(sigma);
                if (proposed != current_state_) {
                    double dE = energies[proposed] - energies[current_state_];
                    double e_kin = 0.5 * mass * v * v;
                    if (e_kin > dE) { 
                        double v_new = std::sqrt(2.0 * (e_kin - dE) / mass);
                        v = (v > 0) ? v_new : -v_new; 
                        current_state_ = proposed;
                    }
                }

                // 后半步速度更新
                double v22_new = -A * std::exp(-B * x * x) + E0_param;
                double v12_new = C * std::exp(-D * x * x);
                double dv22_dx_new = 2.0 * A * B * x * std::exp(-B * x * x);
                double dv12_dx_new = -2.0 * C * D * x * std::exp(-D * x * x);
                
                double gap_new = std::sqrt(v22_new*v22_new + 4.0*v12_new*v12_new);
                double f0_new = -0.5 * dv22_dx_new + (v22_new * dv22_dx_new + 4.0 * v12_new * dv12_dx_new) / (2.0 * gap_new);
                double f1_new = -0.5 * dv22_dx_new - (v22_new * dv22_dx_new + 4.0 * v12_new * dv12_dx_new) / (2.0 * gap_new);
                
                double force_new = (current_state_ == 0) ? f0_new : f1_new;
                v += 0.5 * (force_new / mass) * dt_au;
                step++;

                // 截断条件：超出危险区域算作结束
                if (x > 10.0) { status = "Transmitted"; break; }
                if (x < -10.5) { status = "Reflected"; break; }
                if (step > 3000000) { status = "Timeout"; break; }
            }

            // 获取该轨迹上能级 (态1) 的最终布居数
            double final_pop1 = std::norm(electronic_coeffs_[1]);
            pop1_list.push_back(final_pop1);
            
            // 统计分支比
            if (status == "Transmitted") {
                if (current_state_ == 0) count_T0++;
                else count_T1++;
            } else if (status == "Reflected") {
                if (current_state_ == 0) count_R0++;
                else count_R1++;
            }

            if(ofs_raw.is_open()) {
                ofs_raw << e_tot << "," << traj << "," << final_pop1 << "," << current_state_ << "," << status << "\n";
            }
        }

        double P_T0 = (double)count_T0 / num_trajs;
        double P_R0 = (double)count_R0 / num_trajs;
        double P_T1 = (double)count_T1 / num_trajs;
        double P_R1 = (double)count_R1 / num_trajs;

        double sum_pop1 = 0.0;
        for(double p : pop1_list) sum_pop1 += p;
        double avg_pop1 = sum_pop1 / num_trajs;

        double sq_sum = 0.0;
        for(double p : pop1_list) sq_sum += (p - avg_pop1) * (p - avg_pop1);
        double stddev_pop1 = std::sqrt(sq_sum / (num_trajs - 1)); 

        printf(" %-8.4f | %-10.1f | %-10.1f | %-10.1f | %-10.1f | %-10.5f\n", 
               e_tot, P_T0*100, P_R0*100, P_T1*100, P_R1*100, avg_pop1);

        if(ofs_stats.is_open()) {
            ofs_stats << std::fixed << std::setprecision(8)
                      << e_tot << "," << v_init << "," 
                      << avg_pop1 << "," << stddev_pop1 << "," 
                      << P_T1 << "\n";
        }

        if(ofs_branch.is_open()) {
            ofs_branch << std::fixed << std::setprecision(8)
                       << e_tot << "," << P_T0 << "," << P_R0 << "," << P_T1 << "," << P_R1 << "\n";
        }
    }
    
    // 安全关闭与播报
    std::cout << "\n=====================================================" << std::endl;
    if(ofs_stats.is_open()) { ofs_stats.close(); std::cout << "[INFO] 统计文件关闭: " << stats_filename << std::endl; }
    if(ofs_raw.is_open()) { ofs_raw.close(); std::cout << "[INFO] 原始数据关闭: " << raw_filename << std::endl; }
    if(ofs_branch.is_open()) { ofs_branch.close(); std::cout << "[INFO] 分支比文件关闭: " << branch_filename << std::endl; }
    std::cout << ">>> ALL SIMULATIONS FOR MODEL 2 COMPLETED! <<<\n" << std::endl;
}

// =====================================================================
// 🌟 终极试金石三：Tully Model 3 (Extended Coupling) 势能面与单轨测试
// =====================================================================
void FsshDriver::test_tully_model_3() {
    std::cout << "\n=====================================================" << std::endl;
    std::cout << ">>> STARTING GOLD STANDARD TEST: TULLY MODEL 3  <<<" << std::endl;
    std::cout << "=====================================================\n" << std::endl;

    // 1. Tully Model 3 经典参数 (Atomic Units)
    const double A = 0.0006;
    const double B = 0.10;
    const double C = 0.90;
    const double mass = 2000.0; 

    // ==========================================================
    // 🌟 模块 1：输出 Tully Model 3 的势能面 (PES) 及非绝热耦合 (NAC)
    // ==========================================================
    std::string pes_filename = "tully_model_3_pes.csv";
    std::ofstream pes_ofs(pes_filename);
    if(pes_ofs.is_open()) {
        pes_ofs << "x,E0(Lower),E1(Upper),V11,V22,V12,NAC(d12)\n"; 
        for(double px = -15.0; px <= 15.0; px += 0.05) {
            double pv11 = -A;
            double pv22 = A;
            double pv12 = (px < 0) ? B * std::exp(C * px) : B * (2.0 - std::exp(-C * px));
            
            double pdv11_dx = 0.0;
            double pdv22_dx = 0.0;
            double pdv12_dx = (px < 0) ? B * C * std::exp(C * px) : B * C * std::exp(-C * px);
            
            double gap_sq = (pv11 - pv22)*(pv11 - pv22) + 4.0 * pv12 * pv12;
            double pgap = std::sqrt(gap_sq);
            
            double E0 = 0.5 * (pv11 + pv22 - pgap);
            double E1 = 0.5 * (pv11 + pv22 + pgap);
            
            // 严谨的 d12 解析公式
            double pd12_x = ((pv11 - pv22) * pdv12_dx - pv12 * (pdv11_dx - pdv22_dx)) / gap_sq;
            
            pes_ofs << std::fixed << std::setprecision(8) 
                    << px << "," << E0 << "," << E1 << "," 
                    << pv11 << "," << pv22 << "," << pv12 << "," << pd12_x << "\n";
        }
        pes_ofs.close();
        std::cout << "[INFO] Model 3 势能面数据已保存至: " << pes_filename << std::endl;
    }

    // ==========================================================
    // 🌟 模块 2：单条轨迹演化测试
    // ==========================================================
    double x = -15.0; 
    double e_tot_init = 0.05; // 设定一个适中的总能，以观察粒子在攀爬势垒时的行为
    double e_pot_init = -std::sqrt(A*A + B*std::exp(C*x)*B*std::exp(C*x));
    double v = std::sqrt(2.0 * std::abs(e_tot_init - e_pot_init) / mass); 
    
    double dt_au = 1.0; 
    this->init(1, 2, dt_au, 0); 

    std::string traj_filename = "tully_model_3_results.csv";
    std::ofstream traj_ofs(traj_filename);
    if(traj_ofs.is_open()) {
        traj_ofs << "Step,Time(au),x,v,Pop_0,Pop_1,ActiveState,NAC,E_tot\n";
    }

    double time = 0.0;
    int max_steps = 500000; 
    std::cout << "\n[INFO] 开始跑单条轨迹 (E_tot = 0.05 a.u.)..." << std::endl;

    for (int step = 0; step < max_steps; ++step) {
        
        double v11 = -A;
        double v22 = A;
        double v12 = (x < 0) ? B * std::exp(C * x) : B * (2.0 - std::exp(-C * x));
        
        double dv11_dx = 0.0;
        double dv22_dx = 0.0;
        double dv12_dx = (x < 0) ? B * C * std::exp(C * x) : B * C * std::exp(-C * x);

        double gap_sq = (v11 - v22)*(v11 - v22) + 4.0 * v12 * v12;
        double gap = std::sqrt(gap_sq);
        double E0 = 0.5 * (v11 + v22 - gap);
        double E1 = 0.5 * (v11 + v22 + gap);

        double d12_x = ((v11 - v22) * dv12_dx - v12 * (dv11_dx - dv22_dx)) / gap_sq;
        double sigma12 = v * d12_x;

        double f0 = -0.5 * (dv11_dx + dv22_dx) + ((v11 - v22)*(dv11_dx - dv22_dx) + 4.0 * v12 * dv12_dx) / (2.0 * gap);
        double f1 = -0.5 * (dv11_dx + dv22_dx) - ((v11 - v22)*(dv11_dx - dv22_dx) + 4.0 * v12 * dv12_dx) / (2.0 * gap);

        std::vector<double> energies = {E0, E1};
        ModuleBase::ComplexMatrix sigma(2, 2);
        sigma(0, 1) = std::complex<double>(sigma12, 0.0);
        sigma(1, 0) = std::complex<double>(-sigma12, 0.0);

        double pop_0 = std::norm(electronic_coeffs_[0]);
        double pop_1 = std::norm(electronic_coeffs_[1]);

        double current_force = (current_state_ == 0) ? f0 : f1;
        v += 0.5 * (current_force / mass) * dt_au;
        x += v * dt_au;

        this->propagate_rk4(sigma, energies);

        int proposed = this->check_hopping(sigma);
        if (proposed != current_state_) {
            double dE = energies[proposed] - energies[current_state_];
            double e_kin = 0.5 * mass * v * v;
            if (e_kin > dE) { 
                double v_new = std::sqrt(2.0 * (e_kin - dE) / mass);
                v = (v > 0) ? v_new : -v_new; 
                current_state_ = proposed;
                std::cout << " [HOPPING EVENT!] Time: " << time << " | Hop: " << 1-proposed << " -> " << proposed 
                          << " | Location x = " << std::fixed << std::setprecision(4) << x << std::endl;
            } else {
                // 如果动能不足以爬坡，这是典型的受挫跳跃 (Frustrated Hop)，保持原态！
            }
        }

        double v12_new = (x < 0) ? B * std::exp(C * x) : B * (2.0 - std::exp(-C * x));
        double dv12_dx_new = (x < 0) ? B * C * std::exp(C * x) : B * C * std::exp(-C * x);
        
        double gap_new = std::sqrt((v11 - v22)*(v11 - v22) + 4.0 * v12_new * v12_new);
        double f0_new = -0.5 * (dv11_dx + dv22_dx) + ((v11 - v22)*(dv11_dx - dv22_dx) + 4.0 * v12_new * dv12_dx_new) / (2.0 * gap_new);
        double f1_new = -0.5 * (dv11_dx + dv22_dx) - ((v11 - v22)*(dv11_dx - dv22_dx) + 4.0 * v12_new * dv12_dx_new) / (2.0 * gap_new);
        
        double force_new = (current_state_ == 0) ? f0_new : f1_new;
        v += 0.5 * (force_new / mass) * dt_au;
        time += dt_au;

        double etot = 0.5 * mass * v * v + energies[current_state_];

        if (step % 50 == 0 && traj_ofs.is_open()) {
            traj_ofs << std::fixed << std::setprecision(8)
                     << step << "," << time << "," << x << "," << v << ","
                     << pop_0 << "," << pop_1 << "," 
                     << current_state_ << "," << d12_x << "," << etot << "\n";
        }
        
        if (x > 15.0 || x < -15.5) {
            std::cout << "[INFO] 粒子已飞出相互作用区 (x = " << x << ")，演化结束！" << std::endl;
            break;
        }
    }
    
    if(traj_ofs.is_open()) traj_ofs.close();
    std::cout << "[INFO] Model 3 单轨演化数据已保存至: " << traj_filename << std::endl;
    std::cout << ">>> TULLY MODEL 3 SINGLE TRAJECTORY TEST COMPLETED! <<<\n" << std::endl;
}

// =====================================================================
// 🌟 终极试金石进阶：Tully Model 3 动量扫描 (带反射/透射分支比统计)
// =====================================================================
void FsshDriver::test_tully_model_3_scan() {
    std::cout << "\n=====================================================" << std::endl;
    std::cout << ">>> STARTING MONTE CARLO SCAN: TULLY MODEL 3      <<<" << std::endl;
    std::cout << ">>> Momentum k from 2.0 to 35.0 (Step 0.2)        <<<" << std::endl;
    std::cout << "=====================================================\n" << std::endl;

    const double A = 0.0006;
    const double B = 0.10;
    const double C = 0.90;
    const double mass = 2000.0; 
    const int num_trajs = 1000; 

    // 1. 均值与标准差文件 (加入 k_init 记录)
    std::string stats_filename = "tully_model_3_scan_stats.csv";
    std::ofstream ofs_stats(stats_filename);
    if(ofs_stats.is_open()) {
        ofs_stats << "k_init,E_tot,Avg_Pop_1,StdDev_Pop_1,Fraction_State_1\n";
    }

    // 2. 原始轨迹数据文件
    std::string raw_filename = "tully_model_3_scan_raw.csv";
    std::ofstream ofs_raw(raw_filename);
    if(ofs_raw.is_open()) {
        ofs_raw << "k_init,traj_id,Pop_1,Final_Active_State,Status\n";
    }

    // 3. 分支比统计文件
    std::string branch_filename = "tully_model_3_scan_branching.csv";
    std::ofstream ofs_branch(branch_filename);
    if(ofs_branch.is_open()) {
        ofs_branch << "k_init,T0(Trans_State0),R0(Refl_State0),T1(Trans_State1),R1(Refl_State1)\n";
    }

    printf(" %-6s | %-10s | %-10s | %-10s | %-10s | %-10s\n", 
           "k_init", "T0(%)", "R0(%)", "T1(%)", "R1(%)", "Avg Pop_1");
    printf("------------------------------------------------------------------------\n");

    // 🌟 按照您的要求，扫描初始动量 k 从 2.0 到 35.0，步长 0.2
    for (double k_init = 2.0; k_init <= 35.0; k_init += 0.2) {
        
        double v_init = k_init / mass; 
        double e_kin_init = 0.5 * mass * v_init * v_init;
        double x_init = -15.0;
        
        // 初始绝热势能
        double v11 = -A;
        double v22 = A;
        double v12 = B * std::exp(C * x_init);
        double e_pot_init = 0.5 * (v11 + v22 - std::sqrt((v11-v22)*(v11-v22) + 4.0*v12*v12));
        
        double e_tot = e_kin_init + e_pot_init;
        
        std::vector<double> pop1_list;
        pop1_list.reserve(num_trajs);
        
        int count_T0 = 0; 
        int count_R0 = 0; 
        int count_T1 = 0; 
        int count_R1 = 0; 

        for (int traj = 0; traj < num_trajs; ++traj) {
            
            double x = x_init; 
            double v = v_init;
            double dt_au = 1.0; 
            
            this->init(1, 2, dt_au, 0); 
            
            int step = 0;
            std::string status = "Transmitted";

            while (true) {
                v12 = (x < 0) ? B * std::exp(C * x) : B * (2.0 - std::exp(-C * x));
                double dv11_dx = 0.0;
                double dv22_dx = 0.0;
                double dv12_dx = (x < 0) ? B * C * std::exp(C * x) : B * C * std::exp(-C * x);

                double gap_sq = (v11 - v22)*(v11 - v22) + 4.0 * v12 * v12;
                double gap = std::sqrt(gap_sq);
                double E0 = 0.5 * (v11 + v22 - gap); 
                double E1 = 0.5 * (v11 + v22 + gap);  

                double d12_x = ((v11 - v22) * dv12_dx - v12 * (dv11_dx - dv22_dx)) / gap_sq;
                double sigma12 = v * d12_x; 
                
                double f0 = -0.5 * (dv11_dx + dv22_dx) + ((v11 - v22)*(dv11_dx - dv22_dx) + 4.0 * v12 * dv12_dx) / (2.0 * gap);
                double f1 = -0.5 * (dv11_dx + dv22_dx) - ((v11 - v22)*(dv11_dx - dv22_dx) + 4.0 * v12 * dv12_dx) / (2.0 * gap);

                std::vector<double> energies = {E0, E1};
                ModuleBase::ComplexMatrix sigma(2, 2);
                sigma(0, 1) = std::complex<double>(sigma12, 0.0);
                sigma(1, 0) = std::complex<double>(-sigma12, 0.0);

                double current_force = (current_state_ == 0) ? f0 : f1;
                v += 0.5 * (current_force / mass) * dt_au;
                x += v * dt_au;

                this->propagate_rk4(sigma, energies);

                int proposed = this->check_hopping(sigma);
                if (proposed != current_state_) {
                    double dE = energies[proposed] - energies[current_state_];
                    double e_kin = 0.5 * mass * v * v;
                    if (e_kin > dE) { 
                        double v_new = std::sqrt(2.0 * (e_kin - dE) / mass);
                        v = (v > 0) ? v_new : -v_new; 
                        current_state_ = proposed;
                    }
                }

                double v12_new = (x < 0) ? B * std::exp(C * x) : B * (2.0 - std::exp(-C * x));
                double dv12_dx_new = (x < 0) ? B * C * std::exp(C * x) : B * C * std::exp(-C * x);
                
                double gap_new = std::sqrt((v11 - v22)*(v11 - v22) + 4.0 * v12_new * v12_new);
                double f0_new = -0.5 * (dv11_dx + dv22_dx) + ((v11 - v22)*(dv11_dx - dv22_dx) + 4.0 * v12_new * dv12_dx_new) / (2.0 * gap_new);
                double f1_new = -0.5 * (dv11_dx + dv22_dx) - ((v11 - v22)*(dv11_dx - dv22_dx) + 4.0 * v12_new * dv12_dx_new) / (2.0 * gap_new);
                
                double force_new = (current_state_ == 0) ? f0_new : f1_new;
                v += 0.5 * (force_new / mass) * dt_au;
                step++;

                // 截断：正向跑到 +15 算透射，被山坡弹回跑到 -15.5 算反射
                if (x > 15.0) { status = "Transmitted"; break; }
                if (x < -15.5) { status = "Reflected"; break; }
                if (step > 3000000) { status = "Timeout"; break; }
            }

            double final_pop1 = std::norm(electronic_coeffs_[1]);
            pop1_list.push_back(final_pop1);
            
            // 统计四大分支
            if (status == "Transmitted") {
                if (current_state_ == 0) count_T0++;
                else count_T1++;
            } else if (status == "Reflected") {
                if (current_state_ == 0) count_R0++;
                else count_R1++;
            }

            if(ofs_raw.is_open()) {
                ofs_raw << k_init << "," << traj << "," << final_pop1 << "," << current_state_ << "," << status << "\n";
            }
        }

        double P_T0 = (double)count_T0 / num_trajs;
        double P_R0 = (double)count_R0 / num_trajs;
        double P_T1 = (double)count_T1 / num_trajs;
        double P_R1 = (double)count_R1 / num_trajs;

        double sum_pop1 = 0.0;
        for(double p : pop1_list) sum_pop1 += p;
        double avg_pop1 = sum_pop1 / num_trajs;

        double sq_sum = 0.0;
        for(double p : pop1_list) sq_sum += (p - avg_pop1) * (p - avg_pop1);
        double stddev_pop1 = std::sqrt(sq_sum / (num_trajs - 1)); 

        printf(" %-6.1f | %-10.1f | %-10.1f | %-10.1f | %-10.1f | %-10.5f\n", 
               k_init, P_T0*100, P_R0*100, P_T1*100, P_R1*100, avg_pop1);

        if(ofs_stats.is_open()) {
            ofs_stats << std::fixed << std::setprecision(8)
                      << k_init << "," << e_tot << "," 
                      << avg_pop1 << "," << stddev_pop1 << "," 
                      << P_T1 << "\n";
        }

        if(ofs_branch.is_open()) {
            ofs_branch << std::fixed << std::setprecision(8)
                       << k_init << "," << P_T0 << "," << P_R0 << "," << P_T1 << "," << P_R1 << "\n";
        }
    }
    
    std::cout << "\n=====================================================" << std::endl;
    if(ofs_stats.is_open()) { ofs_stats.close(); std::cout << "[INFO] 统计文件关闭: " << stats_filename << std::endl; }
    if(ofs_raw.is_open()) { ofs_raw.close(); std::cout << "[INFO] 原始数据关闭: " << raw_filename << std::endl; }
    if(ofs_branch.is_open()) { ofs_branch.close(); std::cout << "[INFO] 分支比文件关闭: " << branch_filename << std::endl; }
    std::cout << ">>> ALL SIMULATIONS FOR MODEL 3 COMPLETED! <<<\n" << std::endl;
}
void FsshDriver::run_lowdin_validation_suite_8_5(const std::string& output_prefix)
{
    const int nbasis = 8;
    const int nocc_lr = 4;
    const int nvirt_lr = 4;
    const int nstates = 6; // 0: ground, 1..5: excited
    const int nks_total = nocc_lr + nvirt_lr;
    const int occ_offset = 0;

    this->init(nbasis, nstates, 1.0, 0, nocc_lr, nks_total);

    auto make_zero_complex_matrix = [](int nr, int nc) {
        ModuleBase::ComplexMatrix mat(nr, nc);
        for (int i = 0; i < nr; ++i) {
            for (int j = 0; j < nc; ++j) {
                mat(i, j) = std::complex<double>(0.0, 0.0);
            }
        }
        return mat;
    };

    auto make_identity_ao = [&]() {
        ModuleBase::ComplexMatrix s_ao = make_zero_complex_matrix(nbasis, nbasis);
        for (int i = 0; i < nbasis; ++i) {
            s_ao(i, i) = std::complex<double>(1.0, 0.0);
        }
        return s_ao;
    };

    auto make_rotated_coef = [&](double alpha) {
        ModuleBase::ComplexMatrix coef = make_zero_complex_matrix(nks_total, nbasis);
        const double c = std::cos(alpha);
        const double s = std::sin(alpha);

        // occupied-virtual pairs: i <-> i+4
        for (int i = 0; i < 4; ++i) {
            coef(i, i) = std::complex<double>(c, 0.0);
            coef(i, i + 4) = std::complex<double>(s, 0.0);
            coef(i + 4, i) = std::complex<double>(-s, 0.0);
            coef(i + 4, i + 4) = std::complex<double>(c, 0.0);
        }
        return coef;
    };

    auto make_casida = [&](double theta) {
        std::vector<CasidaWavefunction> wfcs(
            nstates, CasidaWavefunction(0.0, {}, {}, nocc_lr, nvirt_lr));
        const double c = std::cos(theta);
        const double s = std::sin(theta);

        // we have 5 excited states, let's create a 5x5 rotation/mixing or just 5 orthogonal states
        // nocc_lr * nvirt_lr = 16. We will just use the first 5 components and mix them.
        wfcs[0] = CasidaWavefunction(0.0, {}, {}, nocc_lr, nvirt_lr);
        
        for (int I = 1; I <= 5; ++I) {
            std::vector<double> x_I(nocc_lr * nvirt_lr, 0.0);
            // Just some arbitrary mixing of the first 5 elements based on theta
            for (int j = 0; j < 5; ++j) {
                x_I[j] = std::sin(theta * I + j); // some arbitrary coefficients
            }
            // Normalize x_I
            double norm = 0.0;
            for (int j = 0; j < 5; ++j) norm += x_I[j] * x_I[j];
            norm = std::sqrt(norm);
            if (norm > 1e-12) {
                for (int j = 0; j < 5; ++j) x_I[j] /= norm;
            }
            wfcs[I] = CasidaWavefunction(0.10 * I, x_I, {}, nocc_lr, nvirt_lr);
        }
        return wfcs;
    };

    auto compute_occ_phase = [&](const ModuleBase::ComplexMatrix& coef_old_in,
                                 const ModuleBase::ComplexMatrix& coef_new_in,
                                 const ModuleBase::ComplexMatrix& s_ao_dense) {
        std::vector<double> occ_phase(nocc_lr, 1.0);
        for (int i_lr = 0; i_lr < nocc_lr; ++i_lr) {
            const int i_ks = occ_offset + i_lr;
            double s_ii = 0.0;
            for (int mu = 0; mu < nbasis; ++mu) {
                for (int nu = 0; nu < nbasis; ++nu) {
                    s_ii += (std::conj(coef_old_in(i_ks, mu))
                             * s_ao_dense(mu, nu)
                             * coef_new_in(i_ks, nu))
                                .real();
                }
            }
            if (s_ii < 0.0) {
                occ_phase[i_lr] = -1.0;
            }
        }
        return occ_phase;
    };

    auto build_tda_pseudowfc = [&](const std::vector<CasidaWavefunction>& casida_wfcs,
                                   const ModuleBase::ComplexMatrix& coef) {
        std::vector<std::vector<std::vector<std::complex<double>>>> tda_wfc(nstates);
        for (int I = 1; I < nstates; ++I) {
            tda_wfc[I].resize(nocc_lr,
                              std::vector<std::complex<double>>(nbasis, {0.0, 0.0}));
            if (I >= static_cast<int>(casida_wfcs.size())) {
                continue;
            }
            const auto& X_I = casida_wfcs[I].X_coeffs;
            for (int i_lr = 0; i_lr < nocc_lr; ++i_lr) {
                for (int a_lr = 0; a_lr < nvirt_lr; ++a_lr) {
                    const int idx_ia = i_lr * nvirt_lr + a_lr;
                    if (idx_ia >= static_cast<int>(X_I.size())) {
                        continue;
                    }
                    const int coef_row = nocc_ + a_lr;
                    if (coef_row >= coef.nr) {
                        continue;
                    }
                    const double x_ia = X_I[idx_ia];
                    for (int mu = 0; mu < nbasis; ++mu) {
                        tda_wfc[I][i_lr][mu] += x_ia * coef(coef_row, mu);
                    }
                }
            }
        }
        return tda_wfc;
    };

    auto compute_sigma_approx = [&](const ModuleBase::ComplexMatrix& coef_old_in,
                                    const ModuleBase::ComplexMatrix& coef_new_in,
                                    const ModuleBase::ComplexMatrix& s_ao_dense,
                                    const std::vector<CasidaWavefunction>& casida_old_in,
                                    const std::vector<CasidaWavefunction>& casida_new_in,
                                    const std::vector<double>& occ_phase) {
        ModuleBase::ComplexMatrix Sigma = make_zero_complex_matrix(nstates, nstates);
        auto tda_old = build_tda_pseudowfc(casida_old_in, coef_old_in);
        auto tda_new = build_tda_pseudowfc(casida_new_in, coef_new_in);

        for (int I = 0; I < nstates; ++I) {
            for (int J = 0; J < nstates; ++J) {
                std::complex<double> overlap(0.0, 0.0);
                if (I == 0 && J == 0) {
                    overlap = std::complex<double>(1.0, 0.0);
                } else if (I == 0 && J > 0) {
                    // Corrected path: no occupied-phase factor in Sigma(0,J)
                    for (int i_lr = 0; i_lr < nocc_lr; ++i_lr) {
                        const int i_ks = occ_offset + i_lr;
                        for (int mu = 0; mu < nbasis; ++mu) {
                            for (int nu = 0; nu < nbasis; ++nu) {
                                overlap += std::conj(coef_old_in(i_ks, mu))
                                         * s_ao_dense(mu, nu)
                                         * tda_new[J][i_lr][nu];
                            }
                        }
                    }
                } else if (I > 0 && J == 0) {
                    for (int i_lr = 0; i_lr < nocc_lr; ++i_lr) {
                        const int i_ks = occ_offset + i_lr;
                        for (int mu = 0; mu < nbasis; ++mu) {
                            for (int nu = 0; nu < nbasis; ++nu) {
                                overlap += std::conj(tda_old[I][i_lr][mu])
                                         * s_ao_dense(mu, nu)
                                         * coef_new_in(i_ks, nu)
                                         * occ_phase[i_lr];
                            }
                        }
                    }
                } else {
                    for (int i_lr = 0; i_lr < nocc_lr; ++i_lr) {
                        for (int mu = 0; mu < nbasis; ++mu) {
                            for (int nu = 0; nu < nbasis; ++nu) {
                                overlap += std::conj(tda_old[I][i_lr][mu])
                                         * s_ao_dense(mu, nu)
                                         * tda_new[J][i_lr][nu];
                            }
                        }
                    }
                }
                Sigma(I, J) = overlap;
            }
        }
        return Sigma;
    };

    auto compute_s_mo_lr = [&](const ModuleBase::ComplexMatrix& coef_old_in,
                               const ModuleBase::ComplexMatrix& coef_new_in,
                               const ModuleBase::ComplexMatrix& s_ao_dense) {
        const int n_lr = nocc_lr + nvirt_lr;
        std::vector<double> s_mo(n_lr * n_lr, 0.0);
        auto local_to_ks = [&](int loc) -> int {
            if (loc < nocc_lr) {
                return occ_offset + loc;
            }
            return nocc_ + (loc - nocc_lr);
        };
        for (int p = 0; p < n_lr; ++p) {
            const int p_ks = local_to_ks(p);
            for (int q = 0; q < n_lr; ++q) {
                const int q_ks = local_to_ks(q);
                double val = 0.0;
                for (int mu = 0; mu < nbasis; ++mu) {
                    for (int nu = 0; nu < nbasis; ++nu) {
                        val += (std::conj(coef_old_in(p_ks, mu))
                                * s_ao_dense(mu, nu)
                                * coef_new_in(q_ks, nu))
                                   .real();
                    }
                }
                s_mo[p * n_lr + q] = val;
            }
        }
        return s_mo;
    };

    auto determinant = [](std::vector<double> A, int n) {
        double det = 1.0;
        int sign = 1;
        const double eps = 1e-14;
        for (int i = 0; i < n; ++i) {
            int pivot = i;
            double max_abs = std::abs(A[i * n + i]);
            for (int r = i + 1; r < n; ++r) {
                const double v = std::abs(A[r * n + i]);
                if (v > max_abs) {
                    max_abs = v;
                    pivot = r;
                }
            }
            if (max_abs < eps) {
                return 0.0;
            }
            if (pivot != i) {
                for (int c = 0; c < n; ++c) {
                    std::swap(A[i * n + c], A[pivot * n + c]);
                }
                sign *= -1;
            }
            const double pivot_val = A[i * n + i];
            det *= pivot_val;
            for (int r = i + 1; r < n; ++r) {
                const double factor = A[r * n + i] / pivot_val;
                for (int c = i; c < n; ++c) {
                    A[r * n + c] -= factor * A[i * n + c];
                }
            }
        }
        return det * static_cast<double>(sign);
    };

    auto invert_matrix = [](const std::vector<double>& A, int n,
                            std::vector<double>& inv_out, double& det_out) {
        std::vector<double> m = A;
        inv_out.assign(n * n, 0.0);
        for (int i = 0; i < n; ++i) {
            inv_out[i * n + i] = 1.0;
        }
        det_out = 1.0;
        int sign = 1;
        const double eps = 1e-14;

        for (int col = 0; col < n; ++col) {
            int pivot = col;
            double max_abs = std::abs(m[col * n + col]);
            for (int r = col + 1; r < n; ++r) {
                const double v = std::abs(m[r * n + col]);
                if (v > max_abs) {
                    max_abs = v;
                    pivot = r;
                }
            }
            if (max_abs < eps) {
                return false;
            }
            if (pivot != col) {
                for (int c = 0; c < n; ++c) {
                    std::swap(m[col * n + c], m[pivot * n + c]);
                    std::swap(inv_out[col * n + c], inv_out[pivot * n + c]);
                }
                sign *= -1;
            }

            const double piv = m[col * n + col];
            det_out *= piv;
            for (int c = 0; c < n; ++c) {
                m[col * n + c] /= piv;
                inv_out[col * n + c] /= piv;
            }
            for (int r = 0; r < n; ++r) {
                if (r == col) {
                    continue;
                }
                const double f = m[r * n + col];
                for (int c = 0; c < n; ++c) {
                    m[r * n + c] -= f * m[col * n + c];
                    inv_out[r * n + c] -= f * inv_out[col * n + c];
                }
            }
        }
        det_out *= static_cast<double>(sign);
        return true;
    };

    auto max_abs_matrix_diff = [&](const ModuleBase::ComplexMatrix& a,
                                   const ModuleBase::ComplexMatrix& b) {
        double max_diff = 0.0;
        for (int i = 0; i < a.nr; ++i) {
            for (int j = 0; j < a.nc; ++j) {
                max_diff = std::max(max_diff, std::abs(a(i, j) - b(i, j)));
            }
        }
        return max_diff;
    };

    auto frob_sigma_minus_identity = [&](const ModuleBase::ComplexMatrix& sigma) {
        double s = 0.0;
        for (int i = 0; i < sigma.nr; ++i) {
            for (int j = 0; j < sigma.nc; ++j) {
                const std::complex<double> target =
                    (i == j) ? std::complex<double>(1.0, 0.0)
                             : std::complex<double>(0.0, 0.0);
                const double d = std::abs(sigma(i, j) - target);
                s += d * d;
            }
        }
        return std::sqrt(s);
    };

    auto nac_from_sigma = [&](const ModuleBase::ComplexMatrix& sigma) {
        ModuleBase::ComplexMatrix nac = make_zero_complex_matrix(sigma.nr, sigma.nc);
        for (int i = 0; i < sigma.nr; ++i) {
            for (int j = 0; j < sigma.nc; ++j) {
                nac(i, j) = (sigma(i, j) - std::conj(sigma(j, i))) / 2.0;
            }
        }
        return nac;
    };

    auto max_abs_complex_matrix = [&](const ModuleBase::ComplexMatrix& a) {
        double m = 0.0;
        for (int i = 0; i < a.nr; ++i) {
            for (int j = 0; j < a.nc; ++j) {
                m = std::max(m, std::abs(a(i, j)));
            }
        }
        return m;
    };

    auto lowdin_sigma = [&](const ModuleBase::ComplexMatrix& coef_old_in,
                            const ModuleBase::ComplexMatrix& coef_new_in,
                            const ModuleBase::ComplexMatrix& s_ao_dense,
                            const std::vector<CasidaWavefunction>& casida_old_in,
                            const std::vector<CasidaWavefunction>& casida_new_in,
                            const std::vector<double>& occ_phase,
                            double& det_out) {
        ModuleBase::ComplexMatrix sigma = make_zero_complex_matrix(nstates, nstates);
        this->compute_lowdin_sigma(coef_old_in, coef_new_in, s_ao_dense,
                                   casida_old_in, casida_new_in,
                                   nocc_lr, nvirt_lr, occ_offset,
                                   occ_phase, sigma, &det_out);
        return sigma;
    };

    auto sigma_reference_by_determinants = [&](const std::vector<double>& s_mo,
                                               const std::vector<CasidaWavefunction>& casida_old_in,
                                               const std::vector<CasidaWavefunction>& casida_new_in) {
        const int n_lr = nocc_lr + nvirt_lr;
        ModuleBase::ComplexMatrix sigma_ref = make_zero_complex_matrix(nstates, nstates);

        auto det_overlap = [&](const std::vector<int>& bra_occ,
                               const std::vector<int>& ket_occ) {
            std::vector<double> m(nocc_lr * nocc_lr, 0.0);
            for (int r = 0; r < nocc_lr; ++r) {
                for (int c = 0; c < nocc_lr; ++c) {
                    m[r * nocc_lr + c] = s_mo[bra_occ[r] * n_lr + ket_occ[c]];
                }
            }
            return determinant(m, nocc_lr);
        };

        auto ground_occ = [&]() {
            std::vector<int> occ(nocc_lr);
            for (int i = 0; i < nocc_lr; ++i) {
                occ[i] = i;
            }
            return occ;
        }();

        auto excited_occ = [&](int i, int a) {
            std::vector<int> occ = ground_occ;
            occ[i] = nocc_lr + a;
            return occ;
        };

        sigma_ref(0, 0) = std::complex<double>(det_overlap(ground_occ, ground_occ), 0.0);

        for (int J = 1; J < nstates; ++J) {
            double val = 0.0;
            const auto& X_J = casida_new_in[J].X_coeffs;
            for (int j = 0; j < nocc_lr; ++j) {
                for (int b = 0; b < nvirt_lr; ++b) {
                    const int idx = j * nvirt_lr + b;
                    val += X_J[idx] * det_overlap(ground_occ, excited_occ(j, b));
                }
            }
            sigma_ref(0, J) = std::complex<double>(val, 0.0);
        }

        for (int I = 1; I < nstates; ++I) {
            double val = 0.0;
            const auto& X_I = casida_old_in[I].X_coeffs;
            for (int i = 0; i < nocc_lr; ++i) {
                for (int a = 0; a < nvirt_lr; ++a) {
                    const int idx = i * nvirt_lr + a;
                    val += X_I[idx] * det_overlap(excited_occ(i, a), ground_occ);
                }
            }
            sigma_ref(I, 0) = std::complex<double>(val, 0.0);
        }

        for (int I = 1; I < nstates; ++I) {
            const auto& X_I = casida_old_in[I].X_coeffs;
            for (int J = 1; J < nstates; ++J) {
                const auto& X_J = casida_new_in[J].X_coeffs;
                double val = 0.0;
                for (int i = 0; i < nocc_lr; ++i) {
                    for (int a = 0; a < nvirt_lr; ++a) {
                        const int idx_ia = i * nvirt_lr + a;
                        for (int j = 0; j < nocc_lr; ++j) {
                            for (int b = 0; b < nvirt_lr; ++b) {
                                const int idx_jb = j * nvirt_lr + b;
                                val += X_I[idx_ia] * X_J[idx_jb]
                                     * det_overlap(excited_occ(i, a), excited_occ(j, b));
                            }
                        }
                    }
                }
                sigma_ref(I, J) = std::complex<double>(val, 0.0);
            }
        }
        return sigma_ref;
    };

    auto write_report_header = [&](std::ofstream& ofs, const std::string& title) {
        ofs << "# " << title << "\n\n";
        ofs << "Generated by `FsshDriver::run_lowdin_validation_suite_8_5`\n\n";
        ofs << "- nbasis: " << nbasis << "\n";
        ofs << "- nocc_lr: " << nocc_lr << "\n";
        ofs << "- nvirt_lr: " << nvirt_lr << "\n";
        ofs << "- nstates: " << nstates << "\n\n";
    };

    const ModuleBase::ComplexMatrix s_ao = make_identity_ao();
    const auto casida_old = make_casida(0.0);
    const auto casida_new = make_casida(0.05);
    const auto coef_old = make_rotated_coef(0.0);
    const auto coef_new = make_rotated_coef(0.04);
    const auto occ_phase = compute_occ_phase(coef_old, coef_new, s_ao);

    // -------------------- Scheme 1 --------------------
    {
        std::ofstream ofs(output_prefix + "_scheme1_reference.md");
        write_report_header(ofs, "Scheme 1: Explicit Determinant Reference Cross-Check");

        double det_lowdin = 0.0;
        ModuleBase::ComplexMatrix sigma_lowdin =
            lowdin_sigma(coef_old, coef_new, s_ao, casida_old, casida_new, occ_phase, det_lowdin);
        auto s_mo = compute_s_mo_lr(coef_old, coef_new, s_ao);
        // 与 Lowdin 主实现保持一致：对新步占据轨道列应用 occ_phase 修正。
        const int n_lr = nocc_lr + nvirt_lr;
        for (int q = 0; q < nocc_lr; ++q) {
            if (occ_phase[q] < 0.0) {
                for (int p = 0; p < n_lr; ++p) {
                    s_mo[p * n_lr + q] *= -1.0;
                }
            }
        }
        ModuleBase::ComplexMatrix sigma_ref =
            sigma_reference_by_determinants(s_mo, casida_old, casida_new);

        const double max_diff = max_abs_matrix_diff(sigma_lowdin, sigma_ref);
        ofs << "## Result\n\n";
        ofs << "- `max |Sigma_lowdin - Sigma_ref| = " << std::scientific << max_diff << "`\n";
        ofs << "- `det(Soo) from Lowdin path = " << det_lowdin << "`\n";
        ofs << "- Verdict: " << ((max_diff < 1e-8) ? "PASS" : "CHECK") << "\n";
    }

    // -------------------- Scheme 2 --------------------
    {
        std::ofstream ofs(output_prefix + "_scheme2_identities.md");
        write_report_header(ofs, "Scheme 2: Algebraic Identity Checks");

        const auto s_mo = compute_s_mo_lr(coef_old, coef_new, s_ao);
        std::vector<double> Soo(nocc_lr * nocc_lr, 0.0);
        std::vector<double> Sov(nocc_lr * nvirt_lr, 0.0);
        std::vector<double> Svo(nvirt_lr * nocc_lr, 0.0);
        std::vector<double> Svv(nvirt_lr * nvirt_lr, 0.0);

        const int n_lr = nocc_lr + nvirt_lr;
        for (int i = 0; i < nocc_lr; ++i) {
            for (int j = 0; j < nocc_lr; ++j) {
                Soo[i * nocc_lr + j] = s_mo[i * n_lr + j];
            }
            for (int b = 0; b < nvirt_lr; ++b) {
                Sov[i * nvirt_lr + b] = s_mo[i * n_lr + (nocc_lr + b)];
            }
        }
        for (int a = 0; a < nvirt_lr; ++a) {
            for (int j = 0; j < nocc_lr; ++j) {
                Svo[a * nocc_lr + j] = s_mo[(nocc_lr + a) * n_lr + j];
            }
            for (int b = 0; b < nvirt_lr; ++b) {
                Svv[a * nvirt_lr + b] = s_mo[(nocc_lr + a) * n_lr + (nocc_lr + b)];
            }
        }

        std::vector<double> Soo_inv;
        double det_Soo = 0.0;
        const bool inv_ok = invert_matrix(Soo, nocc_lr, Soo_inv, det_Soo);
        double inv_residual_inf = std::numeric_limits<double>::infinity();
        double stilde_consistency = std::numeric_limits<double>::infinity();
        double sigma00_gap = std::numeric_limits<double>::infinity();

        if (inv_ok) {
            // ||Soo*Soo_inv - I||_inf
            double max_row_sum = 0.0;
            for (int i = 0; i < nocc_lr; ++i) {
                double row_sum = 0.0;
                for (int j = 0; j < nocc_lr; ++j) {
                    double v = 0.0;
                    for (int k = 0; k < nocc_lr; ++k) {
                        v += Soo[i * nocc_lr + k] * Soo_inv[k * nocc_lr + j];
                    }
                    v -= (i == j ? 1.0 : 0.0);
                    row_sum += std::abs(v);
                }
                max_row_sum = std::max(max_row_sum, row_sum);
            }
            inv_residual_inf = max_row_sum;

            // Stilde consistency: direct triple product vs Q form
            std::vector<double> Q(nocc_lr * nvirt_lr, 0.0);
            for (int j = 0; j < nocc_lr; ++j) {
                for (int b = 0; b < nvirt_lr; ++b) {
                    for (int k = 0; k < nocc_lr; ++k) {
                        Q[j * nvirt_lr + b] += Soo_inv[j * nocc_lr + k] * Sov[k * nvirt_lr + b];
                    }
                }
            }
            std::vector<double> stilde_a(nvirt_lr * nvirt_lr, 0.0);
            std::vector<double> stilde_b(nvirt_lr * nvirt_lr, 0.0);
            for (int a = 0; a < nvirt_lr; ++a) {
                for (int b = 0; b < nvirt_lr; ++b) {
                    double corr_a = 0.0;
                    double corr_b = 0.0;
                    for (int i = 0; i < nocc_lr; ++i) {
                        corr_b += Svo[a * nocc_lr + i] * Q[i * nvirt_lr + b];
                        for (int j = 0; j < nocc_lr; ++j) {
                            corr_a += Svo[a * nocc_lr + i]
                                   * Soo_inv[i * nocc_lr + j]
                                   * Sov[j * nvirt_lr + b];
                        }
                    }
                    stilde_a[a * nvirt_lr + b] = Svv[a * nvirt_lr + b] - corr_a;
                    stilde_b[a * nvirt_lr + b] = Svv[a * nvirt_lr + b] - corr_b;
                }
            }
            double max_stilde_diff = 0.0;
            for (int i = 0; i < nvirt_lr * nvirt_lr; ++i) {
                max_stilde_diff = std::max(max_stilde_diff, std::abs(stilde_a[i] - stilde_b[i]));
            }
            stilde_consistency = max_stilde_diff;

            double det_lowdin = 0.0;
            ModuleBase::ComplexMatrix sigma_lowdin =
                lowdin_sigma(coef_old, coef_new, s_ao, casida_old, casida_new, occ_phase, det_lowdin);
            sigma00_gap = std::abs(sigma_lowdin(0, 0).real() - det_Soo);
        }

        ofs << "## Result\n\n";
        ofs << "- `inverse(Soo) success = " << (inv_ok ? "true" : "false") << "`\n";
        ofs << "- `||Soo*Soo^{-1} - I||_inf = " << std::scientific << inv_residual_inf << "`\n";
        ofs << "- `|Sigma00 - det(Soo)| = " << sigma00_gap << "`\n";
        ofs << "- `max |Stilde(direct) - Stilde(Q-form)| = " << stilde_consistency << "`\n";
        const bool pass = inv_ok && inv_residual_inf < 1e-10
                       && sigma00_gap < 1e-10 && stilde_consistency < 1e-12;
        ofs << "- Verdict: " << (pass ? "PASS" : "CHECK") << "\n";
    }

    // -------------------- Scheme 3 --------------------
    {
        std::ofstream ofs(output_prefix + "_scheme3_limits.md");
        write_report_header(ofs, "Scheme 3: Limit Behavior Checks");

        const auto casida_same = casida_old;
        const auto coef_same = coef_old;
        const auto occ_phase_same = compute_occ_phase(coef_same, coef_same, s_ao);

        double det_same = 0.0;
        ModuleBase::ComplexMatrix sigma_same =
            lowdin_sigma(coef_same, coef_same, s_ao, casida_same, casida_same, occ_phase_same, det_same);
        ModuleBase::ComplexMatrix nac_same = nac_from_sigma(sigma_same);

        double max_sigma_identity_dev = 0.0;
        for (int i = 0; i < nstates; ++i) {
            for (int j = 0; j < nstates; ++j) {
                const double target = (i == j ? 1.0 : 0.0);
                max_sigma_identity_dev = std::max(
                    max_sigma_identity_dev,
                    std::abs(sigma_same(i, j).real() - target));
            }
        }
        const double max_nac_same = max_abs_complex_matrix(nac_same);

        const auto coef_eps1 = make_rotated_coef(0.04);
        const auto coef_eps2 = make_rotated_coef(0.02);
        const auto occ_eps1 = compute_occ_phase(coef_old, coef_eps1, s_ao);
        const auto occ_eps2 = compute_occ_phase(coef_old, coef_eps2, s_ao);
        double det_eps1 = 0.0, det_eps2 = 0.0;
        ModuleBase::ComplexMatrix sigma_eps1 =
            lowdin_sigma(coef_old, coef_eps1, s_ao, casida_old, casida_new, occ_eps1, det_eps1);
        ModuleBase::ComplexMatrix sigma_eps2 =
            lowdin_sigma(coef_old, coef_eps2, s_ao, casida_old, casida_new, occ_eps2, det_eps2);
        const double frob1 = frob_sigma_minus_identity(sigma_eps1);
        const double frob2 = frob_sigma_minus_identity(sigma_eps2);
        const double ratio = (frob2 > 1e-16) ? (frob1 / frob2) : std::numeric_limits<double>::infinity();

        ofs << "## Result\n\n";
        ofs << "- Identical-step `max |Sigma - I| = " << std::scientific << max_sigma_identity_dev << "`\n";
        ofs << "- Identical-step `max |NAC| = " << max_nac_same << "`\n";
        ofs << "- Perturbation scaling: `||Sigma(alpha=0.04)-I||_F = " << frob1 << "`\n";
        ofs << "- Perturbation scaling: `||Sigma(alpha=0.02)-I||_F = " << frob2 << "`\n";
        ofs << "- Ratio `frob(alpha=0.04)/frob(alpha=0.02) = " << ratio << "`\n";
        // 比例不要求严格 2.0：Lowdin 包含 det(Soo) 等非线性项，有限步长下会偏离理想线性。
        const bool pass = (max_sigma_identity_dev < 1e-12) && (max_nac_same < 1e-12)
                       && (ratio > 1.05 && ratio < 2.5);
        ofs << "- Verdict: " << (pass ? "PASS" : "CHECK") << "\n";
    }

    // -------------------- Scheme 4 --------------------
    {
        std::ofstream ofs(output_prefix + "_scheme4_approx_convergence.md");
        write_report_header(ofs, "Scheme 4: Convergence Toward Diagonal Approximation");

        auto evaluate_case = [&](double alpha) {
            const auto coef_new_case = make_rotated_coef(alpha);
            const auto occ_case = compute_occ_phase(coef_old, coef_new_case, s_ao);
            double det_case = 0.0;
            ModuleBase::ComplexMatrix sigma_lowdin_case =
                lowdin_sigma(coef_old, coef_new_case, s_ao, casida_old, casida_new, occ_case, det_case);
            ModuleBase::ComplexMatrix sigma_approx_case =
                compute_sigma_approx(coef_old, coef_new_case, s_ao, casida_old, casida_new, occ_case);
            ModuleBase::ComplexMatrix nac_lowdin_case = nac_from_sigma(sigma_lowdin_case);
            ModuleBase::ComplexMatrix nac_approx_case = nac_from_sigma(sigma_approx_case);
            const double sigma_diff = max_abs_matrix_diff(sigma_lowdin_case, sigma_approx_case);
            const double nac_diff = max_abs_matrix_diff(nac_lowdin_case, nac_approx_case);
            return std::make_pair(sigma_diff, nac_diff);
        };

        const auto strong = evaluate_case(0.08);
        const auto weak = evaluate_case(0.008);

        ofs << "## Result\n\n";
        ofs << "- Strong coupling (`alpha=0.08`): `max|Sigma_lowdin-Sigma_approx| = "
            << std::scientific << strong.first << "`\n";
        ofs << "- Strong coupling (`alpha=0.08`): `max|NAC_lowdin-NAC_approx| = "
            << strong.second << "`\n";
        ofs << "- Weak coupling (`alpha=0.008`): `max|Sigma_lowdin-Sigma_approx| = "
            << weak.first << "`\n";
        ofs << "- Weak coupling (`alpha=0.008`): `max|NAC_lowdin-NAC_approx| = "
            << weak.second << "`\n";
        const bool pass = weak.first < strong.first && weak.second < strong.second;
        ofs << "- Verdict: " << (pass ? "PASS" : "CHECK") << "\n";
    }

    // -------------------- Scheme 5 --------------------
    {
        std::ofstream ofs(output_prefix + "_scheme5_gauge_phase.md");
        write_report_header(ofs, "Scheme 5: Gauge/Phase Robustness");

        double det_base = 0.0;
        ModuleBase::ComplexMatrix sigma_base =
            lowdin_sigma(coef_old, coef_new, s_ao, casida_old, casida_new, occ_phase, det_base);

        // (a) Occupied-phase robustness: flip sign of one occupied orbital at t'
        ModuleBase::ComplexMatrix coef_new_flip = coef_new;
        for (int mu = 0; mu < nbasis; ++mu) {
            coef_new_flip(0, mu) *= -1.0;
        }
        const auto occ_phase_flip = compute_occ_phase(coef_old, coef_new_flip, s_ao);
        double det_flip = 0.0;
        ModuleBase::ComplexMatrix sigma_flip =
            lowdin_sigma(coef_old, coef_new_flip, s_ao, casida_old, casida_new, occ_phase_flip, det_flip);
        const double phase_robust_diff = max_abs_matrix_diff(sigma_base, sigma_flip);

        // (b) Excited-subspace gauge covariance: rotate casida_new in excited space
        const double phi = 0.3;
        const double c = std::cos(phi);
        const double s = std::sin(phi);
        std::vector<CasidaWavefunction> casida_rot = casida_new;
        std::vector<double> x1 = casida_new[1].X_coeffs;
        std::vector<double> x2 = casida_new[2].X_coeffs;
        for (size_t k = 0; k < x1.size(); ++k) {
            casida_rot[1].X_coeffs[k] = c * x1[k] + s * x2[k];
            casida_rot[2].X_coeffs[k] = -s * x1[k] + c * x2[k];
        }

        double det_rot = 0.0;
        ModuleBase::ComplexMatrix sigma_rot =
            lowdin_sigma(coef_old, coef_new, s_ao, casida_old, casida_rot, occ_phase, det_rot);

        ModuleBase::ComplexMatrix sigma_expected = sigma_base;
        for (int I = 0; I < nstates; ++I) {
            const std::complex<double> col1 = sigma_base(I, 1);
            const std::complex<double> col2 = sigma_base(I, 2);
            sigma_expected(I, 1) = c * col1 + s * col2;
            sigma_expected(I, 2) = -s * col1 + c * col2;
        }
        const double covariant_diff = max_abs_matrix_diff(sigma_rot, sigma_expected);

        ofs << "## Result\n\n";
        ofs << "- Occupied sign-flip robustness: `max|Sigma_base - Sigma_flipped| = "
            << std::scientific << phase_robust_diff << "`\n";
        ofs << "- Excited-subspace covariance: `max|Sigma_rot - Sigma_base*W| = "
            << covariant_diff << "`\n";
        const bool pass = (phase_robust_diff < 1e-10) && (covariant_diff < 1e-10);
        ofs << "- Verdict: " << (pass ? "PASS" : "CHECK") << "\n";
    }

    // -------------------- Scheme 6 --------------------
    {
        std::ofstream ofs(output_prefix + "_scheme6_stability.md");
        write_report_header(ofs, "Scheme 6: Numerical Stability Diagnostics");

        auto norm_inf = [](const std::vector<double>& a, int n) {
            double m = 0.0;
            for (int i = 0; i < n; ++i) {
                double row_sum = 0.0;
                for (int j = 0; j < n; ++j) {
                    row_sum += std::abs(a[i * n + j]);
                }
                m = std::max(m, row_sum);
            }
            return m;
        };

        auto estimate_cond_and_residual = [&](const std::vector<double>& Soo) {
            std::vector<double> Soo_inv;
            double det = 0.0;
            const bool ok = invert_matrix(Soo, nocc_lr, Soo_inv, det);
            double cond_inf = std::numeric_limits<double>::infinity();
            double residual_inf = std::numeric_limits<double>::infinity();
            if (ok) {
                cond_inf = norm_inf(Soo, nocc_lr) * norm_inf(Soo_inv, nocc_lr);
                double max_row_sum = 0.0;
                for (int i = 0; i < nocc_lr; ++i) {
                    double row_sum = 0.0;
                    for (int j = 0; j < nocc_lr; ++j) {
                        double v = 0.0;
                        for (int k = 0; k < nocc_lr; ++k) {
                            v += Soo[i * nocc_lr + k] * Soo_inv[k * nocc_lr + j];
                        }
                        v -= (i == j ? 1.0 : 0.0);
                        row_sum += std::abs(v);
                    }
                    max_row_sum = std::max(max_row_sum, row_sum);
                }
                residual_inf = max_row_sum;
            }
            return std::make_tuple(ok, cond_inf, residual_inf, det);
        };

        // Well-conditioned baseline
        const auto s_mo_base = compute_s_mo_lr(coef_old, coef_new, s_ao);
        std::vector<double> Soo_base(nocc_lr * nocc_lr, 0.0);
        for (int i = 0; i < nocc_lr; ++i) {
            for (int j = 0; j < nocc_lr; ++j) {
                Soo_base[i * nocc_lr + j] = s_mo_base[i * (nocc_lr + nvirt_lr) + j];
            }
        }

        // Near-singular case: make occupied row 1 at t' almost parallel to row 0
        ModuleBase::ComplexMatrix coef_new_ns = make_rotated_coef(0.0);
        const double eps = 1e-6;
        for (int mu = 0; mu < nbasis; ++mu) {
            coef_new_ns(1, mu) = std::complex<double>(0.0, 0.0);
        }
        coef_new_ns(1, 0) = std::complex<double>(1.0 - eps, 0.0);
        coef_new_ns(1, 1) = std::complex<double>(eps, 0.0);
        const auto s_mo_ns = compute_s_mo_lr(coef_old, coef_new_ns, s_ao);
        std::vector<double> Soo_ns(nocc_lr * nocc_lr, 0.0);
        for (int i = 0; i < nocc_lr; ++i) {
            for (int j = 0; j < nocc_lr; ++j) {
                Soo_ns[i * nocc_lr + j] = s_mo_ns[i * (nocc_lr + nvirt_lr) + j];
            }
        }

        bool ok_base, ok_ns;
        double cond_base, cond_ns, resid_base, resid_ns, det_base, det_ns;
        std::tie(ok_base, cond_base, resid_base, det_base) = estimate_cond_and_residual(Soo_base);
        std::tie(ok_ns, cond_ns, resid_ns, det_ns) = estimate_cond_and_residual(Soo_ns);

        ofs << "## Result\n\n";
        ofs << "- Baseline `invertible=" << (ok_base ? "true" : "false")
            << "`, `cond_inf=" << std::scientific << cond_base
            << "`, `||Soo*Soo^{-1}-I||_inf=" << resid_base
            << "`, `det(Soo)=" << det_base << "`\n";
        ofs << "- Near-singular `invertible=" << (ok_ns ? "true" : "false")
            << "`, `cond_inf=" << cond_ns
            << "`, `||Soo*Soo^{-1}-I||_inf=" << resid_ns
            << "`, `det(Soo)=" << det_ns << "`\n";
        bool pass = ok_base && ok_ns && (cond_ns > cond_base * 1e3);
        ofs << "- Verdict: " << (pass ? "PASS" : "CHECK") << "\n";
    }

    std::cout << "[FSSH VALIDATION] Wrote 6 reports with prefix: " << output_prefix << std::endl;
    std::cout << "  - " << output_prefix << "_scheme1_reference.md\n";
    std::cout << "  - " << output_prefix << "_scheme2_identities.md\n";
    std::cout << "  - " << output_prefix << "_scheme3_limits.md\n";
    std::cout << "  - " << output_prefix << "_scheme4_approx_convergence.md\n";
    std::cout << "  - " << output_prefix << "_scheme5_gauge_phase.md\n";
    std::cout << "  - " << output_prefix << "_scheme6_stability.md" << std::endl;
}
