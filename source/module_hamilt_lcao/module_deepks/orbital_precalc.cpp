#ifdef __DEEPKS

/// cal_orbital_precalc : orbital_precalc is usted for training with orbital label,
///                          which equals gvdm * orbital_pdm_shell,
///                          orbital_pdm_shell[Inl,nm*nm] = dm_hl * overlap * overlap

#include "LCAO_deepks.h"
#include "LCAO_deepks_io.h" // mohan add 2024-07-22
#include "module_base/blas_connector.h"
#include "module_base/constants.h"
#include "module_base/libm/libm.h"
#include "module_base/parallel_reduce.h"
#include "module_hamilt_lcao/module_hcontainer/atom_pair.h"
#include "module_parameter/parameter.h"

// calculates orbital_precalc[1,NAt,NDscrpt] = gvdm * orbital_pdm_shell;
// orbital_pdm_shell[Inl,nm*nm] = dm_hl * overlap * overlap;
template <typename TK, typename TH>
void LCAO_Deepks::cal_orbital_precalc(const std::vector<TH>& dm_hl,
                                      const int lmaxd,
                                      const int inlmax,
                                      const int nat,
                                      const int nks,
                                      const int* inl_l,
                                      const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                                      const std::vector<hamilt::HContainer<double>*> phialpha,
                                      const ModuleBase::IntArray* inl_index,
                                      const UnitCell& ucell,
                                      const LCAO_Orbitals& orb,
                                      const Parallel_Orbitals& pv,
                                      const Grid_Driver& GridD)
{
    ModuleBase::TITLE("LCAO_Deepks", "cal_orbital_precalc");
    ModuleBase::timer::tick("LCAO_Deepks", "calc_orbital_precalc");

    this->cal_gvdm(nat);

    const double Rcut_Alpha = orb.Alpha[0].getRcut();

    this->init_orbital_pdm_shell(nks);

    for (int T0 = 0; T0 < ucell.ntype; T0++)
    {
        Atom* atom0 = &ucell.atoms[T0];

        for (int I0 = 0; I0 < atom0->na; I0++)
        {
            const int iat = ucell.itia2iat(T0, I0);
            const ModuleBase::Vector3<double> tau0 = atom0->tau[I0];
            GridD.Find_atom(ucell, atom0->tau[I0], T0, I0);

            // trace alpha orbital
            std::vector<int> trace_alpha_row;
            std::vector<int> trace_alpha_col;
            int ib = 0;
            for (int L0 = 0; L0 <= orb.Alpha[0].getLmax(); ++L0)
            {
                for (int N0 = 0; N0 < orb.Alpha[0].getNchi(L0); ++N0)
                {
                    const int inl = inl_index[T0](I0, L0, N0);
                    const int nm = 2 * L0 + 1;

                    for (int m1 = 0; m1 < nm; ++m1) // m1 = 1 for s, 3 for p, 5 for d
                    {
                        for (int m2 = 0; m2 < nm; ++m2) // m1 = 1 for s, 3 for p, 5 for d
                        {
                            trace_alpha_row.push_back(ib + m1);
                            trace_alpha_col.push_back(ib + m2);
                        }
                    }
                    ib += nm;
                }
            }
            const int trace_alpha_size = trace_alpha_row.size();

            for (int ad1 = 0; ad1 < GridD.getAdjacentNum() + 1; ++ad1)
            {
                const int T1 = GridD.getType(ad1);
                const int I1 = GridD.getNatom(ad1);
                const int ibt1 = ucell.itia2iat(T1, I1);
                const ModuleBase::Vector3<double> tau1 = GridD.getAdjacentTau(ad1);
                const Atom* atom1 = &ucell.atoms[T1];
                const int nw1_tot = atom1->nw * PARAM.globalv.npol;
                const double Rcut_AO1 = orb.Phi[T1].getRcut();

                const double dist1 = (tau1 - tau0).norm() * ucell.lat0;
                if (dist1 >= Rcut_Alpha + Rcut_AO1)
                {
                    continue;
                }

                ModuleBase::Vector3<int> dR1(GridD.getBox(ad1).x, GridD.getBox(ad1).y, GridD.getBox(ad1).z);

                if constexpr (std::is_same<TK, std::complex<double>>::value)
                {
                    if (phialpha[0]->find_matrix(iat, ibt1, dR1.x, dR1.y, dR1.z) == nullptr)
                    {
                        continue;
                    }
                }

                auto row_indexes = pv.get_indexes_row(ibt1);

                const int row_size = row_indexes.size();

                if (row_size == 0)
                {
                    continue;
                }

                std::vector<double> s_1t(trace_alpha_size * row_size);

                std::vector<double> g_1dmt(nks * trace_alpha_size * row_size, 0.0);

                for (int irow = 0; irow < row_size; irow++)
                {
                    hamilt::BaseMatrix<double>* row_ptr = phialpha[0]->find_matrix(iat, ibt1, dR1);
                    for (int i = 0; i < trace_alpha_size; i++)
                    {
                        s_1t[i * row_size + irow] = row_ptr->get_value(row_indexes[irow], trace_alpha_row[i]);
                    }
                }

                for (int ad2 = 0; ad2 < GridD.getAdjacentNum() + 1; ad2++)
                {
                    const int T2 = GridD.getType(ad2);
                    const int I2 = GridD.getNatom(ad2);
                    const int ibt2 = ucell.itia2iat(T2, I2);
                    if constexpr (std::is_same<TK, std::complex<double>>::value) // Why only for multi-k?
                    {
                        if (ibt1 > ibt2)
                        {
                            continue;
                        }
                    }
                    const ModuleBase::Vector3<double> tau2 = GridD.getAdjacentTau(ad2);
                    const Atom* atom2 = &ucell.atoms[T2];
                    const int nw2_tot = atom2->nw * PARAM.globalv.npol;

                    const double Rcut_AO2 = orb.Phi[T2].getRcut();
                    const double dist2 = (tau2 - tau0).norm() * ucell.lat0;

                    if (dist2 >= Rcut_Alpha + Rcut_AO2)
                    {
                        continue;
                    }

                    ModuleBase::Vector3<int> dR2(GridD.getBox(ad2).x, GridD.getBox(ad2).y, GridD.getBox(ad2).z);

                    if constexpr (std::is_same<TK, std::complex<double>>::value)
                    {
                        if (phialpha[0]->find_matrix(iat, ibt2, dR2.x, dR2.y, dR2.z) == nullptr)
                        {
                            continue;
                        }
                    }

                    auto col_indexes = pv.get_indexes_col(ibt2);
                    const int col_size = col_indexes.size();
                    if (col_size == 0)
                    {
                        continue;
                    }

                    std::vector<double> s_2t(trace_alpha_size * col_size);
                    for (int icol = 0; icol < col_size; icol++)
                    {
                        hamilt::BaseMatrix<double>* col_ptr = phialpha[0]->find_matrix(iat, ibt2, dR2);
                        for (int i = 0; i < trace_alpha_size; i++)
                        {
                            s_2t[i * col_size + icol] = col_ptr->get_value(col_indexes[icol], trace_alpha_col[i]);
                        }
                    }

                    std::vector<double> dm_array(row_size * nks * col_size, 0.0);

                    const int row_size_nks = row_size * nks;

                    if constexpr (std::is_same<TK, double>::value)
                    {
                        for (int is = 0; is < PARAM.inp.nspin; is++)
                        {
                            hamilt::AtomPair<double> dm_pair(ibt1, ibt2, 0, 0, 0, &pv);

                            dm_pair.allocate(dm_array.data(), 0);

                            if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER(PARAM.inp.ks_solver))
                            {
                                dm_pair.add_from_matrix(dm_hl[is].c, pv.get_row_size(), 1.0, 1);
                            }
                            else
                            {
                                dm_pair.add_from_matrix(dm_hl[is].c, pv.get_col_size(), 1.0, 0);
                            }
                        } // is
                    }
                    else
                    {
                        for (int ik = 0; ik < nks; ik++)
                        {
                            hamilt::AtomPair<double> dm_pair(ibt1,
                                                             ibt2,
                                                             (dR2 - dR1).x,
                                                             (dR2 - dR1).y,
                                                             (dR2 - dR1).z,
                                                             &pv);

                            dm_pair.allocate(&dm_array[ik * row_size * col_size], 0);

                            const double arg
                                = -(kvec_d[ik] * ModuleBase::Vector3<double>(dR2 - dR1)) * ModuleBase::TWO_PI;

                            double sinp, cosp;

                            ModuleBase::libm::sincos(arg, &sinp, &cosp);

                            const std::complex<double> kphase = std::complex<double>(cosp, sinp);

                            if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER(PARAM.inp.ks_solver))
                            {
                                dm_pair.add_from_matrix(dm_hl[ik].c, pv.get_row_size(), kphase, 1);
                            }
                            else
                            {
                                dm_pair.add_from_matrix(dm_hl[ik].c, pv.get_col_size(), kphase, 0);
                            }
                        } // ik
                    }

                    // dgemm for s_2t and dm_array to get g_1dmt
                    constexpr char transa = 'T', transb = 'N';
                    double gemm_alpha = 1.0, gemm_beta = 1.0;

                    if constexpr (std::is_same<TK, std::complex<double>>::value)
                    {
                        if (ibt1 < ibt2)
                        {
                            gemm_alpha = 2.0;
                        }
                    }

                    dgemm_(&transa,
                           &transb,
                           &row_size_nks,
                           &trace_alpha_size,
                           &col_size,
                           &gemm_alpha,
                           dm_array.data(),
                           &col_size,
                           s_2t.data(),
                           &col_size,
                           &gemm_beta,
                           g_1dmt.data(),
                           &row_size_nks);
                } // ad2

                for (int ik = 0; ik < nks; ik++)
                {
                    // do dot of g_1dmt and s_1t to get orbital_pdm_shell
                    const double* p_g1dmt = g_1dmt.data() + ik * row_size;

                    int ib = 0, index = 0, inc = 1;

                    for (int L0 = 0; L0 <= orb.Alpha[0].getLmax(); ++L0)
                    {
                        for (int N0 = 0; N0 < orb.Alpha[0].getNchi(L0); ++N0)
                        {
                            const int inl = inl_index[T0](I0, L0, N0);
                            const int nm = 2 * L0 + 1;

                            for (int m1 = 0; m1 < nm; ++m1) // m1 = 1 for s, 3 for p, 5 for d
                            {
                                for (int m2 = 0; m2 < nm; ++m2) // m1 = 1 for s, 3 for p, 5 for d
                                {
                                    orbital_pdm_shell[ik][inl][m1 * nm + m2]
                                        += ddot_(&row_size,
                                                 p_g1dmt + index * row_size * nks,
                                                 &inc,
                                                 s_1t.data() + index * row_size,
                                                 &inc);
                                    index++;
                                }
                            }
                            ib += nm;
                        }
                    }
                }
            } // ad1
        }
    }
#ifdef __MPI
    for (int iks = 0; iks < nks; iks++)
    {
        for (int inl = 0; inl < inlmax; inl++)
        {
            Parallel_Reduce::reduce_all(this->orbital_pdm_shell[iks][inl],
                                        (2 * lmaxd + 1) * (2 * lmaxd + 1));
        }
    }
#endif

    // transfer orbital_pdm_shell to orbital_pdm_shell_vector

    int nlmax = inlmax / nat;

    std::vector<torch::Tensor> orbital_pdm_shell_vector;
    for (int nl = 0; nl < nlmax; ++nl)
    {
        std::vector<torch::Tensor> kammv;
        for (int iks = 0; iks < nks; ++iks)
        {
            std::vector<torch::Tensor> ammv;
            for (int iat = 0; iat < nat; ++iat)
            {
                int inl = iat * nlmax + nl;
                int nm = 2 * inl_l[inl] + 1;
                std::vector<double> mmv;

                for (int m1 = 0; m1 < nm; ++m1) // m1 = 1 for s, 3 for p, 5 for d
                {
                    for (int m2 = 0; m2 < nm; ++m2) // m1 = 1 for s, 3 for p, 5 for d
                    {
                        mmv.push_back(this->orbital_pdm_shell[iks][inl][m1 * nm + m2]);
                    }
                }
                torch::Tensor mm
                    = torch::tensor(mmv, torch::TensorOptions().dtype(torch::kFloat64)).reshape({nm, nm}); // nm*nm

                ammv.push_back(mm);
            }
            torch::Tensor amm = torch::stack(ammv, 0);
            kammv.push_back(amm);
        }
        torch::Tensor kamm = torch::stack(kammv, 0);
        orbital_pdm_shell_vector.push_back(kamm);
    }

    assert(orbital_pdm_shell_vector.size() == nlmax);

    // einsum for each nl:
    std::vector<torch::Tensor> orbital_precalc_vector;
    for (int nl = 0; nl < nlmax; ++nl)
    {
        orbital_precalc_vector.push_back(
            at::einsum("kamn, avmn->kav", {orbital_pdm_shell_vector[nl], this->gevdm_vector[nl]}));
    }

    this->orbital_precalc_tensor = torch::cat(orbital_precalc_vector, -1);
    this->del_orbital_pdm_shell(nks);
    ModuleBase::timer::tick("LCAO_Deepks", "calc_orbital_precalc");
    return;
}

template void LCAO_Deepks::cal_orbital_precalc<double, ModuleBase::matrix>(
    const std::vector<ModuleBase::matrix>& dm_hl,
    const int lmaxd,
    const int inlmax,
    const int nat,
    const int nks,
    const int* inl_l,
    const std::vector<ModuleBase::Vector3<double>>& kvec_d,
    const std::vector<hamilt::HContainer<double>*> phialpha,
    const ModuleBase::IntArray* inl_index,
    const UnitCell& ucell,
    const LCAO_Orbitals& orb,
    const Parallel_Orbitals& pv,
    const Grid_Driver& GridD);

template void LCAO_Deepks::cal_orbital_precalc<std::complex<double>, ModuleBase::ComplexMatrix>(
    const std::vector<ModuleBase::ComplexMatrix>& dm_hl,
    const int lmaxd,
    const int inlmax,
    const int nat,
    const int nks,
    const int* inl_l,
    const std::vector<ModuleBase::Vector3<double>>& kvec_d,
    const std::vector<hamilt::HContainer<double>*> phialpha,
    const ModuleBase::IntArray* inl_index,
    const UnitCell& ucell,
    const LCAO_Orbitals& orb,
    const Parallel_Orbitals& pv,
    const Grid_Driver& GridD);
#endif
