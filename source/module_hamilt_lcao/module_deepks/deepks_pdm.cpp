// wenfei 2022-1-11
// This file contains subroutines for calculating pdm,
#include "module_parameter/parameter.h"
// which is defind as sum_mu,nu rho_mu,nu (<chi_mu|alpha><alpha|chi_nu>);
// as well as gdmx, which is the gradient of pdm, defined as
// sum_mu,nu rho_mu,nu d/dX(<chi_mu|alpha><alpha|chi_nu>)

// It also contains subroutines for printing pdm and gdmx
// for checking purpose

// There are 3 subroutines in this file:
// 1. read_pdm, which reads pdm from file
// 2. cal_pdm, which is used for calculating pdm
// 3. check_pdm, which prints pdm to descriptor.dat

#ifdef __DEEPKS

#include "deepks_iterate.h"
#include "deepks_pdm.h"
#include "module_base/constants.h"
#include "module_base/libm/libm.h"
#include "module_base/timer.h"
#include "module_hamilt_lcao/module_hcontainer/atom_pair.h"
#ifdef __MPI
#include "module_base/parallel_reduce.h"
#endif

void DeePKS_domain::read_pdm(bool read_pdm_file,
                             bool is_equiv,
                             bool& init_pdm,
                             const int nat,
                             const int inlmax,
                             const int lmaxd,
                             const std::vector<int>& inl2l,
                             const Numerical_Orbital& alpha,
                             std::vector<torch::Tensor>& pdm)
{
    if (read_pdm_file && !init_pdm) // for DeePKS NSCF calculation
    {
        const std::string file_projdm = PARAM.globalv.global_out_dir + "deepks_projdm.dat";
        std::ifstream ifs(file_projdm.c_str());

        if (!ifs)
        {
            ModuleBase::WARNING_QUIT("DeePKS_domain::read_pdm", "Cannot find the file deepks_projdm.dat");
        }
        if (!is_equiv)
        {
            for (int inl = 0; inl < inlmax; inl++)
            {
                int nm = 2 * inl2l[inl] + 1;
                auto accessor = pdm[inl].accessor<double, 2>();
                for (int m1 = 0; m1 < nm; m1++)
                {
                    for (int m2 = 0; m2 < nm; m2++)
                    {
                        double c;
                        ifs >> c;
                        accessor[m1][m2] = c;
                    }
                }
            }
        }
        else
        {
            int pdm_size = 0;
            int nproj = 0;
            for (int il = 0; il < lmaxd + 1; il++)
            {
                nproj += (2 * il + 1) * alpha.getNchi(il);
            }
            pdm_size = nproj * nproj;
            for (int iat = 0; iat < nat; iat++)
            {
                auto accessor = pdm[iat].accessor<double, 1>();
                for (int ind = 0; ind < pdm_size; ind++)
                {
                    double c;
                    ifs >> c;
                    accessor[ind] = c;
                }
            }
        }

        init_pdm = true;
    }
}

template <typename TK>
void DeePKS_domain::update_dmr(const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                               const std::vector<std::vector<TK>>& dmk,
                               const UnitCell& ucell,
                               const LCAO_Orbitals& orb,
                               const Parallel_Orbitals& pv,
                               const Grid_Driver& GridD,
                               hamilt::HContainer<double>* dmr_deepks)
{
    dmr_deepks->set_zero();
    // save whether the pair with R has been calculated
    std::vector<std::tuple<int, int, int, int, int>> calculated_pairs(0);

    DeePKS_domain::iterate_ad2(
        ucell,
        GridD,
        orb,
        false, // no trace_alpha
        [&](const int iat,
            const ModuleBase::Vector3<double>& tau0,
            const int ibt1,
            const ModuleBase::Vector3<double>& tau1,
            const int start1,
            const int nw1_tot,
            ModuleBase::Vector3<int> dR1,
            const int ibt2,
            const ModuleBase::Vector3<double>& tau2,
            const int start2,
            const int nw2_tot,
            ModuleBase::Vector3<int> dR2) 
        {
            auto row_indexes = pv.get_indexes_row(ibt1);
            auto col_indexes = pv.get_indexes_col(ibt2);
            if (row_indexes.size() * col_indexes.size() == 0)
            {
                return; // to next loop
            }

            hamilt::AtomPair<double> dm_pair = dmr_deepks->get_atom_pair(ibt1, ibt2);

            int dRx = 0;
            int dRy = 0;
            int dRz = 0;
            if (std::is_same<TK, std::complex<double>>::value)
            {
                dRx = (dR1 - dR2).x;
                dRy = (dR1 - dR2).y;
                dRz = (dR1 - dR2).z;
            }
            ModuleBase::Vector3<int> dR(dRx, dRy, dRz);

            // avoid duplicate calculation
            if (std::find(calculated_pairs.begin(), calculated_pairs.end(),
                          std::make_tuple(ibt1, ibt2, dR.x, dR.y, dR.z))
                != calculated_pairs.end())
            {
                return;
            }
            calculated_pairs.push_back(std::make_tuple(ibt1, ibt2, dR.x, dR.y, dR.z));

            dm_pair.find_R(dR);
            hamilt::BaseMatrix<double>* dmr_ptr = dm_pair.find_matrix(dR);
            dmr_ptr->set_zero(); // must reset to zero to avoid accumulation!

            for (int ik = 0; ik < dmk.size(); ik++)
            {
                std::complex<double> kphase = std::complex<double>(1, 0);
                if (std::is_same<TK, std::complex<double>>::value)
                {
                    const double arg = -(kvec_d[ik] * ModuleBase::Vector3<double>(dR)) * ModuleBase::TWO_PI;
                    kphase = std::complex<double>(cos(arg), sin(arg));
                }
                TK* kphase_ptr = reinterpret_cast<TK*>(&kphase);
                if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER(PARAM.inp.ks_solver))
                {
                    dm_pair.add_from_matrix(dmk[ik].data(), pv.get_row_size(), *kphase_ptr, 1);
                }
                else
                {
                    dm_pair.add_from_matrix(dmk[ik].data(), pv.get_col_size(), *kphase_ptr, 0);
                }
            }
        }
    );
    return;
}

// this subroutine performs the calculation of projected density matrices
// pdm_m,m'=\sum_{mu,nu} rho_{mu,nu} <chi_mu|alpha_m><alpha_m'|chi_nu>
template <typename TK>
void DeePKS_domain::cal_pdm(bool& init_pdm,
                            const int inlmax,
                            const int lmaxd,
                            const std::vector<int>& inl2l,
                            const ModuleBase::IntArray* inl_index,
                            const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                            const hamilt::HContainer<double>* dmr,
                            const std::vector<hamilt::HContainer<double>*> phialpha,
                            const UnitCell& ucell,
                            const LCAO_Orbitals& orb,
                            const Grid_Driver& GridD,
                            const Parallel_Orbitals& pv,
                            std::vector<torch::Tensor>& pdm)

{
    ModuleBase::TITLE("DeePKS_domain", "cal_pdm");
    ModuleBase::timer::tick("DeePKS_domain", "cal_pdm");

    // if pdm has been initialized, skip the calculation
    if (init_pdm)
    {
        init_pdm = false;
        return;
    }

    if (!PARAM.inp.deepks_equiv)
    {
        for (int inl = 0; inl < inlmax; inl++)
        {
            int nm = 2 * inl2l[inl] + 1;
            pdm[inl] = torch::zeros({nm, nm}, torch::kFloat64);
        }
    }
    else
    {
        int pdm_size = 0;
        int nproj = 0;
        for (int il = 0; il < lmaxd + 1; il++)
        {
            nproj += (2 * il + 1) * orb.Alpha[0].getNchi(il);
        }
        pdm_size = nproj * nproj;
        for (int iat = 0; iat < ucell.nat; iat++)
        {
            pdm[iat] = torch::zeros({pdm_size}, torch::kFloat64);
        }
    }

    const double Rcut_Alpha = orb.Alpha[0].getRcut();

#pragma omp parallel for schedule(dynamic)
    for (int iat = 0; iat < ucell.nat; iat++)
    {
        const int T0 = ucell.iat2it[iat];
        const int I0 = ucell.iat2ia[iat];
        Atom* atom0 = &ucell.atoms[T0];
        const ModuleBase::Vector3<double> tau0 = atom0->tau[I0];
        AdjacentAtomInfo adjs;
        GridD.Find_atom(ucell, tau0, T0, I0, &adjs);

        // trace alpha orbital
        std::vector<int> trace_alpha_row;
        std::vector<int> trace_alpha_col;
        if (!PARAM.inp.deepks_equiv)
        {
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
        }
        else
        {
            int nproj = 0;
            for (int il = 0; il < lmaxd + 1; il++)
            {
                nproj += (2 * il + 1) * orb.Alpha[0].getNchi(il);
            }
            for (int iproj = 0; iproj < nproj; iproj++)
            {
                for (int jproj = 0; jproj < nproj; jproj++)
                {
                    trace_alpha_row.push_back(iproj);
                    trace_alpha_col.push_back(jproj);
                }
            }
        }
        const int trace_alpha_size = trace_alpha_row.size();

        for (int ad1 = 0; ad1 < adjs.adj_num + 1; ++ad1)
        {
            const int T1 = adjs.ntype[ad1];
            const int I1 = adjs.natom[ad1];
            const int ibt1 = ucell.itia2iat(T1, I1);
            const ModuleBase::Vector3<double> tau1 = adjs.adjacent_tau[ad1];
            const Atom* atom1 = &ucell.atoms[T1];
            const int nw1_tot = atom1->nw * PARAM.globalv.npol;
            const double Rcut_AO1 = orb.Phi[T1].getRcut();
            const double dist1 = (tau1 - tau0).norm() * ucell.lat0;
            if (dist1 >= Rcut_Alpha + Rcut_AO1)
            {
                continue;
            }

            ModuleBase::Vector3<int> dR1(adjs.box[ad1].x, adjs.box[ad1].y, adjs.box[ad1].z);
            if (phialpha[0]->find_matrix(iat, ibt1, dR1.x, dR1.y, dR1.z) == nullptr)
            {
                continue;
            }

            auto row_indexes = pv.get_indexes_row(ibt1);
            const int row_size = row_indexes.size();
            if (row_size == 0)
            {
                continue;
            }

            // no possible to unexist key
            std::vector<double> s_1t(trace_alpha_size * row_size);
            std::vector<double> g_1dmt(trace_alpha_size * row_size, 0.0);
            for (int irow = 0; irow < row_size; irow++)
            {
                hamilt::BaseMatrix<double>* row_ptr = phialpha[0]->find_matrix(iat, ibt1, dR1);

                for (int i = 0; i < trace_alpha_size; i++)
                {
                    s_1t[i * row_size + irow] = row_ptr->get_value(row_indexes[irow], trace_alpha_row[i]);
                }
            }

            for (int ad2 = 0; ad2 < adjs.adj_num + 1; ad2++)
            {
                const int T2 = adjs.ntype[ad2];
                const int I2 = adjs.natom[ad2];
                const int ibt2 = ucell.itia2iat(T2, I2);
                const ModuleBase::Vector3<double> tau2 = adjs.adjacent_tau[ad2];
                const Atom* atom2 = &ucell.atoms[T2];
                const int nw2_tot = atom2->nw * PARAM.globalv.npol;

                ModuleBase::Vector3<int> dR2(adjs.box[ad2].x, adjs.box[ad2].y, adjs.box[ad2].z);
                if (phialpha[0]->find_matrix(iat, ibt2, dR2.x, dR2.y, dR2.z) == nullptr)
                {
                    continue;
                }

                const double Rcut_AO2 = orb.Phi[T2].getRcut();
                const double dist2 = (tau2 - tau0).norm() * ucell.lat0;

                if (dist2 >= Rcut_Alpha + Rcut_AO2)
                {
                    continue;
                }

                auto col_indexes = pv.get_indexes_col(ibt2);
                const int col_size = col_indexes.size();
                if (col_size == 0)
                {
                    continue;
                }

                std::vector<double> s_2t(trace_alpha_size * col_size);
                // no possible to unexist key
                for (int icol = 0; icol < col_size; icol++)
                {
                    hamilt::BaseMatrix<double>* col_ptr = phialpha[0]->find_matrix(iat, ibt2, dR2);
                    for (int i = 0; i < trace_alpha_size; i++)
                    {
                        s_2t[i * col_size + icol] = col_ptr->get_value(col_indexes[icol], trace_alpha_col[i]);
                    }
                }
                // prepare DM from DMR
                int dRx = 0, dRy = 0, dRz = 0;
                if (std::is_same<TK, std::complex<double>>::value)
                {
                    dRx = dR1.x - dR2.x;
                    dRy = dR1.y - dR2.y;
                    dRz = dR1.z - dR2.z;
                }
                ModuleBase::Vector3<double> dR(dRx, dRy, dRz);

                const double* dm_current = dmr->find_matrix(ibt1, ibt2, dR.x, dR.y, dR.z)->get_pointer();

                // use s_2t and dm_current to get g_1dmt
                // dgemm_: C = alpha * A * B + beta * C
                // C = g_1dmt, A = dm_current, B = s_2t
                // all the input should be data pointer
                constexpr char transa = 'T', transb = 'N';
                const double gemm_alpha = 1.0, gemm_beta = 1.0;
                dgemm_(&transa,
                       &transb,
                       &row_size,
                       &trace_alpha_size,
                       &col_size,
                       &gemm_alpha,
                       dm_current,
                       &col_size,
                       s_2t.data(),
                       &col_size,
                       &gemm_beta,
                       g_1dmt.data(),
                       &row_size);
            } // ad2
            if (!PARAM.inp.deepks_equiv)
            {
                int index = 0, inc = 1;
                for (int L0 = 0; L0 <= orb.Alpha[0].getLmax(); ++L0)
                {
                    for (int N0 = 0; N0 < orb.Alpha[0].getNchi(L0); ++N0)
                    {
                        const int inl = inl_index[T0](I0, L0, N0);
                        const int nm = 2 * L0 + 1;

                        auto accessor = pdm[inl].accessor<double, 2>();
                        for (int m1 = 0; m1 < nm; ++m1) // m1 = 1 for s, 3 for p, 5 for d
                        {
                            for (int m2 = 0; m2 < nm; ++m2) // m1 = 1 for s, 3 for p, 5 for d
                            {
                                accessor[m1][m2] += ddot_(&row_size,
                                                          g_1dmt.data() + index * row_size,
                                                          &inc,
                                                          s_1t.data() + index * row_size,
                                                          &inc);
                                index++;
                            }
                        }
                    }
                }
            }
            else
            {
                auto accessor = pdm[iat].accessor<double, 1>();
                int index = 0, inc = 1;
                int nproj = 0;
                for (int il = 0; il < lmaxd + 1; il++)
                {
                    nproj += (2 * il + 1) * orb.Alpha[0].getNchi(il);
                }
                for (int iproj = 0; iproj < nproj; iproj++)
                {
                    for (int jproj = 0; jproj < nproj; jproj++)
                    {
                        // ddot_: dot product of two vectors
                        // inc means the increment of the index
                        accessor[iproj * nproj + jproj] += ddot_(&row_size,
                                                                 g_1dmt.data() + index * row_size,
                                                                 &inc,
                                                                 s_1t.data() + index * row_size,
                                                                 &inc);
                        index++;
                    }
                }
            }
        } // ad1
    }     // iat

#ifdef __MPI
    for (int inl = 0; inl < inlmax; inl++)
    {
        int pdm_size = (2 * inl2l[inl] + 1) * (2 * inl2l[inl] + 1);
        Parallel_Reduce::reduce_all(pdm[inl].data_ptr<double>(), pdm_size);
    }
#endif
    ModuleBase::timer::tick("DeePKS_domain", "cal_pdm");
    return;
}

void DeePKS_domain::check_pdm(const int inlmax, const std::vector<int>& inl2l, const std::vector<torch::Tensor>& pdm)
{
    const std::string file_projdm = PARAM.globalv.global_out_dir + "pdm.dat";
    std::ofstream ofs(file_projdm.c_str());

    ofs << std::setprecision(10);
    for (int inl = 0; inl < inlmax; inl++)
    {
        const int nm = 2 * inl2l[inl] + 1;
        auto accessor = pdm[inl].accessor<double, 2>();
        for (int m1 = 0; m1 < nm; m1++)
        {
            for (int m2 = 0; m2 < nm; m2++)
            {
                ofs << accessor[m1][m2] << " ";
            }
        }
        ofs << std::endl;
    }
}

template void DeePKS_domain::update_dmr<double>(const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                                                const std::vector<std::vector<double>>& dmk,
                                                const UnitCell& ucell,
                                                const LCAO_Orbitals& orb,
                                                const Parallel_Orbitals& pv,
                                                const Grid_Driver& GridD,
                                                hamilt::HContainer<double>* dmr_deepks);

template void DeePKS_domain::update_dmr<std::complex<double>>(const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                                                              const std::vector<std::vector<std::complex<double>>>& dmk,
                                                              const UnitCell& ucell,
                                                              const LCAO_Orbitals& orb,
                                                              const Parallel_Orbitals& pv,
                                                              const Grid_Driver& GridD,
                                                              hamilt::HContainer<double>* dmr_deepks);

template void DeePKS_domain::cal_pdm<double>(bool& init_pdm,
                                             const int inlmax,
                                             const int lmaxd,
                                             const std::vector<int>& inl2l,
                                             const ModuleBase::IntArray* inl_index,
                                             const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                                             const hamilt::HContainer<double>* dmr,
                                             const std::vector<hamilt::HContainer<double>*> phialpha,
                                             const UnitCell& ucell,
                                             const LCAO_Orbitals& orb,
                                             const Grid_Driver& GridD,
                                             const Parallel_Orbitals& pv,
                                             std::vector<torch::Tensor>& pdm);

template void DeePKS_domain::cal_pdm<std::complex<double>>(bool& init_pdm,
                                                           const int inlmax,
                                                           const int lmaxd,
                                                           const std::vector<int>& inl2l,
                                                           const ModuleBase::IntArray* inl_index,
                                                           const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                                                           const hamilt::HContainer<double>* dmr,
                                                           const std::vector<hamilt::HContainer<double>*> phialpha,
                                                           const UnitCell& ucell,
                                                           const LCAO_Orbitals& orb,
                                                           const Grid_Driver& GridD,
                                                           const Parallel_Orbitals& pv,
                                                           std::vector<torch::Tensor>& pdm);

#endif
