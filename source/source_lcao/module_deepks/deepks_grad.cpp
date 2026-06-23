// deepks_grad.cpp  --  Descriptor-gradient label building blocks
// See deepks_grad.h for the full derivation and Python usage.

#ifdef __MLALGO

#include "deepks_grad.h"

#include "deepks_iterate.h"
#include "deepks_vdrpre.h"
#include "source_base/parallel_reduce.h"

// ---------------------------------------------------------------------------
// Step 1: dot_phialpha_hamilt[inl, m1, m2]
//         = <phialpha^I_{nl,m1} | H | phialpha^I_{nl,m2}>
// ---------------------------------------------------------------------------
template <typename TR>
void DeePKS_domain::cal_phialpha_hamilt_proj(const int nlocal,
                                             const int nat,
                                             const DeePKS_Param& deepks_param,
                                             const std::vector<hamilt::HContainer<double>*>& phialpha,
                                             const hamilt::HContainer<TR>& hR,
                                             const UnitCell& ucell,
                                             const LCAO_Orbitals& orb,
                                             const Parallel_Orbitals& pv,
                                             const Grid_Driver& GridD,
                                             torch::Tensor& dot_phialpha_hamilt)
{
    ModuleBase::TITLE("DeePKS_domain", "cal_phialpha_hamilt_proj");
    ModuleBase::timer::start("DeePKS_domain", "cal_phialpha_hamilt_proj");

    using TorchScalar = typename std::conditional<std::is_same<TR, double>::value, double, c10::complex<double>>::type;

    const int inlmax = deepks_param.inlmax;
    const int lmaxd = deepks_param.lmaxd;
    const int nm_max = 2 * lmaxd + 1;

    const torch::Dtype dtype = std::is_same<TR, double>::value ? torch::kFloat64 : torch::kComplexDouble;
    torch::Tensor dot_ph = torch::zeros({inlmax, nm_max, nm_max}, torch::TensorOptions().dtype(dtype));
    auto dot_ph_acc = dot_ph.accessor<TorchScalar, 3>();

    DeePKS_domain::iterate_ad2(
        ucell,
        GridD,
        orb,
        false,
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
            ModuleBase::Vector3<int> dR2) {
            if (phialpha[0]->find_matrix(iat, ibt1, dR1.x, dR1.y, dR1.z) == nullptr
                || phialpha[0]->find_matrix(iat, ibt2, dR2.x, dR2.y, dR2.z) == nullptr)
            {
                return;
            }

            const ModuleBase::Vector3<int> R_H = dR2 - dR1;

            const auto row_indexes = pv.get_indexes_row(ibt1);
            const auto col_indexes = pv.get_indexes_col(ibt2);
            const int row_size = static_cast<int>(row_indexes.size());
            const int col_size = static_cast<int>(col_indexes.size());
            if (row_size == 0 || col_size == 0)
            {
                return;
            }

            const hamilt::BaseMatrix<TR>* h_mat = hR.find_matrix(ibt1, ibt2, R_H.x, R_H.y, R_H.z);
            if (h_mat == nullptr)
            {
                return;
            }
            const TR* h_data = h_mat->get_pointer();

            hamilt::BaseMatrix<double>* pa1 = phialpha[0]->find_matrix(iat, ibt1, dR1);
            hamilt::BaseMatrix<double>* pa2 = phialpha[0]->find_matrix(iat, ibt2, dR2);

            const int T0 = ucell.iat2it[iat];
            const int I0 = ucell.iat2ia[iat];

            int ib = 0;
            for (int L0 = 0; L0 <= orb.Alpha[0].getLmax(); ++L0)
            {
                for (int N0 = 0; N0 < orb.Alpha[0].getNchi(L0); ++N0)
                {
                    const int inl = deepks_param.inl_index[T0](I0, L0, N0);
                    const int nm = 2 * L0 + 1;

                    for (int irow = 0; irow < row_size; ++irow)
                    {
                        for (int icol = 0; icol < col_size; ++icol)
                        {
                            const TR h_val = h_data[irow * col_size + icol];
                            if (h_val == TR{})
                            {
                                continue;
                            }
                            for (int m1 = 0; m1 < nm; ++m1)
                            {
                                const double pa1_val = pa1->get_value(row_indexes[irow], ib + m1);
                                for (int m2 = 0; m2 < nm; ++m2)
                                {
                                    dot_ph_acc[inl][m1][m2] += static_cast<TorchScalar>(
                                        pa1_val * h_val * pa2->get_value(col_indexes[icol], ib + m2));
                                }
                            }
                        }
                    }
                    ib += nm;
                }
            }
        });

#ifdef __MPI
    if (std::is_same<TR, double>::value)
    {
        Parallel_Reduce::reduce_all(dot_ph.data_ptr<double>(), dot_ph.numel());
    }
    else
    {
        Parallel_Reduce::reduce_all(reinterpret_cast<double*>(dot_ph.data_ptr<c10::complex<double>>()),
                                    2 * dot_ph.numel());
    }
#endif

    dot_phialpha_hamilt = dot_ph;

    ModuleBase::timer::end("DeePKS_domain", "cal_phialpha_hamilt_proj");
}

// ---------------------------------------------------------------------------
void DeePKS_domain::cal_vdrpre_square(const int nlocal,
                                      const int nat,
                                      const int R_size,
                                      const DeePKS_Param& deepks_param,
                                      const std::vector<hamilt::HContainer<double>*>& phialpha,
                                      const std::vector<torch::Tensor>& gevdm,
                                      const UnitCell& ucell,
                                      const LCAO_Orbitals& orb,
                                      const Parallel_Orbitals& pv,
                                      const Grid_Driver& GridD,
                                      torch::Tensor& vdrpre_square)
{
    ModuleBase::TITLE("DeePKS_domain", "cal_vdrpre_square");
    ModuleBase::timer::start("DeePKS_domain", "cal_vdrpre_square");

    torch::Tensor vdr_precalc;
    const std::vector<ModuleBase::Vector3<double>> kvec_d_dummy;
    DeePKS_domain::cal_vdr_precalc(nlocal,
                                   nat,
                                   1,
                                   R_size,
                                   deepks_param,
                                   kvec_d_dummy,
                                   phialpha,
                                   gevdm,
                                   ucell,
                                   orb,
                                   pv,
                                   GridD,
                                   vdr_precalc);

    const int N_alpha = nat * deepks_param.des_per_atom;
    const torch::Tensor flat = vdr_precalc.reshape({-1, N_alpha});
    const int K = static_cast<int>(flat.size(0));
    const double* flat_data = flat.data_ptr<double>();

    // vdr_precalc is sparse: only (R, mu, nu) rows with actual phialpha overlaps
    // are non-zero.  Collect them into a compact matrix M and compute M^T @ M
    // instead of flat^T @ flat, skipping the work for the many zero rows.
    std::vector<int> touched;
    touched.reserve(K);
    for (int k = 0; k < K; ++k)
    {
        const double* row = flat_data + k * N_alpha;
        for (int j = 0; j < N_alpha; ++j)
        {
            if (row[j] != 0.0)
            {
                touched.push_back(k);
                break;
            }
        }
    }

    const int ntouched = static_cast<int>(touched.size());
    torch::Tensor M = torch::empty({ntouched, N_alpha}, torch::TensorOptions().dtype(torch::kFloat64));
    double* M_data = M.data_ptr<double>();
    for (int t = 0; t < ntouched; ++t)
    {
        const double* src = flat_data + touched[t] * N_alpha;
        double* dst = M_data + t * N_alpha;
        for (int j = 0; j < N_alpha; ++j)
        {
            dst[j] = src[j];
        }
    }

    vdrpre_square = M.t().mm(M);

    ModuleBase::timer::end("DeePKS_domain", "cal_vdrpre_square");
}

template void DeePKS_domain::cal_phialpha_hamilt_proj<double>(const int,
                                                              const int,
                                                              const DeePKS_Param&,
                                                              const std::vector<hamilt::HContainer<double>*>&,
                                                              const hamilt::HContainer<double>&,
                                                              const UnitCell&,
                                                              const LCAO_Orbitals&,
                                                              const Parallel_Orbitals&,
                                                              const Grid_Driver&,
                                                              torch::Tensor&);

template void DeePKS_domain::cal_phialpha_hamilt_proj<std::complex<double>>(
    const int,
    const int,
    const DeePKS_Param&,
    const std::vector<hamilt::HContainer<double>*>&,
    const hamilt::HContainer<std::complex<double>>&,
    const UnitCell&,
    const LCAO_Orbitals&,
    const Parallel_Orbitals&,
    const Grid_Driver&,
    torch::Tensor&);

#endif
