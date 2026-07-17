//  prepare_phialpha_iRmat : prepare phialpha_r and iR_mat for outputting npy file

#ifdef __MLALGO

#include "deepks_vdrpre.h"

#include "LCAO_deepks_io.h" // mohan add 2024-07-22
#include "deepks_iterate.h"
#include "source_base/constants.h"
#include "source_base/libm/libm.h"
#include "source_base/module_external/blas_connector.h"
#include "source_base/parallel_reduce.h"
#include "source_io/module_parameter/parameter.h"
#include "source_lcao/module_hcontainer/atom_pair.h"

void DeePKS_domain::prepare_phialpha_iRmat(const int nlocal,
                                           const int R_size,
                                           const DeePKS_Param& deepks_param,
                                           const std::vector<hamilt::HContainer<double>*> phialpha,
                                           const UnitCell& ucell,
                                           const LCAO_Orbitals& orb,
                                           const Parallel_Orbitals& pv,
                                           const Grid_Driver& GridD,
                                           torch::Tensor& overlap,
                                           torch::Tensor& iRmat)
{
    ModuleBase::TITLE("DeePKS_domain", "prepare_phialpha_iRmat");
    ModuleBase::timer::start("DeePKS_domain", "prepare_phialpha_iRmat");
    constexpr torch::Dtype dtype = torch::kFloat64;

    // get the maximum nnmax
    std::vector<int> nnmax_vec(ucell.nat, 0);
    DeePKS_domain::iterate_ad1(ucell,
                               GridD,
                               orb,
                               false, // no trace_alpha
                               [&](const int iat,
                                   const ModuleBase::Vector3<double>& tau0,
                                   const int ibt,
                                   const ModuleBase::Vector3<double>& tau1,
                                   const int start,
                                   const int nw_tot,
                                   ModuleBase::Vector3<int> dR) {
                                   if (phialpha[0]->find_matrix(iat, ibt, dR.x, dR.y, dR.z) == nullptr)
                                   {
                                       return; // to next loop
                                   }
                                   nnmax_vec[iat]++;
                               });

    int nnmax = *std::max_element(nnmax_vec.begin(), nnmax_vec.end());
    overlap = torch::zeros({ucell.nat, nnmax, nlocal, deepks_param.des_per_atom}, dtype);
    torch::Tensor dRmat_tmp = torch::zeros({ucell.nat, nnmax, 3}, torch::kInt32);
    auto overlap_accessor = overlap.accessor<double, 4>();
    auto dRmat_accessor = dRmat_tmp.accessor<int, 3>();

    std::fill(nnmax_vec.begin(), nnmax_vec.end(), 0);
    DeePKS_domain::iterate_ad1(ucell,
                               GridD,
                               orb,
                               false, // no trace_alpha
                               [&](const int iat,
                                   const ModuleBase::Vector3<double>& tau0,
                                   const int ibt,
                                   const ModuleBase::Vector3<double>& tau1,
                                   const int start,
                                   const int nw_tot,
                                   ModuleBase::Vector3<int> dR) {
                                   hamilt::BaseMatrix<double>* overlap_mat = phialpha[0]->find_matrix(iat, ibt, dR);
                                   if (overlap_mat == nullptr)
                                   {
                                       return; // to next loop
                                   }
                                   dRmat_accessor[iat][nnmax_vec[iat]][0] = dR.x;
                                   dRmat_accessor[iat][nnmax_vec[iat]][1] = dR.y;
                                   dRmat_accessor[iat][nnmax_vec[iat]][2] = dR.z;

                                   for (int ix = 0; ix < nw_tot; ix++)
                                   {
                                       if (pv.global2local_row(start + ix) < 0 || pv.global2local_col(start + ix) < 0)
                                       {
                                           continue;
                                       }
                                       for (int iy = 0; iy < deepks_param.des_per_atom; iy++)
                                       {
                                           overlap_accessor[iat][nnmax_vec[iat]][start + ix][iy]
                                               = overlap_mat->get_value(ix, iy);
                                       }
                                   }
                                   nnmax_vec[iat]++;
                               });
#ifdef __MPI
    Parallel_Reduce::reduce_all(overlap.data_ptr<double>(), overlap.numel());
#endif
    iRmat = mapping_R(dRmat_tmp.unsqueeze(1) - dRmat_tmp.unsqueeze(2));
    ModuleBase::timer::end("DeePKS_domain", "prepare_phialpha_iRmat");
    return;
}

void DeePKS_domain::cal_vdr_precalc(const int nlocal,
                                    const int nat,
                                    const int nks,
                                    const int R_size,
                                    const DeePKS_Param& deepks_param,
                                    const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                                    const std::vector<hamilt::HContainer<double>*> phialpha,
                                    const std::vector<torch::Tensor> gevdm,
                                    const UnitCell& ucell,
                                    const LCAO_Orbitals& orb,
                                    const Parallel_Orbitals& pv,
                                    const Grid_Driver& GridD,
                                    torch::Tensor& vdr_precalc)
{
    ModuleBase::TITLE("DeePKS_domain", "calc_vdr_precalc");
    ModuleBase::timer::start("DeePKS_domain", "calc_vdr_precalc");

    const int des_per_atom = deepks_param.des_per_atom;

    // Shape: (R_size, R_size, R_size, nlocal, nlocal, nat, des_per_atom)
    torch::Tensor vdr_prec = torch::zeros({R_size, R_size, R_size, nlocal, nlocal, nat, des_per_atom},
                                          torch::TensorOptions().dtype(torch::kFloat64));
    auto accessor = vdr_prec.accessor<double, 7>();

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

            hamilt::BaseMatrix<double>* overlap_1 = phialpha[0]->find_matrix(iat, ibt1, dR1);
            hamilt::BaseMatrix<double>* overlap_2 = phialpha[0]->find_matrix(iat, ibt2, dR2);
            assert(overlap_1->get_col_size() == overlap_2->get_col_size());

            const ModuleBase::Vector3<int> dR = dR2 - dR1;
            const int iRx = DeePKS_domain::mapping_R(dR.x);
            const int iRy = DeePKS_domain::mapping_R(dR.y);
            const int iRz = DeePKS_domain::mapping_R(dR.z);
            if (iRx >= R_size || iRy >= R_size || iRz >= R_size)
            {
                return;
            }

            // Collect MPI-local orbital indices for row (ibt1) and col (ibt2).
            std::vector<int> row_iw, row_iw_all, col_iw, col_iw_all;
            row_iw.reserve(nw1_tot);
            row_iw_all.reserve(nw1_tot);
            col_iw.reserve(nw2_tot);
            col_iw_all.reserve(nw2_tot);
            for (int iw1 = 0; iw1 < nw1_tot; ++iw1)
            {
                if (pv.global2local_row(start1 + iw1) >= 0)
                {
                    row_iw.push_back(iw1);
                    row_iw_all.push_back(start1 + iw1);
                }
            }
            for (int iw2 = 0; iw2 < nw2_tot; ++iw2)
            {
                if (pv.global2local_col(start2 + iw2) >= 0)
                {
                    col_iw.push_back(iw2);
                    col_iw_all.push_back(start2 + iw2);
                }
            }
            const int nrow = static_cast<int>(row_iw.size());
            const int ncol = static_cast<int>(col_iw.size());
            if (nrow == 0 || ncol == 0)
            {
                return;
            }
            const int npairs = nrow * ncol;

            int ib = 0;
            int nl = 0;
            for (int L0 = 0; L0 <= orb.Alpha[0].getLmax(); ++L0)
            {
                for (int N0 = 0; N0 < orb.Alpha[0].getNchi(L0); ++N0, ++nl)
                {
                    const int nm = 2 * L0 + 1;
                    const int nm2 = nm * nm;

                    // T_flat[k1*ncol+k2, m1*nm+m2] = pa1[iw1,ib+m1] * pa2[iw2,ib+m2]
                    std::vector<double> T_flat(npairs * nm2);
                    for (int k1 = 0; k1 < nrow; ++k1)
                    {
                        const int iw1 = row_iw[k1];
                        for (int k2 = 0; k2 < ncol; ++k2)
                        {
                            const int iw2 = col_iw[k2];
                            const int k = k1 * ncol + k2;
                            for (int m1 = 0; m1 < nm; ++m1)
                            {
                                const double a = overlap_1->get_value(iw1, ib + m1);
                                for (int m2 = 0; m2 < nm; ++m2)
                                {
                                    T_flat[k * nm2 + m1 * nm + m2] = a * overlap_2->get_value(iw2, ib + m2);
                                }
                            }
                        }
                    }

                    // gdata[v, m1*nm+m2] = gevdm[nl][iat, v, m1, m2]
                    // gevdm[nl] shape: (nat, nm, nm, nm), contiguous.
                    const double* gdata = gevdm[nl].data_ptr<double>() + iat * nm * nm2;

                    // result[k, v] = sum_m T_flat[k, m] * gdata[v, m]
                    std::vector<double> result(npairs * nm, 0.0);
                    for (int k = 0; k < npairs; ++k)
                    {
                        const double* T_row = T_flat.data() + k * nm2;
                        for (int v = 0; v < nm; ++v)
                        {
                            const double* g_row = gdata + v * nm2;
                            double sum = 0.0;
                            for (int m = 0; m < nm2; ++m)
                            {
                                sum += T_row[m] * g_row[m];
                            }
                            result[k * nm + v] = sum;
                        }
                    }

                    for (int k1 = 0; k1 < nrow; ++k1)
                    {
                        const int iw1_all = row_iw_all[k1];
                        for (int k2 = 0; k2 < ncol; ++k2)
                        {
                            const int iw2_all = col_iw_all[k2];
                            const double* res_row = result.data() + (k1 * ncol + k2) * nm;
                            for (int v = 0; v < nm; ++v)
                            {
                                accessor[iRx][iRy][iRz][iw1_all][iw2_all][iat][ib + v] += res_row[v];
                            }
                        }
                    }

                    ib += nm;
                }
            }
        });

#ifdef __MPI
    Parallel_Reduce::reduce_all(vdr_prec.data_ptr<double>(), vdr_prec.numel());
#endif

    vdr_precalc = vdr_prec;

    ModuleBase::timer::end("DeePKS_domain", "calc_vdr_precalc");
    return;
}

int DeePKS_domain::mapping_R(int R)
{
    // R_index mapping: index(R) = 2R-1 if R > 0, index(R) = -2R if R <= 0
    // after mapping, the new index [0,1,2,3,4,...] -> old index [0,1,-1,2,-2,...]
    // This manipulation makes sure that the new index is natural number
    // which makes it available to be used as index in torch::Tensor
    int R_index = 0;
    if (R > 0)
    {
        R_index = 2 * R - 1;
    }
    else
    {
        R_index = -2 * R;
    }
    return R_index;
}

torch::Tensor DeePKS_domain::mapping_R(const torch::Tensor& R_tensor)
{
    auto R = R_tensor.to(torch::kInt32);
    auto pos = R > 0;
    auto twoR_minus1 = R * 2 - 1;
    auto neg_minus2R = -2 * R;
    return at::where(pos, twoR_minus1, neg_minus2R);
}

template <typename T>
int DeePKS_domain::get_R_size(const hamilt::HContainer<T>& hcontainer)
{
    // get R_size from hcontainer
    int R_size = 0;
    if (hcontainer.size_R_loop() > 0)
    {
        for (int iR = 0; iR < hcontainer.size_R_loop(); ++iR)
        {
            ModuleBase::Vector3<int> R_vec;
            hcontainer.loop_R(iR, R_vec.x, R_vec.y, R_vec.z);
            int R_min = std::min({R_vec.x, R_vec.y, R_vec.z});
            int R_max = std::max({R_vec.x, R_vec.y, R_vec.z});
            int tmp_R_size = std::max(DeePKS_domain::mapping_R(R_min), DeePKS_domain::mapping_R(R_max)) + 1;
            if (tmp_R_size > R_size)
            {
                R_size = tmp_R_size;
            }
        }
    }
    assert(R_size > 0);
    return R_size;
}

template int DeePKS_domain::get_R_size<double>(const hamilt::HContainer<double>& hcontainer);
template int DeePKS_domain::get_R_size<std::complex<double>>(
    const hamilt::HContainer<std::complex<double>>& hcontainer);
#endif
