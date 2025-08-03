//  prepare_phialpha_r : prepare phialpha_r for outputting npy file

#ifdef __MLALGO

#include "deepks_vdrpre.h"

#include "LCAO_deepks_io.h" // mohan add 2024-07-22
#include "deepks_iterate.h"
#include "source_base/module_external/blas_connector.h"
#include "source_base/constants.h"
#include "source_base/libm/libm.h"
#include "source_base/parallel_reduce.h"
#include "source_lcao/module_hcontainer/atom_pair.h"
#include "source_io/module_parameter/parameter.h"

void DeePKS_domain::prepare_phialpha_r(const int nlocal,
                                       const int lmaxd,
                                       const int inlmax,
                                       const int nat,
                                       const int R_size,
                                       const std::vector<hamilt::HContainer<double>*> phialpha,
                                       const UnitCell& ucell,
                                       const LCAO_Orbitals& orb,
                                       const Parallel_Orbitals& pv,
                                       const Grid_Driver& GridD,
                                       torch::Tensor& phialpha_r_out)
{
    ModuleBase::TITLE("DeePKS_domain", "prepare_phialpha_r");
    ModuleBase::timer::tick("DeePKS_domain", "prepare_phialpha_r");
    constexpr torch::Dtype dtype = torch::kFloat64;
    int nlmax = inlmax / nat;
    int mmax = 2 * lmaxd + 1;

    phialpha_r_out = torch::zeros({R_size, R_size, R_size, nat, nlmax, nlocal, mmax}, dtype);
    auto accessor = phialpha_r_out.accessor<double, 7>();

    DeePKS_domain::iterate_ad1(
        ucell,
        GridD,
        orb,
        false, // no trace_alpha
        [&](const int iat,
            const ModuleBase::Vector3<double>& tau0,
            const int ibt,
            const ModuleBase::Vector3<double>& tau,
            const int start,
            const int nw_tot,
            ModuleBase::Vector3<int> dR)
        {
            if (phialpha[0]->find_matrix(iat, ibt, dR.x, dR.y, dR.z) == nullptr)
            {
                return; // to next loop
            }

            // middle loop : all atomic basis on the adjacent atom ad
            for (int iw1 = 0; iw1 < nw_tot; ++iw1)
            {
                const int iw1_all = start + iw1;
                const int iw1_local = pv.global2local_row(iw1_all);
                const int iw2_local = pv.global2local_col(iw1_all);
                if (iw1_local < 0 || iw2_local < 0)
                {
                    continue;
                }
                hamilt::BaseMatrix<double>* overlap = phialpha[0]->find_matrix(iat, ibt, dR);
                const int iR = phialpha[0]->find_R(dR);

                int ib = 0;
                int nl = 0;
                for (int L0 = 0; L0 <= orb.Alpha[0].getLmax(); ++L0)
                {
                    for (int N0 = 0; N0 < orb.Alpha[0].getNchi(L0); ++N0)
                    {
                        const int nm = 2 * L0 + 1;
                        for (int m1 = 0; m1 < nm; ++m1) // nm = 1 for s, 3 for p, 5 for d
                        {
                            int iRx = DeePKS_domain::mapping_R(dR.x);
                            int iRy = DeePKS_domain::mapping_R(dR.y);
                            int iRz = DeePKS_domain::mapping_R(dR.z);
                            accessor[iRx][iRy][iRz][iat][nl][iw1_all][m1]
                                += overlap->get_value(iw1, ib + m1);
                        }
                        ib += nm;
                        nl++;
                    }
                }
            } // end iw
        }
    );

#ifdef __MPI
    int size = R_size * R_size * R_size * nat * nlmax * nlocal * mmax;
    double* data_ptr = phialpha_r_out.data_ptr<double>();
    Parallel_Reduce::reduce_all(data_ptr, size);

#endif

    ModuleBase::timer::tick("DeePKS_domain", "prepare_phialpha_r");
    return;
}

void DeePKS_domain::cal_vdr_precalc(const int nlocal,
                                    const int lmaxd,
                                    const int inlmax,
                                    const int nat,
                                    const int nks,
                                    const int R_size,
                                    const std::vector<int>& inl2l,
                                    const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                                    const std::vector<hamilt::HContainer<double>*> phialpha,
                                    const std::vector<torch::Tensor> gevdm,
                                    const ModuleBase::IntArray* inl_index,
                                    const UnitCell& ucell,
                                    const LCAO_Orbitals& orb,
                                    const Parallel_Orbitals& pv,
                                    const Grid_Driver& GridD,
                                    torch::Tensor& vdr_precalc)
{
    ModuleBase::TITLE("DeePKS_domain", "calc_vdr_precalc");
    ModuleBase::timer::tick("DeePKS_domain", "calc_vdr_precalc");

    torch::Tensor vdr_pdm
        = torch::zeros({R_size, R_size, R_size, nlocal, nlocal, inlmax, (2 * lmaxd + 1), (2 * lmaxd + 1)},
                       torch::TensorOptions().dtype(torch::kFloat64));
    auto accessor = vdr_pdm.accessor<double, 8>();

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
            const int T0 = ucell.iat2it[iat];
            const int I0 = ucell.iat2ia[iat];
            if (phialpha[0]->find_matrix(iat, ibt1, dR1.x, dR1.y, dR1.z) == nullptr
                || phialpha[0]->find_matrix(iat, ibt2, dR2.x, dR2.y, dR2.z) == nullptr)
            {
                return; // to next loop
            }

            hamilt::BaseMatrix<double>* overlap_1 = phialpha[0]->find_matrix(iat, ibt1, dR1);
            hamilt::BaseMatrix<double>* overlap_2 = phialpha[0]->find_matrix(iat, ibt2, dR2);
            assert(overlap_1->get_col_size() == overlap_2->get_col_size());
            ModuleBase::Vector3<int> dR = dR1 - dR2;
            int iRx = DeePKS_domain::mapping_R(dR.x);
            int iRy = DeePKS_domain::mapping_R(dR.y);
            int iRz = DeePKS_domain::mapping_R(dR.z);

            for (int iw1 = 0; iw1 < nw1_tot; ++iw1)
            {
                const int iw1_all = start1 + iw1; // this is \mu
                const int iw1_local = pv.global2local_row(iw1_all);
                if (iw1_local < 0)
                {
                    continue;
                }
                for (int iw2 = 0; iw2 < nw2_tot; ++iw2)
                {
                    const int iw2_all = start2 + iw2; // this is \nu
                    const int iw2_local = pv.global2local_col(iw2_all);
                    if (iw2_local < 0)
                    {
                        continue;
                    }

                    int ib = 0;
                    for (int L0 = 0; L0 <= orb.Alpha[0].getLmax(); ++L0)
                    {
                        for (int N0 = 0; N0 < orb.Alpha[0].getNchi(L0); ++N0)
                        {
                            const int inl = inl_index[T0](I0, L0, N0);
                            const int nm = 2 * L0 + 1;

                            for (int m1 = 0; m1 < nm; ++m1) // nm = 1 for s, 3 for p, 5 for d
                            {
                                for (int m2 = 0; m2 < nm; ++m2) // nm = 1 for s, 3 for p, 5 for d
                                {
                                    double tmp = overlap_1->get_value(iw1, ib + m1)
                                                * overlap_2->get_value(iw2, ib + m2);
                                    accessor[iRx][iRy][iRz][iw1_all][iw2_all][inl][m1][m2]
                                        += tmp;
                                }
                            }
                            ib += nm;
                        }
                    }
                } // iw2
            }     // iw1
        }
    );

#ifdef __MPI
    const int size = R_size * R_size * R_size * nlocal * nlocal * inlmax * (2 * lmaxd + 1) * (2 * lmaxd + 1);
    double* data_ptr = vdr_pdm.data_ptr<double>();
    Parallel_Reduce::reduce_all(data_ptr, size);
#endif

    // transfer v_delta_pdm to v_delta_pdm_vector
    int nlmax = inlmax / nat;
    std::vector<torch::Tensor> vdr_pdm_vector;
    for (int nl = 0; nl < nlmax; ++nl)
    {
        int nm = 2 * inl2l[nl] + 1;
        torch::Tensor vdr_pdm_sliced = vdr_pdm.slice(5, nl, inlmax, nlmax).slice(6, 0, nm, 1).slice(7, 0, nm, 1);
        vdr_pdm_vector.push_back(vdr_pdm_sliced);
    }

    assert(vdr_pdm_vector.size() == nlmax);

    // einsum for each nl:
    std::vector<torch::Tensor> vdr_vector;
    for (int nl = 0; nl < nlmax; ++nl)
    {
        vdr_vector.push_back(at::einsum("pqrxyamn, avmn->pqrxyav", {vdr_pdm_vector[nl], gevdm[nl]}));
    }

    vdr_precalc = torch::cat(vdr_vector, -1);

    ModuleBase::timer::tick("DeePKS_domain", "calc_vdr_precalc");
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
