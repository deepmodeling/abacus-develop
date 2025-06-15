//  prepare_phialpha_r : prepare phialpha_r for outputting npy file

#ifdef __MLALGO

#include "deepks_vdrpre.h"

#include "LCAO_deepks_io.h" // mohan add 2024-07-22
#include "deepks_iterate.h"
#include "module_base/blas_connector.h"
#include "module_base/constants.h"
#include "module_base/libm/libm.h"
#include "module_base/parallel_reduce.h"
#include "module_hamilt_lcao/module_hcontainer/atom_pair.h"
#include "module_parameter/parameter.h"

void DeePKS_domain::prepare_phialpha_r(const int nlocal,
                                     const int lmaxd,
                                     const int inlmax,
                                     const int nat,
                                     const std::vector<hamilt::HContainer<double>*> phialpha,
                                     const UnitCell& ucell,
                                     const LCAO_Orbitals& orb,
                                     const Parallel_Orbitals& pv,
                                     const Grid_Driver& GridD,
                                     torch::Tensor& phialpha_r_out,
                                     torch::Tensor& R_query)
{
    ModuleBase::TITLE("DeePKS_domain", "prepare_phialpha_r");
    ModuleBase::timer::tick("DeePKS_domain", "prepare_phialpha_r");
    constexpr torch::Dtype dtype = torch::kFloat64;
    int nlmax = inlmax / nat;
    int mmax = 2 * lmaxd + 1;
    auto size_R = static_cast<long>(phialpha[0]->size_R_loop());
    phialpha_r_out = torch::zeros({size_R, nat, nlmax, nlocal, mmax}, dtype);
    R_query = torch::zeros({size_R, 3}, torch::kInt32);
    auto accessor = phialpha_r_out.accessor<double, 5>();
    auto R_accessor = R_query.accessor<int, 2>();

    for (int iR = 0; iR < size_R; ++iR)
    {
        phialpha[0]->loop_R(iR, R_accessor[iR][0], R_accessor[iR][1], R_accessor[iR][2]);
    }

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
                            accessor[iR][iat][nl][iw1_all][m1] += overlap->get_value(iw1, ib + m1);
                        }
                        ib += nm;
                        nl++;
                    }
                }
            }     // end iw
        }
    );

#ifdef __MPI
    int size = size_R * nat * nlmax * nlocal * mmax;
    double* data_ptr = phialpha_r_out.data_ptr<double>();
    Parallel_Reduce::reduce_all(data_ptr, size);

#endif

    ModuleBase::timer::tick("DeePKS_domain", "prepare_phialpha_r");
    return;
}
#endif
