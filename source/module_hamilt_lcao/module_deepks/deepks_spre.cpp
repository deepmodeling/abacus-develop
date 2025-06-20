#ifdef __MLALGO

#include "deepks_spre.h"

#include "deepks_iterate.h"
#include "source_base/constants.h"
#include "source_base/libm/libm.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_base/vector3.h"
#include "module_hamilt_lcao/module_hcontainer/atom_pair.h"
#include "module_parameter/parameter.h"

/// this subroutine calculates the gradient of PDM wrt strain tensor:
/// gdmepsl = d/d\epsilon_{ab} *
///           sum_{mu,nu} rho_{mu,nu} <chi_mu|alpha_m><alpha_m'|chi_nu>
template <typename TK>
void DeePKS_domain::cal_gdmepsl(const int lmaxd,
                                const int inlmax,
                                const int nks,
                                const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                                std::vector<hamilt::HContainer<double>*> phialpha,
                                const ModuleBase::IntArray* inl_index,
                                const hamilt::HContainer<double>* dmr,
                                const UnitCell& ucell,
                                const LCAO_Orbitals& orb,
                                const Parallel_Orbitals& pv,
                                const Grid_Driver& GridD,
                                torch::Tensor& gdmepsl)
{
    ModuleBase::TITLE("DeePKS_domain", "cal_gdmepsl");
    ModuleBase::timer::tick("DeePKS_domain", "cal_gdmepsl");
    // get DS_alpha_mu and S_nu_beta

    int nrow = pv.nrow;
    const int nm = 2 * lmaxd + 1;
    // gdmepsl: dD/d\epsilon_{\alpha\beta}
    // size: [6][tot_Inl][2l+1][2l+1]
    gdmepsl = torch::zeros({6, inlmax, nm, nm}, torch::dtype(torch::kFloat64));
    auto accessor = gdmepsl.accessor<double, 4>();

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
            double r1[3] = {0, 0, 0};
            double r2[3] = {0, 0, 0};
            r1[0] = (tau1.x - tau0.x);
            r1[1] = (tau1.y - tau0.y);
            r1[2] = (tau1.z - tau0.z);
            r2[0] = (tau2.x - tau0.x);
            r2[1] = (tau2.y - tau0.y);
            r2[2] = (tau2.z - tau0.z);
            auto row_indexes = pv.get_indexes_row(ibt1);
            auto col_indexes = pv.get_indexes_col(ibt2);
            if (row_indexes.size() * col_indexes.size() == 0)
            {
                return; // to next loop
            }

            int dRx = 0;
            int dRy = 0;
            int dRz = 0;
            if (std::is_same<TK, std::complex<double>>::value)
            {
                dRx = (dR1 - dR2).x;
                dRy = (dR1 - dR2).y;
                dRz = (dR1 - dR2).z;
            }
            ModuleBase::Vector3<double> dR(dRx, dRy, dRz);

            const double* dm_current = dmr->find_matrix(ibt1, ibt2, dR.x, dR.y, dR.z)->get_pointer();

            hamilt::BaseMatrix<double>* overlap_1 = phialpha[0]->find_matrix(iat, ibt1, dR1);
            hamilt::BaseMatrix<double>* overlap_2 = phialpha[0]->find_matrix(iat, ibt2, dR2);
            if (overlap_1 == nullptr || overlap_2 == nullptr)
            {
                return; // to next loop
            }
            std::vector<hamilt::BaseMatrix<double>*> grad_overlap_1(3);
            std::vector<hamilt::BaseMatrix<double>*> grad_overlap_2(3);

            assert(overlap_1->get_col_size() == overlap_2->get_col_size());

            for (int i = 0; i < 3; ++i)
            {
                grad_overlap_1[i] = phialpha[i + 1]->find_matrix(iat, ibt1, dR1);
                grad_overlap_2[i] = phialpha[i + 1]->find_matrix(iat, ibt2, dR2);
            }

            for (int iw1 = 0; iw1 < row_indexes.size(); ++iw1)
            {
                for (int iw2 = 0; iw2 < col_indexes.size(); ++iw2)
                {
                    int ib = 0;
                    for (int L0 = 0; L0 <= orb.Alpha[0].getLmax(); ++L0)
                    {
                        for (int N0 = 0; N0 < orb.Alpha[0].getNchi(L0); ++N0)
                        {
                            const int inl = inl_index[ucell.iat2it[iat]](ucell.iat2ia[iat], L0, N0);
                            const int nm = 2 * L0 + 1;
                            for (int m1 = 0; m1 < nm; ++m1)
                            {
                                for (int m2 = 0; m2 < nm; ++m2)
                                {
                                    int mm = 0;
                                    for (int ipol = 0; ipol < 3; ipol++)
                                    {
                                        for (int jpol = ipol; jpol < 3; jpol++)
                                        {
                                            accessor[mm][inl][m2][m1]
                                                += ucell.lat0 * *dm_current
                                                   * (grad_overlap_2[jpol]->get_value(col_indexes[iw2], ib + m2)
                                                      * overlap_1->get_value(row_indexes[iw1], ib + m1) * r2[ipol]);
                                            accessor[mm][inl][m2][m1]
                                                += ucell.lat0 * *dm_current
                                                   * (overlap_2->get_value(col_indexes[iw2], ib + m1)
                                                      * grad_overlap_1[jpol]->get_value(row_indexes[iw1], ib + m2)
                                                      * r1[ipol]);
                                            mm++;
                                        }
                                    }
                                }
                            }
                            ib += nm;
                        }
                    }
                    assert(ib == overlap_1->get_col_size());
                    dm_current++;
                } // iw2
            }     // iw1
        }
    );

#ifdef __MPI
    Parallel_Reduce::reduce_all(gdmepsl.data_ptr<double>(), 6 * inlmax * nm * nm);
#endif
    ModuleBase::timer::tick("DeePKS_domain", "cal_gdmepsl");
    return;
}

// calculates stress of descriptors from gradient of projected density matrices
// gv_epsl:d(d)/d\epsilon_{\alpha\beta}, [natom][6][des_per_atom]
void DeePKS_domain::cal_gvepsl(const int nat,
                               const int inlmax,
                               const int des_per_atom,
                               const std::vector<int>& inl2l,
                               const std::vector<torch::Tensor>& gevdm,
                               const torch::Tensor& gdmepsl,
                               torch::Tensor& gvepsl,
                               const int rank)
{
    ModuleBase::TITLE("DeePKS_domain", "cal_gvepsl");
    ModuleBase::timer::tick("DeePKS_domain", "cal_gvepsl");
    // dD/d\epsilon_{\alpha\beta}, tensor vector form of gdmepsl
    std::vector<torch::Tensor> gdmepsl_vector;
    auto accessor = gdmepsl.accessor<double, 4>();
    if (rank == 0)
    {
        // make gdmepsl as tensor
        int nlmax = inlmax / nat;
        for (int nl = 0; nl < nlmax; ++nl)
        {
            int nm = 2 * inl2l[nl] + 1;
            torch::Tensor gdmepsl_sliced = gdmepsl.slice(1, nl, inlmax, nlmax).slice(2, 0, nm, 1).slice(3, 0, nm, 1);
            gdmepsl_vector.push_back(gdmepsl_sliced);
        }
        assert(gdmepsl_vector.size() == nlmax);

        // einsum for each inl:
        // gdmepsl_vector : b:npol * a:inl(projector) * m:nm * n:nm
        // gevdm : a:inl * v:nm (descriptor) * m:nm (pdm, dim1) * n:nm
        // (pdm, dim2) gvepsl_vector : b:npol * a:inl(projector) *
        // m:nm(descriptor)
        std::vector<torch::Tensor> gvepsl_vector;
        for (int nl = 0; nl < nlmax; ++nl)
        {
            gvepsl_vector.push_back(at::einsum("bamn, avmn->bav", {gdmepsl_vector[nl], gevdm[nl]}));
        }

        // cat nv-> \sum_nl(nv) = \sum_nl(nm_nl)=des_per_atom
        // concatenate index a(inl) and m(nm)
        gvepsl = torch::cat(gvepsl_vector, -1);
        assert(gvepsl.size(0) == 6);
        assert(gvepsl.size(1) == nat);
        assert(gvepsl.size(2) == des_per_atom);
    }

    ModuleBase::timer::tick("DeePKS_domain", "cal_gvepsl");
    return;
}

template void DeePKS_domain::cal_gdmepsl<double>(const int lmaxd,
                                                 const int inlmax,
                                                 const int nks,
                                                 const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                                                 std::vector<hamilt::HContainer<double>*> phialpha,
                                                 const ModuleBase::IntArray* inl_index,
                                                 const hamilt::HContainer<double>* dmr,
                                                 const UnitCell& ucell,
                                                 const LCAO_Orbitals& orb,
                                                 const Parallel_Orbitals& pv,
                                                 const Grid_Driver& GridD,
                                                 torch::Tensor& gdmepsl);

template void DeePKS_domain::cal_gdmepsl<std::complex<double>>(const int lmaxd,
                                                               const int inlmax,
                                                               const int nks,
                                                               const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                                                               std::vector<hamilt::HContainer<double>*> phialpha,
                                                               const ModuleBase::IntArray* inl_index,
                                                               const hamilt::HContainer<double>* dmr,
                                                               const UnitCell& ucell,
                                                               const LCAO_Orbitals& orb,
                                                               const Parallel_Orbitals& pv,
                                                               const Grid_Driver& GridD,
                                                               torch::Tensor& gdmepsl);

#endif
