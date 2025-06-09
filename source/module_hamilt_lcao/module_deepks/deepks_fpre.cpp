#ifdef __MLALGO

#include "deepks_fpre.h"

#include "deepks_iterate.h"
#include "module_base/constants.h"
#include "module_base/libm/libm.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
#include "module_base/vector3.h"
#include "module_hamilt_lcao/module_hcontainer/atom_pair.h"
#include "module_parameter/parameter.h"

/// this subroutine calculates the gradient of projected density matrices
/// gdmx_m,m = d/dX sum_{mu,nu} rho_{mu,nu} <chi_mu|alpha_m><alpha_m'|chi_nu>
template <typename TK>
void DeePKS_domain::cal_gdmx(const int lmaxd,
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
                             torch::Tensor& gdmx)
{
    ModuleBase::TITLE("DeePKS_domain", "cal_gdmx");
    ModuleBase::timer::tick("DeePKS_domain", "cal_gdmx");
    // get DS_alpha_mu and S_nu_beta

    int nrow = pv.nrow;
    const int nm = 2 * lmaxd + 1;
    // gdmx: dD/dX
    // \sum_{mu,nu} 2*c_mu*c_nu * <dphi_mu/dx|alpha_m><alpha_m'|phi_nu>
    // size: [3][natom][tot_Inl][2l+1][2l+1]
    gdmx = torch::zeros({3, ucell.nat, inlmax, nm, nm}, torch::dtype(torch::kFloat64));
    auto accessor = gdmx.accessor<double, 5>();

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
            if (overlap_1 == nullptr)
            {
                return; // to next loop
            }

            std::vector<hamilt::BaseMatrix<double>*> grad_overlap_2(3);
            for (int i = 0; i < 3; ++i)
            {
                grad_overlap_2[i] = phialpha[i + 1]->find_matrix(iat, ibt2, dR2);
            }

            assert(overlap_1->get_col_size() == grad_overlap_2[0]->get_col_size());

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
                                    for (int i = 0; i < 3; i++)
                                    {
                                        double value = grad_overlap_2[i]->get_value(col_indexes[iw2], ib + m2)
                                                       * overlap_1->get_value(row_indexes[iw1], ib + m1) * *dm_current;
                                        //(<d/dX chi_mu|alpha_m>)<chi_nu|alpha_m'>
                                        accessor[i][iat][inl][m1][m2] += value;

                                        //(<d/dX chi_nu|alpha_m'>)<chi_mu|alpha_m>
                                        accessor[i][iat][inl][m2][m1] += value;

                                        // (<chi_mu|d/dX alpha_m>)<chi_nu|alpha_m'> = -(<d/dX
                                        // chi_mu|alpha_m>)<chi_nu|alpha_m'>
                                        accessor[i][ibt2][inl][m1][m2] -= value;

                                        //(<chi_nu|d/dX alpha_m'>)<chi_mu|alpha_m> = -(<d/dX
                                        // chi_nu|alpha_m'>)<chi_mu|alpha_m>
                                        accessor[i][ibt2][inl][m2][m1] -= value;
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
    Parallel_Reduce::reduce_all(gdmx.data_ptr<double>(), 3 * ucell.nat * inlmax * nm * nm);
#endif
    ModuleBase::timer::tick("DeePKS_domain", "cal_gdmx");
    return;
}

void DeePKS_domain::check_gdmx(const torch::Tensor& gdmx)
{
    std::stringstream ss;
    std::ofstream ofs_x;
    std::ofstream ofs_y;
    std::ofstream ofs_z;

    ofs_x << std::setprecision(10);
    ofs_y << std::setprecision(10);
    ofs_z << std::setprecision(10);

    // size: [3][natom][inlmax][nm][nm]
    auto size = gdmx.sizes();
    auto accessor = gdmx.accessor<double, 5>();
    for (int ia = 0; ia < size[1]; ia++)
    {
        ss.str("");
        ss << "gdmx_" << ia << ".dat";
        ofs_x.open(ss.str().c_str());
        ss.str("");
        ss << "gdmy_" << ia << ".dat";
        ofs_y.open(ss.str().c_str());
        ss.str("");
        ss << "gdmz_" << ia << ".dat";
        ofs_z.open(ss.str().c_str());

        for (int inl = 0; inl < size[2]; inl++)
        {
            for (int m1 = 0; m1 < size[3]; m1++)
            {
                for (int m2 = 0; m2 < size[4]; m2++)
                {
                    ofs_x << accessor[0][ia][inl][m1][m2] << " ";
                    ofs_y << accessor[1][ia][inl][m1][m2] << " ";
                    ofs_z << accessor[2][ia][inl][m1][m2] << " ";
                }
            }
            ofs_x << std::endl;
            ofs_y << std::endl;
            ofs_z << std::endl;
        }
        ofs_x.close();
        ofs_y.close();
        ofs_z.close();
    }
}

// calculates gradient of descriptors from gradient of projected density matrices
void DeePKS_domain::cal_gvx(const int nat,
                            const int inlmax,
                            const int des_per_atom,
                            const std::vector<int>& inl2l,
                            const std::vector<torch::Tensor>& gevdm,
                            const torch::Tensor& gdmx,
                            torch::Tensor& gvx,
                            const int rank)
{
    ModuleBase::TITLE("DeePKS_domain", "cal_gvx");
    ModuleBase::timer::tick("DeePKS_domain", "cal_gvx");
    // gdmr : nat(derivative) * 3 * inl(projector) * nm * nm
    std::vector<torch::Tensor> gdmr;
    auto accessor = gdmx.accessor<double, 5>();

    if (rank == 0)
    {
        // make gdmx as tensor
        int nlmax = inlmax / nat;
        for (int nl = 0; nl < nlmax; ++nl)
        {
            int nm = 2 * inl2l[nl] + 1;
            torch::Tensor gdmx_sliced
                = gdmx.slice(2, nl, inlmax, nlmax).slice(3, 0, nm, 1).slice(4, 0, nm, 1).permute({1, 0, 2, 3, 4});
            gdmr.push_back(gdmx_sliced);
        }

        assert(gdmr.size() == nlmax);

        // einsum for each inl:
        // gdmr : b:nat(derivative) * x:3 * a:inl(projector) * m:nm * n:nm
        // gevdm : a:inl * v:nm (descriptor) * m:nm (pdm, dim1) * n:nm (pdm, dim2)
        // gvx_vector : b:nat(derivative) * x:3 * a:inl(projector) * m:nm(descriptor)
        std::vector<torch::Tensor> gvx_vector;
        for (int nl = 0; nl < nlmax; ++nl)
        {
            gvx_vector.push_back(at::einsum("bxamn, avmn->bxav", {gdmr[nl], gevdm[nl]}));
        }

        // cat nv-> \sum_nl(nv) = \sum_nl(nm_nl)=des_per_atom
        // concatenate index a(inl) and m(nm)
        // gvx:d(d)/dX, size: [natom][3][natom][des_per_atom]
        gvx = torch::cat(gvx_vector, -1);

        assert(gvx.size(0) == nat);
        assert(gvx.size(1) == 3);
        assert(gvx.size(2) == nat);
        assert(gvx.size(3) == des_per_atom);
    }
    ModuleBase::timer::tick("DeePKS_domain", "cal_gvx");
    return;
}

void DeePKS_domain::check_gvx(const torch::Tensor& gvx, const int rank)
{
    std::stringstream ss;
    std::ofstream ofs_x;
    std::ofstream ofs_y;
    std::ofstream ofs_z;

    if (rank != 0)
    {
        return;
    }

    auto size = gvx.sizes();
    auto accessor = gvx.accessor<double, 4>();

    for (int ia = 0; ia < size[0]; ia++)
    {
        ss.str("");
        ss << "gvx_" << ia << ".dat";
        ofs_x.open(ss.str().c_str());
        ss.str("");
        ss << "gvy_" << ia << ".dat";
        ofs_y.open(ss.str().c_str());
        ss.str("");
        ss << "gvz_" << ia << ".dat";
        ofs_z.open(ss.str().c_str());

        ofs_x << std::setprecision(10);
        ofs_y << std::setprecision(10);
        ofs_z << std::setprecision(10);

        for (int ib = 0; ib < size[2]; ib++)
        {
            for (int nlm = 0; nlm < size[3]; nlm++)
            {
                ofs_x << accessor[ia][0][ib][nlm] << " ";
                ofs_y << accessor[ia][1][ib][nlm] << " ";
                ofs_z << accessor[ia][2][ib][nlm] << " ";
            }
            ofs_x << std::endl;
            ofs_y << std::endl;
            ofs_z << std::endl;
        }
        ofs_x.close();
        ofs_y.close();
        ofs_z.close();
    }
}

template void DeePKS_domain::cal_gdmx<double>(const int lmaxd,
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
                                              torch::Tensor& gdmx);

template void DeePKS_domain::cal_gdmx<std::complex<double>>(const int lmaxd,
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
                                                            torch::Tensor& gdmx);

#endif
