#ifdef __DEEPKS

#include "deepks_spre.h"

#include "module_base/constants.h"
#include "module_base/libm/libm.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
#include "module_base/vector3.h"
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
                                const std::vector<std::vector<TK>>& dm,
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

    const double Rcut_Alpha = orb.Alpha[0].getRcut();

    for (int T0 = 0; T0 < ucell.ntype; T0++)
    {
        Atom* atom0 = &ucell.atoms[T0];
        for (int I0 = 0; I0 < atom0->na; I0++)
        {
            const int iat = ucell.itia2iat(T0, I0); // on which alpha is located
            const ModuleBase::Vector3<double> tau0 = atom0->tau[I0];
            GridD.Find_atom(ucell, atom0->tau[I0], T0, I0);

            for (int ad1 = 0; ad1 < GridD.getAdjacentNum() + 1; ++ad1)
            {
                const int T1 = GridD.getType(ad1);
                const int I1 = GridD.getNatom(ad1);
                const int ibt1 = ucell.itia2iat(T1, I1); // on which chi_mu is located
                const int start1 = ucell.itiaiw2iwt(T1, I1, 0);

                const ModuleBase::Vector3<double> tau1 = GridD.getAdjacentTau(ad1);
                const Atom* atom1 = &ucell.atoms[T1];
                const int nw1_tot = atom1->nw * PARAM.globalv.npol;
                const double Rcut_AO1 = orb.Phi[T1].getRcut();

                ModuleBase::Vector3<int> dR1(GridD.getBox(ad1).x, GridD.getBox(ad1).y, GridD.getBox(ad1).z);

                for (int ad2 = 0; ad2 < GridD.getAdjacentNum() + 1; ad2++)
                {
                    const int T2 = GridD.getType(ad2);
                    const int I2 = GridD.getNatom(ad2);
                    const int start2 = ucell.itiaiw2iwt(T2, I2, 0);
                    const int ibt2 = ucell.itia2iat(T2, I2);
                    const ModuleBase::Vector3<double> tau2 = GridD.getAdjacentTau(ad2);
                    const Atom* atom2 = &ucell.atoms[T2];
                    const int nw2_tot = atom2->nw * PARAM.globalv.npol;
                    ModuleBase::Vector3<int> dR2(GridD.getBox(ad2).x, GridD.getBox(ad2).y, GridD.getBox(ad2).z);

                    const double Rcut_AO2 = orb.Phi[T2].getRcut();
                    const double dist1 = (tau1 - tau0).norm() * ucell.lat0;
                    const double dist2 = (tau2 - tau0).norm() * ucell.lat0;

                    if (dist1 > Rcut_Alpha + Rcut_AO1 || dist2 > Rcut_Alpha + Rcut_AO2)
                    {
                        continue;
                    }

                    double r0[3];
                    double r1[3];
                    r1[0] = (tau1.x - tau0.x);
                    r1[1] = (tau1.y - tau0.y);
                    r1[2] = (tau1.z - tau0.z);
                    r0[0] = (tau2.x - tau0.x);
                    r0[1] = (tau2.y - tau0.y);
                    r0[2] = (tau2.z - tau0.z);
                    auto row_indexes = pv.get_indexes_row(ibt1);
                    auto col_indexes = pv.get_indexes_col(ibt2);
                    if (row_indexes.size() * col_indexes.size() == 0)
                    {
                        continue;
                    }

                    double* dm_current;
                    int dRx;
                    int dRy;
                    int dRz;
                    if constexpr (std::is_same<TK, double>::value)
                    {
                        dRx = 0;
                        dRy = 0;
                        dRz = 0;
                    }
                    else
                    {
                        dRx = (dR2 - dR1).x;
                        dRy = (dR2 - dR1).y;
                        dRz = (dR2 - dR1).z;
                    }
                    ModuleBase::Vector3<double> dR(dRx, dRy, dRz);

                    hamilt::AtomPair<double> dm_pair(ibt1, ibt2, dRx, dRy, dRz, &pv);
                    dm_pair.allocate(nullptr, 1);
                    for (int ik = 0; ik < nks; ik++)
                    {
                        TK kphase;
                        if constexpr (std::is_same<TK, double>::value)
                        {
                            kphase = 1.0;
                        }
                        else
                        {
                            const double arg = -(kvec_d[ik] * dR) * ModuleBase::TWO_PI;
                            double sinp, cosp;
                            ModuleBase::libm::sincos(arg, &sinp, &cosp);
                            kphase = TK(cosp, sinp);
                        }
                        if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER(PARAM.inp.ks_solver))
                        {
                            dm_pair.add_from_matrix(dm[ik].data(), pv.get_row_size(), kphase, 1);
                        }
                        else
                        {
                            dm_pair.add_from_matrix(dm[ik].data(), pv.get_col_size(), kphase, 0);
                        }
                    }

                    dm_current = dm_pair.get_pointer();

                    for (int iw1 = 0; iw1 < row_indexes.size(); ++iw1)
                    {
                        for (int iw2 = 0; iw2 < col_indexes.size(); ++iw2)
                        {
                            hamilt::BaseMatrix<double>* overlap_1 = phialpha[0]->find_matrix(iat, ibt1, dR1);
                            hamilt::BaseMatrix<double>* overlap_2 = phialpha[0]->find_matrix(iat, ibt2, dR2);
                            std::vector<hamilt::BaseMatrix<double>*> grad_overlap_1(3);
                            std::vector<hamilt::BaseMatrix<double>*> grad_overlap_2(3);

                            assert(overlap_1->get_col_size() == overlap_2->get_col_size());

                            for (int i = 0; i < 3; ++i)
                            {
                                grad_overlap_1[i] = phialpha[i + 1]->find_matrix(iat, ibt1, dR1);
                                grad_overlap_2[i] = phialpha[i + 1]->find_matrix(iat, ibt2, dR2);
                            }

                            int ib = 0;
                            for (int L0 = 0; L0 <= orb.Alpha[0].getLmax(); ++L0)
                            {
                                for (int N0 = 0; N0 < orb.Alpha[0].getNchi(L0); ++N0)
                                {
                                    const int inl = inl_index[T0](I0, L0, N0);
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
                                                              * overlap_1->get_value(row_indexes[iw1], ib + m1)
                                                              * r0[ipol]);
                                                    accessor[mm][inl][m2][m1]
                                                        += ucell.lat0 * *dm_current
                                                           * (overlap_2->get_value(col_indexes[iw2], ib + m1)
                                                              * grad_overlap_1[jpol]->get_value(row_indexes[iw1],
                                                                                                ib + m2)
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
                }         // ad2
            }             // ad1
        }                 // I0
    }                     // T0

#ifdef __MPI
    Parallel_Reduce::reduce_all(gdmepsl.data_ptr<double>(), 6 * inlmax * nm * nm);
#endif
    ModuleBase::timer::tick("DeePKS_domain", "cal_gdmepsl");
    return;
}

void DeePKS_domain::check_gdmepsl(const torch::Tensor& gdmepsl)
{
    std::stringstream ss;
    std::ofstream ofs;

    ofs << std::setprecision(10);

    // size: [6][inlmax][nm][nm]
    auto size = gdmepsl.sizes();
    auto accessor = gdmepsl.accessor<double, 4>();
    for (int i = 0; i < 6; i++)
    {
        ss.str("");
        ss << "gdmepsl_" << i << ".dat";
        ofs.open(ss.str().c_str());

        for (int inl = 0; inl < size[1]; inl++)
        {
            for (int m1 = 0; m1 < size[2]; m1++)
            {
                for (int m2 = 0; m2 < size[3]; m2++)
                {
                    ofs << accessor[i][inl][m1][m2] << " ";
                }
            }
            ofs << std::endl;
        }
        ofs.close();
    }
}

// calculates stress of descriptors from gradient of projected density matrices
// gv_epsl:d(d)/d\epsilon_{\alpha\beta}, [natom][6][des_per_atom]
void DeePKS_domain::cal_gvepsl(const int nat,
                               const int inlmax,
                               const int des_per_atom,
                               const int* inl_l,
                               const std::vector<torch::Tensor>& gevdm,
                               const torch::Tensor& gdmepsl,
                               torch::Tensor& gvepsl)
{
    ModuleBase::TITLE("DeePKS_domain", "cal_gvepsl");
    // dD/d\epsilon_{\alpha\beta}, tensor vector form of gdmepsl
    std::vector<torch::Tensor> gdmepsl_vector;
    auto accessor = gdmepsl.accessor<double, 4>();
    if (GlobalV::MY_RANK == 0)
    {
        // make gdmx as tensor
        int nlmax = inlmax / nat;
        for (int nl = 0; nl < nlmax; ++nl)
        {
            std::vector<torch::Tensor> bmmv;
            for (int i = 0; i < 6; ++i)
            {
                std::vector<torch::Tensor> ammv;
                for (int iat = 0; iat < nat; ++iat)
                {
                    int inl = iat * nlmax + nl;
                    int nm = 2 * inl_l[inl] + 1;
                    std::vector<double> mmv;
                    for (int m1 = 0; m1 < nm; ++m1)
                    {
                        for (int m2 = 0; m2 < nm; ++m2)
                        {
                            mmv.push_back(accessor[i][inl][m1][m2]);
                        }
                    } // nm^2
                    torch::Tensor mm
                        = torch::tensor(mmv, torch::TensorOptions().dtype(torch::kFloat64)).reshape({nm, nm}); // nm*nm
                    ammv.push_back(mm);
                }
                torch::Tensor bmm = torch::stack(ammv, 0); // nat*nm*nm
                bmmv.push_back(bmm);
            }
            gdmepsl_vector.push_back(torch::stack(bmmv, 0)); // nbt*3*nat*nm*nm
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

    return;
}

void DeePKS_domain::check_gvepsl(const torch::Tensor& gvepsl)
{
    std::stringstream ss;
    std::ofstream ofs;

    if (GlobalV::MY_RANK != 0)
    {
        return;
    }

    auto size = gvepsl.sizes();
    auto accessor = gvepsl.accessor<double, 3>();

    for (int i = 0; i < 6; i++)
    {
        ss.str("");
        ss << "gvepsl_" << i << ".dat";
        ofs.open(ss.str().c_str());

        ofs << std::setprecision(10);

        for (int ia = 0; ia < size[1]; ia++)
        {
            for (int nlm = 0; nlm < size[2]; nlm++)
            {
                ofs << accessor[i][ia][nlm] << " ";
            }
            ofs << std::endl;
        }
        ofs.close();
    }
}

template void DeePKS_domain::cal_gdmepsl<double>(const int lmaxd,
                                                 const int inlmax,
                                                 const int nks,
                                                 const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                                                 std::vector<hamilt::HContainer<double>*> phialpha,
                                                 const ModuleBase::IntArray* inl_index,
                                                 const std::vector<std::vector<double>>& dm,
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
                                                               const std::vector<std::vector<std::complex<double>>>& dm,
                                                               const UnitCell& ucell,
                                                               const LCAO_Orbitals& orb,
                                                               const Parallel_Orbitals& pv,
                                                               const Grid_Driver& GridD,
                                                               torch::Tensor& gdmepsl);

#endif
