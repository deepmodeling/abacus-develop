#include "dftu.h"
#include "module_base/timer.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"

namespace ModuleDFTU
{

void DFTU::fold_dSR_gamma(const int dim1, const int dim2, double* dSR_gamma)
{
    ModuleBase::TITLE("DFTU", "fold_dSR_gamma");

    ModuleBase::GlobalFunc::ZEROS(dSR_gamma, this->LM->ParaV->nloc);

    double* dS_ptr;
    if      (dim1 == 0) dS_ptr = this->LM->DSloc_x;
    else if (dim1 == 1) dS_ptr = this->LM->DSloc_y;
    else if (dim1 == 2) dS_ptr = this->LM->DSloc_z;

    int nnr = 0;
    ModuleBase::Vector3<double> tau1, tau2, dtau;
    ModuleBase::Vector3<double> dtau1, dtau2, tau0;
    for (int T1 = 0; T1 < GlobalC::ucell.ntype; ++T1)
    {
        Atom* atom1 = &GlobalC::ucell.atoms[T1];
        for (int I1 = 0; I1 < atom1->na; ++I1)
        {
            tau1 = atom1->tau[I1];
            const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
            GlobalC::GridD.Find_atom(GlobalC::ucell, tau1, T1, I1);
            for (int ad = 0; ad < GlobalC::GridD.getAdjacentNum() + 1; ++ad)
            {
                const int T2 = GlobalC::GridD.getType(ad);
                const int I2 = GlobalC::GridD.getNatom(ad);
                const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);
                Atom* atom2 = &GlobalC::ucell.atoms[T2];
                tau2 = GlobalC::GridD.getAdjacentTau(ad);
                dtau = tau2 - tau1;
                double distance = dtau.norm() * GlobalC::ucell.lat0;
                double rcut = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut();
                bool adj = false;
                if (distance < rcut)
                    adj = true;
                else if (distance >= rcut)
                {
                    for (int ad0 = 0; ad0 < GlobalC::GridD.getAdjacentNum() + 1; ++ad0)
                    {
                        const int T0 = GlobalC::GridD.getType(ad0);
                        const int I0 = GlobalC::GridD.getNatom(ad0);
                        const int iat0 = GlobalC::ucell.itia2iat(T0, I0);
                        const int start0 = GlobalC::ucell.itiaiw2iwt(T0, I0, 0);
                        tau0 = GlobalC::GridD.getAdjacentTau(ad0);
                        dtau1 = tau0 - tau1;
                        dtau2 = tau0 - tau2;
                        double distance1 = dtau1.norm() * GlobalC::ucell.lat0;
                        double distance2 = dtau2.norm() * GlobalC::ucell.lat0;
                        double rcut1 = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();
                        double rcut2 = GlobalC::ORB.Phi[T2].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();
                        if (distance1 < rcut1 && distance2 < rcut2)
                        {
                            adj = true;
                            break;
                        }
                    }
                }

                if (adj)
                {
                    for (int jj = 0; jj < atom1->nw * GlobalV::NPOL; ++jj)
                    {
                        const int jj0 = jj / GlobalV::NPOL;
                        const int iw1_all = start1 + jj0;
                        const int mu = this->LM->ParaV->trace_loc_row[iw1_all];
                        if (mu < 0) continue;

                        for (int kk = 0; kk < atom2->nw * GlobalV::NPOL; ++kk)
                        {
                            const int kk0 = kk / GlobalV::NPOL;
                            const int iw2_all = start2 + kk0;
                            const int nu = this->LM->ParaV->trace_loc_col[iw2_all];
                            if (nu < 0) continue;

                            dSR_gamma[nu * this->LM->ParaV->nrow + mu] += dS_ptr[nnr] * this->LM->DH_r[nnr * 3 + dim2];

                            ++nnr;
                        } // kk
                    } // jj
                } // adj
            } // ad
        } // I1
    } // T1

    return;
}

void DFTU::folding_matrix_k(const int ik, const int dim1, const int dim2, std::complex<double>* mat_k, std::vector<ModuleBase::Vector3<double>> kvec_d)
{
    ModuleBase::TITLE("DFTU", "folding_matrix_k");
    ModuleBase::timer::tick("DFTU", "folding_matrix_k");
    ModuleBase::GlobalFunc::ZEROS(mat_k, this->LM->ParaV->nloc);

    double* mat_ptr;
    if      (dim1 == 1 || dim1 == 4) mat_ptr = this->LM->DSloc_Rx;
    else if (dim1 == 2 || dim1 == 5) mat_ptr = this->LM->DSloc_Ry;
    else if (dim1 == 3 || dim1 == 6) mat_ptr = this->LM->DSloc_Rz;

    int nnr = 0;
    ModuleBase::Vector3<double> dtau;
    ModuleBase::Vector3<double> tau1;
    ModuleBase::Vector3<double> tau2;

    ModuleBase::Vector3<double> dtau1;
    ModuleBase::Vector3<double> dtau2;
    ModuleBase::Vector3<double> tau0;

    for (int T1 = 0; T1 < GlobalC::ucell.ntype; ++T1)
    {
        Atom* atom1 = &GlobalC::ucell.atoms[T1];
        for (int I1 = 0; I1 < atom1->na; ++I1)
        {
            tau1 = atom1->tau[I1];
            GlobalC::GridD.Find_atom(GlobalC::ucell, tau1, T1, I1);
            Atom* atom1 = &GlobalC::ucell.atoms[T1];
            const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);

            // (2) search among all adjacent atoms.
            for (int ad = 0; ad < GlobalC::GridD.getAdjacentNum() + 1; ++ad)
            {
                const int T2 = GlobalC::GridD.getType(ad);
                const int I2 = GlobalC::GridD.getNatom(ad);
                Atom* atom2 = &GlobalC::ucell.atoms[T2];

                tau2 = GlobalC::GridD.getAdjacentTau(ad);
                dtau = tau2 - tau1;
                double distance = dtau.norm() * GlobalC::ucell.lat0;
                double rcut = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut();

                bool adj = false;

                if (distance < rcut)
                {
                    adj = true;
                }
                else if (distance >= rcut)
                {
                    for (int ad0 = 0; ad0 < GlobalC::GridD.getAdjacentNum() + 1; ++ad0)
                    {
                        const int T0 = GlobalC::GridD.getType(ad0);
                        const int I0 = GlobalC::GridD.getNatom(ad0);

                        tau0 = GlobalC::GridD.getAdjacentTau(ad0);
                        dtau1 = tau0 - tau1;
                        dtau2 = tau0 - tau2;

                        double distance1 = dtau1.norm() * GlobalC::ucell.lat0;
                        double distance2 = dtau2.norm() * GlobalC::ucell.lat0;

                        double rcut1 = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();
                        double rcut2 = GlobalC::ORB.Phi[T2].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();

                        if (distance1 < rcut1 && distance2 < rcut2)
                        {
                            adj = true;
                            break;
                        }
                    }
                }

                if (adj)
                {
                    // (3) calculate the nu of atom (T2, I2)
                    const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);
                    //------------------------------------------------
                    // exp(k dot dR)
                    // dR is the index of box in Crystal coordinates
                    //------------------------------------------------
                    ModuleBase::Vector3<double> dR(GlobalC::GridD.getBox(ad).x,
                                                   GlobalC::GridD.getBox(ad).y,
                                                   GlobalC::GridD.getBox(ad).z);
                    const double arg = (kvec_d[ik] * dR) * ModuleBase::TWO_PI;
                    const std::complex<double> kphase = std::complex<double>(cos(arg), sin(arg));

                    //--------------------------------------------------
                    // calculate how many matrix elements are in
                    // this processor.
                    //--------------------------------------------------
                    for (int ii = 0; ii < atom1->nw * GlobalV::NPOL; ii++)
                    {
                        // the index of orbitals in this processor
                        const int iw1_all = start1 + ii;
                        const int mu = this->LM->ParaV->trace_loc_row[iw1_all];
                        if (mu < 0) continue;

                        for (int jj = 0; jj < atom2->nw * GlobalV::NPOL; jj++)
                        {
                            int iw2_all = start2 + jj;
                            const int nu = this->LM->ParaV->trace_loc_col[iw2_all];
                            if (nu < 0) continue;

                            int iic;
                            if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER())
                            {
                                iic = mu + nu * this->LM->ParaV->nrow;
                            }
                            else
                            {
                                iic = mu * this->LM->ParaV->ncol + nu;
                            }

                            if (dim1 == 0)
                            {
                                if (GlobalV::NSPIN != 4)
                                {
                                    mat_k[iic] += this->LM->SlocR[nnr] * kphase;
                                }
                                else
                                {
                                    mat_k[iic] += this->LM->SlocR_soc[nnr] * kphase;
                                }
                            }
                            else if (dim1 <= 3)
                            {
                                mat_k[iic] += mat_ptr[nnr] * kphase;
                            }
                            else
                            {
                                mat_k[iic] += mat_ptr[nnr] * this->LM->DH_r[nnr * 3 + dim2] * kphase;
                            }

                            ++nnr;
                        } // kk
                    } // jj
                } // adj

            } // ad
        } // I1
    } // T1
    ModuleBase::timer::tick("DFTU", "folding_matrix_k");

    return;
}
} // namespace ModuleDFTU