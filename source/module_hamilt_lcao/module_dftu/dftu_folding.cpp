#ifdef __LCAO
#include "dftu.h"
#include "module_base/timer.h"
#include "module_parameter/parameter.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer_funcs.h"

namespace ModuleDFTU
{

void DFTU::fold_dSR_gamma(const UnitCell& ucell,
                          const Parallel_Orbitals& pv,
                          const Grid_Driver* gd,
                          double* dsloc_x,
                          double* dsloc_y,
                          double* dsloc_z,
                          double* dh_r,
                          const int dim1,
                          const int dim2,
                          double* dSR_gamma)
{
    ModuleBase::TITLE("DFTU", "fold_dSR_gamma");

    ModuleBase::GlobalFunc::ZEROS(dSR_gamma, pv.nloc);

    double* dS_ptr = nullptr;
	if(dim1 == 0) 
	{
		dS_ptr = dsloc_x;
	}
	else if(dim1 == 1) 
	{
		dS_ptr = dsloc_y;
	}
	else if (dim1 == 2) 
	{
		dS_ptr = dsloc_z;
	}

    int nnr = 0;
    ModuleBase::Vector3<double> tau1, tau2, dtau;
    ModuleBase::Vector3<double> dtau1, dtau2, tau0;

    for (int T1 = 0; T1 < ucell.ntype; ++T1)
    {
        Atom* atom1 = &ucell.atoms[T1];
        for (int I1 = 0; I1 < atom1->na; ++I1)
        {
            tau1 = atom1->tau[I1];
            const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
            gd->Find_atom(ucell, tau1, T1, I1);
            for (int ad = 0; ad < gd->getAdjacentNum() + 1; ++ad)
            {
                const int T2 = gd->getType(ad);
                const int I2 = gd->getNatom(ad);
                const int start2 = ucell.itiaiw2iwt(T2, I2, 0);
                Atom* atom2 = &ucell.atoms[T2];
                tau2 = gd->getAdjacentTau(ad);
                dtau = tau2 - tau1;
                double distance = dtau.norm() * ucell.lat0;
                double rcut = orb_cutoff_[T1] + orb_cutoff_[T2];
                bool adj = false;
				if (distance < rcut)
				{
					adj = true;
				}
                else if (distance >= rcut)
                {
                    for (int ad0 = 0; ad0 < gd->getAdjacentNum() + 1; ++ad0)
                    {
                        const int T0 = gd->getType(ad0);
                        const int I0 = gd->getNatom(ad0);
                        const int iat0 = ucell.itia2iat(T0, I0);
                        const int start0 = ucell.itiaiw2iwt(T0, I0, 0);
                        tau0 = gd->getAdjacentTau(ad0);
                        dtau1 = tau0 - tau1;
                        dtau2 = tau0 - tau2;
                        double distance1 = dtau1.norm() * ucell.lat0;
                        double distance2 = dtau2.norm() * ucell.lat0;
                        double rcut1 = orb_cutoff_[T1] + ucell.infoNL.Beta[T0].get_rcut_max();
                        double rcut2 = orb_cutoff_[T2] + ucell.infoNL.Beta[T0].get_rcut_max();
                        if (distance1 < rcut1 && distance2 < rcut2)
                        {
                            adj = true;
                            break;
                        }
                    }
                }

                if (adj)
                {
                    for (int jj = 0; jj < atom1->nw * PARAM.globalv.npol; ++jj)
                    {
                        const int jj0 = jj / PARAM.globalv.npol;
                        const int iw1_all = start1 + jj0;
                        const int mu = pv.global2local_row(iw1_all);
						if (mu < 0) 
						{
							continue;
						}

                        for (int kk = 0; kk < atom2->nw * PARAM.globalv.npol; ++kk)
                        {
                            const int kk0 = kk / PARAM.globalv.npol;
                            const int iw2_all = start2 + kk0;
                            const int nu = pv.global2local_col(iw2_all);
							if (nu < 0) 
							{
								continue;
							}

                            dSR_gamma[nu * pv.nrow + mu] += dS_ptr[nnr] * dh_r[nnr * 3 + dim2];

                            ++nnr;
                        } // kk
                    } // jj
                } // adj
            } // ad
        } // I1
    } // T1

    return;
}

void DFTU::folding_matrix_k(const UnitCell& ucell,
                            const Grid_Driver& gd,
                            ForceStressArrays& fsr,
                            const Parallel_Orbitals& pv,
                            const int ik,
                            const int dim1,
                            const int dim2,
                            std::complex<double>* mat_k,
                            const std::vector<ModuleBase::Vector3<double>>& kvec_d)
{
    ModuleBase::TITLE("DFTU", "folding_matrix_k");
    ModuleBase::timer::tick("DFTU", "folding_matrix_k");
    ModuleBase::GlobalFunc::ZEROS(mat_k, pv.nloc);

    double* mat_ptr;
    if (dim1 == 1 || dim1 == 4)
    {
        mat_ptr = fsr.DSloc_Rx;
    }
    else if (dim1 == 2 || dim1 == 5)
    {
        mat_ptr = fsr.DSloc_Ry;
    }
    else if (dim1 == 3 || dim1 == 6)
    {
        mat_ptr = fsr.DSloc_Rz;
    }

    int nnr = 0;
    ModuleBase::Vector3<double> dtau;
    ModuleBase::Vector3<double> tau1;
    ModuleBase::Vector3<double> tau2;

    ModuleBase::Vector3<double> dtau1;
    ModuleBase::Vector3<double> dtau2;
    ModuleBase::Vector3<double> tau0;

    for (int T1 = 0; T1 < ucell.ntype; ++T1)
    {
        Atom* atom1 = &ucell.atoms[T1];
        for (int I1 = 0; I1 < atom1->na; ++I1)
        {
            tau1 = atom1->tau[I1];
            gd.Find_atom(ucell, tau1, T1, I1);
            Atom* atom1 = &ucell.atoms[T1];
            const int start1 = ucell.itiaiw2iwt(T1, I1, 0);

            // (2) search among all adjacent atoms.
            for (int ad = 0; ad < gd.getAdjacentNum() + 1; ++ad)
            {
                const int T2 = gd.getType(ad);
                const int I2 = gd.getNatom(ad);
                Atom* atom2 = &ucell.atoms[T2];

                tau2 = gd.getAdjacentTau(ad);
                dtau = tau2 - tau1;
                double distance = dtau.norm() * ucell.lat0;
                double rcut = orb_cutoff_[T1] + orb_cutoff_[T2];

                bool adj = false;

                if (distance < rcut)
                {
                    adj = true;
                }
                else if (distance >= rcut)
                {
                    for (int ad0 = 0; ad0 < gd.getAdjacentNum() + 1; ++ad0)
                    {
                        const int T0 = gd.getType(ad0);
                        const int I0 = gd.getNatom(ad0);

                        tau0 = gd.getAdjacentTau(ad0);
                        dtau1 = tau0 - tau1;
                        dtau2 = tau0 - tau2;

                        double distance1 = dtau1.norm() * ucell.lat0;
                        double distance2 = dtau2.norm() * ucell.lat0;

                        double rcut1 = orb_cutoff_[T1] + ucell.infoNL.Beta[T0].get_rcut_max();
                        double rcut2 = orb_cutoff_[T2] + ucell.infoNL.Beta[T0].get_rcut_max();

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
                    const int start2 = ucell.itiaiw2iwt(T2, I2, 0);
                    //------------------------------------------------
                    // exp(k dot dR)
                    // dR is the index of box in Crystal coordinates
                    //------------------------------------------------
                    ModuleBase::Vector3<double> dR(gd.getBox(ad).x, gd.getBox(ad).y, gd.getBox(ad).z);
                    const double arg = (kvec_d[ik] * dR) * ModuleBase::TWO_PI;
                    const std::complex<double> kphase = std::complex<double>(cos(arg), sin(arg));

                    //--------------------------------------------------
                    // calculate how many matrix elements are in
                    // this processor.
                    //--------------------------------------------------
                    for (int ii = 0; ii < atom1->nw * PARAM.globalv.npol; ii++)
                    {
                        // the index of orbitals in this processor
                        const int iw1_all = start1 + ii;
                        const int mu = pv.global2local_row(iw1_all);
                        if (mu < 0) { continue;
}

                        for (int jj = 0; jj < atom2->nw * PARAM.globalv.npol; jj++)
                        {
                            int iw2_all = start2 + jj;
                            const int nu = pv.global2local_col(iw2_all);
                            if (nu < 0) { continue;
}

                            int iic;
                            if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER(PARAM.inp.ks_solver))
                            {
                                iic = mu + nu * pv.nrow;
                            }
                            else
                            {
                                iic = mu * pv.ncol + nu;
                            }

                            if (dim1 <= 3)
                            {
                                mat_k[iic] += mat_ptr[nnr] * kphase;
                            }
                            else
                            {
                                mat_k[iic] += mat_ptr[nnr] * fsr.DH_r[nnr * 3 + dim2] * kphase;
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

void DFTU::folding_matrix_k_new(const int ik,
    hamilt::Hamilt<std::complex<double>>* p_ham)
{
    ModuleBase::TITLE("DFTU", "folding_matrix_k_new");
    ModuleBase::timer::tick("DFTU", "folding_matrix_k_new");

    int hk_type = 0;
    if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER(PARAM.inp.ks_solver))
    {
        hk_type = 1;
    }

    // get SR and fold to mat_k
    if(PARAM.globalv.gamma_only_local)
    {
        dynamic_cast<hamilt::HamiltLCAO<double, double>*>(p_ham)->updateSk(ik, hk_type);
    }
    else
    {
        if(PARAM.inp.nspin != 4)
        {
            dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, double>*>(p_ham)->updateSk(ik, hk_type);
        }
        else
        {
            dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, std::complex<double>>*>(p_ham)->updateSk(ik, hk_type);
        }
    }
}

    

} // namespace ModuleDFTU
#endif // __LCAO