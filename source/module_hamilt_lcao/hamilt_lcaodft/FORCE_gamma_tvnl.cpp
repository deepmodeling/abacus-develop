#include "FORCE_gamma.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include <unordered_map>
#include "module_base/timer.h"

void Force_LCAO_gamma::cal_fvnl_dbeta_new(
	const std::vector<ModuleBase::matrix> &dm2d, 
	const bool isforce, 
	const bool isstress, 
	ModuleBase::matrix& fvnl_dbeta, 
	ModuleBase::matrix& svnl_dbeta)
{
    ModuleBase::TITLE("Force_LCAO_gamma","cal_fvnl_dbeta_new");
    ModuleBase::timer::tick("Force_LCAO_gamma","cal_fvnl_dbeta_new");

    double r0[3];
	double r1[3];

    for(int iat=0; iat<GlobalC::ucell.nat; iat++)
    {
        const int it = GlobalC::ucell.iat2it[iat];
        const int ia = GlobalC::ucell.iat2ia[iat];
        const int T0 = it;
        const int I0 = ia;
        
        const ModuleBase::Vector3<double> tau0 = GlobalC::ucell.atoms[it].tau[ia];
        //find ajacent atom of atom ia
        //GlobalC::GridD.Find_atom( GlobalC::ucell.atoms[it].tau[ia] );
		GlobalC::GridD.Find_atom(GlobalC::ucell, tau0 ,it, ia);
		const double Rcut_Beta = GlobalC::ucell.infoNL.Beta[it].get_rcut_max();

        std::vector<std::unordered_map<int,std::vector<std::vector<double>>>> nlm_tot;
		nlm_tot.resize(GlobalC::GridD.getAdjacentNum()+1); //this saves <psi_i|beta> and <psi_i|\nabla|beta>

//Step 1 : Calculate and save <psi_i|beta> and <psi_i|\nabla|beta>
        for (int ad1 =0 ; ad1 < GlobalC::GridD.getAdjacentNum()+1; ad1++)
        {
            const int T1 = GlobalC::GridD.getType (ad1);
            const Atom* atom1 = &GlobalC::ucell.atoms[T1];
            const int I1 = GlobalC::GridD.getNatom (ad1);
            const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
			const ModuleBase::Vector3<double> tau1 = GlobalC::GridD.getAdjacentTau (ad1);
			const double Rcut_AO1 = GlobalC::ORB.Phi[T1].getRcut();

            nlm_tot[ad1].clear();

            const double dist1 = (tau1-tau0).norm() * GlobalC::ucell.lat0;
            if (dist1 > Rcut_Beta + Rcut_AO1)
            {
                continue;
            }

			for (int iw1=0; iw1<GlobalC::ucell.atoms[T1].nw; ++iw1)
			{
				const int iw1_all = start1 + iw1;
                const int iw1_local = this->ParaV->global2local_row[iw1_all];
                const int iw2_local = this->ParaV->global2local_col[iw1_all];
				if(iw1_local < 0 && iw2_local < 0) continue;
                
                std::vector<std::vector<double>> nlm;

                GlobalC::UOT.snap_psibeta_half(
                    GlobalC::ORB,
                    GlobalC::ucell.infoNL,
                    nlm, tau1, T1,
                    atom1->iw2l[ iw1 ], // L1
                    atom1->iw2m[ iw1 ], // m1
                    atom1->iw2n[ iw1 ], // N1
                    GlobalC::ucell.atoms[T0].tau[I0], T0, 1); //R0,T0

                assert(nlm.size()==4);
                nlm_tot[ad1].insert({iw1,nlm});
            }
        }//ad

//Step 2 : sum to get beta<psi_i|beta><beta|\nabla|psi_j>
		for (int ad1=0; ad1<GlobalC::GridD.getAdjacentNum()+1 ; ++ad1)
        {
            const int T1 = GlobalC::GridD.getType (ad1);
            const Atom* atom1 = &GlobalC::ucell.atoms[T1];
            const int I1 = GlobalC::GridD.getNatom (ad1);
            const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
			const ModuleBase::Vector3<double> tau1 = GlobalC::GridD.getAdjacentTau (ad1);
			const double Rcut_AO1 = GlobalC::ORB.Phi[T1].getRcut();

            for (int ad2=0; ad2 < GlobalC::GridD.getAdjacentNum()+1 ; ad2++)
            {
                const int T2 = GlobalC::GridD.getType (ad2);
                const Atom* atom2 = &GlobalC::ucell.atoms[T2];
                const int I2 = GlobalC::GridD.getNatom (ad2);
                const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);
                const ModuleBase::Vector3<double> tau2 = GlobalC::GridD.getAdjacentTau (ad2);
                const double Rcut_AO2 = GlobalC::ORB.Phi[T2].getRcut();

                const double dist1 = (tau1-tau0).norm() * GlobalC::ucell.lat0;
                const double dist2 = (tau2-tau0).norm() * GlobalC::ucell.lat0;
				if(isstress)
                {
                    r1[0] = ( tau1.x - tau0.x) ;
                    r1[1] = ( tau1.y - tau0.y) ;
                    r1[2] = ( tau1.z - tau0.z) ;
                    r0[0] = ( tau2.x - tau0.x) ;
                    r0[1] = ( tau2.y - tau0.y) ;
                    r0[2] = ( tau2.z - tau0.z) ;
                }

                if (dist1 > Rcut_Beta + Rcut_AO1
                        || dist2 > Rcut_Beta + Rcut_AO2)
                {
                    continue;
                }
                
                for (int iw1=0; iw1<GlobalC::ucell.atoms[T1].nw; ++iw1)
                {
                    const int iw1_all = start1 + iw1;
                    const int iw1_local = this->ParaV->global2local_row[iw1_all];
                    if(iw1_local < 0)continue;
                    for (int iw2=0; iw2<GlobalC::ucell.atoms[T2].nw; ++iw2)
                    {
                        const int iw2_all = start2 + iw2;
                        const int iw2_local = this->ParaV->global2local_col[iw2_all];
                        if(iw2_local < 0)continue;

                        double nlm[3] = {0,0,0};
                        double nlm_t[3] = {0,0,0}; //transpose

                        std::vector<double> nlm1 = nlm_tot[ad1][iw1][0];
                        std::vector<std::vector<double>> nlm2;
                        nlm2.resize(3);
                        for(int i=0;i<3;i++)
                        {
                            nlm2[i] = nlm_tot[ad2][iw2][i+1];
                        }

                        assert(nlm1.size()==nlm2[0].size());

						const int nproj = GlobalC::ucell.infoNL.nproj[T0];
						int ib = 0;
						for (int nb = 0; nb < nproj; nb++)
						{
							const int L0 = GlobalC::ucell.infoNL.Beta[T0].Proj[nb].getL();
							for(int m=0;m<2*L0+1;m++)
							{
                                for(int ir=0;ir<3;ir++)
                                {
                                    nlm[ir] += nlm2[ir][ib]*nlm1[ib]*GlobalC::ucell.atoms[T0].ncpp.dion(nb,nb);
                                }
								ib+=1;
							}
						}
						assert(ib==nlm1.size());

                        if(isstress)
                        {
                            nlm1 = nlm_tot[ad2][iw2][0];
                            for(int i=0;i<3;i++)
                            {
                                nlm2[i] = nlm_tot[ad1][iw1][i+1];
                            }

                            assert(nlm1.size()==nlm2[0].size());
                            
    						const int nproj = GlobalC::ucell.infoNL.nproj[T0];
	    					int ib = 0;
		    				for (int nb = 0; nb < nproj; nb++)
				    		{
		      					const int L0 = GlobalC::ucell.infoNL.Beta[T0].Proj[nb].getL();
		    					for(int m=0;m<2*L0+1;m++)
				    			{
                                   for(int ir=0;ir<3;ir++)
                                    {
                                        nlm_t[ir] += nlm2[ir][ib]*nlm1[ib]*GlobalC::ucell.atoms[T0].ncpp.dion(nb,nb);
                                    }
					    			ib+=1;
					    		}
					    	}
						    assert(ib==nlm1.size());
                        }
                        // dbeta is minus, that's consistent.
                        // only one projector for each atom force.

                        double sum = 0.0;
                        for(int is=0; is<GlobalV::NSPIN; ++is)
                        {
                            //sum += dm2d[is][index];
                            sum += dm2d[is](iw2_local, iw1_local);
                        }
                        sum *= 2.0;

                        if(isforce)
						{
							fvnl_dbeta(iat,0) -= sum * nlm[0];
							fvnl_dbeta(iat,1) -= sum * nlm[1];
							fvnl_dbeta(iat,2) -= sum * nlm[2];
						}

                        if(isstress) 
                        {
                            for(int ipol=0;ipol<3;ipol++)
							{
                                for(int jpol=ipol;jpol<3;jpol++)
                                {
                                    svnl_dbeta(ipol, jpol) += sum/2.0 * (nlm[ipol] * r0[jpol] + nlm_t[ipol] * r1[jpol]);
                                }
                            }
                        }
                    }//iw2
                }//iw1
            }//ad2
        }//ad1
    }//iat
    if(isstress)
    {
        StressTools::stress_fill(GlobalC::ucell.lat0, GlobalC::ucell.omega, svnl_dbeta);
    }
    ModuleBase::timer::tick("Force_LCAO_gamma","cal_fvnl_dbeta_new");
}

void Force_LCAO_gamma::cal_ftvnl_dphi(
	const std::vector<ModuleBase::matrix> &dm2d, 
	const bool isforce, 
	const bool isstress, 
	ModuleBase::matrix& ftvnl_dphi, 
	ModuleBase::matrix& stvnl_dphi)
{
    ModuleBase::TITLE("Force_LCAO_gamma","cal_ftvnl_dphi");
    ModuleBase::timer::tick("Force_LCAO_gamma","cal_ftvnl_dphi");

    for(int i=0; i<GlobalV::NLOCAL; i++)
    {
        const int iat = GlobalC::ucell.iwt2iat[i];
        for(int j=0; j<GlobalV::NLOCAL; j++)
        {
            const int mu = this->ParaV->global2local_row[j];
            const int nu = this->ParaV->global2local_col[i];

            if (mu >= 0 && nu >= 0 )
            {
                const int index = mu * this->ParaV->ncol + nu;
                //contribution from deriv of AO's in T+VNL term

                double sum = 0.0;
                for(int is=0; is<GlobalV::NSPIN; ++is)
                {
                    sum += dm2d[is](nu, mu);
                }
                sum *= 2.0;

                if(isforce)
				{
					ftvnl_dphi(iat,0) += sum * this->UHM->LM->DHloc_fixed_x[index];
					ftvnl_dphi(iat,1) += sum * this->UHM->LM->DHloc_fixed_y[index];
					ftvnl_dphi(iat,2) += sum * this->UHM->LM->DHloc_fixed_z[index];
				}
                if(isstress)
                {
                    stvnl_dphi(0,0) += sum/2.0 * this->UHM->LM->DHloc_fixed_11[index];
                    stvnl_dphi(0,1) += sum/2.0 * this->UHM->LM->DHloc_fixed_12[index];
                    stvnl_dphi(0,2) += sum/2.0 * this->UHM->LM->DHloc_fixed_13[index];
                    stvnl_dphi(1,1) += sum/2.0 * this->UHM->LM->DHloc_fixed_22[index];
                    stvnl_dphi(1,2) += sum/2.0 * this->UHM->LM->DHloc_fixed_23[index];
                    stvnl_dphi(2,2) += sum/2.0 * this->UHM->LM->DHloc_fixed_33[index];   
                }
            }
        }
    }
    if(isstress){
        StressTools::stress_fill(GlobalC::ucell.lat0, GlobalC::ucell.omega, stvnl_dphi);
    }
    ModuleBase::timer::tick("Force_LCAO_gamma","cal_ftvnl_dphi");
    return;
}