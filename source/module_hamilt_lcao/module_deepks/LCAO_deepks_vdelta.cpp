//This file contains subroutines related to V_delta, which is the deepks contribution to Hamiltonian
//defined as |alpha>V(D)<alpha|
//as well as subroutines for printing them for checking
//It also contains subroutine related to calculating e_delta_bands, which is basically
//tr (rho * V_delta)

//Four subroutines are contained in the file:
//1. add_v_delta : adds deepks contribution to hamiltonian, for gamma only
//2. add_v_delta_k : counterpart of 1, for multi-k
//3. check_v_delta : prints H_V_delta for checking
//4. check_v_delta_k : prints H_V_deltaR for checking
//5. cal_e_delta_band : calculates e_delta_bands for gamma only
//6. cal_e_delta_band_k : counterpart of 4, for multi-k

#ifdef __DEEPKS

#include "LCAO_deepks.h"
#include "module_base/vector3.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"

//this subroutine adds dV to the Kohn-Sham Hamiltonian
//for gamma_only calculations
void LCAO_Deepks::add_v_delta(const UnitCell &ucell,
    const LCAO_Orbitals &orb,
    Grid_Driver& GridD)
{
    ModuleBase::TITLE("LCAO_DESCRIPTOR", "add_v_delta");
    ModuleBase::timer::tick ("LCAO_gen_fixedH","add_v_delta");
    ModuleBase::GlobalFunc::ZEROS(this->H_V_delta.data(), pv->nrow * pv->ncol); //init before calculate

    const double Rcut_Alpha = orb.Alpha[0].getRcut();

    for (int T0 = 0; T0 < ucell.ntype; T0++)
    {
		Atom* atom0 = &ucell.atoms[T0]; 
        for (int I0 =0; I0< atom0->na; I0++)
        {
            const int iat = ucell.itia2iat(T0,I0);
            const ModuleBase::Vector3<double> tau0 = atom0->tau[I0];
            GridD.Find_atom(ucell, atom0->tau[I0] ,T0, I0);

            for (int ad1=0; ad1<GridD.getAdjacentNum()+1 ; ++ad1)
            {
                const int T1 = GridD.getType(ad1);
                const int I1 = GridD.getNatom(ad1);
                const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
                const ModuleBase::Vector3<double> tau1 = GridD.getAdjacentTau(ad1);
				const Atom* atom1 = &ucell.atoms[T1];
				const int nw1_tot = atom1->nw*GlobalV::NPOL;
				const double Rcut_AO1 = orb.Phi[T1].getRcut(); 

				for (int ad2=0; ad2 < GridD.getAdjacentNum()+1 ; ad2++)
				{
					const int T2 = GridD.getType(ad2);
					const int I2 = GridD.getNatom(ad2);
					const int start2 = ucell.itiaiw2iwt(T2, I2, 0);
					const ModuleBase::Vector3<double> tau2 = GridD.getAdjacentTau(ad2);
					const Atom* atom2 = &ucell.atoms[T2];
					const int nw2_tot = atom2->nw*GlobalV::NPOL;
					
					const double Rcut_AO2 = orb.Phi[T2].getRcut();
                	const double dist1 = (tau1-tau0).norm() * ucell.lat0;
                	const double dist2 = (tau2-tau0).norm() * ucell.lat0;

					if (dist1 > Rcut_Alpha + Rcut_AO1
							|| dist2 > Rcut_Alpha + Rcut_AO2)
					{
						continue;
					}

					for (int iw1=0; iw1<nw1_tot; ++iw1)
					{
						const int iw1_all = start1 + iw1;
                        const int iw1_local = pv->global2local_row(iw1_all);
						if(iw1_local < 0)continue;
						const int iw1_0 = iw1/GlobalV::NPOL;
						for (int iw2=0; iw2<nw2_tot; ++iw2)
						{
							const int iw2_all = start2 + iw2;
                            const int iw2_local = pv->global2local_col(iw2_all);
							if(iw2_local < 0)continue;
							const int iw2_0 = iw2/GlobalV::NPOL;

                            double nlm=0.0;
                            std::vector<double> nlm1 = this->nlm_save[iat][ad1][iw1_all][0];
                            std::vector<double> nlm2 = this->nlm_save[iat][ad2][iw2_all][0];

                            assert(nlm1.size()==nlm2.size());
                            int ib=0;
                            for (int L0 = 0; L0 <= orb.Alpha[0].getLmax();++L0)
                            {
                                for (int N0 = 0;N0 < orb.Alpha[0].getNchi(L0);++N0)
                                {
                                    const int inl = this->inl_index[T0](I0, L0, N0);
                                    const int nm = 2*L0+1;
                                    for (int m1 = 0;m1 < 2 * L0 + 1;++m1)
                                    {
                                        for (int m2 = 0; m2 < 2 * L0 + 1; ++m2)
                                        {
                                            nlm += this->gedm[inl][m1*nm+m2]*nlm1[ib+m1]*nlm2[ib+m2];
                                        }
                                    }
                                    ib+=nm;
                                }
                            }

                            assert(ib==nlm1.size());
                            
                            int iic;

                            if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER())
                            {
                                iic = iw1_local + iw2_local * pv->nrow;
                            }
                            else
                            {
                                iic = iw1_local * pv->ncol + iw2_local;
                            }
                            this->H_V_delta[iic] += nlm;
						}//iw2
					}//iw1
				}//ad2
			}//ad1
        }
    }

    ModuleBase::timer::tick ("LCAO_DESCRIPTOR","add_v_delta");
	return;
}

void LCAO_Deepks::check_v_delta()
{
	std::ofstream ofs("H_V_delta.dat");
	ofs << std::setprecision(10);
    for (int irow = 0;irow < pv->nrow;irow++)
	{
        for (int icol = 0;icol < pv->ncol;icol++)
		{
			int iic;
            if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER())
			{
                iic = irow + icol * pv->nrow;
			}
			else
			{
                iic = irow * pv->ncol + icol;
			}
			ofs<<H_V_delta[iic]<< " ";
		}
		ofs << std::endl;
	}
}

//this subroutine calculates H_V_deltaR
//used in multi-k calculations
void LCAO_Deepks::add_v_delta_k(const UnitCell &ucell,
    const LCAO_Orbitals &orb,
    Grid_Driver& GridD,
    const int nnr_in)
{
    ModuleBase::TITLE("LCAO_DESCRIPTOR", "add_v_delta_k");
    ModuleBase::timer::tick ("LCAO_DESCRIPTOR","add_v_delta_k");
    ModuleBase::GlobalFunc::ZEROS(this->H_V_deltaR, nnr_in);

    const double Rcut_Alpha = orb.Alpha[0].getRcut();

	int nnr = 0;
	ModuleBase::Vector3<double> tau1, tau2, dtau;
	ModuleBase::Vector3<double> dtau1, dtau2, tau0;
	ModuleBase::Vector3<float> dtau1_f, dtau2_f;
	double distance = 0.0;
	double rcut = 0.0;
	double rcut1, rcut2;
		
	//	Record_adj RA;
	//	RA.for_2d();

	// psi1
	for (int T1 = 0; T1 < ucell.ntype; ++T1)
	{
		const Atom* atom1 = &ucell.atoms[T1];
		for (int I1 =0; I1< atom1->na; ++I1)
		{

			GridD.Find_atom(ucell, atom1->tau[I1] ,T1, I1);
			const int iat1 = ucell.itia2iat(T1, I1);
			const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
			tau1 = atom1->tau[I1];

			// psi2
			for (int ad2=0; ad2<GridD.getAdjacentNum()+1; ++ad2)
			{
				const int T2 = GridD.getType(ad2);
				const Atom* atom2 = &ucell.atoms[T2];
				
				const int I2 = GridD.getNatom(ad2);
				const int iat2 = ucell.itia2iat(T2, I2);
				const int start2 = ucell.itiaiw2iwt(T2, I2, 0);
				tau2 = GridD.getAdjacentTau(ad2);

				bool is_adj = false;

				const int rx2=GridD.getBox(ad2).x;
				const int ry2=GridD.getBox(ad2).y;
				const int rz2=GridD.getBox(ad2).z;

					
				dtau = tau2 - tau1;
				distance = dtau.norm2() * pow(ucell.lat0,2);
				// this rcut is in order to make nnr consistent 
				// with other matrix.
				rcut = pow(orb.Phi[T1].getRcut() + orb.Phi[T2].getRcut(),2);
				if(distance < rcut) is_adj = true;
				else if(distance >= rcut)
				{
					for (int ad0 = 0; ad0 < GridD.getAdjacentNum()+1; ++ad0)
					{
						const int T0 = GridD.getType(ad0);
						//const int I0 = GridD.getNatom(ad0);
						//const int T0 = RA.info[iat1][ad0][3];
						//const int I0 = RA.info[iat1][ad0][4];
						//const int iat0 = ucell.itia2iat(T0, I0);
						//const int start0 = ucell.itiaiw2iwt(T0, I0, 0);

						tau0 = GridD.getAdjacentTau(ad0);
						dtau1 = tau0 - tau1;
						dtau2 = tau0 - tau2;

						const double distance1 = dtau1.norm2() * pow(ucell.lat0,2);
						const double distance2 = dtau2.norm2() * pow(ucell.lat0,2);

						rcut1 = pow(orb.Phi[T1].getRcut() + ucell.infoNL.Beta[T0].get_rcut_max(),2);
						rcut2 = pow(orb.Phi[T2].getRcut() + ucell.infoNL.Beta[T0].get_rcut_max(),2);

						if( distance1 < rcut1 && distance2 < rcut2 )
						{
							is_adj = true;
							break;
						}
					}
				}


				if(is_adj)
				{
					// < psi1 | all projectors | psi2 >
					// ----------------------------- enter the nnr increaing zone -------------------------
					for (int ad0=0; ad0 < GridD.getAdjacentNum()+1 ; ++ad0)
					{
						const int T0 = GridD.getType(ad0);
						const int I0 = GridD.getNatom(ad0);
						const int iat = ucell.itia2iat(T0,I0);

						// mohan add 2010-12-19
						if( ucell.infoNL.nproj[T0] == 0) continue;

						//const int I0 = GridD.getNatom(ad0);
						//const int start0 = ucell.itiaiw2iwt(T0, I0, 0);
						tau0 = GridD.getAdjacentTau(ad0);

						dtau1 = tau0 - tau1;
						dtau2 = tau0 - tau2;
						const double distance1 = dtau1.norm2() * pow(ucell.lat0,2);
						const double distance2 = dtau2.norm2() * pow(ucell.lat0,2);

						// seems a bug here!! mohan 2011-06-17
						rcut1 = pow(orb.Phi[T1].getRcut() + Rcut_Alpha,2);
						rcut2 = pow(orb.Phi[T2].getRcut() + Rcut_Alpha,2);

						if(distance1 >= rcut1 || distance2 >= rcut2)
						{
							continue;
						}
						//const Atom* atom0 = &ucell.atoms[T0];
						const int rx0=GridD.getBox(ad0).x;
						const int ry0=GridD.getBox(ad0).y;
						const int rz0=GridD.getBox(ad0).z;
						key_tuple key_1(iat1,-rx0,-ry0,-rz0);
						key_tuple key_2(iat2,rx2-rx0,ry2-ry0,rz2-rz0);

                        int nnr_inner = 0;

						for (int iw1=0; iw1<atom1->nw*GlobalV::NPOL; iw1++)
						{
							const int iw1_all = start1 + iw1;
                            const int mu = pv->global2local_row(iw1_all);
							if(mu < 0)continue; 

							// fix a serious bug: atom2[T2] -> atom2
							// mohan 2010-12-20
							for (int iw2=0; iw2<atom2->nw*GlobalV::NPOL; iw2++)
							{
								const int iw2_all = start2 + iw2;
                                const int nu = pv->global2local_col(iw2_all);
								if(nu < 0)continue;
  
                                std::vector<double> nlm1 = this->nlm_save_k[iat][key_1][iw1_all][0];
                                std::vector<double> nlm2 = this->nlm_save_k[iat][key_2][iw2_all][0];
                                assert(nlm1.size()==nlm2.size());

                                int ib=0;
                                double nlm = 0.0;
                                for (int L0 = 0; L0 <= orb.Alpha[0].getLmax();++L0)
                                {
                                    for (int N0 = 0;N0 < orb.Alpha[0].getNchi(L0);++N0)
                                    {
                                        const int inl = this->inl_index[T0](I0, L0, N0);
                                        const int nm = 2*L0+1;
                                        for (int m1 = 0;m1 < 2 * L0 + 1;++m1)
                                        {
                                            for (int m2 = 0; m2 < 2 * L0 + 1; ++m2)
                                            {
                                                nlm += this->gedm[inl][m1*nm+m2]*nlm1[ib+m1]*nlm2[ib+m2];
                                            }
                                        }
                                        ib+=(2*L0+1);
                                    }
                                }
                                assert(ib==nlm1.size());

                                this->H_V_deltaR[nnr+nnr_inner] += nlm;
                                nnr_inner++;
                            }//iw2
                        }//iw1
                    }//ad0

					//outer circle : accumulate nnr
					for (int j=0; j<atom1->nw*GlobalV::NPOL; j++)
					{
						const int j0 = j/GlobalV::NPOL;//added by zhengdy-soc
						const int iw1_all = start1 + j;
                        const int mu = pv->global2local_row(iw1_all);
						if(mu < 0)continue; 

						// fix a serious bug: atom2[T2] -> atom2
						// mohan 2010-12-20
						for (int k=0; k<atom2->nw*GlobalV::NPOL; k++)
						{
							const int k0 = k/GlobalV::NPOL;
							const int iw2_all = start2 + k;
                            const int nu = pv->global2local_col(iw2_all);
							if(nu < 0)continue;

							nnr++;
						}
					}
                }//is_adj
            }//ad2
        }//I1
    }//T1

    if( nnr!=nnr_in)
    {
        ModuleBase::WARNING_QUIT("LCAO_DESCRIPTOR","nnr!=LNNR.nnr");
    }

    ModuleBase::timer::tick ("LCAO_DESCRIPTOR","add_v_delta_k");
    return;
}

void LCAO_Deepks::check_v_delta_k(const int nnr)
{
	std::ofstream ofs("H_V_deltaR.dat");
	ofs<<std::setprecision(10);
	for(int iir=0;iir<nnr;iir++)
	{
		ofs<<H_V_deltaR[iir]<<" ";
		if(iir%6==5) ofs<<std::endl;
	}
}

//calculating sum of correction band energies
//for gamma_only calculations
void LCAO_Deepks::cal_e_delta_band(const std::vector<std::vector<double>>& dm)
{
    ModuleBase::TITLE("LCAO_Deepks", "cal_e_delta_band");
    this->e_delta_band = 0;
    for (int i = 0; i < GlobalV::NLOCAL; ++i)
    {
        for (int j = 0; j < GlobalV::NLOCAL; ++j)
        {
            const int mu = pv->global2local_row(j);
            const int nu = pv->global2local_col(i);
            
            if (mu >= 0 && nu >= 0)
            {                
                const int index = nu * pv->nrow + mu;
                for (int is = 0; is < dm.size(); ++is)  //dm.size() == GlobalV::NSPIN
                {
                    //this->e_delta_band += dm[is](nu, mu) * this->H_V_delta[index];
					this->e_delta_band += dm[is][nu*this->pv->nrow+mu] * this->H_V_delta[index];
                }
            }
        }
    }
#ifdef __MPI
    Parallel_Reduce::reduce_all(this->e_delta_band);
#endif
    return;
}

//calculating sum of correction band energies
//for multi_k calculations
void LCAO_Deepks::cal_e_delta_band_k(const std::vector<std::vector<std::complex<double>>>& dm,
    const int nks)
{
    ModuleBase::TITLE("LCAO_Deepks", "cal_e_delta_band");
    std::complex<double> e_delta_band_k=std::complex<double>(0.0,0.0);
    for (int i = 0; i < GlobalV::NLOCAL; ++i)
    {
        for (int j = 0; j < GlobalV::NLOCAL; ++j)
        {
            const int mu = pv->global2local_row(j);
            const int nu = pv->global2local_col(i);
            
            if (mu >= 0 && nu >= 0)
            {                
                int iic;
                if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER())
                {
                    iic = mu + nu * pv->nrow;
                }
                else
                {
                    iic = mu * pv->ncol + nu;
                }
                for(int ik=0;ik<nks;ik++)
                {
                    //e_delta_band_k += dm[ik](nu, mu) * this->H_V_delta_k[ik][iic];
					e_delta_band_k += dm[ik][nu * this->pv->nrow + mu] * this->H_V_delta_k[ik][iic];
                }
            }
        }
    }

    this->e_delta_band = e_delta_band_k.real();
#ifdef __MPI
    Parallel_Reduce::reduce_all(this->e_delta_band);
#endif
    return;
}

#endif
