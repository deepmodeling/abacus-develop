#include "LCAO_gen_fixedH.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_hamilt_pw/hamilt_pwdft/wavefunc.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include <vector>
#include <unordered_map>
#include <map>
#include "module_base/timer.h"

#ifdef __MKL
#include <mkl_service.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

LCAO_gen_fixedH::LCAO_gen_fixedH()
{}

LCAO_gen_fixedH::~LCAO_gen_fixedH()
{}

void LCAO_gen_fixedH::calculate_NL_no(double* HlocR)
{
    ModuleBase::TITLE("LCAO_gen_fixedH","calculate_NL_no");
	if(GlobalV::GAMMA_ONLY_LOCAL)
	{
	  	//for gamma only.
		this->build_Nonlocal_beta_new(HlocR);
	}
	else
	{
		this->build_Nonlocal_mu_new(HlocR, false);
	}

    return;
}

void LCAO_gen_fixedH::calculate_T_no(double* HlocR)
{
    ModuleBase::TITLE("LCAO_gen_fixedH","calculate_T_no");
    this->build_ST_new('T', false, GlobalC::ucell, HlocR);
    return;
}

void LCAO_gen_fixedH::calculate_S_no(double* SlocR)
{
    ModuleBase::TITLE("LCAO_gen_fixedH", "calculate_S_no");
    ModuleBase::timer::tick("LCAO_gen_fixedH","calculate_S_no");
	this->build_ST_new('S', false, GlobalC::ucell, SlocR);
    ModuleBase::timer::tick("LCAO_gen_fixedH","calculate_S_no");
    return;
}


//liaochen modify interface 2010-3-22
void LCAO_gen_fixedH::build_ST_new(const char& dtype, const bool& calc_deri, const UnitCell &ucell, double* HSloc, bool cal_syns, double dmax)
{
    ModuleBase::TITLE("LCAO_gen_fixedH","build_ST_new");

	int total_nnr = 0;
	const Parallel_Orbitals* pv = this->LM->ParaV;
#ifdef _OPENMP
#pragma omp parallel reduction(+:total_nnr)
{
#endif
    //array to store data
    double olm[3]={0.0,0.0,0.0};

    //\sum{T} e**{ikT} <\phi_{ia}|d\phi_{k\beta}(T)>
	ModuleBase::Vector3<double> tau1, tau2, dtau;
	ModuleBase::Vector3<double> dtau1, dtau2, tau0;
#ifdef _OPENMP
// use schedule(dynamic) for load balancing because adj_num is various
#pragma omp for schedule(dynamic)
#endif
    for (int iat1 = 0; iat1 < GlobalC::ucell.nat; iat1++)
    {
		const int T1 = GlobalC::ucell.iat2it[iat1];
		const Atom* atom1 = &GlobalC::ucell.atoms[T1];
		const int I1 = GlobalC::ucell.iat2ia[iat1];
        {
			tau1 = atom1->tau[I1];

            //GlobalC::GridD.Find_atom(tau1);
			AdjacentAtomInfo adjs;
            GlobalC::GridD.Find_atom(ucell, tau1, T1, I1, &adjs);
			// Record_adj.for_2d() may not called in some case
			int nnr = pv->nlocstart ? pv->nlocstart[iat1] : 0;

			if (cal_syns)
            {
                for (int k = 0; k < 3; k++)
                {
                    tau1[k] = tau1[k] - atom1->vel[I1][k] * INPUT.mdp.md_dt / GlobalC::ucell.lat0 ;
                }
            }

            for (int ad = 0; ad < adjs.adj_num+1; ++ad)
            {
                const int T2 = adjs.ntype[ad];
				const int I2 = adjs.natom[ad];
				Atom* atom2 = &ucell.atoms[T2];
				tau2 = adjs.adjacent_tau[ad];
				dtau = tau2 - tau1;
				double distance = dtau.norm() * ucell.lat0;
				double rcut = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut();

				if(distance < rcut)
				{
					int iw1_all = ucell.itiaiw2iwt( T1, I1, 0) ; //iw1_all = combined index (it, ia, iw)

					for(int jj=0; jj<atom1->nw*GlobalV::NPOL; ++jj)
					{
						const int jj0 = jj/GlobalV::NPOL;
						const int L1 = atom1->iw2l[jj0];
						const int N1 = atom1->iw2n[jj0];
						const int m1 = atom1->iw2m[jj0];

						int iw2_all = ucell.itiaiw2iwt( T2, I2, 0);//zhengdy-soc
						for(int kk=0; kk<atom2->nw*GlobalV::NPOL; ++kk)
						{
							const int kk0 = kk/GlobalV::NPOL;
							const int L2 = atom2->iw2l[kk0];
							const int N2 = atom2->iw2n[kk0];
							const int m2 = atom2->iw2m[kk0];

							// mohan add 2010-06-29
							// this is in fact the same as in build_Nonlocal_mu,
							// the difference is that here we use {L,N,m} for ccycle,
							// build_Nonlocal_mu use atom.nw for cycle.
							// so, here we use ParaO::in_this_processor,
							// in build_Non... use global2local_row
                            // and global2local_col directly,
                            if (!pv->in_this_processor(iw1_all, iw2_all))
							{
								++iw2_all;
								continue;
							}

							olm[0] = olm[1] = olm[2] = 0.0;

							if(!calc_deri)
							{
								// PLEASE use UOT as an input parameter of this subroutine
								// mohan add 2021-03-30
			
								GlobalC::UOT.snap_psipsi( GlobalC::ORB, olm, 0, dtype, 
										tau1, T1, L1, m1, N1,                  // info of atom1
										adjs.adjacent_tau[ad], T2, L2, m2, N2, // info of atom2 
										cal_syns,
										dmax);
								// When NSPIN == 4 , only diagonal term is calculated for T or S Operators
								// use olm1 to store the diagonal term with complex data type.
								std::complex<double> olm1[4];
								if(GlobalV::NSPIN == 4)
								{
									olm1[0] = std::complex<double>(olm[0], 0.0);
									olm1[1] = ModuleBase::ZERO;
									olm1[2] = ModuleBase::ZERO;
									olm1[3] = std::complex<double>(olm[0], 0.0);
								}

								if(GlobalV::GAMMA_ONLY_LOCAL)
								{
									// mohan add 2010-06-29
									// set the value in Hloc and Sloc
									// according to global2local_row and global2local_col
									// the last paramete: 1 for Sloc, 2 for Hloc
									// and 3 for Hloc_fixed.
                                    this->LM->set_HSgamma(iw1_all, iw2_all, olm[0], HSloc);
								}
								else // k point algorithm
								{
									// mohan add 2010-10
									// set the values in SlocR and Hloc_fixedR.
									// which is a 1D array.
									if(dtype=='S')
									{
                                        if (GlobalV::NSPIN != 4) HSloc[nnr] = olm[0];
                                        else
										{//only has diagonal term here.
											int is = (jj-jj0*GlobalV::NPOL) + (kk-kk0*GlobalV::NPOL)*2;
											// SlocR_soc is a temporary array with complex data type, it will be refactor soon.
											this->LM->SlocR_soc[nnr] = olm1[is];
                                        }
                                    }
									else if(dtype=='T')
									{
										if(GlobalV::NSPIN!=4) HSloc[nnr] = olm[0];// <phi|kin|d phi>
										else
										{//only has diagonal term here.
											int is = (jj-jj0*GlobalV::NPOL) + (kk-kk0*GlobalV::NPOL)*2;
											// Hloc_fixedR_soc is a temporary array with complex data type, it will be refactor soon.
											this->LM->Hloc_fixedR_soc[nnr] = olm1[is];
                                        }
                                    }
									++total_nnr;
									++nnr;
								}
							}
							else // calculate the derivative
							{
								GlobalC::UOT.snap_psipsi( GlobalC::ORB, olm, 1, dtype, 
									tau1, T1, L1, m1, N1,
									adjs.adjacent_tau[ad], T2, L2, m2, N2
									);

								if(GlobalV::GAMMA_ONLY_LOCAL)
								{
									this->LM->set_force (iw1_all, iw2_all,	olm[0], olm[1], olm[2], dtype);
									if(GlobalV::CAL_STRESS) this->LM->set_stress (iw1_all, iw2_all, olm[0], olm[1], olm[2], dtype, dtau);
								}
								else // k point algorithm
								{
									if(dtype=='S')
									{
										this->LM->DSloc_Rx[nnr] = olm[0];
										this->LM->DSloc_Ry[nnr] = olm[1];
										this->LM->DSloc_Rz[nnr] = olm[2];
										if(GlobalV::CAL_STRESS)
										{
											this->LM->DH_r[nnr*3] = dtau.x;
											this->LM->DH_r[nnr*3 + 1] = dtau.y;
											this->LM->DH_r[nnr*3 + 2] = dtau.z;
										}
									}
									else if(dtype=='T')
									{
										// notice the 'sign'
										this->LM->DHloc_fixedR_x[nnr] = olm[0];
										this->LM->DHloc_fixedR_y[nnr] = olm[1];
										this->LM->DHloc_fixedR_z[nnr] = olm[2];
										if(GlobalV::CAL_STRESS)
										{
											this->LM->stvnl11[nnr] = olm[0] * dtau.x;
											this->LM->stvnl12[nnr] = olm[0] * dtau.y;
											this->LM->stvnl13[nnr] = olm[0] * dtau.z;
											this->LM->stvnl22[nnr] = olm[1] * dtau.y;
											this->LM->stvnl23[nnr] = olm[1] * dtau.z;
											this->LM->stvnl33[nnr] = olm[2] * dtau.z;
										}
									}
									++total_nnr;
									++nnr;
								}
							}
							++iw2_all;
						}// nw2 
						++iw1_all;
					}// nw1
				}// distance
				else if(distance>=rcut && (!GlobalV::GAMMA_ONLY_LOCAL))
				{
					int start1 = ucell.itiaiw2iwt( T1, I1, 0);
					int start2 = ucell.itiaiw2iwt( T2, I2, 0);

					bool is_adj = false;
					for (int ad0=0; ad0 < adjs.adj_num+1; ++ad0)
					{
						const int T0 = adjs.ntype[ad0];
						//const int I0 = GlobalC::GridD.getNatom(ad0);
						//const int iat0 = ucell.itia2iat(T0, I0);
						//const int start0 = ucell.itiaiw2iwt(T0, I0, 0);
						tau0 = adjs.adjacent_tau[ad0];
						dtau1 = tau0 - tau1;
						double distance1 = dtau1.norm() * ucell.lat0;
						double rcut1 = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();
						dtau2 = tau0 - tau2;
						double distance2 = dtau2.norm() * ucell.lat0;
						double rcut2 = GlobalC::ORB.Phi[T2].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();
						if( distance1 < rcut1 && distance2 < rcut2 )
						{
							is_adj = true;
							break;
						}
					}//ad0


					if( is_adj )
					{
						for(int jj=0; jj<atom1->nw * GlobalV::NPOL; ++jj)
						{
                            const int mu = pv->global2local_row(start1 + jj);
							if(mu<0)continue; 
							for(int kk=0; kk<atom2->nw * GlobalV::NPOL; ++kk)
							{
                                const int nu = pv->global2local_col(start2 + kk);
								if(nu<0)continue;
								++total_nnr;
								++nnr;
							}//kk
						}//jj
					}
				}//distance
			}// ad
		}// I1
	}// T1
#ifdef _OPENMP
}
#endif

	if(!GlobalV::GAMMA_ONLY_LOCAL)
	{
		if(total_nnr != pv->nnr)
		{
			std::cout << " nnr=" << total_nnr << " LNNR.nnr=" << pv->nnr << std::endl;
			GlobalV::ofs_running << " nnr=" << total_nnr << " LNNR.nnr=" << pv->nnr << std::endl;
			ModuleBase::WARNING_QUIT("LCAO_gen_fixedH::build_ST_new","nnr != LNNR.nnr");
		}
	}

    return;
}

typedef std::tuple<int,int,int,int> key_tuple;

#include "record_adj.h" //mohan add 2012-07-06
void LCAO_gen_fixedH::build_Nonlocal_mu_new(double* NLloc, const bool &calc_deri)
{
    ModuleBase::TITLE("LCAO_gen_fixedH","b_NL_mu_new");
    ModuleBase::timer::tick("LCAO_gen_fixedH", "b_NL_mu_new");
    const Parallel_Orbitals* pv = this->LM->ParaV;

	// < phi1 | beta > < beta | phi2 >
	// phi1 is within the unitcell.
	// while beta is in the supercell.
	// while phi2 is in the supercell.


	//Step 1 : generate <psi|beta>

	//This is the data structure for storing <psi|beta>
	//It is a 4 layer data structure
	//The outmost layer is std::vector with size being number of atoms in unit cell
	//The second layer is a map, the key being a combination of 4 number (iat, dRx, dRy, dRz)
	//which identifies a unique adjacent atom of the first atom
	//The third layer is an unordered map, with key being the index of atomic basis |psi>
	//The inner layer is a vector, each element representing a projector |beta>
	//It then either stores the number <psi|beta> (nlm_tot)
	//or a vector of 4, storing additionally <d/dx_i psi|beta> (nlm_tot1) x_i=x,y,z
	std::vector<std::map<key_tuple,std::unordered_map<int,std::vector<double>>>> nlm_tot;
	std::vector<std::map<key_tuple,std::unordered_map<int,std::vector<std::vector<double>>>>> nlm_tot1;

	if(!calc_deri)
	{
		nlm_tot.resize(GlobalC::ucell.nat);
	}
	else
	{
		nlm_tot1.resize(GlobalC::ucell.nat);
	}
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
	for(int iat=0;iat<GlobalC::ucell.nat;iat++)
	{
		const int it = GlobalC::ucell.iat2it[iat];
		const int ia = GlobalC::ucell.iat2ia[iat];

		const double Rcut_Beta = GlobalC::ucell.infoNL.Beta[it].get_rcut_max();
		const ModuleBase::Vector3<double> tau = GlobalC::ucell.atoms[it].tau[ia];
		AdjacentAtomInfo adjs;
        GlobalC::GridD.Find_atom(GlobalC::ucell, tau ,it, ia, &adjs);

		if(!calc_deri)
		{
			nlm_tot[iat].clear();
		}
		else
		{
			nlm_tot1[iat].clear();
		}

		for (int ad=0; ad<adjs.adj_num+1 ; ++ad)
		{
			const int T1 = adjs.ntype[ad];
			const int I1 = adjs.natom[ad];
			const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
			const double Rcut_AO1 = GlobalC::ORB.Phi[T1].getRcut();

			const ModuleBase::Vector3<double> &tau1 = adjs.adjacent_tau[ad];
			const Atom* atom1 = &GlobalC::ucell.atoms[T1];
			const int nw1_tot = atom1->nw*GlobalV::NPOL;

			const ModuleBase::Vector3<double> dtau = tau1-tau;
			const double dist1 = dtau.norm2() * pow(GlobalC::ucell.lat0,2);
			
			if (dist1 > pow(Rcut_Beta + Rcut_AO1,2))
			{
				continue;
			}
			std::unordered_map<int,std::vector<double>> nlm_cur;
			std::unordered_map<int,std::vector<std::vector<double>>> nlm_cur1;
			
			if(!calc_deri)
			{
				nlm_cur.clear();
			}
			else
			{
				nlm_cur1.clear();
			}
			for (int iw1=0; iw1<nw1_tot; ++iw1)
			{
				const int iw1_all = start1 + iw1;
                const int iw1_local = pv->global2local_row(iw1_all);
                const int iw2_local = pv->global2local_col(iw1_all);
				if(iw1_local < 0 && iw2_local < 0)continue;
				const int iw1_0 = iw1/GlobalV::NPOL;
				std::vector<std::vector<double>> nlm;
				//nlm is a vector of vectors, but size of outer vector is only 1 here
				//If we are calculating force, we need also to store the gradient
				//and size of outer vector is then 4
				//inner loop : all projectors (L0,M0)
				GlobalC::UOT.snap_psibeta_half(
					GlobalC::ORB,
					GlobalC::ucell.infoNL,
					nlm, tau1, T1,
					atom1->iw2l[ iw1_0 ], // L1
					atom1->iw2m[ iw1_0 ], // m1
					atom1->iw2n[ iw1_0 ], // N1
					tau, it, calc_deri); //R0,T0
				if(!calc_deri)
				{
					nlm_cur.insert({iw1_all,nlm[0]});
				}
				else
				{
					nlm_cur1.insert({iw1_all,nlm});
				}
			}//end iw

			const int iat1=GlobalC::ucell.itia2iat(T1, I1);
			const int rx1=adjs.box[ad].x;
			const int ry1=adjs.box[ad].y;
			const int rz1=adjs.box[ad].z;
			key_tuple key_1(iat1,rx1,ry1,rz1);

			if(!calc_deri)
			{
				nlm_tot[iat][key_1]=nlm_cur;
			}
			else
			{
				nlm_tot1[iat][key_1]=nlm_cur1;
			}
		}//end ad
	}
	//=======================================================
	//Step2:	
	//calculate sum_(L0,M0) beta<psi_i|beta><beta|psi_j>
	//and accumulate the value to Hloc_fixedR(i,j)
	//=======================================================
	int total_nnr = 0;
#ifdef _OPENMP
#pragma omp parallel reduction(+:total_nnr)
{
#endif
	ModuleBase::Vector3<double> tau1, tau2, dtau;
	ModuleBase::Vector3<double> dtau1, dtau2, tau0;
	double distance = 0.0;
	double rcut = 0.0;
	double rcut1, rcut2;
		
	//	Record_adj RA;
	//	RA.for_2d();

	// psi1
#ifdef _OPENMP
// use schedule(dynamic) for load balancing because adj_num is various
#pragma omp for schedule(dynamic)
#endif
    for (int iat1 = 0; iat1 < GlobalC::ucell.nat; iat1++)
    {
		const int T1 = GlobalC::ucell.iat2it[iat1];
		const Atom* atom1 = &GlobalC::ucell.atoms[T1];
		const int I1 = GlobalC::ucell.iat2ia[iat1];
		{
			//GlobalC::GridD.Find_atom( atom1->tau[I1] );
			AdjacentAtomInfo adjs;
			GlobalC::GridD.Find_atom(GlobalC::ucell, atom1->tau[I1] ,T1, I1, &adjs);
			const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
			// Record_adj.for_2d() may not called in some case
			int nnr = pv->nlocstart ? pv->nlocstart[iat1] : 0;
			tau1 = atom1->tau[I1];

			// psi2
			for (int ad2=0; ad2<adjs.adj_num+1; ++ad2)
			{
				const int T2 = adjs.ntype[ad2];
				const Atom* atom2 = &GlobalC::ucell.atoms[T2];
				
				const int I2 = adjs.natom[ad2];
				const int iat2 = GlobalC::ucell.itia2iat(T2, I2);
				const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);
				tau2 = adjs.adjacent_tau[ad2];

				bool is_adj = false;

				const int rx2=adjs.box[ad2].x;
				const int ry2=adjs.box[ad2].y;
				const int rz2=adjs.box[ad2].z;

					
				dtau = tau2 - tau1;
				distance = dtau.norm2() * pow(GlobalC::ucell.lat0,2);
				// this rcut is in order to make nnr consistent 
				// with other matrix.
				rcut = pow(GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut(),2);
				if(distance < rcut) is_adj = true;
				else if(distance >= rcut)
				{
					for (int ad0 = 0; ad0 < adjs.adj_num+1; ++ad0)
					{
						const int T0 = adjs.ntype[ad0];
						//const int I0 = GlobalC::GridD.getNatom(ad0);
						//const int T0 = RA.info[iat1][ad0][3];
						//const int I0 = RA.info[iat1][ad0][4];
						//const int iat0 = GlobalC::ucell.itia2iat(T0, I0);
						//const int start0 = GlobalC::ucell.itiaiw2iwt(T0, I0, 0);

						tau0 = adjs.adjacent_tau[ad0];
						dtau1 = tau0 - tau1;
						dtau2 = tau0 - tau2;

						const double distance1 = dtau1.norm2() * pow(GlobalC::ucell.lat0,2);
						const double distance2 = dtau2.norm2() * pow(GlobalC::ucell.lat0,2);

						rcut1 = pow(GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max(),2);
						rcut2 = pow(GlobalC::ORB.Phi[T2].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max(),2);

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
					for (int ad0=0; ad0 < adjs.adj_num+1 ; ++ad0)
					{
						const int T0 = adjs.ntype[ad0];
						const int I0 = adjs.natom[ad0];
						const int iat = GlobalC::ucell.itia2iat(T0,I0);

						// mohan add 2010-12-19
						if( GlobalC::ucell.infoNL.nproj[T0] == 0) continue;

						//const int I0 = GlobalC::GridD.getNatom(ad0);
						//const int start0 = GlobalC::ucell.itiaiw2iwt(T0, I0, 0);
						tau0 = adjs.adjacent_tau[ad0];

						dtau1 = tau0 - tau1;
						dtau2 = tau0 - tau2;
						const double distance1 = dtau1.norm2() * pow(GlobalC::ucell.lat0,2);
						const double distance2 = dtau2.norm2() * pow(GlobalC::ucell.lat0,2);

						// seems a bug here!! mohan 2011-06-17
						rcut1 = pow(GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max(),2);
						rcut2 = pow(GlobalC::ORB.Phi[T2].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max(),2);

						if(distance1 >= rcut1 || distance2 >= rcut2)
						{
							continue;
						}
						//const Atom* atom0 = &GlobalC::ucell.atoms[T0];
						const int rx0=adjs.box[ad0].x;
						const int ry0=adjs.box[ad0].y;
						const int rz0=adjs.box[ad0].z;
						key_tuple key1(iat1,-rx0,-ry0,-rz0);
						key_tuple key2(iat2,rx2-rx0,ry2-ry0,rz2-rz0);
						
						std::unordered_map<int,std::vector<double>> *nlm_cur1_e; //left hand side, for energy
						std::unordered_map<int,std::vector<std::vector<double>>> *nlm_cur1_f; //lhs, for force
						std::unordered_map<int,std::vector<double>> *nlm_cur2_e; //rhs, for energy
						std::unordered_map<int,std::vector<std::vector<double>>> *nlm_cur2_f; //rhs, for force

						if(!calc_deri)
						{
							nlm_cur1_e = &nlm_tot[iat][key1];
							nlm_cur2_e = &nlm_tot[iat][key2];
						}
						else
						{
							nlm_cur1_f = &nlm_tot1[iat][key1];
							nlm_cur2_f = &nlm_tot1[iat][key2];
						}
				
						int nnr_inner = 0;
						
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

								if(!calc_deri)
								{
									std::vector<double> nlm_1=(*nlm_cur1_e)[iw1_all];
									std::vector<double> nlm_2=(*nlm_cur2_e)[iw2_all];
									if(GlobalV::NSPIN==4)
									{
										std::complex<double> nlm_tmp = ModuleBase::ZERO;
										int is0 = (j-j0*GlobalV::NPOL) + (k-k0*GlobalV::NPOL)*2;
										for (int no = 0; no < GlobalC::ucell.atoms[T0].ncpp.non_zero_count_soc[is0]; no++)
										{
											const int p1 = GlobalC::ucell.atoms[T0].ncpp.index1_soc[is0][no];
											const int p2 = GlobalC::ucell.atoms[T0].ncpp.index2_soc[is0][no];
											nlm_tmp += nlm_1[p1] * nlm_2[p2] * GlobalC::ucell.atoms[T0].ncpp.d_so(is0, p2, p1);
										}
										this->LM->Hloc_fixedR_soc[nnr+nnr_inner] += nlm_tmp;
									}
									else
									{
										double nlm_tmp = 0.0;
										const int nproj = GlobalC::ucell.infoNL.nproj[T0];
										int ib = 0;
										for (int nb = 0; nb < nproj; nb++)
										{
											const int L0 = GlobalC::ucell.infoNL.Beta[T0].Proj[nb].getL();
											for(int m=0;m<2*L0+1;m++)
											{
												if(nlm_1[ib]!=0.0 && nlm_2[ib]!=0.0)
												{
													nlm_tmp += nlm_1[ib]*nlm_2[ib]*GlobalC::ucell.atoms[T0].ncpp.dion(nb,nb);
												}
												ib+=1;
											}
										}
										assert(ib==nlm_1.size());

										if(GlobalV::GAMMA_ONLY_LOCAL)
										{
											// mohan add 2010-12-20
											if( nlm_tmp!=0.0 )
											{
                                                this->LM->set_HSgamma(iw1_all, iw2_all, nlm_tmp, NLloc);//N stands for nonlocal.
											}
										}
										else
										{
											if( nlm_tmp!=0.0 )
											{
												NLloc[nnr+nnr_inner] += nlm_tmp;
											}
										}
									}
								}// calc_deri
								else // calculate the derivative
								{
									if(GlobalV::GAMMA_ONLY_LOCAL)
									{
										double nlm[3]={0,0,0};

										// sum all projectors for one atom.
										std::vector<double> nlm_1 = (*nlm_cur1_f)[iw1_all][0];
										std::vector<std::vector<double>> nlm_2;
										nlm_2.resize(3);
										for(int i=0;i<3;i++)
										{
											nlm_2[i] = (*nlm_cur2_f)[iw2_all][i+1];
										}

										assert(nlm_1.size()==nlm_2[0].size());

										const int nproj = GlobalC::ucell.infoNL.nproj[T0];
										int ib = 0;
										for (int nb = 0; nb < nproj; nb++)
										{
											const int L0 = GlobalC::ucell.infoNL.Beta[T0].Proj[nb].getL();
											for(int m=0;m<2*L0+1;m++)
											{
												for(int ir=0;ir<3;ir++)
												{
													nlm[ir] += nlm_2[ir][ib]*nlm_1[ib]*GlobalC::ucell.atoms[T0].ncpp.dion(nb,nb);
												}
												ib+=1;
											}
										}
										assert(ib==nlm_1.size());
										this->LM->set_force (iw1_all, iw2_all, nlm[0], nlm[1], nlm[2], 'N');
									}
									else
									{
										// mohan change the order on 2011-06-17
										// origin: < psi1 | beta > < beta | dpsi2/dtau >
										//now: < psi1/dtau | beta > < beta | psi2 >
										double nlm[3]={0,0,0};

										// sum all projectors for one atom.
										std::vector<double> nlm_1 = (*nlm_cur2_f)[iw2_all][0];
										std::vector<std::vector<double>> nlm_2;
										nlm_2.resize(3);
										for(int i=0;i<3;i++)
										{
											nlm_2[i] = (*nlm_cur1_f)[iw1_all][i+1];
										}

										assert(nlm_1.size()==nlm_2[0].size());

										const int nproj = GlobalC::ucell.infoNL.nproj[T0];
										int ib = 0;
										for (int nb = 0; nb < nproj; nb++)
										{
											const int L0 = GlobalC::ucell.infoNL.Beta[T0].Proj[nb].getL();
											for(int m=0;m<2*L0+1;m++)
											{
												for(int ir=0;ir<3;ir++)
												{
													nlm[ir] += nlm_2[ir][ib]*nlm_1[ib]*GlobalC::ucell.atoms[T0].ncpp.dion(nb,nb);
												}
												ib+=1;
											}
										}
										assert(ib==nlm_1.size());

										this->LM->DHloc_fixedR_x[nnr+nnr_inner] += nlm[0];
										this->LM->DHloc_fixedR_y[nnr+nnr_inner] += nlm[1];
										this->LM->DHloc_fixedR_z[nnr+nnr_inner] += nlm[2];
									}
								}//!calc_deri
								nnr_inner++;
							}// k
						} // j 
					} // ad0

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
							total_nnr++;
							nnr++;
						}
					}
				}// end is_adj
			} // ad2
		} // I1
	} // T1
#ifdef _OPENMP
}
#endif
	if(!GlobalV::GAMMA_ONLY_LOCAL)
	{
		if( total_nnr!=pv->nnr)
		{
			GlobalV::ofs_running << " nr="  << total_nnr << std::endl;
			GlobalV::ofs_running << " pv->nnr=" << pv->nnr << std::endl;
			ModuleBase::WARNING_QUIT("LCAO_gen_fixedH::build_Nonlocal_mu_new","nnr!=LNNR.nnr");
		}
	}

	ModuleBase::timer::tick ("LCAO_gen_fixedH","b_NL_mu_new");
	return;
}

void LCAO_gen_fixedH::build_Nonlocal_beta_new(double* HSloc) //update by liuyu 2021-04-07
{
    ModuleBase::TITLE("LCAO_gen_fixedH","b_NL_beta_new");
    ModuleBase::timer::tick ("LCAO_gen_fixedH","b_NL_beta_new");

	const Parallel_Orbitals* pv = this->LM->ParaV;

#ifdef __MKL
    const int mkl_threads = mkl_get_max_threads();
    mkl_set_num_threads(1);
#endif

    const std::vector<AdjacentAtomInfo> adjs_all = GlobalC::GridD.get_adjs(GlobalC::ucell);

#ifdef _OPENMP
    #pragma omp parallel
    {
        double* Nonlocal_thread;
        Nonlocal_thread = new double[pv->nloc];
        ModuleBase::GlobalFunc::ZEROS(Nonlocal_thread, pv->nloc);
        #pragma omp for schedule(dynamic)
#endif
        for(int iat=0; iat<GlobalC::ucell.nat; iat++)
        {
            const int T0 = GlobalC::ucell.iat2it[iat];
            const int I0 = GlobalC::ucell.iat2ia[iat];
            Atom* atom0 = &GlobalC::ucell.atoms[T0];

            //=======================================================
            //Step1:
            //saves <beta|psi>, where beta runs over L0,M0 on atom I0
            //and psi runs over atomic basis sets on the current core
            //=======================================================
            #ifdef _OPENMP
                std::vector<std::unordered_map<int,std::vector<double>>> nlm_tot_thread;
                nlm_tot_thread.resize(adjs_all[iat].adj_num + 1);
            #else 
                std::vector<std::unordered_map<int,std::vector<double>>> nlm_tot;
                nlm_tot.resize(adjs_all[iat].adj_num + 1);
            #endif 

            const ModuleBase::Vector3<double> tau0 = atom0->tau[I0];
            const double Rcut_Beta = GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();

            //outermost loop : all adjacent atoms
            for(int ad_count=0; ad_count < adjs_all[iat].adj_num + 1; ad_count++)
            {
                const int T1 = adjs_all[iat].ntype[ad_count];
                const int I1 = adjs_all[iat].natom[ad_count];
                const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
                const double Rcut_AO1 = GlobalC::ORB.Phi[T1].getRcut();
                const ModuleBase::Vector3<double> tau1 = adjs_all[iat].adjacent_tau[ad_count];
                const Atom* atom1 = &GlobalC::ucell.atoms[T1];
                const int nw1_tot = atom1->nw*GlobalV::NPOL;

                #ifdef _OPENMP
                    nlm_tot_thread[ad_count].clear();
                #else 
                    nlm_tot[ad_count].clear();
                #endif 

                //middle loop : atomic basis on current processor (either row or column)
                const double dist1 = (tau1-tau0).norm() * GlobalC::ucell.lat0;
                if (dist1 > Rcut_Beta + Rcut_AO1)
                {
                    continue;
                }

                for(int iw1=0; iw1<nw1_tot; ++iw1)
                {
                    const int iw1_all = start1 + iw1;
                    const int iw1_local = pv->global2local_row(iw1_all);
                    const int iw2_local = pv->global2local_col(iw1_all);

                    if(iw1_local < 0 && iw2_local < 0) continue;

                    const int iw1_0 = iw1/GlobalV::NPOL;
                    std::vector<std::vector<double>> nlm;
                    //2D, but first dimension is only 1 here
                    //for force, the right hand side is the gradient
                    //and the first dimension is then 3
                    //inner loop : all projectors (L0,M0)
                    GlobalC::UOT.snap_psibeta_half(
                        GlobalC::ORB,
                        GlobalC::ucell.infoNL,
                        nlm, tau1, T1,
                        atom1->iw2l[ iw1_0 ], // L1
                        atom1->iw2m[ iw1_0 ], // m1
                        atom1->iw2n[ iw1_0 ], // N1
                        GlobalC::ucell.atoms[T0].tau[I0], T0, 0); //R0,T0

                    #ifdef _OPENMP
                        nlm_tot_thread[ad_count].insert({iw1_all,nlm[0]});
                    #else 
                        nlm_tot[ad_count].insert({iw1_all,nlm[0]});
                    #endif 
                }//end iw
            }//end ad

            //=======================================================
            //Step2:
            //calculate sum_(L0,M0) beta<psi_i|beta><beta|psi_j>
            //and accumulate the value to Hloc_fixed(i,j)
            //=======================================================
            for(int ad1_count=0; ad1_count < adjs_all[iat].adj_num + 1; ad1_count++)
            {
                const int T1 = adjs_all[iat].ntype[ad1_count];
                const int I1 = adjs_all[iat].natom[ad1_count];
                const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
                const ModuleBase::Vector3<double> tau1 = adjs_all[iat].adjacent_tau[ad1_count];
                const Atom* atom1 = &GlobalC::ucell.atoms[T1];
                const int nw1_tot = atom1->nw*GlobalV::NPOL;
                const double Rcut_AO1 = GlobalC::ORB.Phi[T1].getRcut();

                for (int ad2_count=0; ad2_count < adjs_all[iat].adj_num + 1; ad2_count++)
                {
                    const int T2 = adjs_all[iat].ntype[ad2_count];
                    const int I2 = adjs_all[iat].natom[ad2_count];
                    const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);
                    const ModuleBase::Vector3<double> tau2 = adjs_all[iat].adjacent_tau[ad2_count];
                    const Atom* atom2 = &GlobalC::ucell.atoms[T2];
                    const int nw2_tot = atom2->nw*GlobalV::NPOL;
                    const double Rcut_AO2 = GlobalC::ORB.Phi[T2].getRcut();
                    const double dist1 = (tau1-tau0).norm() * GlobalC::ucell.lat0;
                    const double dist2 = (tau2-tau0).norm() * GlobalC::ucell.lat0;

                    if (dist1 > Rcut_Beta + Rcut_AO1
                            || dist2 > Rcut_Beta + Rcut_AO2)
                    {
                        continue;
                    }

                    for(int iw1=0; iw1<nw1_tot; ++iw1)
                    {
                        const int iw1_all = start1 + iw1;
                        const int iw1_local = pv->global2local_row(iw1_all);
                        if(iw1_local < 0) continue;
                        const int iw1_0 = iw1/GlobalV::NPOL;
                        for(int iw2=0; iw2<nw2_tot; ++iw2)
                        {
                            const int iw2_all = start2 + iw2;
                            const int iw2_local = pv->global2local_col(iw2_all);
                            if(iw2_local < 0) continue;
                            const int iw2_0 = iw2/GlobalV::NPOL;
                            #ifdef _OPENMP
                                std::vector<double> nlm1 = nlm_tot_thread[ad1_count][iw1_all];
                                std::vector<double> nlm2 = nlm_tot_thread[ad2_count][iw2_all];
                            #else 
                                std::vector<double> nlm1 = nlm_tot[ad1_count][iw1_all];
                                std::vector<double> nlm2 = nlm_tot[ad2_count][iw2_all];
                            #endif 

                            assert(nlm1.size()==nlm2.size());
                            #ifdef _OPENMP
                                double nlm_thread=0.0;
                            #else 
                                double nlm=0.0;
                            #endif
                            const int nproj = GlobalC::ucell.infoNL.nproj[T0];
                            int ib = 0;
                            for(int nb = 0; nb < nproj; nb++)
                            {
                                const int L0 = GlobalC::ucell.infoNL.Beta[T0].Proj[nb].getL();
                                for(int m=0;m<2*L0+1;m++)
                                {
                                    #ifdef _OPENMP
                                        nlm_thread += nlm1[ib]*nlm2[ib]*GlobalC::ucell.atoms[T0].ncpp.dion(nb,nb);
                                    #else 
                                        nlm += nlm1[ib]*nlm2[ib]*GlobalC::ucell.atoms[T0].ncpp.dion(nb,nb);
                                    #endif
                                    ib+=1;
                                }
                            }
                            assert(ib==nlm1.size());

                            const int ir = pv->global2local_row(iw1_all);
                            const int ic = pv->global2local_col(iw2_all);
                            long index=0;
                            if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER())
                            {
                                index=ic*pv->nrow+ir;
                            }
                            else
                            {
                                index=ir*pv->ncol+ic;
                            }
#ifdef _OPENMP
                            Nonlocal_thread[index] += nlm_thread;
#else
                            this->LM->set_HSgamma(iw1_all, iw2_all, nlm, HSloc);
#endif
                        }//iw2
                    }//iw1
                }//ad2
            }//ad1
        }//end iat

        #ifdef _OPENMP
            #pragma omp critical(cal_nonlocal)
            {
                for(int i=0; i<pv->nloc; i++)
                {
                    this->LM->Hloc_fixed[i] += Nonlocal_thread[i];
                }
            }
            delete[] Nonlocal_thread;
        #endif
#ifdef _OPENMP
    }
#endif

#ifdef __MKL
    mkl_set_num_threads(mkl_threads);
#endif
	
    ModuleBase::timer::tick ("LCAO_gen_fixedH","b_NL_beta_new");
	return;
}
