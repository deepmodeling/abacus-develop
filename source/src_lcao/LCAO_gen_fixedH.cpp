#include "LCAO_gen_fixedH.h"
#include "../src_pw/global.h"
#include "../src_pw/wavefunc.h"
#include "global_fp.h"
#include <vector>
#include <unordered_map>
#include <map>
#include "../module_base/timer.h"

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


void LCAO_gen_fixedH::calculate_NL_no(void)
{
    ModuleBase::TITLE("LCAO_gen_fixedH","calculate_NL_no");

	// PLEASE rebuild the following two functions,
	// 'build_Nonlocal_beta' and  'build_Nonlocal_mu',
	// because the two functions are extremely time consuming
	// for small systems, especially for multiple-k points
	// mohan note 2021-03-23

	if(GlobalV::GAMMA_ONLY_LOCAL)
	{
	  	//for gamma only.
		if(GlobalV::NSPIN!=4 && GlobalV::vnl_method == 1)
		{
  			this->build_Nonlocal_beta_new();
		}
		else
		{
  			this->build_Nonlocal_beta(false);
		}
	}
	else
	{
		// can work for gamma also,
		// only if search_radius is 
		// (Phi.rcutmax + Beta.rcutmax)*2.
		// check in sltk_atom_arrange.
    	if(GlobalV::NSPIN!=4 && GlobalV::vnl_method == 1)
		{
			this->build_Nonlocal_mu_new(false);
		}
		else
		{
			this->build_Nonlocal_mu(false);
		}

//		this->test_Nonlocal();
	}

    return;
}


void LCAO_gen_fixedH::calculate_T_no(void)
{
    ModuleBase::TITLE("LCAO_gen_fixedH","calculate_T_no");
    this->build_ST_new('T', false, GlobalC::ucell);
    return;
}

void LCAO_gen_fixedH::calculate_S_no(void)
{
    ModuleBase::TITLE("LCAO_gen_fixedH", "calculate_S_no");
    ModuleBase::timer::tick("LCAO_gen_fixedH","calculate_S_no");
	this->build_ST_new('S', false, GlobalC::ucell);
    ModuleBase::timer::tick("LCAO_gen_fixedH","calculate_S_no");
    return;
}


//liaochen modify interface 2010-3-22
void LCAO_gen_fixedH::build_ST_new(const char& dtype, const bool& calc_deri, const UnitCell_pseudo &ucell)
{
    ModuleBase::TITLE("LCAO_gen_fixedH","build_ST_new");

    //array to store data
    double olm[3]={0.0,0.0,0.0};
    int nnr = 0; // used onlyh for k points.

    const Parallel_Orbitals* pv = this->LM->ParaV;

    //\sum{T} e**{ikT} <\phi_{ia}|d\phi_{k\beta}(T)>
	ModuleBase::Vector3<double> tau1, tau2, dtau;
	ModuleBase::Vector3<double> dtau1, dtau2, tau0;
    for (int T1=0; T1<ucell.ntype; ++T1)
    {
		Atom* atom1 = &ucell.atoms[T1];
        for (int I1=0; I1<atom1->na; ++I1)
        {
			tau1 = atom1->tau[I1];
            //GlobalC::GridD.Find_atom(tau1);
            GlobalC::GridD.Find_atom(ucell, tau1, T1, I1);

            for (int ad = 0; ad < GlobalC::GridD.getAdjacentNum()+1; ++ad)
            {
                const int T2 = GlobalC::GridD.getType(ad);
				const int I2 = GlobalC::GridD.getNatom(ad);
				Atom* atom2 = &ucell.atoms[T2];
				tau2 = GlobalC::GridD.getAdjacentTau(ad);
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
							// in build_Non... use trace_loc_row
                            // and trace_loc_col directly,
                            if (!pv->in_this_processor(iw1_all, iw2_all))
							{
								++iw2_all;
								continue;
							}

							olm[0] = olm[1] = olm[2] = 0.0;

							std::complex<double> olm1[4]={ModuleBase::ZERO, ModuleBase::ZERO, ModuleBase::ZERO, ModuleBase::ZERO};
							std::complex<double> *olm2 = &olm1[0];
							if(!calc_deri)
							{
								// PLEASE use UOT as an input parameter of this subroutine
								// mohan add 2021-03-30
								GlobalC::UOT.snap_psipsi( GlobalC::ORB, olm, 0, dtype, tau1, 
										T1, L1, m1, N1, GlobalC::GridD.getAdjacentTau(ad), 
										T2, L2, m2, N2, GlobalV::NSPIN,
										olm2//for soc
										);

								if(GlobalV::GAMMA_ONLY_LOCAL)
								{
									// mohan add 2010-06-29
									// set the value in Hloc and Sloc
									// according to trace_loc_row and trace_loc_col
									// the last paramete: 1 for Sloc, 2 for Hloc
									// and 3 for Hloc_fixed.
									this->LM->set_HSgamma(iw1_all, iw2_all, olm[0], dtype);
								}
								else // k point algorithm
								{
									// mohan add 2010-10
									// set the values in SlocR and Hloc_fixedR.
									// which is a 1D array.
									if(dtype=='S')
									{
                                        if (GlobalV::NSPIN != 4) this->LM->SlocR[nnr] = olm[0];
                                        else
										{//only has diagonal term here.
												int is = (jj-jj0*GlobalV::NPOL) + (kk-kk0*GlobalV::NPOL)*2;
											this->LM->SlocR_soc[nnr] = olm1[is];
                                        }
                                    }
									else if(dtype=='T')
									{
										if(GlobalV::NSPIN!=4) this->LM->Hloc_fixedR[nnr] = olm[0];// <phi|kin|d phi>
										else
										{//only has diagonal term here.
												int is = (jj-jj0*GlobalV::NPOL) + (kk-kk0*GlobalV::NPOL)*2;
											this->LM->Hloc_fixedR_soc[nnr] = olm1[is];
                                        }
                                    }
									++nnr;
								}
							}
							else // calculate the derivative
							{
								GlobalC::UOT.snap_psipsi( GlobalC::ORB, olm, 1, dtype, 
									tau1, T1, L1, m1, N1,
									GlobalC::GridD.getAdjacentTau(ad), T2, L2, m2, N2, GlobalV::NSPIN
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
					for (int ad0=0; ad0 < GlobalC::GridD.getAdjacentNum()+1; ++ad0)
					{
						const int T0 = GlobalC::GridD.getType(ad0);
						//const int I0 = GlobalC::GridD.getNatom(ad0);
						//const int iat0 = ucell.itia2iat(T0, I0);
						//const int start0 = ucell.itiaiw2iwt(T0, I0, 0);
						tau0 = GlobalC::GridD.getAdjacentTau(ad0);
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
							const int mu = pv->trace_loc_row[start1+jj];
							if(mu<0)continue; 
							for(int kk=0; kk<atom2->nw * GlobalV::NPOL; ++kk)
							{
								const int nu = pv->trace_loc_col[start2+kk];
								if(nu<0)continue;
								++nnr;
							}//kk
						}//jj
					}
				}//distance
			}// ad
		}// I1
	}// T1

	if(!GlobalV::GAMMA_ONLY_LOCAL)
	{
		if(nnr != pv->nnr)
		{
			std::cout << " nnr=" << nnr << " LNNR.nnr=" << pv->nnr << std::endl;
			GlobalV::ofs_running << " nnr=" << nnr << " LNNR.nnr=" << pv->nnr << std::endl;
			ModuleBase::WARNING_QUIT("LCAO_gen_fixedH::build_ST_new","nnr != LNNR.nnr");
		}
	}

    return;
}

void LCAO_gen_fixedH::test_Nonlocal()
{
    const Parallel_Orbitals* pv = this->LM->ParaV;
    int nnr = 0;
	ModuleBase::Vector3<double> tau1, tau2, dtau_12, tau0, dtau_10, dtau_20;
	double distance = 0.0;
	double rcut = 0.0;

//	double* vnltest = new double[pv->nloc];	
//	ModuleBase::GlobalFunc::ZEROS(vnltest, pv->nloc);

	// psi1
	double sum = 0.0;
	int count = 0;
    for (int T1 = 0; T1 < GlobalC::ucell.ntype; T1++)
    {
		const Atom* atom1 = &GlobalC::ucell.atoms[T1];
        for (int I1 =0; I1< atom1->na; I1++)
        {
            //GlobalC::GridD.Find_atom( atom1->tau[I1] );
            GlobalC::GridD.Find_atom(GlobalC::ucell, atom1->tau[I1] ,T1, I1);
			//const int iat1 = GlobalC::ucell.itia2iat(T1, I1);
			const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
            tau1 = atom1->tau[I1];

			// psi2
            for (int ad2=0; ad2<GlobalC::GridD.getAdjacentNum()+1 ; ad2++)
			{
				const int T2 = GlobalC::GridD.getType(ad2);
				const Atom* atom2 = &GlobalC::ucell.atoms[T2];
                
				const int I2 = GlobalC::GridD.getNatom(ad2);
				//const int iat2 = GlobalC::ucell.itia2iat(T2, I2);
                const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);
                tau2 = GlobalC::GridD.getAdjacentTau(ad2);

				dtau_12 = tau2 - tau1;
				distance = dtau_12.norm() * GlobalC::ucell.lat0;
				rcut = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut();

				if(distance >= rcut)
				{
					for (int j=0; j<atom1->nw*GlobalV::NPOL; j++)
					{
						int j0 = j/GlobalV::NPOL;
						const int iw1_all = start1 + j;
						const int mu = pv->trace_loc_row[iw1_all];
						if(mu < 0)continue; 
						for (int k=0; k<atom2->nw*GlobalV::NPOL; k++)
						{
							int k0 = k/GlobalV::NPOL;
							const int iw2_all = start2 + k;
							const int nu = pv->trace_loc_col[iw2_all];						
							if(nu < 0)continue;
							for (int ad0=0; ad0 < GlobalC::GridD.getAdjacentNum()+1 ; ad0++)
							{
								const int T0 = GlobalC::GridD.getType(ad0);
								if( GlobalC::ucell.infoNL.nproj[T0] == 0) continue; 
								//const int I0 = GlobalC::GridD.getNatom(ad0);
								//const int start0 = GlobalC::ucell.itiaiw2iwt(T0, I0, 0);
								tau0 = GlobalC::GridD.getAdjacentTau(ad0);

								dtau_10 = tau0 - tau1;
								dtau_20 = tau0 - tau2;
								double distance1 = dtau_10.norm() * GlobalC::ucell.lat0;
								double distance2 = dtau_20.norm() * GlobalC::ucell.lat0;

								double rcut_10 = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T0].getRcut();
								double rcut_20 = GlobalC::ORB.Phi[T2].getRcut() + GlobalC::ORB.Phi[T0].getRcut();

								if(distance1 < rcut_10 && distance2 < rcut_20)
								{
									++count;
									double nlm[3]={0,0,0};
									GlobalC::UOT.snap_psibeta(
											GlobalC::ORB,
											GlobalC::ucell.infoNL,
											nlm, 0, tau1, T1,
											atom1->iw2l[ j0 ], // L1
											atom1->iw2m[ j0 ], // m1
											atom1->iw2n[ j0 ], // N1
											tau2, T2,
											atom2->iw2l[ k0 ], // L2
											atom2->iw2m[ k0 ], // m2
											atom2->iw2n[ k0 ], // n2
											tau0, T0, GlobalC::ucell.atoms[T0].dion, GlobalV::NSPIN,
											GlobalC::ucell.atoms[T0].d_so, // mohan  add 2021-05-07
											GlobalC::ucell.atoms[T0].non_zero_count_soc[0], // 0 stands for spin
											GlobalC::ucell.atoms[T0].index1_soc[0],
											GlobalC::ucell.atoms[T0].index2_soc[0],
											GlobalC::ucell.atoms[T0].nproj_soc
											);
									
									//vnltest[ mu * pv->ncol + nu ] += nlm[0];
									sum += abs( nlm[0] );
								}// distance
							} // ad0
							++nnr;
						}// k
					} // j 
				}// end distance
				//----------------------------------------------------------------------------------
			} // ad2
		} // I1
	} // T1

	GlobalV::ofs_running << " not included Vnl std::pair = " << count << std::endl;
	GlobalV::ofs_running << " sum for the correction = " << sum << std::endl;

	/*
	std::cout << " correction for Vnl matrix " << std::endl;
	for(int i=0; i<GlobalV::NLOCAL; i++)
	{
		for(int j=0; j<GlobalV::NLOCAL; j++)
		{
			double a = vnltest[i*pv->ncol+j];
			if( abs(a) > 1.0e-6 )
			{
				std::cout << std::setw(15) << vnltest[i*pv->ncol+j];
			}
			else
			{
				std::cout << std::setw(15) << "0";
			}
		}
		std::cout << std::endl;
	}
	delete[] vnltest;
	*/
		
	return;
}

typedef std::tuple<int,int,int,int> key_tuple;

#include "record_adj.h" //mohan add 2012-07-06
void LCAO_gen_fixedH::build_Nonlocal_mu_new(const bool &calc_deri)
{
    ModuleBase::TITLE("LCAO_gen_fixedH","build_Nonlocal_mu_new");
    ModuleBase::timer::tick("LCAO_gen_fixedH", "build_Nonlocal_mu_new");
    const Parallel_Orbitals* pv = this->LM->ParaV;

	// < phi1 | beta > < beta | phi2 >
	// phi1 is within the unitcell.
	// while beta is in the supercell.
	// while phi2 is in the supercell.


	//Step 1 : generate <psi|beta>
	//type of atom; distance; atomic basis; projectors
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
	for(int iat=0;iat<GlobalC::ucell.nat;iat++)
	{
		const int it = GlobalC::ucell.iat2it[iat];
		const int ia = GlobalC::ucell.iat2ia[iat];

		const double Rcut_Beta = GlobalC::ucell.infoNL.Beta[it].get_rcut_max();
		const ModuleBase::Vector3<double> tau = GlobalC::ucell.atoms[it].tau[ia];
        GlobalC::GridD.Find_atom(GlobalC::ucell, tau ,it, ia);

		if(!calc_deri)
		{
			nlm_tot[iat].clear();
		}
		else
		{
			nlm_tot1[iat].clear();
		}

		for (int ad=0; ad<GlobalC::GridD.getAdjacentNum()+1 ; ++ad)
		{
			const int T1 = GlobalC::GridD.getType(ad);
			const int I1 = GlobalC::GridD.getNatom(ad);
			const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
			const double Rcut_AO1 = GlobalC::ORB.Phi[T1].getRcut();

			const ModuleBase::Vector3<double> tau1 = GlobalC::GridD.getAdjacentTau(ad);
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
				const int iw1_local = pv->trace_loc_row[iw1_all];
				const int iw2_local = pv->trace_loc_col[iw1_all];
				if(iw1_local < 0 && iw2_local < 0)continue;
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
			const int rx1=GlobalC::GridD.getBox(ad).x;
			const int ry1=GlobalC::GridD.getBox(ad).y;
			const int rz1=GlobalC::GridD.getBox(ad).z;
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
	for (int T1 = 0; T1 < GlobalC::ucell.ntype; ++T1)
	{
		const Atom* atom1 = &GlobalC::ucell.atoms[T1];
		for (int I1 =0; I1< atom1->na; ++I1)
		{
			//GlobalC::GridD.Find_atom( atom1->tau[I1] );
			GlobalC::GridD.Find_atom(GlobalC::ucell, atom1->tau[I1] ,T1, I1);
			const int iat1 = GlobalC::ucell.itia2iat(T1, I1);
			const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
			tau1 = atom1->tau[I1];

			// psi2
			for (int ad2=0; ad2<GlobalC::GridD.getAdjacentNum()+1; ++ad2)
			{
				const int T2 = GlobalC::GridD.getType(ad2);
				const Atom* atom2 = &GlobalC::ucell.atoms[T2];
				
				const int I2 = GlobalC::GridD.getNatom(ad2);
				const int iat2 = GlobalC::ucell.itia2iat(T2, I2);
				const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);
				tau2 = GlobalC::GridD.getAdjacentTau(ad2);

				bool is_adj = false;

				const int rx2=GlobalC::GridD.getBox(ad2).x;
				const int ry2=GlobalC::GridD.getBox(ad2).y;
				const int rz2=GlobalC::GridD.getBox(ad2).z;

					
				dtau = tau2 - tau1;
				distance = dtau.norm2() * pow(GlobalC::ucell.lat0,2);
				// this rcut is in order to make nnr consistent 
				// with other matrix.
				rcut = pow(GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut(),2);
				if(distance < rcut) is_adj = true;
				else if(distance >= rcut)
				{
					for (int ad0 = 0; ad0 < GlobalC::GridD.getAdjacentNum()+1; ++ad0)
					{
						const int T0 = GlobalC::GridD.getType(ad0);
						//const int I0 = GlobalC::GridD.getNatom(ad0);
						//const int T0 = RA.info[iat1][ad0][3];
						//const int I0 = RA.info[iat1][ad0][4];
						//const int iat0 = GlobalC::ucell.itia2iat(T0, I0);
						//const int start0 = GlobalC::ucell.itiaiw2iwt(T0, I0, 0);

						tau0 = GlobalC::GridD.getAdjacentTau(ad0);
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
					for (int ad0=0; ad0 < GlobalC::GridD.getAdjacentNum()+1 ; ++ad0)
					{
						const int T0 = GlobalC::GridD.getType(ad0);
						const int I0 = GlobalC::GridD.getNatom(ad0);
						const int iat = GlobalC::ucell.itia2iat(T0,I0);

						// mohan add 2010-12-19
						if( GlobalC::ucell.infoNL.nproj[T0] == 0) continue;

						//const int I0 = GlobalC::GridD.getNatom(ad0);
						//const int start0 = GlobalC::ucell.itiaiw2iwt(T0, I0, 0);
						tau0 = GlobalC::GridD.getAdjacentTau(ad0);

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
						const int rx0=GlobalC::GridD.getBox(ad0).x;
						const int ry0=GlobalC::GridD.getBox(ad0).y;
						const int rz0=GlobalC::GridD.getBox(ad0).z;
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
							const int mu = pv->trace_loc_row[iw1_all];
							if(mu < 0)continue; 

							// fix a serious bug: atom2[T2] -> atom2
							// mohan 2010-12-20
							for (int k=0; k<atom2->nw*GlobalV::NPOL; k++)
							{
								const int k0 = k/GlobalV::NPOL;
								const int iw2_all = start2 + k;
								const int nu = pv->trace_loc_col[iw2_all];						
								if(nu < 0)continue;

								if(!calc_deri)
								{
									std::vector<double> nlm_1=(*nlm_cur1_e)[iw1_all];
									std::vector<double> nlm_2=(*nlm_cur2_e)[iw2_all];
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
												nlm_tmp += nlm_1[ib]*nlm_2[ib]*GlobalC::ucell.atoms[T0].dion(nb,nb);
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
											// GlobalV::ofs_running << std::setw(10) << iw1_all << std::setw(10) 
											// << iw2_all << std::setw(20) << nlm[0] << std::endl; 
											this->LM->set_HSgamma(iw1_all,iw2_all,nlm_tmp,'N');//N stands for nonlocal.
										}
									}
									else
									{
										if( nlm_tmp!=0.0 )
										{
											this->LM->Hloc_fixedR[nnr+nnr_inner] += nlm_tmp;
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
													nlm[ir] += nlm_2[ir][ib]*nlm_1[ib]*GlobalC::ucell.atoms[T0].dion(nb,nb);
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
													nlm[ir] += nlm_2[ir][ib]*nlm_1[ib]*GlobalC::ucell.atoms[T0].dion(nb,nb);
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
						const int mu = pv->trace_loc_row[iw1_all];
						if(mu < 0)continue; 

						// fix a serious bug: atom2[T2] -> atom2
						// mohan 2010-12-20
						for (int k=0; k<atom2->nw*GlobalV::NPOL; k++)
						{
							const int k0 = k/GlobalV::NPOL;
							const int iw2_all = start2 + k;
							const int nu = pv->trace_loc_col[iw2_all];						
							if(nu < 0)continue;

							nnr++;
						}
					}
				}// end is_adj
			} // ad2
		} // I1
	} // T1

	if(!GlobalV::GAMMA_ONLY_LOCAL)
	{
	//		std::cout << " nr="  << nnr << std::endl;
	//		std::cout << " pv->nnr=" << pv->nnr << std::endl;
	//		GlobalV::ofs_running << " nr="  << nnr << std::endl;
	//		GlobalV::ofs_running << " pv->nnr=" << pv->nnr << std::endl;
		if( nnr!=pv->nnr)
		{
			ModuleBase::WARNING_QUIT("LCAO_gen_fixedH::build_Nonlocal_mu_new","nnr!=LNNR.nnr");
		}
	}

//	std::cout << " build_Nonlocal_mu done" << std::endl;

	ModuleBase::timer::tick ("LCAO_gen_fixedH","build_Nonlocal_mu_new");
	return;
}

void LCAO_gen_fixedH::build_Nonlocal_mu(const bool &calc_deri)
{
    ModuleBase::TITLE("LCAO_gen_fixedH","build_Nonlocal_mu");
    ModuleBase::timer::tick("LCAO_gen_fixedH", "build_Nonlocal_mu");
    const Parallel_Orbitals* pv = this->LM->ParaV;

	// < phi1 | beta > < beta | phi2 >
	// phi1 is within the unitcell.
	// while beta is in the supercell.
	// while phi2 is in the supercell.

	int nnr = 0;
	ModuleBase::Vector3<double> tau1, tau2, dtau;
	ModuleBase::Vector3<double> dtau1, dtau2, tau0;
	double distance = 0.0;
	double distance1, distance2;
	double rcut = 0.0;
	double rcut1, rcut2;
		
//	Record_adj RA;
//	RA.for_2d();

	// psi1
    for (int T1 = 0; T1 < GlobalC::ucell.ntype; ++T1)
    {
		const Atom* atom1 = &GlobalC::ucell.atoms[T1];
        for (int I1 =0; I1< atom1->na; ++I1)
        {
            //GlobalC::GridD.Find_atom( atom1->tau[I1] );
            GlobalC::GridD.Find_atom(GlobalC::ucell, atom1->tau[I1] ,T1, I1);
			//const int iat1 = GlobalC::ucell.itia2iat(T1, I1);
			const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
            tau1 = atom1->tau[I1];

			// psi2
            for (int ad2=0; ad2<GlobalC::GridD.getAdjacentNum()+1; ++ad2)
			{
				const int T2 = GlobalC::GridD.getType(ad2);
				const Atom* atom2 = &GlobalC::ucell.atoms[T2];
                
				const int I2 = GlobalC::GridD.getNatom(ad2);
				//const int iat2 = GlobalC::ucell.itia2iat(T2, I2);
                const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);
                tau2 = GlobalC::GridD.getAdjacentTau(ad2);

				bool is_adj = false;
					
				dtau = tau2 - tau1;
				distance = dtau.norm() * GlobalC::ucell.lat0;
				// this rcut is in order to make nnr consistent 
				// with other matrix.
				rcut = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut();
				if(distance < rcut) is_adj = true;
				else if(distance >= rcut)
				{
                    for (int ad0 = 0; ad0 < GlobalC::GridD.getAdjacentNum()+1; ++ad0)
                    {
						const int T0 = GlobalC::GridD.getType(ad0);
						//const int I0 = GlobalC::GridD.getNatom(ad0);
						//const int T0 = RA.info[iat1][ad0][3];
						//const int I0 = RA.info[iat1][ad0][4];
                        //const int iat0 = GlobalC::ucell.itia2iat(T0, I0);
                        //const int start0 = GlobalC::ucell.itiaiw2iwt(T0, I0, 0);

                        tau0 = GlobalC::GridD.getAdjacentTau(ad0);
                        dtau1 = tau0 - tau1;
                        dtau2 = tau0 - tau2;

                        double distance1 = dtau1.norm() * GlobalC::ucell.lat0;
                        double distance2 = dtau2.norm() * GlobalC::ucell.lat0;

                        rcut1 = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();
                        rcut2 = GlobalC::ORB.Phi[T2].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();

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
					for (int j=0; j<atom1->nw*GlobalV::NPOL; j++)
					{
						const int j0 = j/GlobalV::NPOL;//added by zhengdy-soc
						const int iw1_all = start1 + j;
						const int mu = pv->trace_loc_row[iw1_all];
						if(mu < 0)continue; 

						// fix a serious bug: atom2[T2] -> atom2
						// mohan 2010-12-20
						for (int k=0; k<atom2->nw*GlobalV::NPOL; k++)
						{
							const int k0 = k/GlobalV::NPOL;
							const int iw2_all = start2 + k;
							const int nu = pv->trace_loc_col[iw2_all];						
							if(nu < 0)continue;


							//(3) run over all projectors in nonlocal pseudopotential.
							for (int ad0=0; ad0 < GlobalC::GridD.getAdjacentNum()+1 ; ++ad0)
							{
								const int T0 = GlobalC::GridD.getType(ad0);

								// mohan add 2010-12-19
								if( GlobalC::ucell.infoNL.nproj[T0] == 0) continue; 

								//const int I0 = GlobalC::GridD.getNatom(ad0);
								//const int start0 = GlobalC::ucell.itiaiw2iwt(T0, I0, 0);
								tau0 = GlobalC::GridD.getAdjacentTau(ad0);

								dtau1 = tau0 - tau1;
								dtau2 = tau0 - tau2;
								distance1 = dtau1.norm() * GlobalC::ucell.lat0;
								distance2 = dtau2.norm() * GlobalC::ucell.lat0;

								// seems a bug here!! mohan 2011-06-17
								rcut1 = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();
								rcut2 = GlobalC::ORB.Phi[T2].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();

								if(distance1 < rcut1 && distance2 < rcut2)
								{
									//const Atom* atom0 = &GlobalC::ucell.atoms[T0];
									double nlm[3]={0,0,0};
									std::complex<double> nlm1[4]={0,0,0,0};//modified by zhengdy-soc
									std::complex<double> *nlm2 = NULL;
									if(GlobalV::NSPIN==4) nlm2 = &nlm1[0];
									if(!calc_deri)
									{
										int is0 = (j-j0*GlobalV::NPOL) + (k-k0*GlobalV::NPOL)*2;
										GlobalC::UOT.snap_psibeta(
												GlobalC::ORB,
												GlobalC::ucell.infoNL,
												nlm, 0, tau1, T1,
												atom1->iw2l[ j0 ], // L1
												atom1->iw2m[ j0 ], // m1
												atom1->iw2n[ j0 ], // N1
												tau2, T2,
												atom2->iw2l[ k0 ], // L2
												atom2->iw2m[ k0 ], // m2
												atom2->iw2n[ k0 ], // n2
												tau0, T0, GlobalC::ucell.atoms[T0].dion, GlobalV::NSPIN,
												GlobalC::ucell.atoms[T0].d_so, // mohan  add 2021-05-07
												GlobalC::ucell.atoms[T0].non_zero_count_soc[is0], // index stands for spin
												GlobalC::ucell.atoms[T0].index1_soc[is0],
												GlobalC::ucell.atoms[T0].index2_soc[is0],
												GlobalC::ucell.atoms[T0].nproj_soc,
												nlm2, is0 //for soc
												);


										if(GlobalV::GAMMA_ONLY_LOCAL)
										{
											// mohan add 2010-12-20
											if( nlm[0]!=0.0 )
											{
												// GlobalV::ofs_running << std::setw(10) << iw1_all << std::setw(10) 
												// << iw2_all << std::setw(20) << nlm[0] << std::endl; 
												this->LM->set_HSgamma(iw1_all,iw2_all,nlm[0],'N');//N stands for nonlocal.
											}
										}
										else
										{
											if(GlobalV::NSPIN!=4) this->LM->Hloc_fixedR[nnr] += nlm[0];
											else
											{
												int is = (j-j0*GlobalV::NPOL) + (k-k0*GlobalV::NPOL)*2;
												this->LM->Hloc_fixedR_soc[nnr] += nlm1[is];
											}
										}
									}// calc_deri
									else // calculate the derivative
									{
										if(GlobalV::GAMMA_ONLY_LOCAL)
										{
											GlobalC::UOT.snap_psibeta(
													GlobalC::ORB,
													GlobalC::ucell.infoNL,
													nlm, 1, 
													tau1, 
													T1,
													atom1->iw2l[ j0 ], // L1
													atom1->iw2m[ j0 ], // m1
													atom1->iw2n[ j0 ], // N1
													tau2, 
													T2,
													atom2->iw2l[ k0 ], // L2
													atom2->iw2m[ k0 ], // m2
													atom2->iw2n[ k0 ], // n2
													tau0, T0, GlobalC::ucell.atoms[T0].dion, GlobalV::NSPIN,
													GlobalC::ucell.atoms[T0].d_so, // mohan  add 2021-05-07
													GlobalC::ucell.atoms[T0].non_zero_count_soc[0], // index stands for spin
													GlobalC::ucell.atoms[T0].index1_soc[0],
													GlobalC::ucell.atoms[T0].index2_soc[0],
													GlobalC::ucell.atoms[T0].nproj_soc
													);

											// sum all projectors for one atom.
											this->LM->set_force (iw1_all, iw2_all,	nlm[0], nlm[1], nlm[2], 'N');
										}
										else
										{
											// mohan change the order on 2011-06-17
											// origin: < psi1 | beta > < beta | dpsi2/dtau >
											//now: < psi1/dtau | beta > < beta | psi2 >
											GlobalC::UOT.snap_psibeta(
													GlobalC::ORB,
													GlobalC::ucell.infoNL,
													nlm, 1, 
													tau2, 
													T2,
													atom2->iw2l[ k0 ], // L2
													atom2->iw2m[ k0 ], // m2
													atom2->iw2n[ k0 ], // n2
													tau1, 
													T1,
													atom1->iw2l[ j0 ], // L1
													atom1->iw2m[ j0 ], // m1
													atom1->iw2n[ j0 ], // N1
													tau0, T0, GlobalC::ucell.atoms[T0].dion, GlobalV::NSPIN,
													GlobalC::ucell.atoms[T0].d_so, // mohan  add 2021-05-07
													GlobalC::ucell.atoms[T0].non_zero_count_soc[0], // index stands for spin
													GlobalC::ucell.atoms[T0].index1_soc[0],
													GlobalC::ucell.atoms[T0].index2_soc[0],
													GlobalC::ucell.atoms[T0].nproj_soc
													);


											this->LM->DHloc_fixedR_x[nnr] += nlm[0];
											this->LM->DHloc_fixedR_y[nnr] += nlm[1];
											this->LM->DHloc_fixedR_z[nnr] += nlm[2];
										}
									}//!calc_deri
								}// distance
							} // ad0
							++nnr;
						}// k
					} // j 
				}// end is_adj
				//----------------------------------------------------------------------------------
			} // ad2
		} // I1
	} // T1


	if(!GlobalV::GAMMA_ONLY_LOCAL)
	{
//		std::cout << " nr="  << nnr << std::endl;
//		std::cout << " pv->nnr=" << pv->nnr << std::endl;
//		GlobalV::ofs_running << " nr="  << nnr << std::endl;
//		GlobalV::ofs_running << " pv->nnr=" << pv->nnr << std::endl;
		if( nnr!=pv->nnr)
		{
			ModuleBase::WARNING_QUIT("LCAO_gen_fixedH::build_Nonlocal_mu","nnr!=LNNR.nnr");
		}
	}

//	std::cout << " build_Nonlocal_mu done" << std::endl;

    ModuleBase::timer::tick ("LCAO_gen_fixedH","build_Nonlocal_mu");
	return;
}

void LCAO_gen_fixedH::build_Nonlocal_beta_new() //update by liuyu 2021-04-07
{
    ModuleBase::TITLE("LCAO_gen_fixedH","build_Nonlocal_beta_new");
    ModuleBase::timer::tick ("LCAO_gen_fixedH","build_Nonlocal_beta_new");

	const Parallel_Orbitals* pv = this->LM->ParaV;

#ifdef __MKL
    const int mkl_threads = mkl_get_max_threads();
    mkl_set_num_threads(1);
#endif

    const std::vector<std::vector<std::tuple<int, int, ModuleBase::Vector3<int>, ModuleBase::Vector3<double>>>> adjs_all = GlobalC::GridD.get_adjs(GlobalC::ucell);

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
                nlm_tot_thread.resize(adjs_all[iat].size());
            #else 
                std::vector<std::unordered_map<int,std::vector<double>>> nlm_tot;
                nlm_tot.resize(GlobalC::GridD.getAdjacentNum()+1);
            #endif 

            const ModuleBase::Vector3<double> tau0 = atom0->tau[I0];
            const double Rcut_Beta = GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();

            //outermost loop : all adjacent atoms
            for(int ad_count=0; ad_count < adjs_all[iat].size(); ad_count++)
            {
                const int T1 = std::get<0>(adjs_all[iat][ad_count]);
                const int I1 = std::get<1>(adjs_all[iat][ad_count]);
                const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
                const double Rcut_AO1 = GlobalC::ORB.Phi[T1].getRcut();
                const ModuleBase::Vector3<double> tau1 = std::get<3>(adjs_all[iat][ad_count]);
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
                    const int iw1_local = pv->trace_loc_row[iw1_all];
                    const int iw2_local = pv->trace_loc_col[iw1_all];

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
            for(int ad1_count=0; ad1_count < adjs_all[iat].size(); ad1_count++)
            {
                const int T1 = std::get<0>(adjs_all[iat][ad1_count]);
                const int I1 = std::get<1>(adjs_all[iat][ad1_count]);
                const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
                const ModuleBase::Vector3<double> tau1 = std::get<3>(adjs_all[iat][ad1_count]);
                const Atom* atom1 = &GlobalC::ucell.atoms[T1];
                const int nw1_tot = atom1->nw*GlobalV::NPOL;
                const double Rcut_AO1 = GlobalC::ORB.Phi[T1].getRcut();

                for (int ad2_count=0; ad2_count < adjs_all[iat].size(); ad2_count++)
                {
                    const int T2 = std::get<0>(adjs_all[iat][ad2_count]);
                    const int I2 = std::get<1>(adjs_all[iat][ad2_count]);
                    const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);
                    const ModuleBase::Vector3<double> tau2 = std::get<3>(adjs_all[iat][ad2_count]);
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
                        const int iw1_local = pv->trace_loc_row[iw1_all];
                        if(iw1_local < 0) continue;
                        const int iw1_0 = iw1/GlobalV::NPOL;
                        for(int iw2=0; iw2<nw2_tot; ++iw2)
                        {
                            const int iw2_all = start2 + iw2;
                            const int iw2_local = pv->trace_loc_col[iw2_all];
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
                                        nlm_thread += nlm1[ib]*nlm2[ib]*GlobalC::ucell.atoms[T0].dion(nb,nb);
                                    #else 
                                        nlm += nlm1[ib]*nlm2[ib]*GlobalC::ucell.atoms[T0].dion(nb,nb);
                                    #endif
                                    ib+=1;
                                }
                            }
                            assert(ib==nlm1.size());

                            const int ir = pv->trace_loc_row[ iw1_all ];
                            const int ic = pv->trace_loc_col[ iw2_all ];
                            long index=0;
                            if(GlobalV::KS_SOLVER=="genelpa" || GlobalV::KS_SOLVER=="scalapack_gvx")
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
                                this->LM->set_HSgamma(iw1_all,iw2_all,nlm,'N');
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
/*
    for (int T0 = 0; T0 < GlobalC::ucell.ntype; T0++)
    {
		Atom* atom0 = &GlobalC::ucell.atoms[T0]; 
        for (int I0 =0; I0< atom0->na; I0++)
        {
			//=======================================================
			//Step1:	
			//saves <beta|psi>, where beta runs over L0,M0 on atom I0
			//and psi runs over atomic basis sets on the current core
			//=======================================================
			std::vector<std::unordered_map<int,std::vector<double>>> nlm_tot;

            //GlobalC::GridD.Find_atom( atom0->tau[I0] );
			const ModuleBase::Vector3<double> tau0 = atom0->tau[I0];
            GlobalC::GridD.Find_atom(GlobalC::ucell, atom0->tau[I0] ,T0, I0);

			//outermost loop : all adjacent atoms
			nlm_tot.resize(GlobalC::GridD.getAdjacentNum()+1);
			const double Rcut_Beta = GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();

            for (int ad=0; ad<GlobalC::GridD.getAdjacentNum()+1 ; ++ad)
            {
                const int T1 = GlobalC::GridD.getType(ad);
                const int I1 = GlobalC::GridD.getNatom(ad);
                const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
				const double Rcut_AO1 = GlobalC::ORB.Phi[T1].getRcut();

                const ModuleBase::Vector3<double> tau1 = GlobalC::GridD.getAdjacentTau(ad);
				const Atom* atom1 = &GlobalC::ucell.atoms[T1];
				const int nw1_tot = atom1->nw*GlobalV::NPOL;

				//middle loop : atomic basis on current processor (either row or column)
				nlm_tot[ad].clear();

				const double dist1 = (tau1-tau0).norm() * GlobalC::ucell.lat0;
				if (dist1 > Rcut_Beta + Rcut_AO1)
				{
					continue;
				}

				for (int iw1=0; iw1<nw1_tot; ++iw1)
				{
					const int iw1_all = start1 + iw1;
					const int iw1_local = pv->trace_loc_row[iw1_all];
					const int iw2_local = pv->trace_loc_col[iw1_all];
					if(iw1_local < 0 && iw2_local < 0)continue;
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

					nlm_tot[ad].insert({iw1_all,nlm[0]});
				}//end iw
			}//end ad

			//=======================================================
			//Step2:	
			//calculate sum_(L0,M0) beta<psi_i|beta><beta|psi_j>
			//and accumulate the value to Hloc_fixed(i,j)
			//=======================================================

			for (int ad1=0; ad1<GlobalC::GridD.getAdjacentNum()+1 ; ++ad1)
            {
                const int T1 = GlobalC::GridD.getType(ad1);
                const int I1 = GlobalC::GridD.getNatom(ad1);
                const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
                const ModuleBase::Vector3<double> tau1 = GlobalC::GridD.getAdjacentTau(ad1);
				const Atom* atom1 = &GlobalC::ucell.atoms[T1];
				const int nw1_tot = atom1->nw*GlobalV::NPOL;
				const double Rcut_AO1 = GlobalC::ORB.Phi[T1].getRcut();

				for (int ad2=0; ad2 < GlobalC::GridD.getAdjacentNum()+1 ; ad2++)
				{
					const int T2 = GlobalC::GridD.getType(ad2);
					const int I2 = GlobalC::GridD.getNatom(ad2);
					const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);
					const ModuleBase::Vector3<double> tau2 = GlobalC::GridD.getAdjacentTau(ad2);
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

					for (int iw1=0; iw1<nw1_tot; ++iw1)
					{
						const int iw1_all = start1 + iw1;
						const int iw1_local = pv->trace_loc_row[iw1_all];
						if(iw1_local < 0)continue;
						const int iw1_0 = iw1/GlobalV::NPOL;
						for (int iw2=0; iw2<nw2_tot; ++iw2)
						{
							const int iw2_all = start2 + iw2;
							const int iw2_local = pv->trace_loc_col[iw2_all];
							if(iw2_local < 0)continue;
							const int iw2_0 = iw2/GlobalV::NPOL;

							std::vector<double> nlm1 = nlm_tot[ad1][iw1_all];
							std::vector<double> nlm2 = nlm_tot[ad2][iw2_all];

							assert(nlm1.size()==nlm2.size());

							double nlm=0.0;
							const int nproj = GlobalC::ucell.infoNL.nproj[T0];
							int ib = 0;
							for (int nb = 0; nb < nproj; nb++)
							{
								const int L0 = GlobalC::ucell.infoNL.Beta[T0].Proj[nb].getL();
								for(int m=0;m<2*L0+1;m++)
								{
									nlm += nlm1[ib]*nlm2[ib]*GlobalC::ucell.atoms[T0].dion(nb,nb);
									ib+=1;
								}
							}
							assert(ib==nlm1.size());
							this->LM->set_HSgamma(iw1_all,iw2_all,nlm,'N');
						}//iw2
					}//iw1
				}//ad2
			}//ad1

		}//end I0
	}//end T0
*/
	
    ModuleBase::timer::tick ("LCAO_gen_fixedH","build_Nonlocal_beta_new");
	return;
}

void LCAO_gen_fixedH::build_Nonlocal_beta(const bool& calc_deri) //update by liuyu 2021-04-07
{
    ModuleBase::TITLE("LCAO_gen_fixedH","build_Nonlocal_beta");
    ModuleBase::timer::tick ("LCAO_gen_fixedH","build_Nonlocal_beta");

	ModuleBase::matrix Rcut(GlobalC::ucell.ntype, GlobalC::ucell.ntype);
	for(int it1=0; it1<GlobalC::ucell.ntype; ++it1)
        for(int it2=0; it2<GlobalC::ucell.ntype; ++it2)
            Rcut(it1,it2) = GlobalC::ORB.Phi[it1].getRcut() + GlobalC::ORB.Phi[it2].getRcut();
	
    for (int T0 = 0; T0 < GlobalC::ucell.ntype; T0++)
    {
		Atom* atom0 = &GlobalC::ucell.atoms[T0]; 
        for (int I0 =0; I0< atom0->na; I0++)
        {
            //GlobalC::GridD.Find_atom( atom0->tau[I0] );
            GlobalC::GridD.Find_atom(GlobalC::ucell, atom0->tau[I0] ,T0, I0);

            //(2)
            //for each projector (T0, I0), one std::pair of ads are used
            for (int ad1=0; ad1<GlobalC::GridD.getAdjacentNum()+1 ; ++ad1)
            {
                const int T1 = GlobalC::GridD.getType(ad1);
                const int I1 = GlobalC::GridD.getNatom(ad1);
				//const int iat1 = GlobalC::ucell.itia2iat(T1, I1);
                const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
                const ModuleBase::Vector3<double> tau1 = GlobalC::GridD.getAdjacentTau(ad1);
				const Atom* atom1 = &GlobalC::ucell.atoms[T1];
				const int nw1_tot = atom1->nw*GlobalV::NPOL;

				// use to label < mu | H | nu(prime) >
				//int nnr = GlobalC::LNNR.nlocstart[iat];
            
				//(3)
				for (int ad2=0; ad2 < GlobalC::GridD.getAdjacentNum()+1 ; ad2++)
				{
					//if(ad2<ad && !calc_deri) continue; //add by liuyu 20210406
					const int T2 = GlobalC::GridD.getType(ad2);
					const int I2 = GlobalC::GridD.getNatom(ad2);
					const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);
					const ModuleBase::Vector3<double> tau2 = GlobalC::GridD.getAdjacentTau(ad2);
					const Atom* atom2 = &GlobalC::ucell.atoms[T2];
					const int nw2_tot = atom2->nw*GlobalV::NPOL;

					ModuleBase::Vector3<double> dtau = tau2 - tau1;
					double distance = dtau.norm() * GlobalC::ucell.lat0;
					double rcut = Rcut(T1,T2);
					//double rcut = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut();
//					if(distance < rcut)
//					{
						// ------------- enter the nnr increaing zone --------------
						//for (int iw1=0; iw1<atom1->nw*GlobalV::NPOL; ++iw1)
						for (int iw1=0; iw1<nw1_tot; ++iw1)
						{
							const int iw1_all = start1 + iw1;
							const int iw1_local = this->LM->ParaV->trace_loc_row[iw1_all];
							if(iw1_local < 0)continue;
							const int iw1_0 = iw1/GlobalV::NPOL;

							std::vector<std::vector<double>> nlm1;

							//if(!calc_deri && GlobalV::NSPIN!=4)
							if(false) //this method turns out to be slower!
							{
								GlobalC::UOT.snap_psibeta_half(
									GlobalC::ORB,
									GlobalC::ucell.infoNL,
									nlm1, tau1, T1,
									atom1->iw2l[ iw1_0 ], // L1
									atom1->iw2m[ iw1_0 ], // m1
									atom1->iw2n[ iw1_0 ], // N1
									GlobalC::ucell.atoms[T0].tau[I0], T0, 0);
							}

							// mohan fix bug 2010-12-20
							// atom2[T2] -> atom2.
							//for (int k=0; k<atom2->nw*GlobalV::NPOL; k++)
							for (int iw2=0; iw2<nw2_tot; ++iw2)
							{
								const int iw2_all = start2 + iw2;
								const int iw2_local = this->LM->ParaV->trace_loc_col[iw2_all];
								if(iw2_local < 0)continue;
								const int iw2_0 = iw2/GlobalV::NPOL;

								double nlm[3];
								nlm[0] = nlm[1] = nlm[2] = 0.0;

								std::vector<std::vector<double>> nlm2;
								if(!calc_deri)
								{
									//if(GlobalV::NSPIN!=4)
									GlobalC::UOT.snap_psibeta(
										GlobalC::ORB,
										GlobalC::ucell.infoNL,
										nlm, 0, tau1, T1,
										atom1->iw2l[ iw1_0 ], // L1
										atom1->iw2m[ iw1_0 ], // m1
										atom1->iw2n[ iw1_0 ], // N1
										tau2, T2,
										atom2->iw2l[ iw2_0 ], // L2
										atom2->iw2m[ iw2_0 ], // m2
										atom2->iw2n[ iw2_0 ], // n2
										GlobalC::ucell.atoms[T0].tau[I0], T0, GlobalC::ucell.atoms[T0].dion, GlobalV::NSPIN,
										GlobalC::ucell.atoms[T0].d_so,
										GlobalC::ucell.atoms[T0].non_zero_count_soc[0], // index stands for spin
										GlobalC::ucell.atoms[T0].index1_soc[0],
										GlobalC::ucell.atoms[T0].index2_soc[0],
										GlobalC::ucell.atoms[T0].nproj_soc
										);
									//if(GlobalV::GAMMA_ONLY_LOCAL)
									//{
									this->LM->set_HSgamma(iw1_all,iw2_all,nlm[0],'N');//N stands for nonlocal.
										//if(ad!=ad2) this->LM->set_HSgamma(iw2_all,iw1_all,nlm[0],'N'); //add by liuyu 20210406
									//}
								//	else
								//	{
								//		ModuleBase::WARNING_QUIT("LCAO_gen_fixedH::build_Nonlocal_beta","not consistent with k point algorithm.");
//										assert( nnr < pv->nnr );
//										this->LM->Hloc_fixedR[ nnr ] += nlm[0];
//										++nnr;
								//	}
								}
								else  // calculate force
								{
									GlobalC::UOT.snap_psibeta(
											GlobalC::ORB,
											GlobalC::ucell.infoNL,
											nlm, 1, tau1, T1,
											atom1->iw2l[ iw1_0 ], // L1
											atom1->iw2m[ iw1_0 ], // m1
											atom1->iw2n[ iw1_0 ], // N1
											tau2, T2,
											atom2->iw2l[ iw2_0 ], // L2
											atom2->iw2m[ iw2_0 ], // m2
											atom2->iw2n[ iw2_0 ], // n2
											GlobalC::ucell.atoms[T0].tau[I0], T0, GlobalC::ucell.atoms[T0].dion, GlobalV::NSPIN,
											GlobalC::ucell.atoms[T0].d_so,
											GlobalC::ucell.atoms[T0].non_zero_count_soc[0], // index stands for spin
											GlobalC::ucell.atoms[T0].index1_soc[0],
											GlobalC::ucell.atoms[T0].index2_soc[0],
											GlobalC::ucell.atoms[T0].nproj_soc
											);

									//if(GlobalV::GAMMA_ONLY_LOCAL)
									//{
										//add part of nonlocal ps derivatives to T matrix
										this->LM->set_force(iw1_all, iw2_all, nlm[0], nlm[1], nlm[2], 'N');
									//}
									//else
									//{
										//ModuleBase::WARNING_QUIT("LCAO_gen_fixedH::build_Nonlocal_beta","not consistent with k point algorithm.");
										//this->LM->DHloc_fixedR_x[ nnr ] += nlm[0];
										//this->LM->DHloc_fixedR_y[ nnr ] += nlm[1];
										//this->LM->DHloc_fixedR_z[ nnr ] += nlm[2];
										//++nnr;
									//}
								}
							}// end iw2
						}// end iw1
//					} // end distance
                }// end ad2
				// mohan add 2011-06-16

				/*if(!GlobalV::GAMMA_ONLY_LOCAL) // mohan fix bug 2011-06-26
				{
					if( iat < GlobalC::ucell.nat-1 )
					{
						if( nnr != GlobalC::LNNR.nlocstart[iat+1] )
						{
							std::cout << " nnr = " << nnr << std::endl;
							std::cout << " nlocstart[iat] = " << GlobalC::LNNR.nlocstart[iat] << std::endl;
							std::cout << " nlocstart[iat+1] = " << GlobalC::LNNR.nlocstart[iat+1] << std::endl;
							ModuleBase::WARNING_QUIT("build_Nonlocal_beta","nnr");
						}
					}
				}*/
            }// end ad1
        }// end I0
    }// end T0

    ModuleBase::timer::tick ("LCAO_gen_fixedH","build_Nonlocal_beta");
    return;
}


