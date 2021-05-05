#include "LCAO_gen_fixedH.h"
#include "../src_pw/global.h"
#include "../src_pw/wavefunc.h"
#include "LCAO_nnr.h"
#include "global_fp.h"

LCAO_gen_fixedH::LCAO_gen_fixedH()
{}

LCAO_gen_fixedH::~LCAO_gen_fixedH()
{}


void LCAO_gen_fixedH::calculate_NL_no(void)
{
    TITLE("LCAO_gen_fixedH","calculate_NL_no");

	// PLEASE rebuild the following two functions,
	// 'build_Nonlocal_beta' and  'build_Nonlocal_mu',
	// because the two functions are extremely time consuming
	// for small systems, especially for multiple-k points
	// mohan note 2021-03-23

	if(GAMMA_ONLY_LOCAL)
	{
	  	//for gamma only.
  		this->build_Nonlocal_beta(false);
	}
	else
	{
		// can work for gamma also,
		// only if search_radius is 
		// (Phi.rcutmax + Beta.rcutmax)*2.
		// check in sltk_atom_arrange.
    	this->build_Nonlocal_mu(false);
//		this->test_Nonlocal();
	}
    return;
}


void LCAO_gen_fixedH::calculate_T_no(void)
{
    TITLE("LCAO_gen_fixedH","calculate_T_no");
    this->build_ST_new('T', false);
    return;
}

void LCAO_gen_fixedH::calculate_S_no(void)
{
    TITLE("LCAO_gen_fixedH", "calculate_S_no");
    timer::tick("LCAO_gen_fixedH","calculate_S_no");
	this->build_ST_new('S', false);
    timer::tick("LCAO_gen_fixedH","calculate_S_no");
    return;
}


//liaochen modify interface 2010-3-22
void LCAO_gen_fixedH::build_ST_new(const char& dtype, const bool& calc_deri)
{
    TITLE("LCAO_gen_fixedH","build_ST_new");

    //array to store data
    double olm[3]={0.0,0.0,0.0};
	int nnr = 0; // used onlyh for k points.

    //\sum{T} e**{ikT} <\phi_{ia}|d\phi_{k\beta}(T)>
	Vector3<double> tau1, tau2, dtau;
	Vector3<double> dtau1, dtau2, tau0;
    for (int T1=0; T1<ucell.ntype; ++T1)
    {
		Atom* atom1 = &ucell.atoms[T1];
        for (int I1=0; I1<atom1->na; ++I1)
        {
			tau1 = atom1->tau[I1];
            //GridD.Find_atom(tau1);
            GridD.Find_atom(tau1, T1, I1);
            for (int ad = 0; ad < GridD.getAdjacentNum()+1; ++ad)
            {
                const int T2 = GridD.getType(ad);
				const int I2 = GridD.getNatom(ad);
				Atom* atom2 = &ucell.atoms[T2];
				tau2 = GridD.getAdjacentTau(ad);
				dtau = tau2 - tau1;
				double distance = dtau.norm() * ucell.lat0;
				double rcut = ORB.Phi[T1].getRcut() + ORB.Phi[T2].getRcut();
				if(distance < rcut)
				{
					int iw1_all = ucell.itiaiw2iwt( T1, I1, 0) ; //iw1_all = combined index (it, ia, iw)

					for(int jj=0; jj<atom1->nw*NPOL; ++jj)
					{
						const int jj0 = jj/NPOL;
						const int L1 = atom1->iw2l[jj0];
						const int N1 = atom1->iw2n[jj0];
						const int m1 = atom1->iw2m[jj0];

						int iw2_all = ucell.itiaiw2iwt( T2, I2, 0);//zhengdy-soc
						for(int kk=0; kk<atom2->nw*NPOL; ++kk)
						{
							const int kk0 = kk/NPOL;
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
							if ( !ParaO.in_this_processor(iw1_all,iw2_all) )
							{
								++iw2_all;
								continue;
							}

							olm[0] = olm[1] = olm[2] = 0.0;

							complex<double> olm1[4]={ZERO, ZERO, ZERO, ZERO};
							complex<double> *olm2 = &olm1[0];
							if(!calc_deri)
							{
								// PLEASE use UOT as an input parameter of this subroutine
								// mohan add 2021-03-30
								UOT.snap_psipsi( olm, 0, dtype, tau1, 
										T1, L1, m1, N1, GridD.getAdjacentTau(ad), 
										T2, L2, m2, N2,
										olm2//for soc
										);

								if(GAMMA_ONLY_LOCAL)
								{
									// mohan add 2010-06-29
									// set the value in Hloc and Sloc
									// according to trace_loc_row and trace_loc_col
									// the last paramete: 1 for Sloc, 2 for Hloc
									// and 3 for Hloc_fixed.
									LM.set_HSgamma(iw1_all, iw2_all, olm[0], dtype);
								}
								else // k point algorithm
								{
									// mohan add 2010-10
									// set the values in SlocR and Hloc_fixedR.
									// which is a 1D array.
									if(dtype=='S')
									{
										if(NSPIN!=4) LM.SlocR[nnr] = olm[0];
										else
										{//only has diagonal term here.
												int is = (jj-jj0*NPOL) + (kk-kk0*NPOL)*2;
											LM.SlocR_soc[nnr] = olm1[is];
										}
									}
									else if(dtype=='T')
									{
										if(NSPIN!=4) LM.Hloc_fixedR[nnr] = olm[0];// <phi|kin|d phi>
										else
										{//only has diagonal term here.
												int is = (jj-jj0*NPOL) + (kk-kk0*NPOL)*2;
											LM.Hloc_fixedR_soc[nnr] = olm1[is];
										}
									}
									++nnr;
								}
							}
							else // calculate the derivative
							{
								UOT.snap_psipsi( olm, 1, dtype, 
									tau1, T1, L1, m1, N1,
									GridD.getAdjacentTau(ad), T2, L2, m2, N2
									);

								if(GAMMA_ONLY_LOCAL)
								{
									LM.set_force (iw1_all, iw2_all,	olm[0], olm[1], olm[2], dtype);
									if(STRESS) LM.set_stress (iw1_all, iw2_all, olm[0], olm[1], olm[2], dtype, dtau);
								}
								else // k point algorithm
								{
									if(dtype=='S')
									{
										LM.DSloc_Rx[nnr] = olm[0];
										LM.DSloc_Ry[nnr] = olm[1];
										LM.DSloc_Rz[nnr] = olm[2];
										if(STRESS)
										{
											LM.DH_r[nnr*3] = dtau.x;
											LM.DH_r[nnr*3 + 1] = dtau.y;
											LM.DH_r[nnr*3 + 2] = dtau.z;
										}
									}
									else if(dtype=='T')
									{
										// notice the 'sign'
										LM.DHloc_fixedR_x[nnr] = olm[0];
										LM.DHloc_fixedR_y[nnr] = olm[1];
										LM.DHloc_fixedR_z[nnr] = olm[2];
										if(STRESS)
										{
											LM.stvnl11[nnr] = olm[0] * dtau.x;
											LM.stvnl12[nnr] = olm[0] * dtau.y;
											LM.stvnl13[nnr] = olm[0] * dtau.z;
											LM.stvnl22[nnr] = olm[1] * dtau.y;
											LM.stvnl23[nnr] = olm[1] * dtau.z;
											LM.stvnl33[nnr] = olm[2] * dtau.z;
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
				else if(distance>=rcut && (!GAMMA_ONLY_LOCAL))
				{
					int start1 = ucell.itiaiw2iwt( T1, I1, 0);
					int start2 = ucell.itiaiw2iwt( T2, I2, 0);

					bool is_adj = false;
					for (int ad0=0; ad0 < GridD.getAdjacentNum()+1; ++ad0)
					{
						const int T0 = GridD.getType(ad0);
						//const int I0 = GridD.getNatom(ad0);
						//const int iat0 = ucell.itia2iat(T0, I0);
						//const int start0 = ucell.itiaiw2iwt(T0, I0, 0);
						tau0 = GridD.getAdjacentTau(ad0);
						dtau1 = tau0 - tau1;
						double distance1 = dtau1.norm() * ucell.lat0;
						double rcut1 = ORB.Phi[T1].getRcut() + ORB.Beta[T0].get_rcut_max();
						dtau2 = tau0 - tau2;
						double distance2 = dtau2.norm() * ucell.lat0;
						double rcut2 = ORB.Phi[T2].getRcut() + ORB.Beta[T0].get_rcut_max();
						if( distance1 < rcut1 && distance2 < rcut2 )
						{
							is_adj = true;
							break;
						}
					}//ad0


					if( is_adj )
					{
						for(int jj=0; jj<atom1->nw * NPOL; ++jj)
						{
							const int mu = ParaO.trace_loc_row[start1+jj];
							if(mu<0)continue; 
							for(int kk=0; kk<atom2->nw * NPOL; ++kk)
							{
								const int nu = ParaO.trace_loc_col[start2+kk];
								if(nu<0)continue;
								++nnr;
							}//kk
						}//jj
					}
				}//distance
			}// ad
		}// I1
	}// T1

	if(!GAMMA_ONLY_LOCAL)
	{
		if(nnr != LNNR.nnr)
		{
			cout << " nnr=" << nnr << " LNNR.nnr=" << LNNR.nnr << endl;
			ofs_running << " nnr=" << nnr << " LNNR.nnr=" << LNNR.nnr << endl;
			WARNING_QUIT("LCAO_gen_fixedH::build_ST_new","nnr != LNNR.nnr");
		}
	}

    return;
}

void LCAO_gen_fixedH::test_Nonlocal()
{
	int nnr = 0;
	Vector3<double> tau1, tau2, dtau_12, tau0, dtau_10, dtau_20;
	double distance = 0.0;
	double rcut = 0.0;

//	double* vnltest = new double[ParaO.nloc];	
//	ZEROS(vnltest, ParaO.nloc);

	// psi1
	double sum = 0.0;
	int count = 0;
    for (int T1 = 0; T1 < ucell.ntype; T1++)
    {
		const Atom* atom1 = &ucell.atoms[T1];
        for (int I1 =0; I1< atom1->na; I1++)
        {
            //GridD.Find_atom( atom1->tau[I1] );
            GridD.Find_atom( atom1->tau[I1] ,T1, I1);
			//const int iat1 = ucell.itia2iat(T1, I1);
			const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
            tau1 = atom1->tau[I1];

			// psi2
            for (int ad2=0; ad2<GridD.getAdjacentNum()+1 ; ad2++)
			{
				const int T2 = GridD.getType(ad2);
				const Atom* atom2 = &ucell.atoms[T2];
                
				const int I2 = GridD.getNatom(ad2);
				//const int iat2 = ucell.itia2iat(T2, I2);
                const int start2 = ucell.itiaiw2iwt(T2, I2, 0);
                tau2 = GridD.getAdjacentTau(ad2);

				dtau_12 = tau2 - tau1;
				distance = dtau_12.norm() * ucell.lat0;
				rcut = ORB.Phi[T1].getRcut() + ORB.Phi[T2].getRcut();

				if(distance >= rcut)
				{
					for (int j=0; j<atom1->nw*NPOL; j++)
					{
						int j0 = j/NPOL;
						const int iw1_all = start1 + j;
						const int mu = ParaO.trace_loc_row[iw1_all];
						if(mu < 0)continue; 
						for (int k=0; k<atom2->nw*NPOL; k++)
						{
							int k0 = k/NPOL;
							const int iw2_all = start2 + k;
							const int nu = ParaO.trace_loc_col[iw2_all];						
							if(nu < 0)continue;
							for (int ad0=0; ad0 < GridD.getAdjacentNum()+1 ; ad0++)
							{
								const int T0 = GridD.getType(ad0);
								if( ORB.nproj[T0] == 0) continue; 
								//const int I0 = GridD.getNatom(ad0);
								//const int start0 = ucell.itiaiw2iwt(T0, I0, 0);
								tau0 = GridD.getAdjacentTau(ad0);

								dtau_10 = tau0 - tau1;
								dtau_20 = tau0 - tau2;
								double distance1 = dtau_10.norm() * ucell.lat0;
								double distance2 = dtau_20.norm() * ucell.lat0;

								double rcut_10 = ORB.Phi[T1].getRcut() + ORB.Phi[T0].getRcut();
								double rcut_20 = ORB.Phi[T2].getRcut() + ORB.Phi[T0].getRcut();

								if(distance1 < rcut_10 && distance2 < rcut_20)
								{
									++count;
									double nlm[3]={0,0,0};
									UOT.snap_psibeta(
											nlm, 0, tau1, T1,
											atom1->iw2l[ j0 ], // L1
											atom1->iw2m[ j0 ], // m1
											atom1->iw2n[ j0 ], // N1
											tau2, T2,
											atom2->iw2l[ k0 ], // L2
											atom2->iw2m[ k0 ], // m2
											atom2->iw2n[ k0 ], // n2
											tau0, T0
											);
									
									//vnltest[ mu * ParaO.ncol + nu ] += nlm[0];
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

	ofs_running << " not included Vnl pair = " << count << endl;
	ofs_running << " sum for the correction = " << sum << endl;

	/*
	cout << " correction for Vnl matrix " << endl;
	for(int i=0; i<NLOCAL; i++)
	{
		for(int j=0; j<NLOCAL; j++)
		{
			double a = vnltest[i*ParaO.ncol+j];
			if( abs(a) > 1.0e-6 )
			{
				cout << setw(15) << vnltest[i*ParaO.ncol+j];
			}
			else
			{
				cout << setw(15) << "0";
			}
		}
		cout << endl;
	}
	delete[] vnltest;
	*/
		
	return;
}


#include "record_adj.h" //mohan add 2012-07-06
void LCAO_gen_fixedH::build_Nonlocal_mu(const bool &calc_deri)
{
    TITLE("LCAO_gen_fixedH","build_Nonlocal_mu");
    timer::tick ("LCAO_gen_fixedH","build_Nonlocal_mu",'G');

	// < phi1 | beta > < beta | phi2 >
	// phi1 is within the unitcell.
	// while beta is in the supercell.
	// while phi2 is in the supercell.

	int nnr = 0;
	Vector3<double> tau1, tau2, dtau;
	Vector3<double> dtau1, dtau2, tau0;
	double distance = 0.0;
	double distance1, distance2;
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
            //GridD.Find_atom( atom1->tau[I1] );
            GridD.Find_atom( atom1->tau[I1] ,T1, I1);
			//const int iat1 = ucell.itia2iat(T1, I1);
			const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
            tau1 = atom1->tau[I1];

			// psi2
            for (int ad2=0; ad2<GridD.getAdjacentNum()+1; ++ad2)
			{
				const int T2 = GridD.getType(ad2);
				const Atom* atom2 = &ucell.atoms[T2];
                
				const int I2 = GridD.getNatom(ad2);
				//const int iat2 = ucell.itia2iat(T2, I2);
                const int start2 = ucell.itiaiw2iwt(T2, I2, 0);
                tau2 = GridD.getAdjacentTau(ad2);

				bool is_adj = false;
					
				dtau = tau2 - tau1;
				distance = dtau.norm() * ucell.lat0;
				// this rcut is in order to make nnr consistent 
				// with other matrix.
				rcut = ORB.Phi[T1].getRcut() + ORB.Phi[T2].getRcut();
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

                        double distance1 = dtau1.norm() * ucell.lat0;
                        double distance2 = dtau2.norm() * ucell.lat0;

                        rcut1 = ORB.Phi[T1].getRcut() + ORB.Beta[T0].get_rcut_max();
                        rcut2 = ORB.Phi[T2].getRcut() + ORB.Beta[T0].get_rcut_max();

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
					for (int j=0; j<atom1->nw*NPOL; j++)
					{
						const int j0 = j/NPOL;//added by zhengdy-soc
						const int iw1_all = start1 + j;
						const int mu = ParaO.trace_loc_row[iw1_all];
						if(mu < 0)continue; 

						// fix a serious bug: atom2[T2] -> atom2
						// mohan 2010-12-20
						for (int k=0; k<atom2->nw*NPOL; k++)
						{
							const int k0 = k/NPOL;
							const int iw2_all = start2 + k;
							const int nu = ParaO.trace_loc_col[iw2_all];						
							if(nu < 0)continue;


							//(3) run over all projectors in nonlocal pseudopotential.
							for (int ad0=0; ad0 < GridD.getAdjacentNum()+1 ; ++ad0)
							{
								const int T0 = GridD.getType(ad0);

								// mohan add 2010-12-19
								if( ORB.nproj[T0] == 0) continue; 

								//const int I0 = GridD.getNatom(ad0);
								//const int start0 = ucell.itiaiw2iwt(T0, I0, 0);
								tau0 = GridD.getAdjacentTau(ad0);

								dtau1 = tau0 - tau1;
								dtau2 = tau0 - tau2;
								distance1 = dtau1.norm() * ucell.lat0;
								distance2 = dtau2.norm() * ucell.lat0;

								// seems a bug here!! mohan 2011-06-17
								rcut1 = ORB.Phi[T1].getRcut() + ORB.Beta[T0].get_rcut_max();
								rcut2 = ORB.Phi[T2].getRcut() + ORB.Beta[T0].get_rcut_max();

								if(distance1 < rcut1 && distance2 < rcut2)
								{
									//const Atom* atom0 = &ucell.atoms[T0];
									double nlm[3]={0,0,0};
									complex<double> nlm1[4]={0,0,0,0};//modified by zhengdy-soc
									complex<double> *nlm2 = NULL;
									if(NSPIN==4) nlm2 = &nlm1[0];
									if(!calc_deri)
									{
										int is0 = (j-j0*NPOL) + (k-k0*NPOL)*2;
										UOT.snap_psibeta(
												nlm, 0, tau1, T1,
												atom1->iw2l[ j0 ], // L1
												atom1->iw2m[ j0 ], // m1
												atom1->iw2n[ j0 ], // N1
												tau2, T2,
												atom2->iw2l[ k0 ], // L2
												atom2->iw2m[ k0 ], // m2
												atom2->iw2n[ k0 ], // n2
												tau0, T0,
												nlm2, is0 //for soc
												);


										if(GAMMA_ONLY_LOCAL)
										{
											// mohan add 2010-12-20
											if( nlm[0]!=0.0 )
											{
												// ofs_running << setw(10) << iw1_all << setw(10) 
												// << iw2_all << setw(20) << nlm[0] << endl; 
												LM.set_HSgamma(iw1_all,iw2_all,nlm[0],'N');//N stands for nonlocal.
											}
										}
										else
										{
											if(NSPIN!=4) LM.Hloc_fixedR[nnr] += nlm[0];
											else
											{
												int is = (j-j0*NPOL) + (k-k0*NPOL)*2;
												LM.Hloc_fixedR_soc[nnr] += nlm1[is];
											}
										}
									}// calc_deri
									else // calculate the derivative
									{
										if(GAMMA_ONLY_LOCAL)
										{
											UOT.snap_psibeta(
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
													tau0, T0
													);

											// sum all projectors for one atom.
											LM.set_force (iw1_all, iw2_all,	nlm[0], nlm[1], nlm[2], 'N');
										}
										else
										{
											// mohan change the order on 2011-06-17
											// origin: < psi1 | beta > < beta | dpsi2/dtau >
											//now: < psi1/dtau | beta > < beta | psi2 >
											UOT.snap_psibeta(
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
													tau0, T0
													);


											LM.DHloc_fixedR_x[nnr] += nlm[0];
											LM.DHloc_fixedR_y[nnr] += nlm[1];
											LM.DHloc_fixedR_z[nnr] += nlm[2];
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


	if(!GAMMA_ONLY_LOCAL)
	{
//		cout << " nr="  << nnr << endl;
//		cout << " LNNR.nnr=" << LNNR.nnr << endl;
//		ofs_running << " nr="  << nnr << endl;
//		ofs_running << " LNNR.nnr=" << LNNR.nnr << endl;
		if( nnr!=LNNR.nnr)
		{
			WARNING_QUIT("LCAO_gen_fixedH::build_Nonlocal_mu","nnr!=LNNR.nnr");
		}
	}

//	cout << " build_Nonlocal_mu done" << endl;

    timer::tick ("LCAO_gen_fixedH","build_Nonlocal_mu",'G');
	return;
}


void LCAO_gen_fixedH::build_Nonlocal_beta(const bool& calc_deri)
{
    TITLE("LCAO_gen_fixedH","build_Nonlocal_beta");
    timer::tick ("LCAO_gen_fixedH","build_Nonlocal_beta",'G');

    for (int T0 = 0; T0 < ucell.ntype; T0++)
    {
		Atom* atom0 = &ucell.atoms[T0]; 
        for (int I0 =0; I0< atom0->na; I0++)
        {
            //GridD.Find_atom( atom0->tau[I0] );
            GridD.Find_atom( atom0->tau[I0] ,T0, I0);

            //(2)
            //for each projector (T0, I0), one pair of ads are used
            for (int ad=0; ad<GridD.getAdjacentNum()+1 ; ad++)
            {
                const int T1 = GridD.getType(ad);
                const int I1 = GridD.getNatom(ad);
				const int iat = ucell.itia2iat(T1, I1);
                const int start = ucell.itiaiw2iwt(T1, I1, 0);
                const Vector3<double> tau1 = GridD.getAdjacentTau(ad);
				const Atom* atom1 = &ucell.atoms[T1];

				// use to label < mu | H | nu(prime) >
				int nnr = LNNR.nlocstart[iat];
            
				//(3)
				for (int ad2=0; ad2 < GridD.getAdjacentNum()+1 ; ad2++)
				{
					const int T2 = GridD.getType(ad2);
					const int I2 = GridD.getNatom(ad2);
					const int start2 = ucell.itiaiw2iwt(T2, I2, 0);
					const Vector3<double> tau2 = GridD.getAdjacentTau(ad2);
					const Atom* atom2 = &ucell.atoms[T2];

					Vector3<double> dtau = tau2 - tau1;
					double distance = dtau.norm() * ucell.lat0;
					double rcut = ORB.Phi[T1].getRcut() + ORB.Phi[T2].getRcut();
					if(distance < rcut)
					{
						// ------------- enter the nnr increaing zone --------------
						for (int j=0; j<atom1->nw*NPOL; j++)
						{
							const int j0 = j/NPOL;
							const int iw1_all = start + j;
							const int mu = ParaO.trace_loc_row[iw1_all];
							if(mu < 0)continue; 

							// mohan fix bug 2010-12-20
							// atom2[T2] -> atom2.
							for (int k=0; k<atom2->nw*NPOL; k++)
							{
								const int k0 = k/NPOL;
								const int iw2_all = start2 + k;
								const int nu = ParaO.trace_loc_col[iw2_all];
								if(nu < 0)continue;

								double nlm[3];
								nlm[0] = nlm[1] = nlm[2] = 0.0;

								if(!calc_deri)
								{
									UOT.snap_psibeta(
											nlm, 0, tau1, T1,
											atom1->iw2l[ j0 ], // L1
											atom1->iw2m[ j0 ], // m1
											atom1->iw2n[ j0 ], // N1
											tau2, T2,
											atom2->iw2l[ k0 ], // L2
											atom2->iw2m[ k0 ], // m2
											atom2->iw2n[ k0 ], // n2
											ucell.atoms[T0].tau[I0], T0
											);

									if(GAMMA_ONLY_LOCAL)
									{
										LM.set_HSgamma(iw1_all,iw2_all,nlm[0],'N');//N stands for nonlocal.
									}
									else
									{
										WARNING_QUIT("LCAO_gen_fixedH::build_Nonlocal_beta","not consistent with k point algorithm.");
//										assert( nnr < LNNR.nnr );
//										LM.Hloc_fixedR[ nnr ] += nlm[0];
//										++nnr;
									}
								}
								else  // calculate force
								{
									UOT.snap_psibeta(
											nlm, 1, tau1, T1,
											atom1->iw2l[ j0 ], // L1
											atom1->iw2m[ j0 ], // m1
											atom1->iw2n[ j0 ], // N1
											tau2, T2,
											atom2->iw2l[ k0 ], // L2
											atom2->iw2m[ k0 ], // m2
											atom2->iw2n[ k0 ], // n2
											ucell.atoms[T0].tau[I0], T0
											);

									if(GAMMA_ONLY_LOCAL)
									{
										//add part of nonlocal ps derivatives to T matrix
										LM.set_force(iw1_all, iw2_all, nlm[0], nlm[1], nlm[2], 'N');
									}
									else
									{
										WARNING_QUIT("LCAO_gen_fixedH::build_Nonlocal_beta","not consistent with k point algorithm.");
										//LM.DHloc_fixedR_x[ nnr ] += nlm[0];
										//LM.DHloc_fixedR_y[ nnr ] += nlm[1];
										//LM.DHloc_fixedR_z[ nnr ] += nlm[2];
										++nnr;
									}
								}
							}// end k
						}// j 
					} // end distance
                }// ad2
				// mohan add 2011-06-16

				if(!GAMMA_ONLY_LOCAL) // mohan fix bug 2011-06-26
				{
					if( iat < ucell.nat-1 )
					{
						if( nnr != LNNR.nlocstart[iat+1] )
						{
							cout << " nnr = " << nnr << endl;
							cout << " nlocstart[iat] = " << LNNR.nlocstart[iat] << endl;
							cout << " nlocstart[iat+1] = " << LNNR.nlocstart[iat+1] << endl;
							WARNING_QUIT("build_Nonlocal_beta","nnr");
						}
					}
				}
            }// ad
        }// end I0
    }// end T0

    timer::tick ("LCAO_gen_fixedH","build_Nonlocal_beta",'G');
    return;
}


