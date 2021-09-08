#include "LCAO_nnr.h"
#include "../src_pw/global.h"
#include "record_adj.h" //mohan add 2012-07-06

//----------------------------
// define a global class obj.
//----------------------------
namespace GlobalC
{
LCAO_nnr LNNR;
}

LCAO_nnr::LCAO_nnr()
{
	nnr = 1;
	nlocdimg = new int[1];	
	nlocstartg = new int[1];
	nlocdim = new int[1];	
	nlocstart = new int[1];
	
	// number of adjacent atoms for each atom.
	nad = new int[1];
	allocate_find_R2 = false;
}

LCAO_nnr::~LCAO_nnr()
{
	delete[] nlocdimg;
	delete[] nlocstartg;
	delete[] nad;
	delete[] nlocdim;
	delete[] nlocstart;

	if(allocate_find_R2)
	{
		for(int iat=0; iat<GlobalC::ucell.nat; iat++)
		{
			delete[] find_R2[iat];
			delete[] find_R2st[iat];
		}
		delete[] find_R2;
		delete[] find_R2st;
	}
}

// be called in LOOP_ions.cpp
void LCAO_nnr::cal_nnr(void)
{
	ModuleBase::TITLE("LCAO_nnr","cal_nnr");

	delete[] nlocdim;
	delete[] nlocstart;
	nlocdim = new int[GlobalC::ucell.nat];	
	nlocstart = new int[GlobalC::ucell.nat];
	ModuleBase::GlobalFunc::ZEROS(nlocdim, GlobalC::ucell.nat);
	ModuleBase::GlobalFunc::ZEROS(nlocstart, GlobalC::ucell.nat);

	this->nnr = 0;
	int start = 0;
	//int ind1 = 0;
	int iat = 0;

	// (1) find the adjacent atoms of atom[T1,I1];
	ModuleBase::Vector3<double> tau1;
	ModuleBase::Vector3<double> tau2;
	ModuleBase::Vector3<double> dtau;
	ModuleBase::Vector3<double> tau0;
	ModuleBase::Vector3<double> dtau1;
	ModuleBase::Vector3<double> dtau2;

	for (int T1 = 0; T1 < GlobalC::ucell.ntype; T1++)
	{
		for (int I1 = 0; I1 < GlobalC::ucell.atoms[T1].na; I1++)
		{
			tau1 = GlobalC::ucell.atoms[T1].tau[I1];
			//GlobalC::GridD.Find_atom( tau1 );
			GlobalC::GridD.Find_atom(GlobalC::ucell,  tau1 ,T1, I1);
			const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
			this->nlocstart[iat] = nnr;
			int nw1 = GlobalC::ucell.atoms[T1].nw * GlobalV::NPOL;

			// (2) search among all adjacent atoms.
			for (int ad = 0; ad < GlobalC::GridD.getAdjacentNum()+1; ad++)
			{
				const int T2 = GlobalC::GridD.getType(ad);
				const int I2 = GlobalC::GridD.getNatom(ad);
				//const int iat2 = GlobalC::ucell.itia2iat(T2, I2);
				const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);
				int nw2 = GlobalC::ucell.atoms[T2].nw * GlobalV::NPOL;

				tau2 = GlobalC::GridD.getAdjacentTau(ad);

				dtau = tau2 - tau1;
				double distance = dtau.norm() * GlobalC::ucell.lat0;
				double rcut = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut();
				
				if(distance < rcut)
				{
					//--------------------------------------------------
					// calculate how many matrix elements are in 
					// this processor.
					for(int ii=0; ii<nw1; ii++)
					{
						// the index of orbitals in this processor
						// according to HPSEPS's division method.
						const int iw1_all = start1 + ii;
						const int mu = GlobalC::ParaO.trace_loc_row[iw1_all];
						if(mu<0)continue;

						for(int jj=0; jj<nw2; jj++)
						{
							const int iw2_all = start2 + jj;
							const int nu = GlobalC::ParaO.trace_loc_col[iw2_all];
							if(nu<0)continue;

							// orbital numbers for this atom (iat),
							// seperated by atoms in different cells.
							this->nlocdim[iat]++; 

							++nnr;
						}// end jj
					}// end ii
				}//end distance
				// there is another possibility that i and j are adjacent atoms.
				// which is that <i|beta> are adjacents while <beta|j> are also
				// adjacents, these considerations are only considered in k-point
				// algorithm, 
				// mohan fix bug 2012-07-03
				else if(distance >= rcut)
				{
					for (int ad0 = 0; ad0 < GlobalC::GridD.getAdjacentNum()+1; ++ad0)
					{
						const int T0 = GlobalC::GridD.getType(ad0);
						const int I0 = GlobalC::GridD.getNatom(ad0);
						//const int iat0 = GlobalC::ucell.itia2iat(T0, I0);
						//const int start0 = GlobalC::ucell.itiaiw2iwt(T0, I0, 0);
					
						tau0 = GlobalC::GridD.getAdjacentTau(ad0);
						dtau1 = tau0 - tau1; 
						double distance1 = dtau1.norm() * GlobalC::ucell.lat0;
						double rcut1 = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();

						dtau2 = tau0 - tau2;
						double distance2 = dtau2.norm() * GlobalC::ucell.lat0;
						double rcut2 = GlobalC::ORB.Phi[T2].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();

						if( distance1 < rcut1 && distance2 < rcut2 )
						{
							for(int ii=0; ii<nw1; ++ii)
							{
								const int iw1_all = start1 + ii;
								const int mu = GlobalC::ParaO.trace_loc_row[iw1_all];
								if(mu<0)continue;

								for(int jj=0; jj<nw2; ++jj)
								{
									const int iw2_all = start2 + jj;
									const int nu = GlobalC::ParaO.trace_loc_col[iw2_all];
									if(nu<0)continue;

									// orbital numbers for this atom (iat),
									// seperated by atoms in different cells.
									this->nlocdim[iat]++;

									++nnr;
								}
							}
							break;	
						} // dis1, dis2
					}//ad0
				}
			}// end ad

			//start position of atom[T1,I1]
			start += nw1;
			++iat;
		}// end I1
	} // end T1

	//xiaohui add 'GlobalV::OUT_LEVEL' line, 2015-09-16
	if(GlobalV::OUT_LEVEL != "m") ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nnr",nnr);
//	for(int iat=0; iat<GlobalC::ucell.nat; iat++)
//	{
//		std::cout << " nlocdim[" << iat << "]=" << nlocdim[iat];
//		std::cout << " nlocstart[" << iat << "]=" << nlocstart[iat] << std::endl;
//	}

	return;
}

// This is for cell R dependent part. 
void LCAO_nnr::cal_nnrg(const Grid_Technique &GT)
{
	ModuleBase::TITLE("LCAO_nnr","cal_nnrg");

	this->cal_max_box_index();

	this->nnrg = 0;

	delete[] nlocdimg;
	delete[] nlocstartg;
	delete[] nad;
	
	this->nad = new int[GlobalC::ucell.nat];
	this->nlocdimg = new int[GlobalC::ucell.nat];	
	this->nlocstartg = new int[GlobalC::ucell.nat];
	
	ModuleBase::GlobalFunc::ZEROS(nad, GlobalC::ucell.nat);
	ModuleBase::GlobalFunc::ZEROS(nlocdimg, GlobalC::ucell.nat);
	ModuleBase::GlobalFunc::ZEROS(nlocstartg, GlobalC::ucell.nat);


	ModuleBase::Vector3<double> tau1, tau2, dtau;
	ModuleBase::Vector3<double> dtau1, dtau2, tau0;
	for (int T1 = 0; T1 < GlobalC::ucell.ntype; ++T1)
	{
		Atom* atom1 = &GlobalC::ucell.atoms[T1];
		for (int I1 = 0; I1 < atom1->na; ++I1)
		{
			tau1 = atom1->tau[I1];

			GlobalC::GridD.Find_atom(GlobalC::ucell, tau1, T1, I1);

			const int iat = GlobalC::ucell.itia2iat(T1,I1);

			// for grid integration (on FFT box),
			// we only need to consider <phi_i | phi_j>,
			// which is different from non-local term,
			// which we need to consdier <phi_i|beta_k><beta_k|phi_j>

			// whether this atom is in this processor.
			if(GT.in_this_processor[iat])
			{
				// starting index of adjacents.
				this->nlocstartg[iat] = this->nnrg;

				// number of adjacent atoms for atom 'iat'
				this->nad[iat] = 0;

				int count = 0;
				for (int ad = 0; ad < GlobalC::GridD.getAdjacentNum()+1; ++ad)
				{
					const int T2 = GlobalC::GridD.getType(ad);
					const int I2 = GlobalC::GridD.getNatom(ad);
					const int iat2 = GlobalC::ucell.itia2iat(T2, I2);
					Atom* atom2 = &GlobalC::ucell.atoms[T2]; 

					// if the adjacent atom is in this processor.
					if(GT.in_this_processor[iat2])
					{
						tau2 = GlobalC::GridD.getAdjacentTau(ad);
						dtau = GlobalC::GridD.getAdjacentTau(ad) - tau1;
						double distance = dtau.norm() * GlobalC::ucell.lat0;
						double rcut = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut();


						//if(distance < rcut)
						// mohan reset this 2013-07-02 in Princeton
						// we should make absolutely sure that the distance is smaller than GlobalC::ORB.Phi[it].getRcut
						// this should be consistant with LCAO_nnr::cal_nnrg function 
						// typical example : 7 Bohr cutoff Si orbital in 14 Bohr length of cell.
						// distance = 7.0000000000000000
						// GlobalC::ORB.Phi[it].getRcut = 7.0000000000000008
						if(distance < rcut - 1.0e-15)
						{
							const int nelement = atom1->nw * atom2->nw;//modified by zhengdy-soc, no need to double
							this->nnrg += nelement;
							this->nlocdimg[iat] += nelement; 
							this->nad[iat]++;
							++count;

							/*
							   ofs << std::setw(10) << iat << std::setw(10) << iat2
							   << std::setw(10) << GlobalC::GridD.getBox(ad).x 
							   << std::setw(10) << GlobalC::GridD.getBox(ad).y 
							   << std::setw(10) << GlobalC::GridD.getBox(ad).z 
							   << std::setw(20) << distance << std::endl;
							 */
						}
						// there is another possibility that i and j are adjacent atoms.
						// which is that <i|beta> are adjacents while <beta|j> are also
						// adjacents, these considerations are only considered in k-point
						// algorithm, 
						// mohan fix bug 2012-07-03
						/*
						else if(distance >= rcut)
						{
							for (int ad0 = 0; ad0 < GlobalC::GridD.getAdjacentNum()+1; ++ad0)
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

								if( distance1 < rcut1 && distance2 < rcut2 )
								{
									const int nelement = atom1->nw * atom2->nw;
									this->nnrg += nelement;
									this->nlocdimg[iat] += nelement; 
									this->nad[iat]++;
									++count;
									break;	
								} // dis1, dis2
							}//ad0
						}//distance
						*/
					}// end iat2
				}// end ad
//				GlobalV::ofs_running << " iat=" << iat << " nlocstartg=" << nlocstartg[iat] << " nad=" << nad[iat] << std::endl;
			}// end iat
		}// end I1
	}// end T1

	if(GlobalV::OUT_LEVEL != "m") ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nnrg",this->nnrg);

	//--------------------------------------------------
	// search again, to order each (iat2, b1, b2, b3)
	// find_R2 is used to target DM_R.
	// because DM_R is allocated with nnrg.
	// So once we had dR = R2 - R1 and iat2,
	// we need to find out the corresponding positions
	// in DM_R
	//--------------------------------------------------
	if(allocate_find_R2)
	{
		for(int iat=0; iat<GlobalC::ucell.nat; iat++)
		{
			delete[] find_R2[iat];
		}
		delete[] find_R2;

		for(int iat=0; iat<GlobalC::ucell.nat; iat++)
		{
			delete[] find_R2st[iat];
		}
		delete[] find_R2st;
		allocate_find_R2 = false;
	}

	this->find_R2 = new int*[GlobalC::ucell.nat];
	for(int iat=0; iat<GlobalC::ucell.nat; iat++)
	{
		// at least nad contains itself, so nad[iat] can not be 0.
		this->find_R2[iat] = new int[nad[iat]];
		ModuleBase::GlobalFunc::ZEROS(find_R2[iat], nad[iat]);
	}

	this->find_R2st = new int*[GlobalC::ucell.nat];
	for(int iat=0; iat<GlobalC::ucell.nat; iat++)
	{
		this->find_R2st[iat] = new int[nad[iat]];
		ModuleBase::GlobalFunc::ZEROS(find_R2st[iat], nad[iat]);
	}
	allocate_find_R2 = true;

//	GlobalV::ofs_running << std::setw(5) << "b1" << std::setw(5) << "b2" << std::setw(5) << "b3"
//	<< std::setw(8) << "iat" << std::setw(8) << "ad" << std::setw(8) << "iat2"
//	<< std::setw(8) << "find_R2" << std::setw(8) << "find_R2st" << std::setw(8) << "dis" << std::endl;
	for (int T1 = 0; T1 < GlobalC::ucell.ntype; T1++)
	{
		for (int I1 = 0; I1 < GlobalC::ucell.atoms[T1].na; I1++)
		{
//			std::cout << " T1=" << T1 << " I1=" << I1 << std::endl; 
			tau1 = GlobalC::ucell.atoms[T1].tau[I1];
			GlobalC::GridD.Find_atom(GlobalC::ucell, tau1, T1, I1);
			const int iat = GlobalC::ucell.itia2iat(T1,I1);

//			std::cout << " Number of adjacent = " << GlobalC::GridD.getAdjacentNum()+1 << std::endl;
			
			int count=0;
			for (int ad = 0; ad < GlobalC::GridD.getAdjacentNum()+1; ad++)
			{
		//		std::cout << " ad=" << ad << std::endl;
				const int T2 = GlobalC::GridD.getType(ad);
				const int I2 = GlobalC::GridD.getNatom(ad);
				const int iat2 = GlobalC::ucell.itia2iat(T2,I2);

				
				// if this atom is in this processor.
				if(GT.in_this_processor[iat])
				{
					if(GT.in_this_processor[iat2])
					{
						dtau = GlobalC::GridD.getAdjacentTau(ad) - tau1;
                        double distance = dtau.norm() * GlobalC::ucell.lat0;
                        double rcut = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut();

						const int b1 = GlobalC::GridD.getBox(ad).x;
						const int b2 = GlobalC::GridD.getBox(ad).y;
						const int b3 = GlobalC::GridD.getBox(ad).z;
                        
						// for test
						/*
						if( this->cal_RindexAtom(b1, b2, b3, iat2) == 232 )
						{
							std::cout << " ====== nnrg =========" << std::endl;
							std::cout << " index=" << cal_RindexAtom(b1, b2, b3, iat2) << std::endl;
							std::cout << " iat=" << iat << " iat2=" << iat2 << std::endl;
							std::cout << " R1 = " << tau1.x << " " << tau1.y << " " << tau1.z << std::endl;
							std::cout << " R2 = " << GlobalC::GridD.getAdjacentTau(ad).x 
							<< " " << GlobalC::GridD.getAdjacentTau(ad).y 
							<< " " << GlobalC::GridD.getAdjacentTau(ad).z << std::endl;
							std::cout << std::setprecision(25);
							std::cout << " distance = " << distance << std::endl;
							std::cout << " box = " << b1 << " " << b2 << " " << b3 << std::endl;
							std::cout << " rcut = " << rcut << std::endl;
						}
						*/
						

						//std::cout << " iat=" << iat << " find_R2=" << this->cal_RindexAtom(b1, b2, b3, iat2) <<
						// " b1=" << b1 << " b2=" << b2 << " b3=" << b3 << " iat2=" << iat2 << " distance=" << distance << std::endl;
					
						// mohan fix bug 2011-06-26, should be '<', not '<='	
						//			if(distance < rcut)

						// mohan reset this 2013-07-02 in Princeton
						// we should make absolutely sure that the distance is smaller than GlobalC::ORB.Phi[it].getRcut
						// this should be consistant with LCAO_nnr::cal_nnrg function 
						// typical example : 7 Bohr cutoff Si orbital in 14 Bohr length of cell.
						// distance = 7.0000000000000000
						// GlobalC::ORB.Phi[it].getRcut = 7.0000000000000008
						if(distance < rcut - 1.0e-15)
						{
						//	assert( count < nad[iat] );
							//--------------------------------------------------------------
							// start positions of adjacent atom of 'iat'
							// note: the first is not zero.
							//--------------------------------------------------------------
							find_R2[iat][count] = this->cal_RindexAtom(b1, b2, b3, iat2);


							/*if(iat==50 && iat2==96)
							{
								GlobalV::ofs_running << " ************** iat=" << iat << " count=" << count << " find_R2=" << find_R2[iat][count] << 
								" b1=" << b1 << " b2=" << b2 << " b3=" << b3 << " iat2=" << iat2 << " distance=" << distance 
								<< " rcut=" << rcut <<std::endl;
							}
							else if(find_R2[iat][count]==10536)
							{
								GlobalV::ofs_running << " ************** iat=" << iat << " count=" << count << " find_R2=" << find_R2[iat][count] << 
								" b1=" << b1 << " b2=" << b2 << " b3=" << b3 << " iat2=" << iat2 << " distance=" << distance 
								<< " rcut=" << rcut <<std::endl;
							}*/

							// find_R2st
							// note: the first must be zero.
							// find_R2st: start position of each adjacen atom.
							if( count + 1 < nad[iat] )
							{
								find_R2st[iat][count+1] = find_R2st[iat][count] 
								+ GlobalC::ucell.atoms[T1].nw * GlobalC::ucell.atoms[T2].nw; //modified by zhengdy-soc
							}
							++count;
						}
					}
				}
			}
		}
	}

	//---------
	// for test
	//---------
	/*
	GlobalV::ofs_running << " print find_R2 " << std::endl;
	for(int i=0; i<GlobalC::ucell.nat; i++)
	{
		for(int j=0; j<nad[i]; j++)
		{
			GlobalV::ofs_running << " i=" << i << " j=" << j << " find_R2=" << find_R2[i][j] << std::endl;
		}
	}
	GlobalV::ofs_running << std::endl;
	*/

	return;
}

void LCAO_nnr::cal_max_box_index(void)
{
	ModuleBase::TITLE("LCAO_nnr","cal_max_box_index");
	this->maxB1 = this->maxB2 = this->maxB3 = -10000;
	this->minB1 = this->minB2 = this->minB3 = 10000;
	for (int T1 = 0; T1 < GlobalC::ucell.ntype; T1++)
	{
		for (int I1 = 0; I1 < GlobalC::ucell.atoms[T1].na; I1++)
		{
			ModuleBase::Vector3<double> tau1 = GlobalC::ucell.atoms[T1].tau[I1];
			//GlobalC::GridD.Find_atom(tau1);
			GlobalC::GridD.Find_atom(GlobalC::ucell, tau1, T1, I1);
			for (int ad = 0; ad < GlobalC::GridD.getAdjacentNum()+1; ad++)
			{
				this->maxB1 = max( GlobalC::GridD.getBox(ad).x, maxB1 ); 
				this->maxB2 = max( GlobalC::GridD.getBox(ad).y, maxB2 ); 
				this->maxB3 = max( GlobalC::GridD.getBox(ad).z, maxB3 ); 

				this->minB1 = min( GlobalC::GridD.getBox(ad).x, minB1 ); 
				this->minB2 = min( GlobalC::GridD.getBox(ad).y, minB2 ); 
				this->minB3 = min( GlobalC::GridD.getBox(ad).z, minB3 ); 
			}
		}
	}

	/*
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"maxB1",maxB1);
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"maxB2",maxB2);
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"maxB3",maxB3);
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"minB1",minB1);
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"minB2",minB2);
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"minB3",minB3);
	*/

	nB1 = maxB1-minB1+1;
	nB2 = maxB2-minB2+1;
	nB3 = maxB3-minB3+1;

	/*
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nB1",nB1);
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nB2",nB2);
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nB3",nB3);
	*/

	nbox = nB1 * nB2 * nB3;
	
	//ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nbox",nbox);

	return;
}

int LCAO_nnr::cal_RindexAtom(const int &u1, const int &u2, const int &u3, const int &iat2)
{
	const int x1 = u1 - this->minB1;
	const int x2 = u2 - this->minB2;
	const int x3 = u3 - this->minB3;
	
	if(x1<0 || x2<0 || x3<0)
	{
		std::cout << " u1=" << u1 << " minB1=" << minB1 << std::endl;
		std::cout << " u2=" << u2 << " minB2=" << minB2 << std::endl;
		std::cout << " u3=" << u3 << " minB3=" << minB3 << std::endl;
		ModuleBase::WARNING_QUIT("LCAO_nnr::cal_Rindex","x1<0 || x2<0 || x3<0 !");
	}

	assert(x1>=0);
	assert(x2>=0);
	assert(x3>=0);

	return (iat2 + (x3 + x2 * this->nB3 + x1 * this->nB2 * this->nB3) * GlobalC::ucell.nat);
}


// be called in LCAO_Hamilt::calculate_Hk.
void LCAO_nnr::folding_fixedH(const int &ik)
{
	ModuleBase::TITLE("LCAO_nnr","folding_fixedH");
	ModuleBase::timer::tick("LCAO_nnr","folding_fixedH");

	int iat = 0;
	int index = 0;
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
			//GlobalC::GridD.Find_atom(tau1);
			GlobalC::GridD.Find_atom(GlobalC::ucell, tau1, T1, I1);
			Atom* atom1 = &GlobalC::ucell.atoms[T1];
			const int start = GlobalC::ucell.itiaiw2iwt(T1,I1,0);

			// (2) search among all adjacent atoms.
			for (int ad = 0; ad < GlobalC::GridD.getAdjacentNum()+1; ++ad)
			{
				const int T2 = GlobalC::GridD.getType(ad);
				const int I2 = GlobalC::GridD.getNatom(ad);
				Atom* atom2 = &GlobalC::ucell.atoms[T2];

				tau2 = GlobalC::GridD.getAdjacentTau(ad);
				dtau = tau2 - tau1;
				double distance = dtau.norm() * GlobalC::ucell.lat0;
				double rcut = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut();

				bool adj = false;

				if(distance < rcut) 
				{
					adj = true;
				}
				else if(distance >= rcut)
				{
					for (int ad0 = 0; ad0 < GlobalC::GridD.getAdjacentNum()+1; ++ad0)
					{
						const int T0 = GlobalC::GridD.getType(ad0); 
						const int I0 = GlobalC::GridD.getNatom(ad0); 
						//const int iat0 = GlobalC::ucell.itia2iat(T0, I0);
						//const int start0 = GlobalC::ucell.itiaiw2iwt(T0, I0, 0);

						tau0 = GlobalC::GridD.getAdjacentTau(ad0);
						dtau1 = tau0 - tau1;
						dtau2 = tau0 - tau2;

						double distance1 = dtau1.norm() * GlobalC::ucell.lat0;
						double distance2 = dtau2.norm() * GlobalC::ucell.lat0;

						double rcut1 = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();
						double rcut2 = GlobalC::ORB.Phi[T2].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();

						if( distance1 < rcut1 && distance2 < rcut2 )
						{
							adj = true;
							break;
						}
					}
				}

				if(adj) // mohan fix bug 2011-06-26, should not be '<='
				{
					// (3) calculate the nu of atom (T2, I2)
					const int start2 = GlobalC::ucell.itiaiw2iwt(T2,I2,0);
					//------------------------------------------------
					// exp(k dot dR)
					// dR is the index of box in Crystal coordinates
					//------------------------------------------------
					ModuleBase::Vector3<double> dR(GlobalC::GridD.getBox(ad).x, GlobalC::GridD.getBox(ad).y, GlobalC::GridD.getBox(ad).z); 
					const double arg = ( GlobalC::kv.kvec_d[ik] * dR ) * ModuleBase::TWO_PI;
					//const double arg = ( GlobalC::kv.kvec_d[ik] * GlobalC::GridD.getBox(ad) ) * ModuleBase::TWO_PI;
					const std::complex<double> kphase = std::complex <double> ( cos(arg),  sin(arg) );

					//--------------------------------------------------
					// calculate how many matrix elements are in 
					// this processor.
					//--------------------------------------------------
					for(int ii=0; ii<atom1->nw*GlobalV::NPOL; ii++)
					{
						// the index of orbitals in this processor
						const int iw1_all = start + ii;
						const int mu = GlobalC::ParaO.trace_loc_row[iw1_all];
						if(mu<0)continue;

						for(int jj=0; jj<atom2->nw*GlobalV::NPOL; jj++)
						{
							int iw2_all = start2 + jj;
							const int nu = GlobalC::ParaO.trace_loc_col[iw2_all];

							if(nu<0)continue;
							//const int iic = mu*GlobalC::ParaO.ncol+nu;
                            int iic;
                            if(GlobalV::KS_SOLVER=="genelpa" || GlobalV::KS_SOLVER=="scalapack_gvx")  // save the matrix as column major format
                            {
                                iic=mu+nu*GlobalC::ParaO.nrow;
                            }
                            else
                            {
                                iic=mu*GlobalC::ParaO.ncol+nu;
                            }

							//########################### EXPLAIN ###############################
							// 1. overlap matrix with k point
							// GlobalC::LM.SlocR = < phi_0i | phi_Rj >, where 0, R are the cell index
							// while i,j are the orbital index.

							// 2. H_fixed=T+Vnl matrix element with k point (if Vna is not used).
							// H_fixed=T+Vnl+Vna matrix element with k point (if Vna is used).
							// GlobalC::LM.Hloc_fixed = < phi_0i | H_fixed | phi_Rj>

							// 3. H(k) |psi(k)> = S(k) | psi(k)> 
							// Sloc2 is used to diagonalize for a give k point.
							// Hloc_fixed2 is used to diagonalize (eliminate index R).
							//###################################################################
							
							if(GlobalV::NSPIN!=4)
							{
								GlobalC::LM.Sloc2[iic] += GlobalC::LM.SlocR[index] * kphase;
								GlobalC::LM.Hloc_fixed2[iic] += GlobalC::LM.Hloc_fixedR[index] * kphase;
							}
							else
							{
								GlobalC::LM.Sloc2[iic] += GlobalC::LM.SlocR_soc[index] * kphase;
								GlobalC::LM.Hloc_fixed2[iic] += GlobalC::LM.Hloc_fixedR_soc[index] * kphase;
							}
							++index;

						}//end jj
					}//end ii
				}
			}// end ad
			++iat;
		}// end I1
	} // end T1

	assert(index==this->nnr);

	ModuleBase::timer::tick("LCAO_nnr","folding_fixedH");
	return;
}
