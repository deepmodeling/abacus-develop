#include "LCAO_nnr.h"
#include "../src_pw/global.h"
#include "record_adj.h" //mohan add 2012-07-06
#include "dftu.h"       //quxin add 2020-10-29

//----------------------------
// define a global class obj.
LCAO_nnr LNNR;
//----------------------------

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
		for(int iat=0; iat<ucell.nat; iat++)
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
	TITLE("LCAO_nnr","cal_nnr");

	delete[] nlocdim;
	delete[] nlocstart;
	nlocdim = new int[ucell.nat];	
	nlocstart = new int[ucell.nat];
	ZEROS(nlocdim, ucell.nat);
	ZEROS(nlocstart, ucell.nat);

	this->nnr = 0;
	int start = 0;
	int ind1 = 0;
	int iat = 0;

	// (1) find the adjacent atoms of atom[T1,I1];
	Vector3<double> tau1, tau2, dtau;
	Vector3<double> tau0, dtau1, dtau2;
	for (int T1 = 0; T1 < ucell.ntype; T1++)
	{
		for (int I1 = 0; I1 < ucell.atoms[T1].na; I1++)
		{
			tau1 = ucell.atoms[T1].tau[I1];
			//GridD.Find_atom( tau1 );
			GridD.Find_atom( tau1 ,T1, I1);
			const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
			this->nlocstart[iat] = nnr;
			int nw1 = ucell.atoms[T1].nw * NPOL;

			// (2) search among all adjacent atoms.
			for (int ad = 0; ad < GridD.getAdjacentNum()+1; ad++)
			{
				const int T2 = GridD.getType(ad);
				const int I2 = GridD.getNatom(ad);
				const int iat2 = ucell.itia2iat(T2, I2);
				const int start2 = ucell.itiaiw2iwt(T2, I2, 0);
				int nw2 = ucell.atoms[T2].nw * NPOL;

				tau2 = GridD.getAdjacentTau(ad);

				dtau = tau2 - tau1;
				double distance = dtau.norm() * ucell.lat0;
				double rcut = ORB.Phi[T1].getRcut() + ORB.Phi[T2].getRcut();
				
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
						const int mu = ParaO.trace_loc_row[iw1_all];
						if(mu<0)continue;

						for(int jj=0; jj<nw2; jj++)
						{
							const int iw2_all = start2 + jj;
							const int nu = ParaO.trace_loc_col[iw2_all];
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
					for (int ad0 = 0; ad0 < GridD.getAdjacentNum()+1; ++ad0)
					{
						const int T0 = GridD.getType(ad0);
						const int I0 = GridD.getNatom(ad0);
						const int iat0 = ucell.itia2iat(T0, I0);
						const int start0 = ucell.itiaiw2iwt(T0, I0, 0);
					
						tau0 = GridD.getAdjacentTau(ad0);
						dtau1 = tau0 - tau1; 
						double distance1 = dtau1.norm() * ucell.lat0;
						double rcut1 = ORB.Phi[T1].getRcut() + ORB.Beta[T0].get_rcut_max();

						dtau2 = tau0 - tau2;
						double distance2 = dtau2.norm() * ucell.lat0;
						double rcut2 = ORB.Phi[T2].getRcut() + ORB.Beta[T0].get_rcut_max();

						if( distance1 < rcut1 && distance2 < rcut2 )
						{
							for(int ii=0; ii<nw1; ++ii)
							{
								const int iw1_all = start1 + ii;
								const int mu = ParaO.trace_loc_row[iw1_all];
								if(mu<0)continue;

								for(int jj=0; jj<nw2; ++jj)
								{
									const int iw2_all = start2 + jj;
									const int nu = ParaO.trace_loc_col[iw2_all];
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

	//xiaohui add 'OUT_LEVEL' line, 2015-09-16
	if(OUT_LEVEL != "m") OUT(ofs_running,"nnr",nnr);
//	for(int iat=0; iat<ucell.nat; iat++)
//	{
//		cout << " nlocdim[" << iat << "]=" << nlocdim[iat];
//		cout << " nlocstart[" << iat << "]=" << nlocstart[iat] << endl;
//	}



	return;
}

// This is for cell R dependent part. 
void LCAO_nnr::cal_nnrg(const Grid_Technique &GT)
{
	TITLE("LCAO_nnr","cal_nnrg");

	this->cal_max_box_index();

	this->nnrg = 0;

	delete[] nlocdimg;
	delete[] nlocstartg;
	delete[] nad;
	
	this->nad = new int[ucell.nat];
	this->nlocdimg = new int[ucell.nat];	
	this->nlocstartg = new int[ucell.nat];
	
	ZEROS(nad, ucell.nat);
	ZEROS(nlocdimg, ucell.nat);
	ZEROS(nlocstartg, ucell.nat);

//	stringstream ss;
//	ss << global_out_dir << "sltk_box.dat";
//	ofstream ofs(ss.str().c_str());

	Vector3<double> tau1, tau2, dtau;
	Vector3<double> dtau1, dtau2, tau0;
	for (int T1 = 0; T1 < ucell.ntype; ++T1)
	{
		Atom* atom1 = &ucell.atoms[T1];
		for (int I1 = 0; I1 < atom1->na; ++I1)
		{
			tau1 = atom1->tau[I1];
			//GridD.Find_atom(tau1);
			GridD.Find_atom(tau1, T1, I1);
			const int iat = ucell.itia2iat(T1,I1);

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
				for (int ad = 0; ad < GridD.getAdjacentNum()+1; ++ad)
				{
					const int T2 = GridD.getType(ad);
					const int I2 = GridD.getNatom(ad);
					const int iat2 = ucell.itia2iat(T2, I2);
					Atom* atom2 = &ucell.atoms[T2]; 

					// if the adjacent atom is in this processor.
					if(GT.in_this_processor[iat2])
					{
						tau2 = GridD.getAdjacentTau(ad);
						dtau = GridD.getAdjacentTau(ad) - tau1;
						double distance = dtau.norm() * ucell.lat0;
						double rcut = ORB.Phi[T1].getRcut() + ORB.Phi[T2].getRcut();


						//if(distance < rcut)
			// mohan reset this 2013-07-02 in Princeton
			// we should make absolutely sure that the distance is smaller than ORB.Phi[it].getRcut
			// this should be consistant with LCAO_nnr::cal_nnrg function 
			// typical example : 7 Bohr cutoff Si orbital in 14 Bohr length of cell.
			// distance = 7.0000000000000000
			// ORB.Phi[it].getRcut = 7.0000000000000008
						if(distance < rcut - 1.0e-15)
						{
							const int nelement = atom1->nw * atom2->nw;//modified by zhengdy-soc, no need to double
							this->nnrg += nelement;
							this->nlocdimg[iat] += nelement; 
							this->nad[iat]++;
							++count;

							/*
							   ofs << setw(10) << iat << setw(10) << iat2
							   << setw(10) << GridD.getBox(ad).x 
							   << setw(10) << GridD.getBox(ad).y 
							   << setw(10) << GridD.getBox(ad).z 
							   << setw(20) << distance << endl;
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
							for (int ad0 = 0; ad0 < GridD.getAdjacentNum()+1; ++ad0)
							{
								const int T0 = GridD.getType(ad0);
								const int I0 = GridD.getNatom(ad0);
								const int iat0 = ucell.itia2iat(T0, I0);
								const int start0 = ucell.itiaiw2iwt(T0, I0, 0);

								tau0 = GridD.getAdjacentTau(ad0);	
								dtau1 = tau0 - tau1; 
								dtau2 = tau0 - tau2;

								double distance1 = dtau1.norm() * ucell.lat0;
								double distance2 = dtau2.norm() * ucell.lat0;

								double rcut1 = ORB.Phi[T1].getRcut() + ORB.Beta[T0].get_rcut_max();
								double rcut2 = ORB.Phi[T2].getRcut() + ORB.Beta[T0].get_rcut_max();

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
//				ofs_running << " iat=" << iat << " nlocstartg=" << nlocstartg[iat] << " nad=" << nad[iat] << endl;
			}// end iat
		}// end I1
	}// end T1
	//ofs.close();
	if(OUT_LEVEL != "m") OUT(ofs_running,"nnrg",this->nnrg);


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
		for(int iat=0; iat<ucell.nat; iat++)delete[] find_R2[iat];
		delete[] find_R2;

		for(int iat=0; iat<ucell.nat; iat++)delete[] find_R2st[iat];
		delete[] find_R2st;
		allocate_find_R2 = false;
	}

	this->find_R2 = new int*[ucell.nat];
	for(int iat=0; iat<ucell.nat; iat++)
	{
		// at least nad contains itself, so nad[iat] can not be 0.
		this->find_R2[iat] = new int[nad[iat]];
		ZEROS(find_R2[iat], nad[iat]);
	}

	this->find_R2st = new int*[ucell.nat];
	for(int iat=0; iat<ucell.nat; iat++)
	{
		this->find_R2st[iat] = new int[nad[iat]];
		ZEROS(find_R2st[iat], nad[iat]);
	}
	allocate_find_R2 = true;

//	ofs_running << setw(5) << "b1" << setw(5) << "b2" << setw(5) << "b3"
//	<< setw(8) << "iat" << setw(8) << "ad" << setw(8) << "iat2"
//	<< setw(8) << "find_R2" << setw(8) << "find_R2st" << setw(8) << "dis" << endl;
	for (int T1 = 0; T1 < ucell.ntype; T1++)
	{
		for (int I1 = 0; I1 < ucell.atoms[T1].na; I1++)
		{
//			cout << " T1=" << T1 << " I1=" << I1 << endl; 
			tau1 = ucell.atoms[T1].tau[I1];
			//GridD.Find_atom(tau1);
			GridD.Find_atom(tau1, T1, I1);
			const int iat = ucell.itia2iat(T1,I1);

//			cout << " Number of adjacent = " << GridD.getAdjacentNum()+1 << endl;
			
			int count=0;
			for (int ad = 0; ad < GridD.getAdjacentNum()+1; ad++)
			{
		//		cout << " ad=" << ad << endl;
				const int T2 = GridD.getType(ad);
				const int I2 = GridD.getNatom(ad);
				const int iat2 = ucell.itia2iat(T2,I2);

				
				// if this atom is in this processor.
				if(GT.in_this_processor[iat])
				{
					if(GT.in_this_processor[iat2])
					{
						dtau = GridD.getAdjacentTau(ad) - tau1;
                        double distance = dtau.norm() * ucell.lat0;
                        double rcut = ORB.Phi[T1].getRcut() + ORB.Phi[T2].getRcut();

						const int b1 = GridD.getBox(ad).x;
						const int b2 = GridD.getBox(ad).y;
						const int b3 = GridD.getBox(ad).z;
                        
						// for test
						/*
						if( this->cal_RindexAtom(b1, b2, b3, iat2) == 232 )
						{
							cout << " ====== nnrg =========" << endl;
							cout << " index=" << cal_RindexAtom(b1, b2, b3, iat2) << endl;
							cout << " iat=" << iat << " iat2=" << iat2 << endl;
							cout << " R1 = " << tau1.x << " " << tau1.y << " " << tau1.z << endl;
							cout << " R2 = " << GridD.getAdjacentTau(ad).x 
							<< " " << GridD.getAdjacentTau(ad).y 
							<< " " << GridD.getAdjacentTau(ad).z << endl;
							cout << setprecision(25);
							cout << " distance = " << distance << endl;
							cout << " box = " << b1 << " " << b2 << " " << b3 << endl;
							cout << " rcut = " << rcut << endl;
						}
						*/
						

						//cout << " iat=" << iat << " find_R2=" << this->cal_RindexAtom(b1, b2, b3, iat2) <<
						// " b1=" << b1 << " b2=" << b2 << " b3=" << b3 << " iat2=" << iat2 << " distance=" << distance << endl;
					
						// mohan fix bug 2011-06-26, should be '<', not '<='	
			//			if(distance < rcut)

			// mohan reset this 2013-07-02 in Princeton
			// we should make absolutely sure that the distance is smaller than ORB.Phi[it].getRcut
			// this should be consistant with LCAO_nnr::cal_nnrg function 
			// typical example : 7 Bohr cutoff Si orbital in 14 Bohr length of cell.
			// distance = 7.0000000000000000
			// ORB.Phi[it].getRcut = 7.0000000000000008
						if(distance < rcut - 1.0e-15)
						{
						//	assert( count < nad[iat] );
							//--------------------------------------------------------------
							// start positions of adjacent atom of 'iat'
							// note: the first is not zero.
							//--------------------------------------------------------------
							find_R2[iat][count] = this->cal_RindexAtom(b1, b2, b3, iat2);


							if(iat==50 && iat2==96)
							{
								ofs_running << " ************** iat=" << iat << " count=" << count << " find_R2=" << find_R2[iat][count] << 
								" b1=" << b1 << " b2=" << b2 << " b3=" << b3 << " iat2=" << iat2 << " distance=" << distance 
								<< " rcut=" << rcut <<endl;
							}
							else if(find_R2[iat][count]==10536)
							{
								ofs_running << " ************** iat=" << iat << " count=" << count << " find_R2=" << find_R2[iat][count] << 
								" b1=" << b1 << " b2=" << b2 << " b3=" << b3 << " iat2=" << iat2 << " distance=" << distance 
								<< " rcut=" << rcut <<endl;
							}

							// find_R2st
							// note: the first must be zero.
							// find_R2st: start position of each adjacen atom.
							if( count + 1 < nad[iat] )
							{
								find_R2st[iat][count+1] = find_R2st[iat][count] + ucell.atoms[T1].nw * ucell.atoms[T2].nw; //modified by zhengdy-soc
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
	ofs_running << " print find_R2 " << endl;
	for(int i=0; i<ucell.nat; i++)
	{
		for(int j=0; j<nad[i]; j++)
		{
			ofs_running << " i=" << i << " j=" << j << " find_R2=" << find_R2[i][j] << endl;
		}
	}
	ofs_running << endl;
	*/



	return;
}

void LCAO_nnr::cal_max_box_index(void)
{
	TITLE("LCAO_nnr","cal_max_box_index");
	this->maxB1 = this->maxB2 = this->maxB3 = -10000;
	this->minB1 = this->minB2 = this->minB3 = 10000;
	for (int T1 = 0; T1 < ucell.ntype; T1++)
	{
		for (int I1 = 0; I1 < ucell.atoms[T1].na; I1++)
		{
			Vector3<double> tau1 = ucell.atoms[T1].tau[I1];
			//GridD.Find_atom(tau1);
			GridD.Find_atom(tau1, T1, I1);
			for (int ad = 0; ad < GridD.getAdjacentNum()+1; ad++)
			{
				this->maxB1 = max( GridD.getBox(ad).x, maxB1 ); 
				this->maxB2 = max( GridD.getBox(ad).y, maxB2 ); 
				this->maxB3 = max( GridD.getBox(ad).z, maxB3 ); 

				this->minB1 = min( GridD.getBox(ad).x, minB1 ); 
				this->minB2 = min( GridD.getBox(ad).y, minB2 ); 
				this->minB3 = min( GridD.getBox(ad).z, minB3 ); 
			}
		}
	}

	/*
	OUT(ofs_running,"maxB1",maxB1);
	OUT(ofs_running,"maxB2",maxB2);
	OUT(ofs_running,"maxB3",maxB3);
	OUT(ofs_running,"minB1",minB1);
	OUT(ofs_running,"minB2",minB2);
	OUT(ofs_running,"minB3",minB3);
	*/

	nB1 = maxB1-minB1+1;
	nB2 = maxB2-minB2+1;
	nB3 = maxB3-minB3+1;

	/*
	OUT(ofs_running,"nB1",nB1);
	OUT(ofs_running,"nB2",nB2);
	OUT(ofs_running,"nB3",nB3);
	*/

	nbox = nB1 * nB2 * nB3;
	
	//OUT(ofs_running,"nbox",nbox);

	return;
}

int LCAO_nnr::cal_RindexAtom(const int &u1, const int &u2, const int &u3, const int &iat2)
{
	const int x1 = u1 - this->minB1;
	const int x2 = u2 - this->minB2;
	const int x3 = u3 - this->minB3;
	
	if(x1<0 || x2<0 || x3<0)
	{
		cout << " u1=" << u1 << " minB1=" << minB1 << endl;
		cout << " u2=" << u2 << " minB2=" << minB2 << endl;
		cout << " u3=" << u3 << " minB3=" << minB3 << endl;
		WARNING_QUIT("LCAO_nnr::cal_Rindex","x1<0 || x2<0 || x3<0 !");
	}

	assert(x1>=0);
	assert(x2>=0);
	assert(x3>=0);

	return (iat2 + (x3 + x2 * this->nB3 + x1 * this->nB2 * this->nB3) * ucell.nat);
}


// be called in LCAO_Hamilt::calculate_Hk.
void LCAO_nnr::folding_fixedH(const int &ik)
{
	TITLE("LCAO_nnr","folding_fixedH");
	timer::tick("LCAO_nnr","folding_fixedH",'G');

//	cout << " kvec_c = " << kv.kvec_c[ik].x << " " << kv.kvec_c[ik].y << " " << kv.kvec_c[ik].z << endl;
//	Record_adj RA;
//	RA.for_2d();

	//Quxin added for DFT+U calculation
 	if(INPUT.dft_plus_u) 
	{
		ZEROS( VECTOR_TO_PTR(dftu.Sm_k.at(ik)), ParaO.nloc);
	}

	int iat = 0;
	int index = 0;
	Vector3<double> dtau, tau1, tau2;
	Vector3<double> dtau1, dtau2, tau0;
	for (int T1 = 0; T1 < ucell.ntype; ++T1)
	{
		Atom* atom1 = &ucell.atoms[T1];
		for (int I1 = 0; I1 < atom1->na; ++I1)
		{
			tau1 = atom1->tau[I1];
			//GridD.Find_atom(tau1);
			GridD.Find_atom(tau1, T1, I1);
			Atom* atom1 = &ucell.atoms[T1];
			const int start = ucell.itiaiw2iwt(T1,I1,0);

			// (2) search among all adjacent atoms.
			for (int ad = 0; ad < GridD.getAdjacentNum()+1; ++ad)
			{
				const int T2 = GridD.getType(ad);
				const int I2 = GridD.getNatom(ad);
				Atom* atom2 = &ucell.atoms[T2];

				tau2 = GridD.getAdjacentTau(ad);
				dtau = tau2 - tau1;
				double distance = dtau.norm() * ucell.lat0;
				double rcut = ORB.Phi[T1].getRcut() + ORB.Phi[T2].getRcut();

				bool adj = false;

				if(distance < rcut) adj = true;
				else if(distance >= rcut)
				{
					for (int ad0 = 0; ad0 < GridD.getAdjacentNum()+1; ++ad0)
					{
						const int T0 = GridD.getType(ad0); 
						const int I0 = GridD.getNatom(ad0); 
						const int iat0 = ucell.itia2iat(T0, I0);
						const int start0 = ucell.itiaiw2iwt(T0, I0, 0);

						tau0 = GridD.getAdjacentTau(ad0);
						dtau1 = tau0 - tau1;
						dtau2 = tau0 - tau2;

						double distance1 = dtau1.norm() * ucell.lat0;
						double distance2 = dtau2.norm() * ucell.lat0;

						double rcut1 = ORB.Phi[T1].getRcut() + ORB.Beta[T0].get_rcut_max();
						double rcut2 = ORB.Phi[T2].getRcut() + ORB.Beta[T0].get_rcut_max();

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
					const int start2 = ucell.itiaiw2iwt(T2,I2,0);
					//------------------------------------------------
					// exp(k dot dR)
					// dR is the index of box in Crystal coordinates
					//------------------------------------------------
					Vector3<double> dR(GridD.getBox(ad).x, GridD.getBox(ad).y, GridD.getBox(ad).z); 
					const double arg = ( kv.kvec_d[ik] * dR ) * TWO_PI;
					//const double arg = ( kv.kvec_d[ik] * GridD.getBox(ad) ) * TWO_PI;
					const complex<double> kphase = complex <double> ( cos(arg),  sin(arg) );

					//--------------------------------------------------
					// calculate how many matrix elements are in 
					// this processor.
					//--------------------------------------------------
					for(int ii=0; ii<atom1->nw*NPOL; ii++)
					{
						// the index of orbitals in this processor
						const int iw1_all = start + ii;
						const int mu = ParaO.trace_loc_row[iw1_all];
						if(mu<0)continue;

						for(int jj=0; jj<atom2->nw*NPOL; jj++)
						{
							int iw2_all = start2 + jj;
							const int nu = ParaO.trace_loc_col[iw2_all];

							if(nu<0)continue;
							//const int iic = mu*ParaO.ncol+nu;
                            int iic;
                            if(KS_SOLVER=="genelpa" || KS_SOLVER=="scalapack_gvx")  // save the matrix as column major format
                            {
                                iic=mu+nu*ParaO.nrow;
                            }
                            else
                            {
                                iic=mu*ParaO.ncol+nu;
                            }

							//########################### EXPLAIN ###############################
							// 1. overlap matrix with k point
							// LM.SlocR = < phi_0i | phi_Rj >, where 0, R are the cell index
							// while i,j are the orbital index.

							// 2. H_fixed=T+Vnl matrix element with k point (if Vna is not used).
							// H_fixed=T+Vnl+Vna matrix element with k point (if Vna is used).
							// LM.Hloc_fixed = < phi_0i | H_fixed | phi_Rj>

							// 3. H(k) |psi(k)> = S(k) | psi(k)> 
							// Sloc2 is used to diagonalize for a give k point.
							// Hloc_fixed2 is used to diagonalize (eliminate index R).
							//###################################################################
							
							if(NSPIN!=4)
							{
								LM.Sloc2[iic] += LM.SlocR[index] * kphase;
								LM.Hloc_fixed2[iic] += LM.Hloc_fixedR[index] * kphase;

								//quxin added for DFT+U calculation
							 	if(INPUT.dft_plus_u) dftu.Sm_k.at(ik).at(iic) += LM.SlocR[index] * kphase;
							}
							else
							{
								LM.Sloc2[iic] += LM.SlocR_soc[index] * kphase;
								LM.Hloc_fixed2[iic] += LM.Hloc_fixedR_soc[index] * kphase;

								//quxin added for DFT+U calculation
							 	if(INPUT.dft_plus_u) 
								{
									dftu.Sm_k.at(ik).at(iic) += LM.SlocR_soc[index] * kphase;
								}
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

//	cout << " folding Hfixed" << endl;

	timer::tick("LCAO_nnr","folding_fixedH",'G');
	return;
}
