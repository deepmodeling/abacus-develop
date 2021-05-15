#include "FORCE_k.h"
#include "../src_pw/global.h"
#include "dftu.h"  //Quxin add for DFT+U on 20201029

Force_LCAO_k::Force_LCAO_k ()
{
}

Force_LCAO_k::~Force_LCAO_k ()
{
}

#include "LCAO_nnr.h"
// be called in Force_LCAO::start_force_calculation
void Force_LCAO_k::ftable_k (
		const bool isforce,
		const bool isstress,
		matrix& foverlap,
		matrix& ftvnl_dphi,
		matrix& fvnl_dbeta,	
		matrix& fvl_dphi,
		matrix& soverlap,
		matrix& stvnl_dphi,
		matrix& svnl_dbeta,
		matrix& svl_dphi
		)
{
    TITLE("Force_LCAO_k", "ftable_k");
	timer::tick("Force_LCAO_k","ftable_k",'F');
	
	this->allocate_k();

	// calculate the energy density matrix
	// and the force related to overlap matrix and energy density matrix.
	this->cal_foverlap_k(isforce, isstress, foverlap, soverlap);

	// calculate the density matrix
	double** dm2d = new double*[NSPIN];
	for(int is=0; is<NSPIN; is++)
	{
		dm2d[is] = new double[LNNR.nnr];
		ZEROS(dm2d[is], LNNR.nnr);
	}
	Memory::record ("Force_LCAO_k", "dm2d", NSPIN*LNNR.nnr, "double");	
	bool with_energy = false;

	
	this->set_EDM_k(dm2d, with_energy);
	
	this->cal_ftvnl_dphi_k(dm2d, isforce, isstress, ftvnl_dphi, stvnl_dphi);

	//Quxin add for DFT+U on 20201029
	if(INPUT.dft_plus_u) 
	{
		dftu.force_stress();
		dftu.erase_force_stress();
	}

	// ---------------------------------------
	// doing on the real space grid.
	// ---------------------------------------
	this->cal_fvl_dphi_k(dm2d, isforce, isstress, fvl_dphi, svl_dphi);

	this->cal_fvnl_dbeta_k(dm2d, isforce, isstress, fvnl_dbeta, svnl_dbeta);

	for(int is=0; is<NSPIN; is++)
	{
		delete[] dm2d[is];
	}
	delete[] dm2d;

	//----------------------------------------------------------------
	// reduce the force according to 2D distribution of H & S matrix.
	//----------------------------------------------------------------
    if(isforce)
	{
        Parallel_Reduce::reduce_double_pool( foverlap.c, foverlap.nr * foverlap.nc);
        Parallel_Reduce::reduce_double_pool( ftvnl_dphi.c, ftvnl_dphi.nr * ftvnl_dphi.nc);
        Parallel_Reduce::reduce_double_pool( fvnl_dbeta.c, fvnl_dbeta.nr * fvnl_dbeta.nc);
        Parallel_Reduce::reduce_double_pool( fvl_dphi.c, fvl_dphi.nr * fvl_dphi.nc);
	}
    if(isstress)
    {
        Parallel_Reduce::reduce_double_pool( soverlap.c, soverlap.nr * soverlap.nc);
        Parallel_Reduce::reduce_double_pool( stvnl_dphi.c, stvnl_dphi.nr * stvnl_dphi.nc);
        Parallel_Reduce::reduce_double_pool( svnl_dbeta.c, svnl_dbeta.nr * svnl_dbeta.nc);
        Parallel_Reduce::reduce_double_pool( svl_dphi.c, svl_dphi.nr * svl_dphi.nc);
    }

	// test the force.
	/*
	cout << " overlap force" << endl;
	for(int iat=0; iat<ucell.nat; ++iat)
	{
		const double fac = Ry_to_eV / 0.529177;
		cout << setw(5) << iat+1 << setw(15) << foverlap[iat][0] *fac<< setw(15) << foverlap[iat][1]*fac << 
		setw(15) << foverlap[iat][2]*fac << endl;
	}
	*/

	this->finish_k();

	timer::tick("Force_LCAO_k","ftable_k",'F');
    return;
}

void Force_LCAO_k::allocate_k(void)
{
	TITLE("Force_LCAO_k","allocate_k");
	timer::tick("Force_LCAO_k","allocate_k");

	const int nnr = LNNR.nnr;
	//--------------------------------
    // (1) allocate for dSx dSy & dSz
	//--------------------------------
	LM.DSloc_Rx = new double [nnr];
    LM.DSloc_Ry = new double [nnr];
    LM.DSloc_Rz = new double [nnr];
    ZEROS(LM.DSloc_Rx, nnr);
    ZEROS(LM.DSloc_Ry, nnr);
    ZEROS(LM.DSloc_Rz, nnr);
	Memory::record("force_lo", "dS", nnr*3, "double");
    
	if(STRESS){
		LM.DH_r = new double [3* nnr];
		ZEROS(LM.DH_r, 3 * nnr);
		LM.stvnl11 = new double [nnr];
		LM.stvnl12 = new double [nnr];
		LM.stvnl13 = new double [nnr];
		LM.stvnl22 = new double [nnr];
		LM.stvnl23 = new double [nnr];
		LM.stvnl33 = new double [nnr];
		ZEROS(LM.stvnl11,  nnr);
		ZEROS(LM.stvnl12,  nnr);
		ZEROS(LM.stvnl13,  nnr);
		ZEROS(LM.stvnl22,  nnr);
		ZEROS(LM.stvnl23,  nnr);
		ZEROS(LM.stvnl33,  nnr);
		Memory::record("stress_lo", "dSR", nnr*6, "double");
	}

	//-----------------------------
	// calculate dS = <phi | dphi> 
	//-----------------------------
    // tips: build_ST_new --> ParaO.set_force 
	bool cal_deri = true;
	UHM.genH.build_ST_new ('S', cal_deri);

	//-----------------------------------------
	// (2) allocate for <phi | T + Vnl | dphi>
	//-----------------------------------------
    LM.DHloc_fixedR_x = new double [nnr];
    LM.DHloc_fixedR_y = new double [nnr];
    LM.DHloc_fixedR_z = new double [nnr];
    ZEROS (LM.DHloc_fixedR_x, nnr);
    ZEROS (LM.DHloc_fixedR_y, nnr);
    ZEROS (LM.DHloc_fixedR_z, nnr);
	Memory::record("force_lo", "dTVNL", nnr*3, "double");
    
    // calculate dT=<phi|kin|dphi> in LCAO
    // calculate T + VNL(P1) in LCAO basis
    UHM.genH.build_ST_new ('T', cal_deri);
	//test(LM.DHloc_fixedR_x,"LM.DHloc_fixedR_x T part");
   
   	// calculate dVnl=<phi|dVnl|dphi> in LCAO 
	UHM.genH.build_Nonlocal_mu (cal_deri);
	//test(LM.DHloc_fixedR_x,"LM.DHloc_fixedR_x Vnl part");

	timer::tick("Force_LCAO_k","allocate");
	return;
}

void Force_LCAO_k::finish_k(void)
{
    delete [] LM.DSloc_Rx;
    delete [] LM.DSloc_Ry;
    delete [] LM.DSloc_Rz;
    delete [] LM.DHloc_fixedR_x;
    delete [] LM.DHloc_fixedR_y;
    delete [] LM.DHloc_fixedR_z;
	if(STRESS)
	{
		delete [] LM.DH_r;
		delete [] LM.stvnl11;
		delete [] LM.stvnl12;
		delete [] LM.stvnl13;
		delete [] LM.stvnl22;
		delete [] LM.stvnl23;
		delete [] LM.stvnl33;
	}
	return;
}

#include "record_adj.h"
#include "LCAO_nnr.h"
void Force_LCAO_k::set_EDM_k(double** dm2d, const bool with_energy)
{
	TITLE("Force_LCAO_k","set_EDM_k");
	timer::tick("Force_LCAO_k","set_EDM_k",'G');

	Vector3<double> tau1, dtau;

	//----------------------------------------------------------
	// RA will set the adjacent information for each atom
	// 2d means this adjacent information is for HPSEPS's kind
	// of division of H matrix.
	//----------------------------------------------------------
 	//xiaohui add "OUT_LEVEL", 2015-09-16
	if(OUT_LEVEL != "m") ofs_running << " Calculate the energy density matrix with k " << endl;
	Record_adj RA;
	RA.for_2d();

	//------------------------
	// circle for each atom
	//------------------------
	for(int T1=0; T1<ucell.ntype; ++T1)
    {
        Atom* atom1 = &ucell.atoms[T1];
		for(int I1=0; I1<atom1->na; ++I1)
		{
			const int iat = ucell.itia2iat(T1,I1);
			const int start1 = ucell.itiaiw2iwt(T1,I1,0);
			const int gstart = LNNR.nlocstart[iat];
			const int irr = LNNR.nlocdim[iat];//number of adjacet orbitals

			complex<double> **vvv = new complex<double>*[NSPIN];
         //xiaohui add 2014-03-17, add "if(irr > 0)", 
         //meaing only allocate memory when number of iat-th atom adjacet orbitals are not zero in this processor
         if(irr > 0)
         {
   			for(int is=0; is<NSPIN; is++)
			   {
				   vvv[is] = new complex<double>[ irr ];
				   ZEROS(vvv[is], irr );
			   }
         }

			int ispin=0;
			complex<double> *dm;
			//----------------
			// circle for k
			//----------------
			for(int ik=0; ik<kv.nks; ++ik)
			{
				// set the spin direction
				if(NSPIN==2)
				{
					ispin = kv.isk[ik];
				}
				dm = vvv[ispin];
				//------------------
				// circle for bands
				//------------------
				for(int ib=0; ib<NBANDS; ++ib)
				{
					const double w1=wf.wg(ik,ib);
					if(w1>0)
					{
						//-----------------------------
						// start the adjacent cycle.
						//-----------------------------
						complex<double> *wfc = LOWF.WFC_K[ik][ib];
						int count = 0;
						for (int cb = 0; cb < RA.na_each[iat]; ++cb)
						{
							const int T2 = RA.info[iat][cb][3];
							const int I2 = RA.info[iat][cb][4];
							Atom* atom2 = &ucell.atoms[T2];

							//-----------------
							// exp[i * R * k]
							//-----------------
							const complex<double> phase = w1 * exp(TWO_PI * IMAG_UNIT * (
										kv.kvec_d[ik].x * RA.info[iat][cb][0] +
										kv.kvec_d[ik].y * RA.info[iat][cb][1] +
										kv.kvec_d[ik].z * RA.info[iat][cb][2]
										) );

							const int start2 = ucell.itiaiw2iwt(T2,I2,0);

							for(int jj=0; jj<atom1->nw; ++jj)
							{
								const int iw1_all = start1 + jj;
								
								// 2D division (HPSEPS)
								const int mu = ParaO.trace_loc_row[iw1_all];
								if(mu<0) continue;
								const int mug = GridT.trace_lo[iw1_all];

								for(int kk=0; kk<atom2->nw; ++kk)
								{
									const int iw2_all = start2 + kk;
									
									// 2D division (HPSEPS)
									const int nu = ParaO.trace_loc_col[iw2_all];
									if(nu<0) continue;
									const int nug = GridT.trace_lo[iw2_all];
	
									if(mug >= 0 && nug >= 0)
									{
										dm[count] += set_EDM_k_element(phase, with_energy, 
										wfc[mug], wfc[nug], wf.ekb[ik][ib]); 
									}
									else if( mug >= 0 && nug <= 0)
									{
										const int a4 = LOWF.trace_aug[iw2_all];
									
										//assert(a4>=0);

										dm[count] += set_EDM_k_element(phase, with_energy, wfc[mug], 
										LOWF.WFC_K_aug[ik][ib][a4], wf.ekb[ik][ib]); 
									}
									else if( mug <= 0 && nug >= 0)
									{
										const int a3 = LOWF.trace_aug[iw1_all]; 

										dm[count] += set_EDM_k_element(phase, with_energy, 
										LOWF.WFC_K_aug[ik][ib][a3], wfc[nug], wf.ekb[ik][ib]); 
									}
									else if( mug <=0 && nug <=0 )
									{
										const int a1 = LOWF.trace_aug[iw1_all];
										const int a2 = LOWF.trace_aug[iw2_all];

										dm[count] += set_EDM_k_element(phase, with_energy, 
										LOWF.WFC_K_aug[ik][ib][a1], LOWF.WFC_K_aug[ik][ib][a2], wf.ekb[ik][ib]); 
									}
									assert(count<irr);
									++ count;
								}//kk
							}//jj
						}// cb
//						ofs_running << " count = " << count << endl;
						assert(count == LNNR.nlocdim[iat]);
					}// w1
				}//ib
			}//ik

			//--------------------------------------
			// get the real value density matrix or
			// energy density matrix 
			//--------------------------------------
			for(int is=0; is<NSPIN; ++is)
			{
				for(int iv=0; iv<irr; ++iv)
				{
					dm2d[is][gstart+iv] = vvv[is][iv].real();
				}
            //xiaohui add 2014-03-17, add "if(irr > 0)"
            //meaning delete memory when number of iat-th atom adjacet orbitals are not zero in this processor 
            if(irr > 0)
            {
				   delete[] vvv[is];
            }
			}
			delete[] vvv;

		}// I1
	}// T1

RA.delete_grid();//xiaohui add 2015-02-04
	timer::tick("Force_LCAO_k","set_EDM_k",'G');
	return;
}


complex<double> Force_LCAO_k::set_EDM_k_element(
	const complex<double> &phase,
	const bool with_energy,
	complex<double> &coef1, complex<double> &coef2,
	const double &ekb)
{
	complex<double> dm = complex<double>(0,0);
	//--------------------------------------
	// for energy density matrix
	// \sum E(i)*exp(iRk)*psi(mu)*psi(nu)
	//--------------------------------------
	if(with_energy)
	{
		dm += phase * ekb * conj(coef1) * coef2;
	}
	//--------------------------------------
	// for density matrix
	// \sum E(i)*psi(mu)*psi(nu)
	//--------------------------------------
	else
	{
		dm += phase * conj(coef1) * coef2 ;
	}
	return dm;
}


void Force_LCAO_k::cal_foverlap_k(
	const bool isforce, 
	const bool isstress, 
	matrix& foverlap, 
	matrix& soverlap)
{
	TITLE("Force_LCAO_k","cal_foverlap_k");
	timer::tick("Force_LCAO_k","cal_foverlap_k",'G');

	//--------------------------------------------
	// (1) allocate energy density matrix (nnr)
	//--------------------------------------------
	double** edm2d = new double*[NSPIN];
	for(int is=0; is<NSPIN; is++)
	{
		edm2d[is] = new double[LNNR.nnr];
		ZEROS(edm2d[is], LNNR.nnr);
	}
	bool with_energy = true;

	//--------------------------------------------	
	// calculate the energy density matrix here.
	//--------------------------------------------	
	this->set_EDM_k(edm2d, with_energy);

	//--------------------------------------------
    //summation \sum_{i,j} E(i,j)*dS(i,j)
    //BEGIN CALCULATION OF FORCE OF EACH ATOM
	//--------------------------------------------
	Vector3<double> tau1, dtau, tau2;

	Record_adj RA;
	RA.for_2d();

	int irr = 0;
	int iat = 0;
    for(int T1=0; T1<ucell.ntype; ++T1)
    {
        Atom* atom1 = &ucell.atoms[T1];
        for(int I1=0; I1<atom1->na; ++I1)
        {
			const int start1 = ucell.itiaiw2iwt(T1,I1,0);
			for (int cb = 0; cb < RA.na_each[iat]; ++cb)
			{
				const int T2 = RA.info[iat][cb][3];
				const int I2 = RA.info[iat][cb][4];
				const int start2 = ucell.itiaiw2iwt(T2, I2, 0);

				Atom* atom2 = &ucell.atoms[T2];

				for(int jj=0; jj<atom1->nw; jj++)
				{
					const int iw1_all = start1 + jj; 

					// HPSEPS
					const int mu = ParaO.trace_loc_row[iw1_all];
					if(mu<0)continue;

					for(int kk=0; kk<atom2->nw; kk++)
					{
						const int iw2_all = start2 + kk;

						// HPSEPS
						const int nu = ParaO.trace_loc_col[iw2_all];
						if(nu<0)continue;
						//==============================================================
						// here we use 'minus', but in GAMMA_ONLY_LOCAL we use 'plus',
						// both are correct because the 'DSloc_Rx' is used in 'row' (-),
						// however, the 'DSloc_x' in GAMMA is used in 'col' (+),
						// mohan update 2011-06-16
						//==============================================================
						for(int is=0; is<NSPIN; ++is)
						{
							double edm2d2 = 2.0 * edm2d[is][irr];
							if(isforce)
							{
								foverlap(iat,0) -= edm2d2 * LM.DSloc_Rx[irr];
								foverlap(iat,1) -= edm2d2 * LM.DSloc_Ry[irr];
								foverlap(iat,2) -= edm2d2 * LM.DSloc_Rz[irr];
							}
							if(isstress)
							{
								for(int ipol = 0;ipol<3;ipol++)
								{
									soverlap(0,ipol) += edm2d[is][irr] * LM.DSloc_Rx[irr] * LM.DH_r[irr * 3 + ipol];
									soverlap(1,ipol) += edm2d[is][irr] * LM.DSloc_Ry[irr] * LM.DH_r[irr * 3 + ipol];
									soverlap(2,ipol) += edm2d[is][irr] * LM.DSloc_Rz[irr] * LM.DH_r[irr * 3 + ipol];
								}
							}
						}
						++irr;
					}// end kk
				}// end jj
			}// end cb
			++iat;
		}
	}

	//-----------------
	// test the force
	//-----------------
	/*
	cout << " overlap force" << endl;
	for(int iat=0; iat<ucell.nat; ++iat)
	{
		const double fac = Ry_to_eV / 0.529177;
		cout << setw(5) << iat+1 << setw(15) << foverlap[iat][0] *fac<< setw(15) << foverlap[iat][1]*fac << 
		setw(15) << foverlap[iat][2]*fac << endl;
	}
	*/
	if(isstress){
		for(int i=0;i<3;i++)
		{
			for(int j=0;j<3;j++)
			{
				soverlap(i,j) *=  ucell.lat0 / ucell.omega;
			}
		}
	}

	if(irr!=LNNR.nnr)
	{
		OUT(ofs_running,"wrong irr",irr);
		OUT(ofs_running,"wrong LNNR.nnr",LNNR.nnr);
		WARNING_QUIT("Force_LCAO_k::cal_foverlap_k","irr!=LNNR.nnr");
	}
	
	for(int is=0; is<NSPIN; is++)
	{
		delete[] edm2d[is];
	}
	delete[] edm2d;

	RA.delete_grid();//xiaohui add 2015-02-04
	timer::tick("Force_LCAO_k","cal_foverlap_k",'G');
	return;
}

void Force_LCAO_k::cal_ftvnl_dphi_k(
	double** dm2d, 
	const bool isforce, 
	const bool isstress, 
	matrix& ftvnl_dphi, 
	matrix& stvnl_dphi)
{	
	TITLE("Force_LCAO_k","cal_ftvnl_dphi");
	timer::tick("Force_LCAO_k","cal_ftvnl_dphi",'G');
	
	// get the adjacent atom's information.

//	ofs_running << " calculate the ftvnl_dphi_k force" << endl;
	Record_adj RA;
	RA.for_2d();

	int irr = 0;
    for(int T1=0; T1<ucell.ntype; ++T1)
    {
        Atom* atom1 = &ucell.atoms[T1];
        for(int I1=0; I1<atom1->na; ++I1)
        {
			const int iat = ucell.itia2iat(T1,I1);
			const int start1 = ucell.itiaiw2iwt(T1,I1,0);
			for (int cb = 0; cb < RA.na_each[iat]; ++cb)
			{
				const int T2 = RA.info[iat][cb][3];
				const int I2 = RA.info[iat][cb][4];
				const int start2 = ucell.itiaiw2iwt(T2,I2,0);
				Atom* atom2 = &ucell.atoms[T2];

				for(int jj=0; jj<atom1->nw; ++jj)
				{
					const int iw1_all = start1 + jj; 
					const int mu = ParaO.trace_loc_row[iw1_all];
					if(mu<0)continue;
					for(int kk=0; kk<atom2->nw; ++kk)
					{
						const int iw2_all = start2 + kk;
						const int nu = ParaO.trace_loc_col[iw2_all];
						if(nu<0)continue;
						//==============================================================
						// here we use 'minus', but in GAMMA_ONLY_LOCAL we use 'plus',
						// both are correct because the 'DSloc_Rx' is used in 'row' (-),
						// however, the 'DSloc_x' is used in 'col' (+),
						// mohan update 2011-06-16
						//==============================================================
						for(int is=0; is<NSPIN; ++is)
						{
							double dm2d2 = 2.0 * dm2d[is][irr];
							if(isforce)
							{
								ftvnl_dphi(iat,0) += dm2d2 * LM.DHloc_fixedR_x[irr];
								ftvnl_dphi(iat,1) += dm2d2 * LM.DHloc_fixedR_y[irr];
								ftvnl_dphi(iat,2) += dm2d2 * LM.DHloc_fixedR_z[irr];
							}
							if(isstress)
							{
								stvnl_dphi(0,0) -= dm2d[is][irr] * LM.stvnl11[irr];
								stvnl_dphi(0,1) -= dm2d[is][irr] * LM.stvnl12[irr];
								stvnl_dphi(0,2) -= dm2d[is][irr] * LM.stvnl13[irr];
								stvnl_dphi(1,1) -= dm2d[is][irr] * LM.stvnl22[irr];
								stvnl_dphi(1,2) -= dm2d[is][irr] * LM.stvnl23[irr];
								stvnl_dphi(2,2) -= dm2d[is][irr] * LM.stvnl33[irr];
							}
						}
						++irr;
					}//end kk
				}//end jj
			}// end cb
		}
	}
	assert(irr==LNNR.nnr);
	
//	test(LM.DSloc_Rx);
//	test(dm2d[0],"dm2d");

	if(isstress){
		for(int i=0;i<3;i++)
		{
			for(int j=0;j<3;j++)
			{
				if(i<j) stvnl_dphi(j,i) = stvnl_dphi(i,j);
				stvnl_dphi(i,j) *=  ucell.lat0 / ucell.omega;
			}
		}
	}

	RA.delete_grid();//xiaohui add 2015-02-04
	timer::tick("Force_LCAO_k","cal_ftvnl_dphi",'G');
	return;
}

	
void Force_LCAO_k::test(double* mmm, const string &name)
{
	if(NPROC!=1)return;
	cout << "test!" << endl;

	int irr = 0;
	int ca = 0;

	ofs_running << " Calculate the test in Force_LCAO_k" << endl;
	Record_adj RA;
	RA.for_2d();
	
	double *test;
	test = new double[NLOCAL * NLOCAL];
	ZEROS(test, NLOCAL *NLOCAL);
	
	for(int T1=0; T1<ucell.ntype; T1++)
    {
        Atom* atom1 = &ucell.atoms[T1];
        for(int I1=0; I1<atom1->na; I1++)
        {
			const int iat = ucell.itia2iat(T1,I1);
			const int start1 = ucell.itiaiw2iwt(T1,I1,0);
			for (int cb = 0; cb < RA.na_each[ca]; cb++ )
			{
				const int T2 = RA.info[ca][cb][3];
				const int I2 = RA.info[ca][cb][4];
				Atom* atom2 = &ucell.atoms[T2];
				const int start2 = ucell.itiaiw2iwt(T2,I2,0);

				for(int jj=0; jj<atom1->nw; jj++)
				{
					const int iw1_all = start1+jj;	
					for(int kk=0; kk<atom2->nw; kk++)
					{
						const int iw2_all = start2+kk;
						assert(irr<LNNR.nnr);
						//test[iw1_all*NLOCAL+iw2_all] += LM.DHloc_fixedR_x[irr];
						test[iw1_all*NLOCAL+iw2_all] += mmm[irr];
						++irr;
					}
				}
			}
			++ca;
		}
	}
		
	cout << "\n " << name << endl;
	cout << setprecision(4);
	for(int i=0; i<NLOCAL; i++)
	{
		for(int j=0; j<NLOCAL; j++)
		{
			if( abs(test[i*NLOCAL+j]) > 1.0e-5)
			cout << setw(12) << test[i*NLOCAL+j];
			else
			cout << setw(12) << "0";
		}
		cout << endl;
	}
	delete[] test;	

RA.delete_grid();//xiaohui add 2015-02-04
	return;
}


// must consider three-center H matrix.
void Force_LCAO_k::cal_fvnl_dbeta_k(
	double** dm2d, 
	const bool isforce, 
	const bool isstress, 
	matrix& fvnl_dbeta, 
	matrix& svnl_dbeta)
{
	TITLE("Force_LCAO_k","cal_fvnl_dbeta_k");
	timer::tick("Force_LCAO_k","cal_fvnl_dbeta_k",'G');
	int iir = 0;
	Vector3<double> tau1;
	Vector3<double> tau2;
	Vector3<double> dtau;
	Vector3<double> tau0;
	Vector3<double>	dtau1;
	Vector3<double> dtau2;

	double rcut;
	double distance;

	double rcut1;
	double rcut2;
	double distance1;
	double distance2;
	
	for (int T1 = 0; T1 < ucell.ntype; ++T1)
	{
		const Atom* atom1 = &ucell.atoms[T1];

		for (int I1 =0; I1< atom1->na; ++I1)
		{
			tau1 = atom1->tau[I1];
			//GridD.Find_atom( tau1 );
			GridD.Find_atom( tau1 ,T1, I1);
			const int iat = ucell.itia2iat(T1, I1);
			const int start1 = ucell.itiaiw2iwt(T1, I1, 0);

			for (int ad2=0; ad2<GridD.getAdjacentNum()+1 ; ++ad2)
			{
				const int T2 = GridD.getType(ad2);
				const Atom* atom2 = &ucell.atoms[T2];
				const int I2 = GridD.getNatom(ad2);
				const int iat2 = ucell.itia2iat(T2, I2);
				const int start2 = ucell.itiaiw2iwt(T2, I2, 0);
				tau2 = GridD.getAdjacentTau(ad2);

				dtau = tau2 - tau1;
				distance = dtau.norm() * ucell.lat0;
				rcut = ORB.Phi[T1].getRcut() + ORB.Phi[T2].getRcut();

				// check if this a adjacent atoms.
				bool is_adj = false;
				if(distance < rcut) is_adj = true;
				else if(distance >= rcut)
				{
					for (int ad0=0; ad0 < GridD.getAdjacentNum()+1 ; ++ad0)
					{
						const int T0 = GridD.getType(ad0);
						if( ORB.nproj[T0] == 0) continue;
						const int I0 = GridD.getNatom(ad0);
						const int iat0 = ucell.itia2iat(T0, I0);
						const int start0 = ucell.itiaiw2iwt(T0, I0, 0);

						tau0 = GridD.getAdjacentTau(ad0);
						dtau1 = tau0 - tau1;
						distance1 = dtau1.norm() * ucell.lat0;
						rcut1 = ORB.Phi[T1].getRcut() + ORB.Beta[T0].get_rcut_max();

						dtau2 = tau0 - tau2;
						distance2 = dtau2.norm() * ucell.lat0;
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
					for (int j=0; j<atom1->nw; ++j)
					{
						const int iw1_all = start1 + j;
						const int mu = ParaO.trace_loc_row[iw1_all];
						if(mu < 0)continue;
						for (int k=0; k<atom2->nw; ++k)
						{
							const int iw2_all = start2 + k;
							const int nu = ParaO.trace_loc_col[iw2_all];
							if(nu < 0)continue;
							
							for (int ad0=0; ad0 < GridD.getAdjacentNum()+1 ; ++ad0)
							{
								const int T0 = GridD.getType(ad0);
								if( ORB.nproj[T0] == 0) continue;
								const int I0 = GridD.getNatom(ad0);
								const int iat0 = ucell.itia2iat(T0, I0);
								const int start0 = ucell.itiaiw2iwt(T0, I0, 0);
								tau0 = GridD.getAdjacentTau(ad0);

								dtau1 = tau0 - tau1;
								distance1 = dtau1.norm() * ucell.lat0;
								rcut1 = ORB.Phi[T1].getRcut() + ORB.Beta[T0].get_rcut_max();

								dtau2 = tau0 - tau2;
								distance2 = dtau2.norm() * ucell.lat0;
								rcut2 = ORB.Phi[T2].getRcut() + ORB.Beta[T0].get_rcut_max();

								double r0[3];
								double r1[3];
								r1[0] = ( tau1.x - tau0.x) ;
								r1[1] = ( tau1.y - tau0.y) ;
								r1[2] = ( tau1.z - tau0.z) ;
								r0[0] = ( tau2.x - tau0.x) ;
								r0[1] = ( tau2.y - tau0.y) ;
								r0[2] = ( tau2.z - tau0.z) ;

								if(distance1 < rcut1 && distance2 < rcut2)
								{
									const Atom* atom0 = &ucell.atoms[T0];
									double nlm[3]={0,0,0};

									UOT.snap_psibeta(
											nlm, 1,
											tau2,
											T2,
											atom2->iw2l[ k ], // L2
											atom2->iw2m[ k ], // m2
											atom2->iw2n[ k ], // n2
											tau1,
											T1,
											atom1->iw2l[ j ], // L1
											atom1->iw2m[ j ], // m1
											atom1->iw2n[ j ], // N1
											tau0, T0, ucell.atoms[T0].dion, NSPIN,
											ucell.atoms[T0].d_so,
											ucell.atoms[T0].non_zero_count_soc[0], // index stands for spin
											ucell.atoms[T0].index1_soc[0],
											ucell.atoms[T0].index2_soc[0],
											ucell.atoms[T0].nproj_soc
											); // mohan  add 2021-05-07

									double nlm1[3]={0,0,0};
									if(isstress)
									{
										UOT.snap_psibeta(
											nlm1, 1,
											tau1,
											T1,
											atom1->iw2l[ j ], // L1
											atom1->iw2m[ j ], // m1
											atom1->iw2n[ j ], // n1
											tau2,
											T2,
											atom2->iw2l[ k ], // L2
											atom2->iw2m[ k ], // m2
											atom2->iw2n[ k ], // N2
											tau0, T0, ucell.atoms[T0].dion, NSPIN,
											ucell.atoms[T0].d_so,
											ucell.atoms[T0].non_zero_count_soc[0], // index stands for spin
											ucell.atoms[T0].index1_soc[0],
											ucell.atoms[T0].index2_soc[0],
											ucell.atoms[T0].nproj_soc
											); // mohan  add 2021-05-07
									}
									/// only one projector for each atom force, but another projector for stress
									for(int is=0; is<NSPIN; ++is)
									{
										double dm2d2 = 2.0 * dm2d[is][iir];
										for(int jpol=0;jpol<3;jpol++)
										{
											if(isforce)
											{
												fvnl_dbeta(iat0, jpol) -= dm2d2 * nlm[jpol];
											}
											if(isstress) 
											{
												for(int ipol=0;ipol<3;ipol++)
												{
													svnl_dbeta(jpol, ipol) += dm2d[is][iir] * 
													(nlm[jpol] * r1[ipol] + nlm1[jpol] * r0[ipol]);
												}
											}
										}
									}

								}// distance
							}// ad0

							++iir;
						}// k
					}// j
				}// distance
			}// ad2
		}// I1
	}// T1

	assert( iir == LNNR.nnr );

	if(isstress)
	{
		for(int i=0;i<3;i++)
		{
			for(int j=0;j<3;j++)
			{
				svnl_dbeta(i,j) *=  ucell.lat0 / ucell.omega;
			}
		}
	}

	timer::tick("Force_LCAO_k","cal_fvnl_dbeta_k",'G');
	return;
}


// calculate the force due to < phi | Vlocal | dphi >
void Force_LCAO_k::cal_fvl_dphi_k(
	double** dm2d, 
	const bool isforce, 
	const bool isstress, 
	matrix& fvl_dphi, 
	matrix& svl_dphi)
{
	TITLE("Force_LCAO_k","cal_fvl_dphi_k");
	timer::tick("Force_LCAO_k","cal_fvl_dphi_k",'G');

	if(!isforce&&!isstress) return;
	assert(LM.DHloc_fixedR_x!=NULL);
	assert(LM.DHloc_fixedR_y!=NULL);
	assert(LM.DHloc_fixedR_z!=NULL);

	int istep = 1;

	// if Vna potential is not used.
	pot.init_pot(istep, pw.strucFac);


	for(int is=0; is<NSPIN; ++is)
	{
		CURRENT_SPIN = is;
//		ZEROS (LM.DHloc_fixedR_x, LNNR.nnr);
//		ZEROS (LM.DHloc_fixedR_y, LNNR.nnr);
//		ZEROS (LM.DHloc_fixedR_z, LNNR.nnr);
//		cout << " CURRENT_SPIN=" << CURRENT_SPIN << endl;

		for(int ir=0; ir<pw.nrxx; ir++)
		{
			pot.vr_eff1[ir] = pot.vr_eff(CURRENT_SPIN, ir);
		}

		//--------------------------------
		// Grid integration here.
		//--------------------------------
		// fvl_dphi can not be set to zero here if Vna is used
		if(isstress&&isforce) 
		{
			UHM.GK.svl_k_RealSpace(fvl_dphi,svl_dphi,pot.vr_eff1);
		}
		else if(isforce) 
		{
			UHM.GK.fvl_k_RealSpace(fvl_dphi,pot.vr_eff1);
		}
	}

	
	if(isstress){
		for(int ipol=0;ipol<3;ipol++){
			for(int jpol=0;jpol<3;jpol++){
				if(ipol < jpol) svl_dphi(jpol, ipol) = svl_dphi(ipol, jpol);
				svl_dphi(ipol, jpol) /= ucell.omega;
			}
		}
	}

	timer::tick("Force_LCAO_k","cal_fvl_dphi_k",'G');
	return;
}


