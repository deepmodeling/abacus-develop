#include "src_pw/tools.h"
#include "gint_k.h"
#include "LCAO_nnr.h"
#include "../module_ORB/ORB_read.h"
#include "grid_technique.h"
#include "src_global/ylm.h"
#include "src_pw/global.h"
#include "global_fp.h" // mohan add 2021-01-30

Gint_k::Gint_k()
{
	ik_now = 0;	
	pvpR_alloc_flag = false;
	spin_now = -1; // for a start value, must not equal 1,2 or 4.
	reduced = true;// use reduced memory for H storage.
}

Gint_k::~Gint_k()
{
	
}

void Gint_k::reset_spin(const int &spin_now_in)
{
	this->spin_now = spin_now_in;
	return;
}




void Gint_k::allocate_pvpR(void)
{
	TITLE("Gint_k","allocate_pvpR");

	if(this->pvpR_alloc_flag)
	{
		return; //Liuxh add, 20181012
		WARNING_QUIT("Gint_k::allocate_pvpR","pvpR has been allocated!");
	}

	//	reduced = NURSE; 
	//xiaohui modify 2015-05-30
	//cout << " reduced algorithm for grid integration = " << reduced << endl;

	if(this->reduced)
	{
		// the number of matrix element <phi_0 | V | phi_R> is LNNR.nnrg.
		this->pvpR_reduced = new double*[NSPIN];
		for(int is =0;is<NSPIN;is++)
		{
			this->pvpR_reduced[is] = new double[LNNR.nnrg];	
			ZEROS( pvpR_reduced[is], LNNR.nnrg);
		}

		double mem = Memory::record("allocate_pvpR", "pvpR_reduced", LNNR.nnrg * NSPIN , "double");

        if(OUT_LEVEL != "m") 
		{
			ofs_running << " Memory of pvpR : " << mem << " MB" << endl;
		}

        if( mem > 800 )
        {
            ofs_warning << " memory for pvpR = " << mem << endl;
            ofs_warning << " which is larger than 800 MB ! " << endl;
			WARNING_QUIT("Gint_k","allocate_pvpR");
		}

	}
	else
	{
		double mem = Memory::record("allocate_pvpR", "pvpR", GridT.lgd * GridT.nutot
				* GridT.lgd * GridT.nutot , "double");

		if(OUT_LEVEL != "m") 
		{
			ofs_running << " Memory of pvpR : " << mem << " MB" << endl;
		}

		if( mem > 800 )
		{
			ofs_warning << " memory for pvpR = " << mem << endl;
			ofs_warning << " which is larger than 800 MB ! " << endl;
			WARNING_QUIT("Gint_k","allocate_pvpR");
		}

		//----------------------------------------------
		// allocate the complex matrix !!
		// nutot : total number of unitcells involved.
		// this may be very large, at least 
		// 3*3*3 = 27.
		//----------------------------------------------
		const int LDIM=GridT.lgd*GridT.nutot;

		this->pvpR_pool = new double[LDIM*LDIM];
		ZEROS(pvpR_pool, LDIM*LDIM);

		this->pvpR = new double*[LDIM];
		for(int i=0; i<LDIM; i++)
		{
			pvpR[i] = &pvpR_pool[i*LDIM];
		}
	}

	this->pvpR_alloc_flag = true;
	return;
}



void Gint_k::destroy_pvpR(void)
{
	TITLE("Gint_k","destroy_pvpR");
	
	if(!pvpR_alloc_flag)
	{
		WARNING_QUIT("Gint_k::destroy_pvpR","<phi_0i | V | phi_Rj> matrix has not been allocated yet!");
	}
	
	if(this->reduced)
	{
		for(int is =0;is<NSPIN;is++) delete[] pvpR_reduced[is];
		delete[] pvpR_reduced;
	}	
	else
	{
		delete[] pvpR;
		delete[] pvpR_pool;
	}

	this->pvpR_alloc_flag = false;
	return;
}





// fold the <phi | vl |dphi(R)> * DM(R) to 
// calculate the force.
void Gint_k::folding_force(
	matrix& fvl_dphi,
	double* pvdpx, 
	double* pvdpy, 
	double* pvdpz)
{
	TITLE("Gint_k","folding_force");
	timer::tick("Gint_k","folding_force");

	//xiaohui modify 2013-12-17, test
//	assert(GridT.lgd > 0); //mohan add 2012-06-10

	// mohan add 2014-01-20
	const int lgd = GridT.lgd;

	double** ppx;
	double** ppy;
	double** ppz;

	if(GridT.lgd>0)
	{
		ppx = new double*[lgd];
	  	ppy = new double*[lgd];
	  	ppz = new double*[lgd];
	  	for(int i=0; i<lgd; i++)
	  	{
			ppx[i] = new double[lgd];
			ppy[i] = new double[lgd];
			ppz[i] = new double[lgd];
			ZEROS( ppx[i], lgd);
			ZEROS( ppy[i], lgd);
			ZEROS( ppz[i], lgd);
		}
	}
	
	Vector3<double> tau1, dtau;
	for(int T1=0; T1<ucell.ntype; ++T1)
	{
		Atom* atom1 = &ucell.atoms[T1];
		for(int I1=0; I1< atom1->na; ++I1)
		{
			const int iat = ucell.itia2iat(T1,I1);
			if(GridT.in_this_processor[iat])
			{
				assert( lgd > 0 );

				const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
				// get the start positions of elements.
				const int DM_start = LNNR.nlocstartg[iat];
				// get the coordinates of adjacent atoms.
				tau1 = atom1->tau[I1];
				//GridD.Find_atom(tau1);
				GridD.Find_atom(tau1, T1, I1);
				// search for the adjacent atoms.
				int nad = 0;
				for (int ad = 0; ad < GridD.getAdjacentNum()+1; ++ad)
				{
					// get iat2
					const int T2 = GridD.getType(ad);
					const int I2 = GridD.getNatom(ad);
					const int iat2 = ucell.itia2iat(T2, I2);
					if(GridT.in_this_processor[iat2])
					{
						Atom* atom2 = &ucell.atoms[T2];
						dtau = GridD.getAdjacentTau(ad) - tau1;
						double distance = dtau.norm() * ucell.lat0;
						double rcut = ORB.Phi[T1].getRcut() + ORB.Phi[T2].getRcut();
						if(distance < rcut)
						{
							const int start2 = ucell.itiaiw2iwt(T2, I2, 0);
							int ixxx = DM_start + LNNR.find_R2st[iat][nad];
							for(int iw=0; iw<atom1->nw; iw++)
							{
								const int iw_all = start1+iw;
								const int iw_local = GridT.trace_lo[iw_all];
								// iw1_lo
								double *vijx = ppx[iw_local];
								double *vijy = ppy[iw_local];
								double *vijz = ppz[iw_local];

								double *vRx = &pvdpx[ixxx]; //just fold R to normal matrix.
								double *vRy = &pvdpy[ixxx];
								double *vRz = &pvdpz[ixxx];

								int* iw2_lo = &GridT.trace_lo[start2];
								int* iw2_end = iw2_lo + atom2->nw;

								for(; iw2_lo<iw2_end; ++iw2_lo, ++vRx, ++vRy, ++vRz)
								{
									vijx[iw2_lo[0]] += vRx[0] ;
									vijy[iw2_lo[0]] += vRy[0] ;
									vijz[iw2_lo[0]] += vRz[0] ;
								}
								ixxx += atom2->nw;
							}
							++nad;
						}//end distance<rcut
					}
				}//end ad
			}
		}//end ia
	}//end it



	double* tmp = new double[NLOCAL*3];
	for(int i=0; i<NLOCAL; ++i)
	{
		ZEROS(tmp, 3*NLOCAL);
		const int mug = GridT.trace_lo[i];
		// if the row element is on this processor
		if(mug>=0)
		{
			//ofs_running << " i=" << i << " mug=" << mug << endl;
			for(int j=0; j<NLOCAL; ++j)
			{
				const int nug = GridT.trace_lo[j];
				// if the col element is on this processor
				if(nug>=0)
				{
	//				if(mug<nug)
	//				{
						const int index = 3*j;
						tmp[index] = ppx[mug][nug];
						tmp[index+1] = ppy[mug][nug];
						tmp[index+2] = ppz[mug][nug];
	//				}
	//				else
	//				{
					//	tmpx[j] = 0.0;
					//	tmpy[j] = 0.0;
					//	tmpz[j] = 0.0;
	//				}
				}
			}
		}
		// collect the matrix after folding.
		Parallel_Reduce::reduce_double_pool( tmp, NLOCAL*3 );
		for (int j=0; j<NLOCAL; j++)
		{
			if (!ParaO.in_this_processor(i,j))
			{
				continue;
			}
			const int iat = ucell.iwt2iat[i];
			const int index = 3*j;
			fvl_dphi(iat,0) += 2.0*tmp[index];	
			fvl_dphi(iat,1) += 2.0*tmp[index+1];	
			fvl_dphi(iat,2) += 2.0*tmp[index+2];	
		}
	}
	delete[] tmp;

	// mohan add 2014-01-20
	if(GridT.lgd > 0)
	{
		//-------------------------
		// delete the tmp matrix.
		//-------------------------
		for(int i=0; i<GridT.lgd; i++)
		{
			delete[] ppx[i];
			delete[] ppy[i];
			delete[] ppz[i];
		}
		delete[] ppx;
		delete[] ppy;
		delete[] ppz;
	}
	timer::tick("Gint_k","folding_force");
	return;
}

// fold the <phi | vl * R_beta|dphi(R_alpha)> * DM(R) to 
// calculate the stress.
void Gint_k::folding_stress(
	matrix& fvl_dphi, 
	matrix& svl_dphi,
	double* pvdpx, 
	double* pvdpy, 
	double* pvdpz,
	double* pvdp11, 
	double* pvdp22, 
	double* pvdp33,
	double* pvdp12, 
	double* pvdp13, 
	double* pvdp23)
{
	TITLE("Gint_k","folding_stress");
	timer::tick("Gint_k","folding_stress");

	const int lgd = GridT.lgd;

	double** ppx;
	double** ppy;
	double** ppz;

	double** pp11;
	double** pp22;
	double** pp33;
	double** pp12;
	double** pp13;
	double** pp23;

	double r0[3],rt[3];

	if(GridT.lgd>0)
	{
		ppx = new double*[lgd];
		ppy = new double*[lgd];
		ppz = new double*[lgd];

		pp11 = new double*[lgd];
		pp22 = new double*[lgd];
		pp33 = new double*[lgd];
		pp12 = new double*[lgd];
		pp13 = new double*[lgd];
		pp23 = new double*[lgd];

		for(int i=0; i<lgd; i++)
		{
			ppx[i] = new double[lgd];
			ppy[i] = new double[lgd];
			ppz[i] = new double[lgd];
			ZEROS( ppx[i], lgd);
			ZEROS( ppy[i], lgd);
			ZEROS( ppz[i], lgd);

			pp11[i] = new double[lgd];
			pp22[i] = new double[lgd];
			pp33[i] = new double[lgd];
			ZEROS( pp11[i], lgd);
			ZEROS( pp22[i], lgd);
			ZEROS( pp33[i], lgd);
			pp12[i] = new double[lgd];
			pp13[i] = new double[lgd];
			pp23[i] = new double[lgd];
			ZEROS( pp12[i], lgd);
			ZEROS( pp13[i], lgd);
			ZEROS( pp23[i], lgd);

		}
	}
	Vector3<double> tau1, dtau;
	for(int T1=0; T1<ucell.ntype; ++T1)
	{
		const Atom* atom1 = &ucell.atoms[T1];
                
		for (int I1 =0; I1< atom1->na; ++I1)
		{
			const int iat = ucell.itia2iat(T1,I1);
			if(GridT.in_this_processor[iat])
			{
				assert( lgd > 0 );

				const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
				// get the start positions of elements.
				const int DM_start = LNNR.nlocstartg[iat];
				// get the coordinates of adjacent atoms.
				tau1 = atom1->tau[I1];
				//GridD.Find_atom(tau1);
				GridD.Find_atom(tau1, T1, I1);
				// search for the adjacent atoms.
				int nad = 0;
				for (int ad = 0; ad < GridD.getAdjacentNum()+1; ++ad)
				{
					// get iat2
					const int T2 = GridD.getType(ad);
					const int I2 = GridD.getNatom(ad);

					const Vector3<double> tau2 = GridD.getAdjacentTau(ad);
					const int iat2 = ucell.itia2iat(T2, I2);
					if(GridT.in_this_processor[iat2])
					{
						Atom* atom2 = &ucell.atoms[T2];
						dtau = GridD.getAdjacentTau(ad) - tau1;
						double distance = dtau.norm() * ucell.lat0;
						double rcut = ORB.Phi[T1].getRcut() + ORB.Phi[T2].getRcut();
						if(distance < rcut)
						{
							const int start2 = ucell.itiaiw2iwt(T2, I2, 0);
							int ixxx = DM_start + LNNR.find_R2st[iat][nad];
							for(int iw=0; iw<atom1->nw; iw++)
							{
								const int iw_all = start1+iw;
								const int iw_local = GridT.trace_lo[iw_all];
								// iw1_lo
								double *vijx = ppx[iw_local];
								double *vijy = ppy[iw_local];
								double *vijz = ppz[iw_local];

								double *vij11 = pp11[iw_local];
								double *vij22 = pp22[iw_local];
								double *vij33 = pp33[iw_local];
								double *vij12 = pp12[iw_local];
								double *vij13 = pp13[iw_local];
								double *vij23 = pp23[iw_local];

								//just fold R to normal matrix.
								double *vRx = &pvdpx[ixxx];
								double *vRy = &pvdpy[ixxx];
								double *vRz = &pvdpz[ixxx];

								double *vR11 = &pvdp11[ixxx]; 
								double *vR22 = &pvdp22[ixxx];
								double *vR33 = &pvdp33[ixxx];
								double *vR12 = &pvdp12[ixxx]; 
								double *vR13 = &pvdp13[ixxx];
								double *vR23 = &pvdp23[ixxx];


								int* iw2_lo = &GridT.trace_lo[start2];
								int* iw2_end = iw2_lo + atom2->nw;

								for(; iw2_lo<iw2_end; ++iw2_lo, ++vRx, ++vRy, ++vRz,++vR11, 
									++vR22,++vR33,++vR12,++vR13,++vR23)
								{
									vijx[iw2_lo[0]] += vRx[0];
									vijy[iw2_lo[0]] += vRy[0];
									vijz[iw2_lo[0]] += vRz[0];

									vij11[iw2_lo[0]] += vR11[0];
									vij22[iw2_lo[0]] += vR22[0];
									vij33[iw2_lo[0]] += vR33[0];
									vij12[iw2_lo[0]] += vR12[0];
									vij13[iw2_lo[0]] += vR13[0];
									vij23[iw2_lo[0]] += vR23[0];

								}
								ixxx += atom2->nw;
							}
							++nad;
						}//end distance<rcut
					}
				}//end ad
			}
		}//end ia
	}//end it

	double* tmp = new double[NLOCAL*3];
	double* tmp1 = new double[NLOCAL*6];
	for(int i=0; i<NLOCAL; ++i)
	{
		ZEROS(tmp, 3*NLOCAL);
		ZEROS(tmp1, 6*NLOCAL);
		const int mug = GridT.trace_lo[i];
		// if the row element is on this processor
		if(mug>=0)
		{
			//ofs_running << " i=" << i << " mug=" << mug << endl;
			for(int j=0; j<NLOCAL; ++j)
			{
				const int nug = GridT.trace_lo[j];
				// if the col element is on this process or
				if(nug>=0)
				{
						const int index = 3*j;
						tmp[index] = ppx[mug][nug];
						tmp[index+1] = ppy[mug][nug];
						tmp[index+2] = ppz[mug][nug];

						const int index1 = 6*j;
						tmp1[index1] = pp11[mug][nug];
						tmp1[index1+1] = pp22[mug][nug];
						tmp1[index1+2] = pp33[mug][nug];
						tmp1[index1+3] = pp12[mug][nug];
						tmp1[index1+4] = pp13[mug][nug];
						tmp1[index1+5] = pp23[mug][nug];

				}
			}
		}
		// collect the matrix after folding.
		Parallel_Reduce::reduce_double_pool( tmp, NLOCAL*3 );
		Parallel_Reduce::reduce_double_pool( tmp1, NLOCAL*6 );
		for (int j=0; j<NLOCAL; j++)
		{
			if (!ParaO.in_this_processor(i,j))
			{
				continue;
			}
			const int iat = ucell.iwt2iat[i];
			const int index = 3*j;
			const int index1 = 6*j;
			fvl_dphi(iat,0) += 2.0*tmp[index];
			fvl_dphi(iat,1) += 2.0*tmp[index+1];
			fvl_dphi(iat,2) += 2.0*tmp[index+2];

			svl_dphi(0,0) -= 2.0*tmp1[index1];
			svl_dphi(1,1) -= 2.0*tmp1[index1+1];
			svl_dphi(2,2) -= 2.0*tmp1[index1+2];
			svl_dphi(0,1) -= 2.0*tmp1[index1+3];
			svl_dphi(0,2) -= 2.0*tmp1[index1+4];
			svl_dphi(1,2) -= 2.0*tmp1[index1+5];
		}
	}
	delete[] tmp;
	delete[] tmp1;
	if(GridT.lgd > 0)
	{
		//-------------------------
		// delete the tmp matrix.
		//-------------------------
		for(int i=0; i<GridT.lgd; i++)
		{
			delete[] ppx[i];
			delete[] ppy[i];
			delete[] ppz[i];

			delete[] pp11[i];
			delete[] pp22[i];
			delete[] pp33[i];
			delete[] pp12[i];
			delete[] pp13[i];
			delete[] pp23[i];
		}
		delete[] ppx;
		delete[] ppy;
		delete[] ppz;

		delete[] pp11;
		delete[] pp22;
		delete[] pp33;
		delete[] pp12;
		delete[] pp13;
		delete[] pp23;
	}
	timer::tick("Gint_k","folding_stress");
	return;
}

// folding the matrix for 'ik' k-point.
// H(k)=\sum{R} H(R)exp(ikR) 
void Gint_k::folding_vl_k(const int &ik)
{
	TITLE("Gint_k","folding_vl_k");
	timer::tick("Gint_k","folding_vl_k",'G');

	if(!pvpR_alloc_flag)
	{
		WARNING_QUIT("Gint_k::destroy_pvpR","pvpR hasnot been allocated yet!");
	}

	//####################### EXPLAIN #################################
	// 1. what is GridT.lgd ?
	// GridT.lgd is the number of orbitals in each processor according
	// to the division of real space FFT grid.
	// 
	// 2. why the folding of vlocal is different from folding of 
	// < phi_0i | T+Vnl | phi_Rj > ?
	// Because the (i,j) is different for T+Vnl and Vlocal
	// The first part is due to 2D division of H and S matrix,
	// The second part is due to real space division. 
	// 
	// here we construct a temporary matrix to store the
	// matrix element < phi_0 | Vlocal | phi_R >
	//#################################################################
	this->ik_now = ik;
	this->pvp = new complex<double>*[GridT.lgd];
	for(int i=0; i<GridT.lgd; i++)
	{
		this->pvp[i] = new complex<double>[GridT.lgd];
		ZEROS( this->pvp[i], GridT.lgd);
	}

	if(!reduced)
	{	
		Vector3<double> dR;
		double arg;
		complex<double> phase;
		complex<double> *pp1;
		double *pp2;
		int count;
		for(int k=0; k<GridT.nutot; k++)
		{
			const int R1x = GridT.ucell_index2x[k];
			const int R1y = GridT.ucell_index2y[k];
			const int R1z = GridT.ucell_index2z[k];

			const int dimk = GridT.lgd*k;
			for(int m=0; m<GridT.nutot; m++)
			{
				//------------------------------------------------
				// exp(k dot dR)
				// dR is the index of box in Crystal coordinates
				//------------------------------------------------
				dR.x = GridT.ucell_index2x[m] - R1x;
				dR.y = GridT.ucell_index2y[m] - R1y;
				dR.z = GridT.ucell_index2z[m] - R1z;

				arg = (kv.kvec_d[ this->ik_now ] * dR) * TWO_PI;
				phase = complex<double>(cos(arg), sin(arg));
				for(int i=0; i<GridT.lgd; i++)
				{
					pp1 = this->pvp[i];
					pp2 = this->pvpR[i+GridT.lgd*m];
					count = dimk;
					for(int j=0; j<GridT.lgd; j++)
					{
						// folding matrix
						pp1[j] += pp2[count] * phase;
						++count;
					}
				}
			}
		}
	}
	else
	{
		int lgd = 0;
		Vector3<double> tau1, dtau, dR;
		for(int T1=0; T1<ucell.ntype; ++T1)
		{
			for(int I1=0; I1<ucell.atoms[T1].na; ++I1)
			{
				// get iat
				const int iat = ucell.itia2iat(T1,I1);
				// atom in this grid piece.
				if(GridT.in_this_processor[iat])
				{
					Atom* atom1 = &ucell.atoms[T1];
					const int start1 = ucell.itiaiw2iwt(T1, I1, 0);

					// get the start positions of elements.
					const int DM_start = LNNR.nlocstartg[iat];

					// get the coordinates of adjacent atoms.
					tau1 = ucell.atoms[T1].tau[I1];
					//GridD.Find_atom(tau1);	
					GridD.Find_atom(tau1, T1, I1);	
					// search for the adjacent atoms.
					int nad = 0;


					for (int ad = 0; ad < GridD.getAdjacentNum()+1; ad++)
					{
						// get iat2
						const int T2 = GridD.getType(ad);
						const int I2 = GridD.getNatom(ad);
						const int iat2 = ucell.itia2iat(T2, I2);


						// adjacent atom is also on the grid.
						if(GridT.in_this_processor[iat2])
						{
							Atom* atom2 = &ucell.atoms[T2];
							dtau = GridD.getAdjacentTau(ad) - tau1;
							double distance = dtau.norm() * ucell.lat0;
							double rcut = ORB.Phi[T1].getRcut() + ORB.Phi[T2].getRcut();

							// for the local part, only need to calculate <phi_i | phi_j> within range
							// mohan note 2012-07-06
							if(distance < rcut)
							{
								const int start2 = ucell.itiaiw2iwt(T2, I2, 0); 

								// calculate the distance between iat1 and iat2.
								// Vector3<double> dR = GridD.getAdjacentTau(ad) - tau1;
								dR.x = GridD.getBox(ad).x;
								dR.y = GridD.getBox(ad).y;
								dR.z = GridD.getBox(ad).z;

								// calculate the phase factor exp(ikR).
								const double arg = (kv.kvec_d[ this->ik_now ] * dR) * TWO_PI;
								complex<double> phase = complex<double>(cos(arg), sin(arg));
								int ixxx = DM_start + LNNR.find_R2st[iat][nad];
								for(int iw=0; iw<atom1->nw; iw++)
								{
									// iw1_lo
									complex<double> *vij = this->pvp[GridT.trace_lo[start1+iw]];


									int* iw2_lo = &GridT.trace_lo[start2];
									int* iw2_end = iw2_lo + atom2->nw;

									// get the <phi | V | phi>(R) Hamiltonian.
									double *vijR = &pvpR_reduced[0][ixxx];
									// complex<double> *vijR_soc = &pvpR_reduced_soc[ixxx];
									for(; iw2_lo<iw2_end; ++iw2_lo, ++vijR)
									{
										vij[iw2_lo[0]] += vijR[0] * phase; 
									}
									ixxx += atom2->nw;
									++lgd;
								}

								++nad;
							}// end distane<rcut
						}
					}// end ad
				}
			}// end ia
		}// end it

		//------------------
		// To test the pvpR
		//------------------
/*
		for(int i=0; i<LNNR.nlocdimg[0]; i++)
		{
			const int DM_start = LNNR.nlocstartg[0];
			const int j = i + DM_start;
			if( abs(pvpR_reduced[j]) > 1.0e-5  )
			{
//				cout << " pvpR_reduced[" << i <<"] = " << pvpR_reduced[j] << endl;
			}
		}
*/

	}

	//----------------------
	// Print the pvp matrix
	//----------------------
/*
	cout << " pvp matrix:" << endl;
	for(int i=0; i<GridT.lgd; i++)
	{
		for(int j=0; j<GridT.lgd; j++)
		{
			cout << setw(15) << pvp[i][j].real();
		}
		cout << endl;
	}
	*/

	// Distribution of data.
	timer::tick("Gint_k","Distri",'G');
	complex<double>* tmp;
	for (int i=0; i<NLOCAL; i++)
	{
		tmp = new complex<double>[NLOCAL];
		ZEROS(tmp, NLOCAL);
		const int mug = GridT.trace_lo[i];
		// if the row element is on this processor.
		if (mug >= 0)
		{
			for (int j=0; j<NLOCAL; j++)
			{
				const int nug = GridT.trace_lo[j];
				// if the col element is on this processor.
				if (nug >=0)
				{
					if (mug <= nug)
					{
						// pvp is symmetric, only half is calculated.
						tmp[j] = this->pvp[mug][nug];
					}
					else
					{
						// need to get elements from the other half.
						// I have question on this! 2011-02-22
						tmp[j] = conj(this->pvp[nug][mug]);
					}
				}
			}
		}
		// collect the matrix after folding.
		Parallel_Reduce::reduce_complex_double_pool( tmp, NLOCAL );

		//-----------------------------------------------------
		// NOW! Redistribute the Hamiltonian matrix elements
		// according to the HPSEPS's 2D distribution methods.
		//-----------------------------------------------------
		for (int j=0; j<NLOCAL; j++)
		{
			if (!ParaO.in_this_processor(i,j))
			{
				continue;
			}
			// set the matrix value.
			LM.set_HSk(i,j,tmp[j],'L');
		}
		delete[] tmp;
	}

	// delete the tmp matrix.
	for(int i=0; i<GridT.lgd; i++)
	{
		delete[] pvp[i];
	}
	delete[] pvp;
	timer::tick("Gint_k","Distri",'G');

	timer::tick("Gint_k","folding_vl_k",'G');
	return;
}

// folding the matrix for 'ik' k-point.
// H(k)=\sum{R} H(R)exp(ikR) 
// for non-collinear case  
void Gint_k::folding_vl_k_nc(const int &ik)
{
	TITLE("Gint_k","folding_vl_k_nc");
	timer::tick("Gint_k","folding_vl_k_nc",'G');

	if(!pvpR_alloc_flag)
	{
		WARNING_QUIT("Gint_k::destroy_pvpR","pvpR hasnot been allocated yet!");
	}

	this->ik_now = ik;
//	complex<double>** pvp_nc[4];
	for(int spin=0;spin<4;spin++)
	{
		pvp_nc[spin] = new complex<double>*[GridT.lgd];
		for(int i=0; i<GridT.lgd; i++)
		{
			pvp_nc[spin][i] = new complex<double>[GridT.lgd];
			ZEROS( this->pvp_nc[spin][i], GridT.lgd);
		}
	}

	if(!reduced)
	{	
		Vector3<double> dR;
		double arg;
		complex<double> phase;
		complex<double> *pp1;
		double *pp2;
		int count;
		for(int k=0; k<GridT.nutot; k++)
		{
			const int R1x = GridT.ucell_index2x[k];
			const int R1y = GridT.ucell_index2y[k];
			const int R1z = GridT.ucell_index2z[k];

			const int dimk = GridT.lgd*k;
			for(int m=0; m<GridT.nutot; m++)
			{
				//------------------------------------------------
				// exp(k dot dR)
				// dR is the index of box in Crystal coordinates
				//------------------------------------------------
				dR.x = GridT.ucell_index2x[m] - R1x;
				dR.y = GridT.ucell_index2y[m] - R1y;
				dR.z = GridT.ucell_index2z[m] - R1z;

				arg = (kv.kvec_d[ this->ik_now ] * dR) * TWO_PI;
				phase = complex<double>(cos(arg), sin(arg));
				for(int i=0; i<GridT.lgd; i++)
				{
					pp1 = this->pvp_nc[0][i];
					pp2 = this->pvpR[i+GridT.lgd*m];
					count = dimk;
					for(int j=0; j<GridT.lgd; j++)
					{
						// folding matrix
						pp1[j] += pp2[count] * phase;
						++count;
					}
				}
			}
		}
	}
	else
	{
		int lgd = 0;
		Vector3<double> tau1, dtau, dR;
		for(int T1=0; T1<ucell.ntype; ++T1)
		{
			for(int I1=0; I1<ucell.atoms[T1].na; ++I1)
			{
				// get iat
				const int iat = ucell.itia2iat(T1,I1);
				// atom in this grid piece.
				if(GridT.in_this_processor[iat])
				{
					Atom* atom1 = &ucell.atoms[T1];
					const int start1 = ucell.itiaiw2iwt(T1, I1, 0);

					// get the start positions of elements.
					const int DM_start = LNNR.nlocstartg[iat];

					// get the coordinates of adjacent atoms.
					tau1 = ucell.atoms[T1].tau[I1];
					//GridD.Find_atom(tau1);	
					GridD.Find_atom(tau1, T1, I1);	
					// search for the adjacent atoms.
					int nad = 0;


					for (int ad = 0; ad < GridD.getAdjacentNum()+1; ad++)
					{
						// get iat2
						const int T2 = GridD.getType(ad);
						const int I2 = GridD.getNatom(ad);
						const int iat2 = ucell.itia2iat(T2, I2);


						// adjacent atom is also on the grid.
						if(GridT.in_this_processor[iat2])
						{
							Atom* atom2 = &ucell.atoms[T2];
							dtau = GridD.getAdjacentTau(ad) - tau1;
							double distance = dtau.norm() * ucell.lat0;
							double rcut = ORB.Phi[T1].getRcut() + ORB.Phi[T2].getRcut();

							// for the local part, only need to calculate <phi_i | phi_j> within range
							// mohan note 2012-07-06
							if(distance < rcut)
							{
								const int start2 = ucell.itiaiw2iwt(T2, I2, 0);

								// calculate the distance between iat1 and iat2.
								// Vector3<double> dR = GridD.getAdjacentTau(ad) - tau1;
								dR.x = GridD.getBox(ad).x;
								dR.y = GridD.getBox(ad).y;
								dR.z = GridD.getBox(ad).z;

								// calculate the phase factor exp(ikR).
								const double arg = (kv.kvec_d[ this->ik_now ] * dR) * TWO_PI;
								complex<double> phase = complex<double>(cos(arg), sin(arg));
								int ixxx = DM_start + LNNR.find_R2st[iat][nad];
								for(int iw=0; iw<atom1->nw; iw++)
								{
									// iw1_lo
									complex<double> *vij[4];
									for(int spin=0;spin<4;spin++)
										vij[spin] = this->pvp_nc[spin][GridT.trace_lo[start1]/NPOL + iw];


									int iw2_lo = GridT.trace_lo[start2]/NPOL;
									int iw2_end = iw2_lo + atom2->nw;

									// get the <phi | V | phi>(R) Hamiltonian.
//									complex<double> *vijR_soc = &pvpR_reduced_soc[ixxx];
									double *vijR[4];
									for(int spin = 0;spin<4;spin++) 
									{
										vijR[spin] = &pvpR_reduced[spin][ixxx];
									}
									for(; iw2_lo<iw2_end; ++iw2_lo, ++vijR[0], ++vijR[1],++vijR[2],++vijR[3])
									{
										for(int spin =0;spin<4;spin++)
										{
											vij[spin][iw2_lo] += vijR[spin][0] * phase; 
										}
//										else vij[iw2_lo[0]] += vijR_soc[0] * phase;
									}
									ixxx += atom2->nw;
									++lgd;
								}

								++nad;
							}// end distane<rcut
						}
					}// end ad
				}
			}// end ia
		}// end it

		//------------------
		// To test the pvpR
		//------------------
/*
		for(int i=0; i<LNNR.nlocdimg[0]; i++)
		{
			const int DM_start = LNNR.nlocstartg[0];
			const int j = i + DM_start;
			if( abs(pvpR_reduced[j]) > 1.0e-5  )
			{
//				cout << " pvpR_reduced[" << i <<"] = " << pvpR_reduced[j] << endl;
			}
		}
*/

	}

	//----------------------
	// Print the pvp matrix
	//----------------------
/*
	cout << " pvp matrix:" << endl;
	for(int i=0; i<GridT.lgd; i++)
	{
		for(int j=0; j<GridT.lgd; j++)
		{
			cout << setw(15) << pvp[i][j].real();
		}
		cout << endl;
	}
	*/

	// Distribution of data.
	timer::tick("Gint_k","Distri",'G');
	complex<double>* tmp;
	for (int i=0; i<NLOCAL; i++)
	{
		tmp = new complex<double>[NLOCAL];
		ZEROS(tmp, NLOCAL);
		const int mug = GridT.trace_lo[i];
		const int mug0 = mug/NPOL;
		// if the row element is on this processor.
		if (mug >= 0)
		{
			for (int j=0; j<NLOCAL; j++)
			{
				const int nug = GridT.trace_lo[j];
				const int nug0 = nug/NPOL;
				// if the col element is on this processor.
				if (nug >=0)
				{
					if (mug <= nug)
					{
						// pvp is symmetric, only half is calculated.
						//int spin=0;
						if(i%2==0&&j%2==0)
						{
							//spin = 0;
							tmp[j] = this->pvp_nc[0][mug0][nug0]+this->pvp_nc[3][mug0][nug0];
						}	
						else if(i%2==1&&j%2==1)
						{
							//spin = 3;
							tmp[j] = this->pvp_nc[0][mug0][nug0]-this->pvp_nc[3][mug0][nug0];
						}
						else if(i%2==0&&j%2==1)
						{
							// spin = 1;
							if(!DOMAG) tmp[j] = 0;
							else tmp[j] = this->pvp_nc[1][mug0][nug0] - complex<double>(0.0,1.0) * this->pvp_nc[2][mug0][nug0];
						}
						else if(i%2==1&&j%2==0) 
						{
							//spin = 2;
							if(!DOMAG) tmp[j] = 0;
							else tmp[j] = this->pvp_nc[1][mug0][nug0] + complex<double>(0.0,1.0) * this->pvp_nc[2][mug0][nug0];
						}
						else
						{
							WARNING_QUIT("Gint_k::folding_vl_k_nc","index is wrong!");
						}
						//tmp[j] = this->pvp[spin][mug][nug];
					}
					else
					{
						// need to get elements from the other half.
						// I have question on this! 2011-02-22
						if(i%2==0&&j%2==0)
						{
							//spin = 0;
							tmp[j] = conj(this->pvp_nc[0][nug0][mug0]+this->pvp_nc[3][nug0][mug0]);
						}	
						else if(i%2==1&&j%2==1)
						{
							//spin = 3;
							tmp[j] = conj(this->pvp_nc[0][nug0][mug0]-this->pvp_nc[3][nug0][mug0]);
						}
						else if(i%2==1&&j%2==0)
						{
							// spin = 1;
							if(!DOMAG) tmp[j] = 0;
							else tmp[j] = conj(this->pvp_nc[1][nug0][mug0] - complex<double>(0.0,1.0) * this->pvp_nc[2][nug0][mug0]);
						}
						else if(i%2==0&&j%2==1) 
						{
							//spin = 2;
							if(!DOMAG) tmp[j] = 0;
							else tmp[j] = conj(this->pvp_nc[1][nug0][mug0] + complex<double>(0.0,1.0) * this->pvp_nc[2][nug0][mug0]);
						}
						else
						{
							WARNING_QUIT("Gint_k::folding_vl_k_nc","index is wrong!");
						}
						//tmp[j] = conj(this->pvp[spin][nug][mug]);
						
					}
				}
			}
		}
		// collect the matrix after folding.
		Parallel_Reduce::reduce_complex_double_pool( tmp, NLOCAL );

		//-----------------------------------------------------
		// NOW! Redistribute the Hamiltonian matrix elements
		// according to the HPSEPS's 2D distribution methods.
		//-----------------------------------------------------
		for (int j=0; j<NLOCAL; j++)
		{
			if (!ParaO.in_this_processor(i,j))
			{
				continue;
			}
			// set the matrix value.
			LM.set_HSk(i,j,tmp[j],'L');
		}
		delete[] tmp;
	}

	// delete the tmp matrix.
	for(int spin =0;spin<4;spin++)
	{
		for(int i=0; i<GridT.lgd; i++)
		{
			delete[] this->pvp_nc[spin][i];
		}
		delete[] this->pvp_nc[spin];
	}
	timer::tick("Gint_k","Distri",'G');

	timer::tick("Gint_k","folding_vl_k_nc",'G');
	return;
}

void Gint_k::set_ijk_atom(
	const int &grid_index, 
	const int &size,
	double*** psir_ylm, 
	double*** dr, 
	bool** cal_flag, 
	double** distance,
	const double &delta_r)
{
	const Numerical_Orbital_Lm* pointer;
	double mt[3];
	for (int id=0; id<size; id++)
	{
		// (2.1) get the atom type and atom index.
		const int mcell_index = GridT.bcell_start[grid_index] + id;	
		const int imcell = GridT.which_bigcell[mcell_index];
		const int iat = GridT.which_atom[mcell_index];
		const int it = ucell.iat2it[ iat ];
		const int ia = ucell.iat2ia[ iat ];

		// (2.2) get the distance between the grid and the atom.
		mt[0] = GridT.meshball_positions[imcell][0] - GridT.tau_in_bigcell[iat][0];
		mt[1] = GridT.meshball_positions[imcell][1] - GridT.tau_in_bigcell[iat][1];
		mt[2] = GridT.meshball_positions[imcell][2] - GridT.tau_in_bigcell[iat][2];

		for(int ib=0; ib<pw.bxyz; ib++)
		{
			// meshcell_pos: z is the fastest
			dr[ib][id][0] = GridT.meshcell_pos[ib][0] + mt[0];
			dr[ib][id][1] = GridT.meshcell_pos[ib][1] + mt[1];
			dr[ib][id][2] = GridT.meshcell_pos[ib][2] + mt[2];

			distance[ib][id] = std::sqrt(dr[ib][id][0]*dr[ib][id][0] 
			+ dr[ib][id][1]*dr[ib][id][1] + dr[ib][id][2]*dr[ib][id][2]);

			//if(distance[ib][id] <= ORB.Phi[it].getRcut())
			// mohan reset this 2012-06-27
			// if(distance[ib][id] < ORB.Phi[it].getRcut() )
			// mohan reset this 2013-07-02 in Princeton
			// we should make absolutely sure that the distance is smaller than ORB.Phi[it].getRcut
			// this should be consistant with LCAO_nnr::cal_nnrg function 
			// typical example : 7 Bohr cutoff Si orbital in 14 Bohr length of cell.
			// distance = 7.0000000000000000
			// ORB.Phi[it].getRcut = 7.0000000000000008
			if(distance[ib][id] < (ORB.Phi[it].getRcut() - 1.0e-15) )
			{
				cal_flag[ib][id]=true;
			}
			else
			{
				cal_flag[ib][id]=false;
				continue;
			}

			if (distance[ib][id] < 1.0E-9) distance[ib][id] += 1.0E-9;

			std::vector<double> ylma;
			Ylm::sph_harm ( ucell.atoms[it].nwl,
					dr[ib][id][0] / distance[ib][id],
					dr[ib][id][1] / distance[ib][id],
					dr[ib][id][2] / distance[ib][id],
					ylma);

			// these parameters are about interpolation
			// because once we know the distance from atom to grid point,
			// we can get the parameters we need to do interpolation and
			// store them first!! these can save a lot of effort.
			const double position = distance[ib][id] / delta_r;

			const int ip = static_cast<int>(position);

			const double dx = position - ip;
			const double dx2 = dx * dx;
			const double dx3 = dx2 * dx;

			const double c3 = 3.0*dx2-2.0*dx3;
			const double c1 = 1.0-c3;
			const double c2 = (dx-2.0*dx2+dx3)*delta_r;
			const double c4 = (dx3-dx2)*delta_r;

			Atom* atom1 = &ucell.atoms[it];
			double tmp=0.0;//mohan fix bug 2011-05-04
			for (int iw=0; iw< atom1->nw; iw++)
			{
				if ( atom1->iw2_new[iw] )
				{
					pointer = &ORB.Phi[it].PhiLN(
							atom1->iw2l[iw],
							atom1->iw2n[iw]);

					// Efficient!! to get the orbital value at this point.
					tmp = c1*pointer->psi_uniform[ip] 
						+ c2*pointer->dpsi_uniform[ip]
						+ c3*pointer->psi_uniform[ip+1] 
						+ c4*pointer->dpsi_uniform[ip+1];
				}
				psir_ylm[ib][id][iw] = tmp * ylma[atom1->iw2_ylm[iw]];
			}// end iw.
		}//end ib
	}// int id
	return;
}




void Gint_k::allocate_pvpR_tr(void)
{
    TITLE("Gint_k","allocate_pvpR_tr");

    int R_x = GridD.getCellX();
    int R_y = GridD.getCellY();
    int R_z = GridD.getCellZ();

//cout<<"R_x: "<<R_x<<endl;
//cout<<"R_y: "<<R_y<<endl;
//cout<<"R_z: "<<R_z<<endl;
//cout<<"GridT.lgd: "<<GridT.lgd<<endl;
    if(NSPIN!=4)
    {
        pvpR_tr = new double****[R_x];
        for(int ix=0; ix<R_x; ix++)
        {
            pvpR_tr[ix] = new double***[R_y];
            for(int iy=0; iy<R_y; iy++)
            {
                pvpR_tr[ix][iy] = new double**[R_z];
                for(int iz=0; iz<R_z; iz++)
                {
                    pvpR_tr[ix][iy][iz] = new double*[GridT.lgd];
                    for(int iw=0; iw<GridT.lgd; iw++)
                    {
                        pvpR_tr[ix][iy][iz][iw] = new double[GridT.lgd];
//do    uble mem = Memory::record("allocate_pvpR_tr", "pvpR_tr[ix][iy][iz][iw]", GridT.lgd , "double");
//co    ut<<" Memory of pvpR_tr[ix][iy][iz][iw]: "<<mem<<" MB"<<endl;
                        ZEROS(pvpR_tr[ix][iy][iz][iw], GridT.lgd);
                    }
                }
            }
        }
    }
    else
    {
        pvpR_tr_soc = new complex<double>****[R_x];
        for(int ix=0; ix<R_x; ix++)
        {
            pvpR_tr_soc[ix] = new complex<double>***[R_y];
            for(int iy=0; iy<R_y; iy++)
            {
                pvpR_tr_soc[ix][iy] = new complex<double>**[R_z];
                for(int iz=0; iz<R_z; iz++)
                {
                    pvpR_tr_soc[ix][iy][iz] = new complex<double>*[GridT.lgd];
                    for(int iw=0; iw<GridT.lgd; iw++)
                    {
                        pvpR_tr_soc[ix][iy][iz][iw] = new complex<double>[GridT.lgd];
//do    uble mem = Memory::record("allocate_pvpR_tr", "pvpR_tr[ix][iy][iz][iw]", GridT.lgd , "double");
//co    ut<<" Memory of pvpR_tr[ix][iy][iz][iw]: "<<mem<<" MB"<<endl;
                        ZEROS(pvpR_tr_soc[ix][iy][iz][iw], GridT.lgd);
                    }
                }
            }
        }
    }

    return;
}

void Gint_k::destroy_pvpR_tr(void)
{
    TITLE("Gint_k","destroy_pvpR_tr");

    int R_x = GridD.getCellX();
    int R_y = GridD.getCellY();
    int R_z = GridD.getCellZ();

    if(NSPIN!=4)
    {
        for(int ix=0; ix<R_x; ix++)
        {
            for(int iy=0; iy<R_y; iy++)
            {
                for(int iz=0; iz<R_z; iz++)
                {
                    for(int iw=0; iw<GridT.lgd; GridT.lgd++)
                    {
                        delete[] pvpR_tr[ix][iy][iz][iw];
                    }
                    delete[] pvpR_tr[ix][iy][iz];
                }
                delete[] pvpR_tr[ix][iy];
            }
            delete[] pvpR_tr[ix];
        }
        delete[] pvpR_tr;
    }
    else
    {
        for(int ix=0; ix<R_x; ix++)
        {
            for(int iy=0; iy<R_y; iy++)
            {
                for(int iz=0; iz<R_z; iz++)
                {
                    for(int iw=0; iw<GridT.lgd; GridT.lgd++)
                    {
                        delete[] pvpR_tr_soc[ix][iy][iz][iw];
                    }
                    delete[] pvpR_tr_soc[ix][iy][iz];
                }
                delete[] pvpR_tr_soc[ix][iy];
            }
            delete[] pvpR_tr_soc[ix];
        }
        delete[] pvpR_tr_soc;
    }

    return;
}

void Gint_k::distribute_pvpR_tr(void)
{
    TITLE("Gint_k","distribute_pvpR_tr");

    int R_x = GridD.getCellX();
    int R_y = GridD.getCellY();
    int R_z = GridD.getCellZ();

    for(int ix=0; ix<R_x; ix++)
    {
        for(int iy=0; iy<R_y; iy++)
        {
            for(int iz=0; iz<R_z; iz++)
            {
                double* tmp;
                complex<double>* tmp_soc;
                for(int i=0; i<NLOCAL; i++)
                {
                    if(NSPIN!=4)
                    {
                        tmp = new double[NLOCAL];
                        ZEROS(tmp, NLOCAL);
                    }
                    else
                    {
                        tmp_soc = new complex<double>[NLOCAL];
                        ZEROS(tmp_soc, NLOCAL);
                    }

                    const int mug = GridT.trace_lo[i];
//cout<<"mug: "<<mug<<endl;
                    if(mug >= 0)
                    {
                        for(int j=0; j<NLOCAL; j++)
                        {
                            const int nug = GridT.trace_lo[j];
//cout<<"nug: "<<nug<<endl;
                            if(nug >= 0)
                            {
                                //if(mug <= nug)
                                //{
//cout<<"ix: "<<ix<<endl;
//cout<<"iy: "<<iy<<endl;
//cout<<"iz: "<<iz<<endl;
//cout<<"mug: "<<mug<<endl;
//cout<<"nug: "<<nug<<endl;
//cout<<"pvpR_tr: "<<pvpR_tr[ix][iy][iz][mug][nug]<<endl;
                                    if(NSPIN!=4) tmp[j] = pvpR_tr[ix][iy][iz][mug][nug];
                                    else tmp_soc[j] = pvpR_tr_soc[ix][iy][iz][mug][nug];
//cout<<"tmp["<<j<<"]: "<<tmp[j]<<endl;
                                //}
                                //else
                                //{
                                    //tmp[j] = pvpR_tr[ix][iy][iz][nug][mug];
//cout<<"tmp["<<j<<"]: "<<tmp[j]<<endl;
                                //}
                            }
                        }
                    }
                    // collect the matrix after folding.
                    if(NSPIN!=4) Parallel_Reduce::reduce_double_pool( tmp, NLOCAL );
                    else Parallel_Reduce::reduce_complex_double_pool( tmp_soc, NLOCAL );
                    for(int j=0; j<NLOCAL; j++)
                    {
                        if(!ParaO.in_this_processor(i,j))
                        {
                            continue;
                        }
                        else
                        {
                            //LM.set_HSk(i,j,tmp[j],'L');
                            if(NSPIN!=4) LM.set_HR_tr(ix,iy,iz,i,j,tmp[j]);
                            else LM.set_HR_tr_soc(ix,iy,iz,i,j,tmp_soc[j]);
                        }
                    }
                    if(NSPIN!=4) delete[] tmp;
                    else delete[] tmp_soc;
                }
            }
        }
    }

    return;
}


void Gint_k::cal_vlocal_R(const int current_spin)
{
    TITLE("Gint_k","cal_vlocal_R");

    allocate_pvpR_tr();

    int lgd = 0;

    double R_minX = GridD.getD_minX();
    double R_minY = GridD.getD_minY();
    double R_minZ = GridD.getD_minZ();

    int R_x;
    int R_y;
    int R_z;

    Vector3<double> tau1, dtau, dR;
    for(int T1=0; T1<ucell.ntype; ++T1)
    {
        for(int I1=0; I1<ucell.atoms[T1].na; ++I1)
        {
            const int iat = ucell.itia2iat(T1,I1);
            if(GridT.in_this_processor[iat])
            {
                Atom* atom1 = &ucell.atoms[T1];
                const int start1 = ucell.itiaiw2iwt(T1, I1, 0);

                const int DM_start = LNNR.nlocstartg[iat];
                tau1 = ucell.atoms[T1].tau[I1];
                //GridD.Find_atom(tau1);        
                GridD.Find_atom(tau1, T1, I1);
                int nad2 = 0;

                for(int ad = 0; ad < GridD.getAdjacentNum()+1; ad++)
                {
                    const int T2 = GridD.getType(ad);
                    const int I2 = GridD.getNatom(ad);
                    const int iat2 = ucell.itia2iat(T2, I2);

                    if(GridT.in_this_processor[iat2])
                    {
                        Atom* atom2 = &ucell.atoms[T2];
                        dtau = GridD.getAdjacentTau(ad) - tau1;
                        double distance = dtau.norm() * ucell.lat0;
                        double rcut = ORB.Phi[T1].getRcut() + ORB.Phi[T2].getRcut();

                        if(distance < rcut)
                        {
                            const int start2 = ucell.itiaiw2iwt(T2, I2, 0);

                            dR.x = GridD.getBox(ad).x;
                            dR.y = GridD.getBox(ad).y;
                            dR.z = GridD.getBox(ad).z;

                            R_x = (int) (dR.x -R_minX);
                            R_y = (int) (dR.y -R_minY);
                            R_z = (int) (dR.z -R_minZ);

                            int ixxx = DM_start + LNNR.find_R2st[iat][nad2];
                            for(int iw=0; iw<atom1->nw * NPOL; iw++)
                            {
                                int* iw2_lo = &GridT.trace_lo[start2];
                                int* iw2_end = iw2_lo + atom2->nw;

                                //double *vijR = &pvpR_reduced[ixxx];
                                double *vijR = &pvpR_reduced[current_spin][ixxx];
								for(int iw2=0;iw2<atom2->nw * NPOL; iw2++)
								{
                                    double *HlocR;
                                    complex<double> *HlocR_soc;
                                    if(NSPIN!=4) HlocR = &pvpR_tr[R_x][R_y][R_z][GridT.trace_lo[start1+iw]][GridT.trace_lo[start2+iw2]];
									else    HlocR_soc = &pvpR_tr_soc[R_x][R_y][R_z][GridT.trace_lo[start1+iw]][GridT.trace_lo[start2+iw2]];
									const int nw = atom2->nw;
									const int mug0 = iw/NPOL;
									const int nug0 = iw2/NPOL;
									const int iw_nowg = ixxx + mug0*nw + nug0;
									const int iw_nowg1 = ixxx + nug0*nw + mug0;
									if(NSPIN==4)
									{
											
											// if the col element is on this processor.
											
												// pvp is symmetric, only half is calculated.
												//int spin=0;
												if(iw%2==0&&iw2%2==0)
												{
													//spin = 0;
													HlocR_soc[0] = complex<double>(1.0,0.0) * pvpR_reduced[0][iw_nowg] + complex<double>(1.0,0.0) * pvpR_reduced[3][iw_nowg];
												}	
												else if(iw%2==1&&iw2%2==1)
												{
													//spin = 3;
													HlocR_soc[0] = complex<double>(1.0,0.0) * pvpR_reduced[0][iw_nowg] - complex<double>(1.0,0.0) * pvpR_reduced[3][iw_nowg];
												}
												else if(iw%2==0&&iw2%2==1)
												{
													// spin = 1;
													if(!DOMAG) HlocR_soc[0] = complex<double>(0.0,0.0);
													else HlocR_soc[0] = pvpR_reduced[1][iw_nowg] - complex<double>(0.0,1.0) * pvpR_reduced[2][iw_nowg];
												}	
												else if(iw%2==1&&iw2%2==0) 
												{
													//spin = 2;
													if(!DOMAG) HlocR_soc[0] = complex<double>(0.0,0.0);
													else HlocR_soc[0] = pvpR_reduced[1][iw_nowg] + complex<double>(0.0,1.0) * pvpR_reduced[2][iw_nowg];
												}
												else
												{
													WARNING_QUIT("Gint_k::folding_vl_k_nc","index is wrong!");
												}
									}//endif NC
									else
									{
										HlocR[0] = pvpR_reduced[current_spin][iw_nowg];
                                            //pvpR_tr[R_x][R_y][R_z][GridT.trace_lo[start1+iw]][iw2_lo[0]] = vijR[0];
									}//endif normal
								}
//                               for(; iw2_lo<iw2_end; ++iw2_lo, ++vijR)
//                               {
//                                   pvpR_tr[R_x][R_y][R_z][GridT.trace_lo[start1+iw]][iw2_lo[0]] = vijR[0];
//                               }
//                               ixxx += atom2->nw;
                                ++lgd;
                            }
                            ++nad2;
                        }
                    }
                }
            }
        }
    }

    return;
}

