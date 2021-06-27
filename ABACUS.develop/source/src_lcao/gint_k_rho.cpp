#include "src_pw/tools.h"
#include "gint_k.h"
#include "LCAO_nnr.h"
#include "../module_ORB/ORB_read.h"
#include "grid_technique.h"
#include "src_global/ylm.h"
#include "src_pw/global.h"
#include "module_base/blas_connector.h"
#include "global_fp.h" // mohan add 2021-01-30

inline void setVindex(const int ncyz, const int ibx, const int jby, const int kbz, int* vindex)
{				
	int bindex = 0;
	// z is the fastest, 
	for(int ii=0; ii<pw.bx; ii++)
	{
		const int ipart = (ibx + ii) * ncyz + kbz;
		for(int jj=0; jj<pw.by; jj++)
		{
			const int jpart = (jby + jj) * pw.nczp + ipart;
			for(int kk=0; kk<pw.bz; kk++)
			{
				vindex[bindex] = kk + jpart; 
				++bindex;
			}
		}
	}
}


inline void cal_psir_ylm(int size, int grid_index, double delta_r,
						int* at, int* uc, int* block_index, int* block_iw, int* block_size, 
						bool** cal_flag, double** psir_ylm)
{
	const Numerical_Orbital_Lm *pointer;
	double mt[3];
	double dr[3];
	double distance;
	block_index[0]=0;
	for (int id=0; id<size; id++)
	{
		// there are two parameters we want to know here:
		// in which bigcell of the meshball the atom in?
		// what's the cartesian coordinate of the bigcell?
		const int mcell_index=GridT.bcell_start[grid_index] + id;
		const int imcell=GridT.which_bigcell[mcell_index];

		const int iat=GridT.which_atom[mcell_index];
		at[id]=iat;
		uc[id]=GridT.which_unitcell[mcell_index];
		
		const int it=ucell.iat2it[iat];
		const int ia=ucell.iat2ia[iat];
		const int start=ucell.itiaiw2iwt(it, ia, 0);
		block_iw[id]=GridT.trace_lo[start]/NPOL;
		Atom* atom=&ucell.atoms[it];
		block_size[id]=atom->nw;
		block_index[id+1]=block_index[id]+atom->nw;
		// meshball_positions should be the bigcell position in meshball
		// to the center of meshball.
		// calculated in cartesian coordinates
		// the vector from the grid which is now being operated to the atom position.
		// in meshball language, is the vector from imcell to the center cel, plus
		// tau_in_bigcell.
		mt[0]=GridT.meshball_positions[imcell][0] - GridT.tau_in_bigcell[iat][0];
		mt[1]=GridT.meshball_positions[imcell][1] - GridT.tau_in_bigcell[iat][1];
		mt[2]=GridT.meshball_positions[imcell][2] - GridT.tau_in_bigcell[iat][2];

		for(int ib=0; ib<pw.bxyz; ib++)
		{
			double *p=&psir_ylm[ib][block_index[id]];
			// meshcell_pos: z is the fastest
			dr[0]=GridT.meshcell_pos[ib][0] + mt[0]; 
			dr[1]=GridT.meshcell_pos[ib][1] + mt[1]; 
			dr[2]=GridT.meshcell_pos[ib][2] + mt[2]; 	

			distance=std::sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
			if(distance > (ORB.Phi[it].getRcut() - 1.0e-15)) 
			{
				cal_flag[ib][id]=false;
				ZEROS(p, block_size[id]);
				continue;
			}

			cal_flag[ib][id]=true;
			
			std::vector<double> ylma;
			//if(distance[id] > GridT.orbital_rmax) continue;
			//	Ylm::get_ylm_real(this->nnn[it], this->dr[id], ylma);
			if (distance < 1.0E-9) distance += 1.0E-9;
			
			Ylm::sph_harm (	ucell.atoms[it].nwl,
					dr[0] / distance,
					dr[1] / distance,
					dr[2] / distance,
					ylma);
			// these parameters are about interpolation
			// because once we know the distance from atom to grid point,
			// we can get the parameters we need to do interpolation and
			// store them first!! these can save a lot of effort.
			const double position=distance / delta_r;
			const int ip=static_cast<int>(position);
			const double dx = position - ip;
			const double dx2 = dx * dx;
			const double dx3 = dx2 * dx;

			const double c3 = 3.0*dx2-2.0*dx3;
			const double c1 = 1.0-c3;
			const double c2 = (dx-2.0*dx2+dx3)*delta_r;
			const double c4 = (dx3-dx2)*delta_r;

			double phi=0;
			for (int iw=0; iw< atom->nw; ++iw, ++p)
			{
				if ( atom->iw2_new[iw] )
				{
					pointer=&ORB.Phi[it].PhiLN(atom->iw2l[iw], atom->iw2n[iw]);
					phi=c1*pointer->psi_uniform[ip]+c2*pointer->dpsi_uniform[ip]
						+ c3*pointer->psi_uniform[ip+1] + c4*pointer->dpsi_uniform[ip+1];
				}
				*p=phi * ylma[atom->iw2_ylm[iw]];
			} // end iw
		}// end ib
	}// end id
}

inline void cal_band_rho(const int size, const int grid_index, const int LD_pool, 
                      int* block_iw, int* block_size, int* at, int *uc, int* block_index,
                      double ** psir_ylm, double **psir_DM, double* psir_DM_pool, 
                      int* vindex, bool** cal_flag)
{
	char trans='N';
	double alpha_diag=1, alpha_nondiag=2, beta=1;
	int inc=1;
	int cal_num;
	
	int iw1_lo;
	for(int is=0; is<NSPIN; ++is)
	{
		ZEROS(psir_DM_pool, pw.bxyz*LD_pool);
		for (int ia1=0; ia1<size; ++ia1)
		{
			const int iw1_lo=block_iw[ia1];
			const int iat1=at[ia1];
			const int id1=uc[ia1];
			const int idx1=block_index[ia1];
			const int R1x = GridT.ucell_index2x[id1];
			const int R1y = GridT.ucell_index2y[id1];
			const int R1z = GridT.ucell_index2z[id1];
			const int T1 = ucell.iat2it[iat1];
			int* find_start = LNNR.find_R2[iat1];
			int* find_end = LNNR.find_R2[iat1] + LNNR.nad[iat1];
			//ia2==ia1
			cal_num=0;
			for(int ib=0; ib<pw.bxyz; ++ib)
			{
    			if(cal_flag[ib][ia1])
    			    ++cal_num;
			}
			if(cal_num>pw.bxyz/4)
			{
				//find offset
				const int dRx=0;
				const int dRy=0;
				const int dRz=0;				
				const int index = LNNR.cal_RindexAtom(dRx, dRy, dRz, iat1);
				int offset = -1;
				for(int* find=find_start; find < find_end; find++)
				{
					//--------------------------------------------------------------
					// start positions of adjacent atom of 'iat'
					//--------------------------------------------------------------
					if( find[0] == index ) 
					{
						offset = find - find_start;
						break;
					}
				}

				if(offset == -1)
				{					
					cout << "== charge ==========================" << endl;
					cout << " grid_index = " << grid_index << endl;
                    cout << " index = " << index << endl;
					cout << " size1=" << ia1 << " size2=" << ia1 << endl;
                    cout << " iat1=" << iat1 << " iat2=" << iat1 << endl;
                    cout << " dR=" << dRx << " " << dRy << " " << dRz << endl;
					WARNING_QUIT("gint_k","evaluate_pDMp wrong");
				}
				//const int offset=AllOffset[ia1][ia2];
				assert(offset < LNNR.nad[iat1]);
				
				const int DM_start = LNNR.nlocstartg[iat1]+ LNNR.find_R2st[iat1][offset];					
				dgemm_(&trans, &trans, &block_size[ia1], &pw.bxyz, &block_size[ia1], &alpha_diag,
					&LOC.DM_R[is][DM_start], &block_size[ia1], 
					&psir_ylm[0][idx1], &LD_pool,  
					&beta, &psir_DM[0][idx1], &LD_pool);
			}
			else if(cal_num>0)
			{
				//find offset
				const int dRx=0;
				const int dRy=0;
				const int dRz=0;				
				const int index = LNNR.cal_RindexAtom(dRx, dRy, dRz, iat1);
				int offset = -1;
				for(int* find=find_start; find < find_end; find++)
				{
					//--------------------------------------------------------------
					// start positions of adjacent atom of 'iat'
					//--------------------------------------------------------------
					if( find[0] == index ) 
					{
						offset = find - find_start;
						break;
					}
				}

				if(offset == -1)
				{					
					cout << "== charge ==========================" << endl;
					cout << " grid_index = " << grid_index << endl;
                    cout << " index = " << index << endl;
					cout << " size1=" << ia1 << " size2=" << ia1 << endl;
                    cout << " iat1=" << iat1 << " iat2=" << iat1 << endl;
                    cout << " dR=" << dRx << " " << dRy << " " << dRz << endl;
					WARNING_QUIT("gint_k","evaluate_pDMp wrong");
				}
				//const int offset=AllOffset[ia1][ia2];
				assert(offset < LNNR.nad[iat1]);
				
				const int DM_start = LNNR.nlocstartg[iat1]+ LNNR.find_R2st[iat1][offset];
				for(int ib=0; ib<pw.bxyz; ++ib					)
				{
    				if(cal_flag[ib][ia1])
    				{
        				dgemv_(&trans, &block_size[ia1], &block_size[ia1], &alpha_diag,
					            &LOC.DM_R[is][DM_start], &block_size[ia1], 
					            &psir_ylm[ib][idx1], &inc,  
					            &beta, &psir_DM[ib][idx1], &inc);
    				}
				}
			}
			//ia2>ia1
			for (int ia2=ia1+1; ia2<size; ++ia2)
			{			
			    cal_num=0;
    			for(int ib=0; ib<pw.bxyz; ++ib)
    			{
        			if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
        			    ++cal_num;
    			}
				if(cal_num>pw.bxyz/4)
				{
    				const int iw2_lo=block_iw[ia2];
    				const int iat2=at[ia2];
    				const int T2 = ucell.iat2it[iat2];
    				
    				// find offset
    				const int id2=uc[ia2];
			        const int idx2=block_index[ia2];
    				const int R2x = GridT.ucell_index2x[id2];
    				const int R2y = GridT.ucell_index2y[id2];
    				const int R2z = GridT.ucell_index2z[id2];
    				const int dRx = R1x - R2x;
    				const int dRy = R1y - R2y;
    				const int dRz = R1z - R2z;
    				const int index = LNNR.cal_RindexAtom(dRx, dRy, dRz, iat2);
    				int offset = -1;
    				for(int* find=find_start; find < find_end; find++)
    				{
    					//--------------------------------------------------------------
    					// start positions of adjacent atom of 'iat'
    					//--------------------------------------------------------------
    					if( find[0] == index ) 
    					{
    						offset = find - find_start;
    						break;
    					}
    				}
    				if(offset == -1)
    				{					
    					cout << "== charge ==========================" << endl;
    					cout << " grid_index = " << grid_index << endl;
                        cout << " index = " << index << endl;
    					cout << " size1=" << ia1 << " size2=" << ia2 << endl;
                        cout << " iat1=" << iat1 << " iat2=" << iat2 << endl;
                        cout << " dR=" << dRx << " " << dRy << " " << dRz << endl;
                        cout << " R1=" << R1x << " " << R1y << " " << R1z << endl;
                        cout << " R2=" << R2x << " " << R2y << " " << R2z << endl;
    					WARNING_QUIT("gint_k","evaluate_pDMp wrong");
    				}				
    				assert(offset < LNNR.nad[iat1]);

    				const int DM_start = LNNR.nlocstartg[iat1]+ LNNR.find_R2st[iat1][offset];
    				dgemm_(&trans, &trans, &block_size[ia2], &pw.bxyz, &block_size[ia1], &alpha_nondiag,
    					&LOC.DM_R[is][DM_start], &block_size[ia2], 
    					&psir_ylm[0][idx1], &LD_pool,
    					&beta, &psir_DM[0][idx2], &LD_pool);
				}
				else if(cal_num>0)
				{
    				const int iw2_lo=block_iw[ia2];
    				const int iat2=at[ia2];
    				const int T2 = ucell.iat2it[iat2];
    				
    				// find offset
    				const int id2=uc[ia2];
    				const int idx2=block_index[ia2];
    				const int R2x = GridT.ucell_index2x[id2];
    				const int R2y = GridT.ucell_index2y[id2];
    				const int R2z = GridT.ucell_index2z[id2];
    				const int dRx = R1x - R2x;
    				const int dRy = R1y - R2y;
    				const int dRz = R1z - R2z;
    				const int index = LNNR.cal_RindexAtom(dRx, dRy, dRz, iat2);
    				int offset = -1;
    				for(int* find=find_start; find < find_end; find++)
    				{
    					//--------------------------------------------------------------
    					// start positions of adjacent atom of 'iat'
    					//--------------------------------------------------------------
    					if( find[0] == index ) 
    					{
    						offset = find - find_start;
    						break;
    					}
    				}
    				if(offset == -1)
    				{					
    					cout << "== charge ==========================" << endl;
    					cout << " grid_index = " << grid_index << endl;
                        cout << " index = " << index << endl;
    					cout << " size1=" << ia1 << " size2=" << ia2 << endl;
                        cout << " iat1=" << iat1 << " iat2=" << iat2 << endl;
                        cout << " dR=" << dRx << " " << dRy << " " << dRz << endl;
                        cout << " R1=" << R1x << " " << R1y << " " << R1z << endl;
                        cout << " R2=" << R2x << " " << R2y << " " << R2z << endl;
    					WARNING_QUIT("gint_k","evaluate_pDMp wrong");
    				}				
    				assert(offset < LNNR.nad[iat1]);

    				const int DM_start = LNNR.nlocstartg[iat1]+ LNNR.find_R2st[iat1][offset];
    				for(int ib=0; ib<pw.bxyz; ++ib					)
    				{
        				if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
        				{
            				dgemv_(&trans, &block_size[ia2], &block_size[ia1], &alpha_nondiag,
            					&LOC.DM_R[is][DM_start], &block_size[ia2], 
            					&psir_ylm[ib][idx1], &inc,
            					&beta, &psir_DM[ib][idx2], &inc);
        				}
    				}
				} // cal_num
			}// ia2
		} // ia1
		
		// calculate rho
		double *rhop = CHR.rho[is];
		for(int ib=0; ib<pw.bxyz; ++ib)
		{
			double r=ddot_(&block_index[size], psir_ylm[ib], &inc, psir_DM[ib], &inc);
			const int grid = vindex[ib];
			rhop[ grid ] += r;
		}
	}
}

void Gint_k::calculate_charge(void)
{
	TITLE("Gint_k","calculate_charge");
	timer::tick("Gint_k","charge",'F');

	const double delta_r = ORB.dr_uniform;
	// it's a uniform grid to save orbital values, so the delta_r is a constant.
	const int max_size = GridT.max_atom;
	
	double *psir_ylm_pool, **psir_ylm;
	double *psir_DM_pool, **psir_DM;
	int *block_iw; // index of wave functions of each block;	
	int *block_size; //band size: number of columns of a band
	int *at;
	int *uc;
	int *block_index;
	int* vindex;	
	bool** cal_flag;
	//int** AllOffset;
	int LD_pool=max_size*ucell.nwmax;
    if(max_size!=0)
    {
		psir_ylm_pool=new double[pw.bxyz*LD_pool];
		ZEROS(psir_ylm_pool, pw.bxyz*LD_pool);
		psir_ylm=new double *[pw.bxyz];
		psir_DM_pool=new double[pw.bxyz*LD_pool];
		ZEROS(psir_DM_pool, pw.bxyz*LD_pool);
		psir_DM=new double *[pw.bxyz];
		block_iw=new int[max_size];
		block_size=new int[max_size];
		at=new int[max_size];
		uc=new int[max_size];
		block_index=new int[max_size+1];		

		// mohan fix bug 2011-05-02
		int nn = 0;
		for(int it=0; it<ucell.ntype; it++)
		{
			nn = max(nn,(ucell.atoms[it].nwl+1)*(ucell.atoms[it].nwl+1));
		}

		vindex = new int[pw.bxyz];
		ZEROS(vindex, pw.bxyz);
		
        cal_flag = new bool*[pw.bxyz];

        for(int i=0; i<pw.bxyz; i++)
        {
			psir_ylm[i]=&psir_ylm_pool[i*LD_pool];
			psir_DM[i]=&psir_DM_pool[i*LD_pool];
            cal_flag[i] = new bool[max_size];
        }
    }
	const int nbx = GridT.nbx;
	const int nby = GridT.nby;
	const int nbz_start = GridT.nbzp_start;
	const int nbz = GridT.nbzp;
	
	const int ncyz = pw.ncy*pw.nczp;
	const int nbyz = nby*nbz;	
	for(int i=0; i<nbx; i++)
	{
		const int ibx = i*pw.bx; // mohan add 2012-03-25
		for(int j=0; j<nby; j++)
		{
			const int jby = j*pw.by; // mohan add 2012-03-25
			for(int k=nbz_start; k<nbz_start+nbz; k++)
			{
				const int kbz = k*pw.bz-pw.nczp_start; //mohan add 2012-03-25
				const int grid_index = (k-nbz_start) + j * nbz + i * nbyz;
				const int size = GridT.how_many_atoms[grid_index];
				if(size==0) continue;
				setVindex(ncyz, ibx, jby, kbz, vindex);
				//timer::tick("Gint_k","cal_psir_ylm",'G');
				cal_psir_ylm(size, grid_index, delta_r,
						at, uc, block_index, block_iw, block_size, 
						cal_flag, psir_ylm);
				//timer::tick("Gint_k","cal_psir_ylm",'G');

				cal_band_rho(size, grid_index, LD_pool, block_iw, block_size, at, uc, block_index, psir_ylm, psir_DM, psir_DM_pool, vindex, cal_flag);
			}// int k
		}// int j
	} // int i

	if(max_size!=0)
    {
        for(int i=0; i<pw.bxyz; i++)
        {
        	delete[] cal_flag[i];
        }
		delete[] psir_ylm;
		delete[] psir_ylm_pool;
		delete[] psir_DM;
		delete[] psir_DM_pool;
		delete[] block_iw;
		delete[] block_size;
		delete[] at;
		delete[] block_index;
        delete[] cal_flag;

		delete[] vindex;
    }	

//	cout << " calculate the charge density from density matrix " << endl;

	timer::tick("Gint_k","charge",'F');
	return;
}

void Gint_k::evaluate_pDMp(const int &grid_index, const int &size,
	bool** cal_flag, double*** psir_ylm, int* vindex)
{
	//-----------------------------------------------------
	// in order to calculate <i,alpha,R1 | DM_R | j,beta,R2>
	//-----------------------------------------------------
	double **tchg = new double*[NSPIN];
	for(int is=0; is<NSPIN; is++)
	{
		tchg[is] = new double[pw.bxyz];
		ZEROS(tchg[is], pw.bxyz);
	}

	bool *all_out_of_range = new bool[size];
	for(int ia=0; ia<size; ia++)
	{
		all_out_of_range[ia] = true;
		for(int ib=0; ib<pw.bxyz; ib++)
		{
			if(cal_flag[ib][ia])
			{
				all_out_of_range[ia] = false;
			}
		}
	}

	double *psi1, *psi2;
	double *iw1p, *iw2p;
	double *end1, *end2;
	double *dm, *tchgs;
	int ixxx, ixxx2;
	int iw1_lo, iw2_lo;
	double psi1_2;			
	
	// get (i,alpha,R1)
	// size: how many atoms in this big cell.
	for (int ia1=0; ia1<size; ia1++)
	{
		if(all_out_of_range[ia1]) continue;
		const int mcell_index1 = GridT.bcell_start[grid_index] + ia1;
		const int iat = GridT.which_atom[mcell_index1];
		const int T1 = ucell.iat2it[iat];
		const int I1 = ucell.iat2ia[iat];
		const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
		Atom *atom1 = &ucell.atoms[T1];
		const int nw1 = atom1->nw;

		//~~~~~~~~~~~~~~~~
		// get cell R1.
		//~~~~~~~~~~~~~~~~
		const int id1 = GridT.which_unitcell[mcell_index1];	
		const int R1x = GridT.ucell_index2x[id1];
		const int R1y = GridT.ucell_index2y[id1];
		const int R1z = GridT.ucell_index2z[id1];
		const int DM_start = LNNR.nlocstartg[iat];

		// get (j,beta,R2)
		for (int ia2=0; ia2<size; ia2++)
		{
			if(all_out_of_range[ia2]) continue;
		
			
			bool same_flag = false;
			for(int ib=0; ib<pw.bxyz; ib++)
			{
				if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
				{
					same_flag = true;
					break;
				}
			}
			
			// only go on if the two atoms are in the same meshcell! 
			if(!same_flag)continue;

			const int mcell_index2 = GridT.bcell_start[grid_index] + ia2;
			const int iat2 = GridT.which_atom[mcell_index2];
			const int T2 = ucell.iat2it[iat2];

			if (T2 >= T1)
			{
				Atom *atom2 = &ucell.atoms[T2];
				const int nw2 = atom2->nw;
				const int I2 = ucell.iat2ia[iat2];
				const int start2 = ucell.itiaiw2iwt(T2, I2, 0);

				//~~~~~~~~~~~~~~~~
				// get cell R2.
				//~~~~~~~~~~~~~~~~
				const int id2 = GridT.which_unitcell[mcell_index2];	
				const int R2x = GridT.ucell_index2x[id2];
				const int R2y = GridT.ucell_index2y[id2];
				const int R2z = GridT.ucell_index2z[id2];

				const int dRx = R1x - R2x;
				const int dRy = R1y - R2y;
				const int dRz = R1z - R2z;

// get the index from dRx, dRy, dRz.				
// in fact I can calculate this once, and update only when meshcell is changed.
				const int index = LNNR.cal_RindexAtom(dRx, dRy, dRz, iat2);

				int offset = -1;
				int* find_start = LNNR.find_R2[iat];
				int* find_end = LNNR.find_R2[iat] + LNNR.nad[iat];

				for(int* find=find_start; find < find_end; find++)
				{
					//--------------------------------------------------------------
					// start positions of adjacent atom of 'iat'
					//--------------------------------------------------------------
					if( find[0] == index ) 
					{
						offset = find - find_start;
						break;
					}
				}

				if(offset == -1)
				{					
					cout << "== charge ==========================" << endl;
					cout << " grid_index = " << grid_index << endl;
                    cout << " index = " << index << endl;
					cout << " size1=" << ia1 << " size2=" << ia2 << endl;
                    cout << " iat=" << iat << " iat2=" << iat2 << endl;
                    cout << " dR=" << dRx << " " << dRy << " " << dRz << endl;
                    cout << " R1=" << R1x << " " << R1y << " " << R1z << endl;
                    cout << " R2=" << R2x << " " << R2y << " " << R2z << endl;
					WARNING_QUIT("gint_k","evaluate_pDMp wrong");
				}
				assert(offset < LNNR.nad[iat]);


				// key variable:
				ixxx = DM_start + LNNR.find_R2st[iat][offset];	
				
				for(int is=0; is<NSPIN; is++)
				{
					dm = LOC.DM_R[is];
					tchgs = tchg[is];
					for(int ib=0; ib<pw.bxyz; ib++)
					{
						if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
						{
							psi1 = psir_ylm[ib][ia1];
							psi2 = psir_ylm[ib][ia2];
							end1 = psi1 + nw1;
							end2 = psi2 + nw2;
							
							iw1_lo = GridT.trace_lo[start1]/NPOL;
							ixxx2 = ixxx;
							//------------------------------------
							// circle for wave functions of atom 1.
							//------------------------------------
							for (iw1p=psi1; iw1p<end1; ++iw1p)
							{
								iw2_lo = GridT.trace_lo[start2]/NPOL;
								// 2.0 counts for the undiagonalized part
								psi1_2 = 2.0 * iw1p[0];
								//------------------------------------
								// circle for wave functions of atom 2.
								// dmt: temperary density matrix
								//------------------------------------
								double* dmt = &dm[ixxx2];
								for (iw2p=psi2; iw2p<end2; ++iw2p)
								{
									if ( iw1_lo > iw2_lo)
									{
										++iw2_lo;
										++dmt;
										continue;
									}
									else if ( iw1_lo < iw2_lo)
									{
										tchgs[ib] += 
											psi1_2
											* iw2p[0] 
											* dmt[0];
									}
									else// means iw1_lo == iw2_lo
									{
										tchgs[ib] += 
											iw1p[0]
											* iw2p[0] 
											* dmt[0];
									}
									++iw2_lo;
									++dmt;
								}// iw2
								++iw1_lo;
								ixxx2 += nw2;
							}// iw
						}//end cal_flag
					}// end ib
				}// end is
			}// T
		}// ia2
	}// ia1

	for(int is=0; is<NSPIN; is++)
	{
		for(int ib=0; ib<pw.bxyz; ib++)
		{
			CHR.rho[is][vindex[ib]] += tchg[is][ib];
		}
	}

	delete[] all_out_of_range;

	for(int is=0; is<NSPIN; is++)
	{
		delete[] tchg[is];
	}
	delete[] tchg;

	return;
}


