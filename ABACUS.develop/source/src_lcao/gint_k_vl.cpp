#include "src_pw/tools.h"
#include "gint_k.h"
#include "LCAO_nnr.h"
#include "../module_ORB/ORB_read.h"
#include "grid_technique.h"
#include "src_global/ylm.h"
#include "src_pw/global.h"
#include "module_base/blas_connector.h"
//#include <mkl_cblas.h>

inline int find_offset(const int size, const int grid_index, 
				const int ibx, const int jby, const int kbz, 
				const int bx, const int by, const int bz, 
				const int ia1, const int iat1, const int id1, const int T1, 
				const int ia2, const int iat2, const int id2, const int T2,
				double** distance, int* find_start, int* find_end)
{
	const int R1x=GridT.ucell_index2x[id1];
	const int R2x=GridT.ucell_index2x[id2];
	const int dRx=R1x-R2x;
	const int R1y=GridT.ucell_index2y[id1];
	const int R2y=GridT.ucell_index2y[id2];
	const int dRy=R1y-R2y;
	const int R1z=GridT.ucell_index2z[id1];
	const int R2z=GridT.ucell_index2z[id2];
	const int dRz=R1z-R2z;

	const int index=LNNR.cal_RindexAtom(dRx, dRy, dRz, iat2);
	
	int offset=-1;
	for(int* find=find_start; find < find_end; ++find)
	{
		if( find[0] == index )
		{
			offset = find - find_start;
			break;
		}
	}

	if(offset == -1 )
	{
		ofs_running << "================ BUG REPORT ===================" << endl;
		ofs_running << " grid_index = " << grid_index << endl;
		ofs_running << " index of adjacent atom according to (dRx, dRy, dRz, iat)= " << index << endl;
    	ofs_running << " find list:"<<endl;
		for(int* find=find_start; find < find_end; ++find)
			ofs_running << *find << endl;
		ofs_running << " id2 = " << id2 << endl;
		ofs_running << " T1=" << ucell.atoms[T1].label << " T2=" << ucell.atoms[T2].label << endl;
		ofs_running << " size (how many atoms on this grid) = " << size << endl;
		ofs_running << " ia1=" << ia1 << " ia2=" << ia2 << endl;
		ofs_running << " iat1=" << iat1 << " iat2=" << iat2 << endl;
		ofs_running << " dR=" << dRx << " " << dRy << " " << dRz << endl;
		ofs_running << " R1=" << R1x << " " << R1y << " " << R1z << endl;
		int bindex = 0;
		// z is the fastest,
		for(int ii=0; ii<bx; ii++)
		{
			//const int iii = ibx + ii;
			for(int jj=0; jj<by; jj++)
			{
				//const int jjj = jby + jj;
				for(int kk=0; kk<bz; kk++)
				{
					//const int kkk = kbz + kk;
					if(distance[bindex][ia1] < ORB.Phi[T1].getRcut() )
					{
						ofs_running << " ib=" << bindex << " dis=" << distance[bindex][ia1] 
//									<< " ncxyz (" << (kkk-	.nczp_start) + jjj*	.nczp + iii*	.ncy*	.nczp 
//									<< " = " << iii << " " << jjj << " " << kkk <<") "
//						<< " nbxyz (" << i << " " << j << " " << k << ") "
						<< " bxyz  (" << ii << " " << jj << " " << kk << ") "
						<< " smaller than cutoff = " << setprecision(20) << distance[bindex][ia1] - ORB.Phi[T1].getRcut()
						<< endl;
					}
					else
					{
						ofs_running << " ib=" << bindex << " dis=" << distance[bindex][ia2]
						<< " rcut=" << ORB.Phi[T1].getRcut() << endl;
					}
					++bindex;
				}
			}
		}

		ofs_running << " R2=" << R2x << " " << R2y << " " << R2z << endl;
		bindex = 0;
		// z is the fastest,
		for(int ii=0; ii<bx; ii++)
		{
			//const int iii = ibx + ii;
			for(int jj=0; jj<by; jj++)
			{
				//const int jjj = jby + jj;
				for(int kk=0; kk<bz; kk++)
				{
					//const int kkk = kbz + kk;
					if(distance[bindex][ia2] < ORB.Phi[T2].getRcut() )//mohan T1->T2
					{
						ofs_running << " ib=" << bindex << " dis=" << distance[bindex][ia2] 
//									<< " ncxyz (" << (kkk-	.nczp_start) + jjj*pw.nczp + iii*	.ncy*pw.nczp 
//									<< " = " << iii << " " << jjj << " " << kkk <<") "
//						<< " nbxyz (" << i << " " << j << " " << k << ") "
						<< " bxyz  (" << ii << " " << jj << " " << kk << ") "
						<< endl;
					}
					else
					{
						ofs_running << " ib=" << bindex << " dis=" << distance[bindex][ia2]
						<< " rcut=" << ORB.Phi[T2].getRcut() << endl;
					}
					++bindex;
				}
			}
		}

		ofs_running << " target index = " << index << endl;
		ofs_running << " iat=" << iat1 << " nad=" << LNNR.nad[iat1] << endl;
		for(int iii = 0; iii < LNNR.nad[iat1]; iii++)
		{
			ofs_running << " ad=" << iii << " find_R2=" << LNNR.find_R2[iat1][iii] << endl;
		}
		ofs_warning << " The adjacent atom found by 	 is not found by SLTK_Adjacent program!" << endl;
		WARNING_QUIT("gint_k","evaluate_pvpR_reduced wrong");
	}
	assert(offset < LNNR.nad[iat1]);
	return offset;
}

inline void cal_psir_ylm(int size, int grid_index, double delta_r,
						double** distance,
						int* at, int* block_index, int* block_iw, int* block_size, 
						bool** cal_flag, double** psir_ylm)
{
	const Numerical_Orbital_Lm *pointer;
	double mt[3];
	double dr[3];
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

			distance[ib][id]=std::sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
			if(distance[ib][id] > (ORB.Phi[it].getRcut()- 1.0e-15)) 
			{
				cal_flag[ib][id]=false;
				ZEROS(p, block_size[id]);
				continue;
			}

			cal_flag[ib][id]=true;
			
			std::vector<double> ylma;
			//if(distance[id] > GridT.orbital_rmax) continue;
			//	Ylm::get_ylm_real(this->nnn[it], this->dr[id], ylma);
			if (distance[ib][id] < 1.0E-9) distance[ib][id] += 1.0E-9;
			
			Ylm::sph_harm (	ucell.atoms[it].nwl,
					dr[0] / distance[ib][id],
					dr[1] / distance[ib][id],
					dr[2] / distance[ib][id],
					ylma);
			// these parameters are about interpolation
			// because once we know the distance from atom to grid point,
			// we can get the parameters we need to do interpolation and
			// store them first!! these can save a lot of effort.
			const double position=distance[ib][id] / delta_r;
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
					pointer=&ORB.Phi[it].PhiLN(
							atom->iw2l[iw],
							atom->iw2n[iw]);
					phi=c1*pointer->psi_uniform[ip]+c2*pointer->dpsi_uniform[ip]
						+ c3*pointer->psi_uniform[ip+1] + c4*pointer->dpsi_uniform[ip+1];
				}
				*p=phi * ylma[atom->iw2_ylm[iw]];
			} // end iw
		}// end ib
	}// end id
}

inline void cal_pvpR_reduced(int size, int LD_pool, int grid_index, 
							const int ibx, const int jby, const int kbz,
							int* block_size, int* at, int* block_index, int* block_iw,
							double* vldr3, double** psir_ylm, double** psir_vlbr3, 
							double** distance, bool** cal_flag, double* pvpR)
{
	char transa='N', transb='T';
	double alpha=1, beta=1;
	int allnw=block_index[size];
	for(int i=0; i<pw.bxyz; ++i)
	{
		for(int j=0; j<allnw; ++j)
		{
			psir_vlbr3[i][j]=psir_ylm[i][j]*vldr3[i];
		}
	}

	int k=pw.bxyz;
	for(int ia1=0; ia1<size; ++ia1)
	{
		//if(all_out_of_range[ia1]) continue;
		//const int iw1_lo=block_iw[ia1];
		const int idx1=block_index[ia1];
		int m=block_size[ia1];
		const int iat1=at[ia1];
		const int T1 = ucell.iat2it[iat1];
		const int mcell_index1 = GridT.bcell_start[grid_index] + ia1;
		const int id1 = GridT.which_unitcell[mcell_index1];
		const int DM_start = LNNR.nlocstartg[iat1];
		// nad : how many adjacent atoms for atom 'iat'
		int* find_start = LNNR.find_R2[iat1];
		int* find_end = LNNR.find_R2[iat1] + LNNR.nad[iat1];
		for(int ia2=0; ia2<size; ++ia2)
		{
			const int iat2=at[ia2];
			const int T2 = ucell.iat2it[iat2];
			if (T1 <= T2)
			{
    			int cal_num=0;
    			for(int ib=0; ib<pw.bxyz; ++ib)
    			{
    				if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
    				    ++cal_num;
    			}

    			if(cal_num==0) continue;
    			
    			//const int iw2_lo=block_iw[ia2];
                const int idx2=block_index[ia2];
        		int n=block_size[ia2];
				//const int I2 = ucell.iat2ia[iat2];
				const int mcell_index2 = GridT.bcell_start[grid_index] + ia2;
				const int id2 = GridT.which_unitcell[mcell_index2];
				int offset;

				offset=find_offset(size, grid_index, 
						ibx, jby, kbz, 
						pw.bx, pw.by, pw.bz, 
						ia1, iat1, id1, T1, 
						ia2, iat2, id2, T2, 
						distance, find_start, find_end);

				const int iatw = DM_start + LNNR.find_R2st[iat1][offset];	

			    if(cal_num>pw.bxyz/4)
			    {
    				//if(iw1_lo<=iw2_lo)
    				//{
    			        k=pw.bxyz;
    					dgemm_(&transa, &transb, &n, &m, &k, &alpha,
    						&psir_vlbr3[0][idx2], &LD_pool, 
    						&psir_ylm[0][idx1], &LD_pool,
    						&beta, &pvpR[iatw], &n);
    				//}
				}
    			else
    			{
        			//if(iw1_lo<=iw2_lo)
    				//{
            			for(int ib=0; ib<pw.bxyz; ++ib)
            			{
                			if(cal_flag[ib][ia1]&&cal_flag[ib][ia2])
                			{
                    			k=1;
            					dgemm_(&transa, &transb, &n, &m, &k, &alpha,
            						&psir_vlbr3[ib][idx2], &LD_pool, 
            						&psir_ylm[ib][idx1], &LD_pool,
            						&beta, &pvpR[iatw], &n);	
                			}
            			}
        			//}
    			}
			}
		}
	}
}

void Gint_k::cal_vlocal_k(const double *vrs1, const Grid_Technique &GridT, const int spin)
{
	TITLE("Gint_k","cal_vlocal_k");

	if(!pvpR_alloc_flag)
	{
		WARNING_QUIT("Gint_k::destroy_pvpR","pvpR has not been allocated yet!");
	}
	else
	{
		// reduce the dimension of array, which
		// is used to save <phi | Vl | phi>
		if(this->reduced)
		{
			ZEROS(this->pvpR_reduced[spin], LNNR.nnrg);
		}
		// else one needs to consdier all cell with a vector R
		// the number of cells is GridT.nutot,
		// and the elments in each processor is GridT.lgd.
		else
		{
			for(int i=0; i<GridT.lgd * GridT.nutot; i++)
			{
				ZEROS(pvpR[i], GridT.lgd * GridT.nutot);
			}
		}
	}

	timer::tick("Gint_k","vlocal", 'K');

	// it's a uniform grid to save orbital values, so the delta_r is a constant.
	double delta_r = ORB.dr_uniform;
	// possible max atom number in real space grid. 
	const int max_size = GridT.max_atom;
	// how many meshcells in bigcell.
	const int bxyz = pw.bxyz;
	const int LD_pool=max_size*ucell.nwmax;
	
	double **distance = nullptr; // distance between atom and grid: [bxyz, maxsize]
	double *psir_ylm_pool = nullptr;
	double **psir_ylm = nullptr;
	double *psir_vlbr3_pool = nullptr;
	double **psir_vlbr3 = nullptr;
	bool** cal_flag = nullptr;
	int *block_iw = nullptr; // index of wave functions of each block;	
	int *block_size = nullptr; //band size: number of columns of a band
	int *at = nullptr;
	int *block_index = nullptr;
	if(max_size!=0)
	{
		// save the small box information for a big box.
		distance = new double*[bxyz];
		psir_ylm_pool=new double[bxyz*LD_pool];
		ZEROS(psir_ylm_pool, bxyz*LD_pool);
		psir_ylm=new double *[bxyz];
		psir_vlbr3_pool=new double[bxyz*LD_pool];
		ZEROS(psir_vlbr3_pool, bxyz*LD_pool);
		psir_vlbr3=new double *[bxyz];
		cal_flag = new bool*[bxyz];
		block_iw=new int[max_size];
		block_size=new int[max_size];
		at=new int[max_size];
		block_index=new int[max_size+1];		

		// mohan fix bug 2011-05-02
		// possible number of atom configureation (l,m)
		int nn = 0;
		for(int it=0; it<ucell.ntype; it++)
		{
			nn = std::max(nn, (ucell.atoms[it].nwl+1)*(ucell.atoms[it].nwl+1));
		}

		for(int i=0; i<bxyz; i++)
		{
			// possible max atom number in a big box.
			psir_ylm[i]=&psir_ylm_pool[i*LD_pool];
			psir_vlbr3[i]=&psir_vlbr3_pool[i*LD_pool];
			distance[i] = new double[max_size];
			cal_flag[i] = new bool[max_size];

			ZEROS(distance[i], max_size);
			ZEROS(cal_flag[i], max_size);
		}
	}
	
	assert(this->ncxyz!=0);
	const double dv = ucell.omega/this->ncxyz;
	int vl_index=0;

	// array to store local potential for each small box in
	// a big box.
	double* vldr3 = new double[bxyz];
	ZEROS(vldr3, bxyz);

	for(int i=0; i<nbx; i++)
	{
		const int ibx=i*pw.bx;
		for(int j=0; j<nby; j++)
		{
			const int jby=j*pw.by;
			// count the z according to big box.
			for(int k=nbz_start; k<nbz_start+nbz; k++)
			{
				const int kbz=k*pw.bz-pw.nczp_start;
				const int grid_index = (k-nbz_start) + j * nbz + i * nby * nbz;
				const int size = GridT.how_many_atoms[ grid_index ];
				if(size==0) continue;

				//--------------------------------- 
				// get the wave functions in this
				// grid.
				//--------------------------------- 
				cal_psir_ylm(size, grid_index, delta_r, distance,
						at, block_index, block_iw, block_size, 
						cal_flag, psir_ylm);

				int bindex = 0;
				// z is the fastest,
				for(int ii=0; ii<pw.bx; ii++)
				{
					const int iii = ibx + ii;
					for(int jj=0; jj<pw.by; jj++)
					{
						const int jjj = jby + jj;
						for(int kk=0; kk<pw.bz; kk++)
						{
							const int kkk = kbz + kk;
							vl_index = kkk + jjj*pw.nczp + iii*pw.ncy*pw.nczp;
							vldr3[bindex] = vrs1[ vl_index ] * dv;
							++bindex;
						}
					}
				}

				if(this->reduced)
				{
					cal_pvpR_reduced(size, LD_pool, grid_index, 
									ibx, jby, kbz, 
									block_size, at, block_index, block_iw, 
									vldr3, psir_ylm, psir_vlbr3, 
									distance, cal_flag, this->pvpR_reduced[spin]);
				}
				else
				{
					//this->evaluate_pvpR_full(grid_index, size, psir_ylm, cal_flag, vldr3);
					cout<<"call pvpR_full"<<endl;
				}
			}// int k
		}// int j
	} // int i

	delete[] vldr3;
	if(max_size!=0)
	{
		for(int i=0; i<bxyz; i++)
		{
			delete[] distance[i];
			delete[] cal_flag[i];
		}
		delete[] distance;
		delete[] psir_ylm;
		delete[] psir_ylm_pool;
		delete[] psir_vlbr3;
		delete[] psir_vlbr3_pool;
		delete[] cal_flag;
		delete[] block_iw;
		delete[] block_size;
		delete[] block_index;
	}	

	timer::tick("Gint_k","vlocal", 'K');
	return;
}

void Gint_k::evaluate_pvpR_reduced(
	double* pvpR, 
	const int &grid_index, 
	const int &size, // Atom number on this grid point r. 
	const int &i, const int &j, const int &k,
	double*** psir_ylm, 
	bool** cal_flag, // Flag to determine if an atom is the adjacent atom. 
	double* vldr3, 
	double** distance, 
	const Grid_Technique &gt)
{
	double *psi1, *psi2;
	double *iw1p, *iw2p;
	double *end1, *end2;
	double *pvp;
	int iw1_lo, iw2_lo;
	int iwi, iww;
	double vpsir1;

	bool *all_out_of_range = new bool[size];
	for(int ia=0; ia<size; ia++)
	{
		all_out_of_range[ia] = true;
		for(int ib=0; ib<gt.bxyz; ib++)
		{
			if(cal_flag[ib][ia])
			{
				all_out_of_range[ia] = false;
			}
		}
	}

	for (int ia1=0; ia1<size; ++ia1)
	{
		if(all_out_of_range[ia1]) continue;

        const int mcell_index1 = gt.bcell_start[grid_index] + ia1;
		const int iat = gt.which_atom[mcell_index1];
        const int T1 = ucell.iat2it[iat];
        const int I1 = ucell.iat2ia[iat];
        const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
		const int iw1_start = gt.trace_lo[start1]/NPOL;
        Atom *atom1 = &ucell.atoms[T1];
	
        //~~~~~~~~~~~~~~~~
        // get cell R1.
        //~~~~~~~~~~~~~~~~
        const int id1 = gt.which_unitcell[mcell_index1];
        const int R1x = gt.ucell_index2x[id1];
        const int R1y = gt.ucell_index2y[id1];
        const int R1z = gt.ucell_index2z[id1];
        const int DM_start = LNNR.nlocstartg[iat];

        // get (j,beta,R2)
        for (int ia2=0; ia2<size; ++ia2)
        {
			if(all_out_of_range[ia2]) continue;

			//---------------------------------------------
			// check if we need to calculate the big cell.
			//---------------------------------------------
			bool same_flag = false;
			for(int ib=0; ib<gt.bxyz; ++ib)
			{
				if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
				{
//					ofs_running << " ia1=" << ia1 << " ia2=" << ia2;
//					ofs_running << " cal_flag1=" << cal_flag[ib][ia1] << " cal_flag2=" << cal_flag[ib][ia2] << endl;
					same_flag = true;
					break;
				}
			} 

			if(!same_flag) continue;

            const int bcell2 = gt.bcell_start[grid_index] + ia2;
			const int iat2 = gt.which_atom[bcell2];
            const int T2 = ucell.iat2it[iat2];

            if (T2 >= T1)
            {
                Atom *atom2 = &ucell.atoms[T2];
                const int I2 = ucell.iat2ia[iat2];
                const int start2 = ucell.itiaiw2iwt(T2, I2, 0);
				const int iw2_start = gt.trace_lo[start2]/NPOL;

	            //~~~~~~~~~~~~~~~~
                // get cell R2.
                //~~~~~~~~~~~~~~~~
                const int id2 = gt.which_unitcell[bcell2];
                const int R2x = gt.ucell_index2x[id2];
                const int R2y = gt.ucell_index2y[id2];
                const int R2z = gt.ucell_index2z[id2];

				//------------------------------------------------
				// calculate the 'offset': R2 position relative
				// to  R1 atom.
				//------------------------------------------------
                const int dRx = R1x - R2x;
                const int dRy = R1y - R2y;
                const int dRz = R1z - R2z;
	
				const int index = LNNR.cal_RindexAtom(dRx, dRy, dRz, iat2);
                int offset = -1;

				// nad : how many adjacent atoms for atom 'iat'
				int* find_start = LNNR.find_R2[iat];
				int* findend = LNNR.find_R2[iat] + LNNR.nad[iat];
				
				// the nad should be a large expense of time.
				for(int* find=find_start; find < findend; ++find)
				{
					if( find[0] == index )
					{
						offset = find - find_start;
						break;
					}
				}

				if(offset == -1 )
                {
					ofs_running << "================ BUG REPORT ===================" << endl;
					ofs_running << " grid_index = " << grid_index << endl;
                    ofs_running << " index of adjacent atom according to (dRx, dRy, dRz, iat)= " << index << endl;
					ofs_running << " id2 = " << id2 << endl;
					ofs_running << " T1=" << ucell.atoms[T1].label << " T2=" << ucell.atoms[T2].label << endl;
					ofs_running << " size (how many atoms on this grid) = " << size << endl;
					ofs_running << " ia1=" << ia1 << " ia2=" << ia2 << endl;
                    ofs_running << " iat=" << iat << " iat2=" << iat2 << endl;
                    ofs_running << " dR=" << dRx << " " << dRy << " " << dRz << endl;

                    ofs_running << " R1=" << R1x << " " << R1y << " " << R1z << endl;
					int bindex = 0;
					// z is the fastest,
					for(int ii=0; ii<gt.bx; ii++)
					{
						for(int jj=0; jj<gt.by; jj++)
						{
							for(int kk=0; kk<gt.bz; kk++)
							{
								//const int iii = i*gt.bx + ii;
								//const int jjj = j*gt.by + jj;
								//const int kkk = k*gt.bz + kk;
								if(distance[bindex][ia1] < ORB.Phi[T1].getRcut() )
								{
									ofs_running << " ib=" << bindex << " dis=" << distance[bindex][ia1] 
//									<< " ncxyz (" << (kkk-gt.nczp_start) + jjj*gt.nczp + iii*gt.ncy*gt.nczp 
//									<< " = " << iii << " " << jjj << " " << kkk <<") "
									<< " nbxyz (" << i << " " << j << " " << k << ") "
									<< " bxyz  (" << ii << " " << jj << " " << kk << ") "
									<< " smaller than cutoff = " << setprecision(20) << distance[bindex][ia1] - ORB.Phi[T1].getRcut()
									<< endl;
								}
								else
								{
									ofs_running << " ib=" << bindex << " dis=" << distance[bindex][ia2]
									<< " rcut=" << ORB.Phi[T1].getRcut() << endl;
								}
								++bindex;
							}
						}
					}


                    ofs_running << " R2=" << R2x << " " << R2y << " " << R2z << endl;
					bindex = 0;
					// z is the fastest,
					for(int ii=0; ii<gt.bx; ii++)
					{
						for(int jj=0; jj<gt.by; jj++)
						{
							for(int kk=0; kk<gt.bz; kk++)
							{
								//const int iii = i*gt.bx + ii;
								//const int jjj = j*gt.by + jj;
								//const int kkk = k*gt.bz + kk;
								if(distance[bindex][ia2] < ORB.Phi[T2].getRcut() )//mohan T1->T2
								{
									ofs_running << " ib=" << bindex << " dis=" << distance[bindex][ia2] 
//									<< " ncxyz (" << (kkk-gt.nczp_start) + jjj*pw.nczp + iii*gt.ncy*pw.nczp 
//									<< " = " << iii << " " << jjj << " " << kkk <<") "
									<< " nbxyz (" << i << " " << j << " " << k << ") "
									<< " bxyz  (" << ii << " " << jj << " " << kk << ") "
									<< endl;
								}
								else
								{
									ofs_running << " ib=" << bindex << " dis=" << distance[bindex][ia2]
									<< " rcut=" << ORB.Phi[T2].getRcut() << endl;
								}
								++bindex;
							}
						}
					}

					ofs_running << " target index = " << index << endl;
					ofs_running << " iat=" << iat << " nad=" << LNNR.nad[iat] << endl;
                    for(int iii = 0; iii < LNNR.nad[iat]; iii++)
                    {
                        ofs_running << " ad=" << iii << " find_R2=" << LNNR.find_R2[iat][iii] << endl;
                    }
					ofs_warning << " The adjacent atom found by gt is not found by SLTK_Adjacent program!" << endl;
                    WARNING_QUIT("gint_k","evaluate_pvpR_reduced wrong");
                }
                assert(offset < LNNR.nad[iat]);

				// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				// what I do above is to get 'offset' for atom pair (iat1, iat2)
				// if I want to simplify this searching for offset,
				// I should take advantage of gt.which_unitcell.
				// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

				const int iatw = DM_start + LNNR.find_R2st[iat][offset];
				
				for(int ib=0; ib<gt.bxyz; ib++)
				{
					if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
					{
						psi1 = psir_ylm[ib][ia1];
						psi2 = psir_ylm[ib][ia2];
						end1 = psi1 + atom1->nw;
						end2 = psi2 + atom2->nw;
						iw1_lo = iw1_start;	
						//------------------------------------
						// circle for wave functions of atom 1.
						//------------------------------------
						iwi = 0;
						for (iw1p=psi1; iw1p < end1; ++ iw1p)
						{
							vpsir1 = iw1p[0] * vldr3[ib];
							iw2_lo = iw2_start;
							iww = iatw + iwi;// -1 because ++iww from below.
							iwi += atom2->nw;
							pvp = &pvpR[iww];

							//------------------------------------
							// circle for wave functions of atom 2.
							//------------------------------------
							for(iw2p=psi2; iw2p < end2; ++ iw2p)
							{
								if( iw1_lo > iw2_lo)
								{
									++iw2_lo;
									++pvp;
									continue;
								}
								pvp[0] += vpsir1 * iw2p[0];
								++iw2_lo;
								++pvp;
							}
							++iw1_lo;
						}// iw
					}//end flag
				}//end ib
            }// T
        }// ia2
	}

	delete[] all_out_of_range;
	return;
}


void Gint_k::evaluate_pvpR_full(const int &grid_index, const int &size, double*** psir_ylm, bool** cal_flag, double* vldr3)
{
	//-----------------------------------------------------
	// in order to calculate <i,alpha,R1 | V | j,beta,R2>
	//-----------------------------------------------------
	// get (i,alpha,R1)
	for (int ia1=0; ia1<size; ia1++)
	{
		const int mcell_index1 = GridT.bcell_start[grid_index] + ia1;
		const int T1 = ucell.iat2it[ GridT.which_atom[mcell_index1] ];
		const int I1 = ucell.iat2ia[ GridT.which_atom[mcell_index1] ];
		const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
		Atom *atom1 = &ucell.atoms[T1];

		//~~~~~~~~~~~~~~~~
		// get cell R1.
		//~~~~~~~~~~~~~~~~
		const int id1 = GridT.which_unitcell[mcell_index1];	
		const int dim1 = id1 * GridT.lgd;

		// get (j,beta,R2)
		for (int ia2=0; ia2<size; ia2++)
		{
			const int mcell_index2 = GridT.bcell_start[grid_index] + ia2;
			const int T2 = ucell.iat2it[ GridT.which_atom[mcell_index2]];

			if (T2 >= T1)
			{
				Atom *atom2 = &ucell.atoms[T2];
				const int I2 = ucell.iat2ia[ GridT.which_atom[mcell_index2]];
				const int start2 = ucell.itiaiw2iwt(T2, I2, 0);

				//~~~~~~~~~~~~~~~~
				// get cell R2.
				//~~~~~~~~~~~~~~~~
				const int id2 = GridT.which_unitcell[mcell_index2];
				const int dim2 = id2 * GridT.lgd;

				// circle for wave functions of atom 1.
				for(int ib=0; ib<pw.bxyz; ib++)
				{
					if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
					{
						double *psi1 = psir_ylm[ib][ia1];
						double *psi2 = psir_ylm[ib][ia2];
						
						int iw1_lo = GridT.trace_lo[start1]/NPOL;
						for (int iw=0; iw< atom1->nw; iw++)
						{
							const double vpsir1 = psi1[iw] * vldr3[ib];
							int iw2_lo = GridT.trace_lo[start2]/NPOL;
							double *ppp = this->pvpR[iw1_lo + dim1];
							double *pppp = &ppp[iw2_lo+dim2];
							double *pppp_end = &ppp[iw2_lo+dim2] + atom2->nw;
							double *psi2_p = psi2;
							
							// circle for wave functions of atom 2.
							for (;pppp < pppp_end; ++iw2_lo, ++pppp, ++psi2_p)
							{
								if ( iw1_lo > iw2_lo)
								{
									continue;
								}
								pppp[0] += vpsir1 * psi2_p[0];
							}// iw2
							++iw1_lo;
						}// iw
					}// T
				}// ia2
			}// ia1
		}// end cal_flag
	}// end ib

	return;
}

