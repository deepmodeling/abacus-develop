#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "gint_k.h"
#include "../module_orbital/ORB_read.h"
#include "grid_technique.h"
#include "../module_base/ylm.h"
#include "../src_pw/global.h"
#include "../module_base/blas_connector.h"
#include "../module_base/timer.h"
//#include <mkl_cblas.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef __MKL
#include <mkl_service.h>
#endif

inline int find_offset(const int size, const int grid_index, 
				const int ibx, const int jby, const int kbz, 
				const int bx, const int by, const int bz, 
				const int ia1, const int iat1, const int id1, const int T1, 
				const int ia2, const int iat2, const int id2, const int T2,
				double** distance, int* find_start, int* find_end)
{
	const int R1x=GlobalC::GridT.ucell_index2x[id1];
	const int R2x=GlobalC::GridT.ucell_index2x[id2];
	const int dRx=R1x-R2x;
	const int R1y=GlobalC::GridT.ucell_index2y[id1];
	const int R2y=GlobalC::GridT.ucell_index2y[id2];
	const int dRy=R1y-R2y;
	const int R1z=GlobalC::GridT.ucell_index2z[id1];
	const int R2z=GlobalC::GridT.ucell_index2z[id2];
	const int dRz=R1z-R2z;

	const int index=GlobalC::GridT.cal_RindexAtom(dRx, dRy, dRz, iat2);
	
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
		GlobalV::ofs_running << "================ BUG REPORT ===================" << std::endl;
		GlobalV::ofs_running << " grid_index = " << grid_index << std::endl;
		GlobalV::ofs_running << " index of adjacent atom according to (dRx, dRy, dRz, iat)= " << index << std::endl;
    	GlobalV::ofs_running << " find list:"<<std::endl;
		for(int* find=find_start; find < find_end; ++find)
			GlobalV::ofs_running << *find << std::endl;
		GlobalV::ofs_running << " id2 = " << id2 << std::endl;
		GlobalV::ofs_running << " T1=" << GlobalC::ucell.atoms[T1].label << " T2=" << GlobalC::ucell.atoms[T2].label << std::endl;
		GlobalV::ofs_running << " size (how many atoms on this grid) = " << size << std::endl;
		GlobalV::ofs_running << " ia1=" << ia1 << " ia2=" << ia2 << std::endl;
		GlobalV::ofs_running << " iat1=" << iat1 << " iat2=" << iat2 << std::endl;
		GlobalV::ofs_running << " dR=" << dRx << " " << dRy << " " << dRz << std::endl;
		GlobalV::ofs_running << " R1=" << R1x << " " << R1y << " " << R1z << std::endl;
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
					if(distance[bindex][ia1] < GlobalC::ORB.Phi[T1].getRcut() )
					{
						GlobalV::ofs_running << " ib=" << bindex << " dis=" << distance[bindex][ia1] 
//									<< " ncxyz (" << (kkk-	.nczp_start) + jjj*	.nczp + iii*	.ncy*	.nczp 
//									<< " = " << iii << " " << jjj << " " << kkk <<") "
//						<< " nbxyz (" << i << " " << j << " " << k << ") "
						<< " bxyz  (" << ii << " " << jj << " " << kk << ") "
						<< " smaller than cutoff = " << std::setprecision(20) << distance[bindex][ia1] - GlobalC::ORB.Phi[T1].getRcut()
						<< std::endl;
					}
					else
					{
						GlobalV::ofs_running << " ib=" << bindex << " dis=" << distance[bindex][ia2]
						<< " rcut=" << GlobalC::ORB.Phi[T1].getRcut() << std::endl;
					}
					++bindex;
				}
			}
		}

		GlobalV::ofs_running << " R2=" << R2x << " " << R2y << " " << R2z << std::endl;
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
					if(distance[bindex][ia2] < GlobalC::ORB.Phi[T2].getRcut() )//mohan T1->T2
					{
						GlobalV::ofs_running << " ib=" << bindex << " dis=" << distance[bindex][ia2] 
//									<< " ncxyz (" << (kkk-	.nczp_start) + jjj*GlobalC::pw.nczp + iii*	.ncy*GlobalC::pw.nczp 
//									<< " = " << iii << " " << jjj << " " << kkk <<") "
//						<< " nbxyz (" << i << " " << j << " " << k << ") "
						<< " bxyz  (" << ii << " " << jj << " " << kk << ") "
						<< std::endl;
					}
					else
					{
						GlobalV::ofs_running << " ib=" << bindex << " dis=" << distance[bindex][ia2]
						<< " rcut=" << GlobalC::ORB.Phi[T2].getRcut() << std::endl;
					}
					++bindex;
				}
			}
		}

		GlobalV::ofs_running << " target index = " << index << std::endl;
		GlobalV::ofs_running << " iat=" << iat1 << " nad=" << GlobalC::GridT.nad[iat1] << std::endl;
		for(int iii = 0; iii < GlobalC::GridT.nad[iat1]; iii++)
		{
			GlobalV::ofs_running << " ad=" << iii << " find_R2=" << GlobalC::GridT.find_R2[iat1][iii] << std::endl;
		}
		GlobalV::ofs_warning << " The adjacent atom found by 	 is not found by SLTK_Adjacent program!" << std::endl;
		ModuleBase::WARNING_QUIT("gint_k","find_offset wrong");
	}
	assert(offset < GlobalC::GridT.nad[iat1]);
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
		const int mcell_index=GlobalC::GridT.bcell_start[grid_index] + id;
		const int imcell=GlobalC::GridT.which_bigcell[mcell_index];

		const int iat=GlobalC::GridT.which_atom[mcell_index];
		at[id]=iat;
		
		const int it=GlobalC::ucell.iat2it[iat];
		const int ia=GlobalC::ucell.iat2ia[iat];
		const int start=GlobalC::ucell.itiaiw2iwt(it, ia, 0);
		block_iw[id]=GlobalC::GridT.trace_lo[start]/GlobalV::NPOL;
		Atom* atom=&GlobalC::ucell.atoms[it];
		block_size[id]=atom->nw;
		block_index[id+1]=block_index[id]+atom->nw;
		// meshball_positions should be the bigcell position in meshball
		// to the center of meshball.
		// calculated in cartesian coordinates
		// the std::vector from the grid which is now being operated to the atom position.
		// in meshball language, is the std::vector from imcell to the center cel, plus
		// tau_in_bigcell.
		mt[0]=GlobalC::GridT.meshball_positions[imcell][0] - GlobalC::GridT.tau_in_bigcell[iat][0];
		mt[1]=GlobalC::GridT.meshball_positions[imcell][1] - GlobalC::GridT.tau_in_bigcell[iat][1];
		mt[2]=GlobalC::GridT.meshball_positions[imcell][2] - GlobalC::GridT.tau_in_bigcell[iat][2];

		for(int ib=0; ib<GlobalC::pw.bxyz; ib++)
		{
			double *p=&psir_ylm[ib][block_index[id]];
			// meshcell_pos: z is the fastest
			dr[0]=GlobalC::GridT.meshcell_pos[ib][0] + mt[0]; 
			dr[1]=GlobalC::GridT.meshcell_pos[ib][1] + mt[1]; 
			dr[2]=GlobalC::GridT.meshcell_pos[ib][2] + mt[2]; 	

			distance[ib][id]=std::sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
			if(distance[ib][id] > (GlobalC::ORB.Phi[it].getRcut()- 1.0e-15)) 
			{
				cal_flag[ib][id]=false;
				ModuleBase::GlobalFunc::ZEROS(p, block_size[id]);
				continue;
			}

			cal_flag[ib][id]=true;
			
			std::vector<double> ylma;
			//if(distance[id] > GlobalC::GridT.orbital_rmax) continue;
			//	Ylm::get_ylm_real(this->nnn[it], this->dr[id], ylma);
			if (distance[ib][id] < 1.0E-9) distance[ib][id] += 1.0E-9;
			
			ModuleBase::Ylm::sph_harm (	GlobalC::ucell.atoms[it].nwl,
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
					pointer=&GlobalC::ORB.Phi[it].PhiLN(
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
	for(int i=0; i<GlobalC::pw.bxyz; ++i)
	{
		for(int j=0; j<allnw; ++j)
		{
			psir_vlbr3[i][j]=psir_ylm[i][j]*vldr3[i];
		}
	}

	int k=GlobalC::pw.bxyz;
	for(int ia1=0; ia1<size; ++ia1)
	{
		//if(all_out_of_range[ia1]) continue;
		//const int iw1_lo=block_iw[ia1];
		const int idx1=block_index[ia1];
		int m=block_size[ia1];
		const int iat1=at[ia1];
		const int T1 = GlobalC::ucell.iat2it[iat1];
		const int mcell_index1 = GlobalC::GridT.bcell_start[grid_index] + ia1;
		const int id1 = GlobalC::GridT.which_unitcell[mcell_index1];
		const int DM_start = GlobalC::GridT.nlocstartg[iat1];
		// nad : how many adjacent atoms for atom 'iat'
		int* find_start = GlobalC::GridT.find_R2[iat1];
		int* find_end = GlobalC::GridT.find_R2[iat1] + GlobalC::GridT.nad[iat1];
		for(int ia2=0; ia2<size; ++ia2)
		{
			const int iat2=at[ia2];
			const int T2 = GlobalC::ucell.iat2it[iat2];
			if (iat1 <= iat2)
			{
    			int cal_num=0;
    			for(int ib=0; ib<GlobalC::pw.bxyz; ++ib)
    			{
    				if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
    				    ++cal_num;
    			}

    			if(cal_num==0) continue;
    			
                const int idx2=block_index[ia2];
        		int n=block_size[ia2];
				//const int I2 = GlobalC::ucell.iat2ia[iat2];
				const int mcell_index2 = GlobalC::GridT.bcell_start[grid_index] + ia2;
				const int id2 = GlobalC::GridT.which_unitcell[mcell_index2];
				int offset;

				offset=find_offset(size, grid_index, 
						ibx, jby, kbz, 
						GlobalC::pw.bx, GlobalC::pw.by, GlobalC::pw.bz, 
						ia1, iat1, id1, T1, 
						ia2, iat2, id2, T2, 
						distance, find_start, find_end);

				const int iatw = DM_start + GlobalC::GridT.find_R2st[iat1][offset];	

			    if(cal_num>GlobalC::pw.bxyz/4)
			    {
					k=GlobalC::pw.bxyz;
					dgemm_(&transa, &transb, &n, &m, &k, &alpha,
						&psir_vlbr3[0][idx2], &LD_pool, 
						&psir_ylm[0][idx1], &LD_pool,
						&beta, &pvpR[iatw], &n);
				}
    			else
    			{
					for(int ib=0; ib<GlobalC::pw.bxyz; ++ib)
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
    			}
			}
		}
	}
}

void Gint_k::cal_vlocal_k(const double *vrs1, const Grid_Technique &GridT, const int spin)
{
	ModuleBase::TITLE("Gint_k","cal_vlocal_k");

	if(!pvpR_alloc_flag)
	{
		ModuleBase::WARNING_QUIT("Gint_k::cal_vlocal_k","pvpR has not been allocated yet!");
	}
	else
	{
		ModuleBase::GlobalFunc::ZEROS(this->pvpR_reduced[spin], GlobalC::GridT.nnrg);
	}

	ModuleBase::timer::tick("Gint_k","vlocal");

#ifdef __MKL
    const int mkl_threads = mkl_get_max_threads();
    mkl_set_num_threads(1);
#endif

#ifdef _OPENMP
    #pragma omp parallel
    {
        double* pvpR_reduced_thread;
        pvpR_reduced_thread = new double[GlobalC::GridT.nnrg];
        ModuleBase::GlobalFunc::ZEROS(pvpR_reduced_thread, GlobalC::GridT.nnrg);
#endif
	// it's a uniform grid to save orbital values, so the delta_r is a constant.
	double delta_r = GlobalC::ORB.dr_uniform;
	// possible max atom number in real space grid. 
	const int max_size = GridT.max_atom;
	// how many meshcells in bigcell.
	const int bxyz = GlobalC::pw.bxyz;
	const int LD_pool=max_size*GlobalC::ucell.nwmax;
	
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
		ModuleBase::GlobalFunc::ZEROS(psir_ylm_pool, bxyz*LD_pool);
		psir_ylm=new double *[bxyz];
		psir_vlbr3_pool=new double[bxyz*LD_pool];
		ModuleBase::GlobalFunc::ZEROS(psir_vlbr3_pool, bxyz*LD_pool);
		psir_vlbr3=new double *[bxyz];
		cal_flag = new bool*[bxyz];
		block_iw=new int[max_size];
		block_size=new int[max_size];
		at=new int[max_size];
		block_index=new int[max_size+1];		

		// mohan fix bug 2011-05-02
		// possible number of atom configureation (l,m)
		int nn = 0;
		for(int it=0; it<GlobalC::ucell.ntype; it++)
		{
			nn = std::max(nn, (GlobalC::ucell.atoms[it].nwl+1)*(GlobalC::ucell.atoms[it].nwl+1));
		}

		for(int i=0; i<bxyz; i++)
		{
			// possible max atom number in a big box.
			psir_ylm[i]=&psir_ylm_pool[i*LD_pool];
			psir_vlbr3[i]=&psir_vlbr3_pool[i*LD_pool];
			distance[i] = new double[max_size];
			cal_flag[i] = new bool[max_size];

			ModuleBase::GlobalFunc::ZEROS(distance[i], max_size);
			ModuleBase::GlobalFunc::ZEROS(cal_flag[i], max_size);
		}
	}
	
	assert(this->ncxyz!=0);
	const double dv = GlobalC::ucell.omega/this->ncxyz;
	int vl_index=0;

	// array to store local potential for each small box in
	// a big box.
	double* vldr3 = new double[bxyz];
	ModuleBase::GlobalFunc::ZEROS(vldr3, bxyz);

#ifdef _OPENMP
        #pragma omp for
#endif
	for(int i=0; i<nbx; i++)
	{
		const int ibx=i*GlobalC::pw.bx;
		for(int j=0; j<nby; j++)
		{
			const int jby=j*GlobalC::pw.by;
			// count the z according to big box.
			for(int k=nbz_start; k<nbz_start+nbz; k++)
			{
				const int kbz=k*GlobalC::pw.bz-GlobalC::pw.nczp_start;
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
				for(int ii=0; ii<GlobalC::pw.bx; ii++)
				{
					const int iii = ibx + ii;
					for(int jj=0; jj<GlobalC::pw.by; jj++)
					{
						const int jjj = jby + jj;
						for(int kk=0; kk<GlobalC::pw.bz; kk++)
						{
							const int kkk = kbz + kk;
							vl_index = kkk + jjj*GlobalC::pw.nczp + iii*GlobalC::pw.ncy*GlobalC::pw.nczp;
							vldr3[bindex] = vrs1[ vl_index ] * dv;
							++bindex;
						}
					}
				}

#ifdef _OPENMP
				cal_pvpR_reduced(size, LD_pool, grid_index, 
								ibx, jby, kbz, 
								block_size, at, block_index, block_iw, 
								vldr3, psir_ylm, psir_vlbr3, 
								distance, cal_flag, pvpR_reduced_thread);
#else
				cal_pvpR_reduced(size, LD_pool, grid_index, 
								ibx, jby, kbz, 
								block_size, at, block_index, block_iw, 
								vldr3, psir_ylm, psir_vlbr3, 
								distance, cal_flag, this->pvpR_reduced[spin]);
#endif
			}// int k
		}// int j
	} // int i

#ifdef _OPENMP
        #pragma omp critical(cal_vl_k)
        for(int innrg=0; innrg<GlobalC::GridT.nnrg; innrg++)
        {
            pvpR_reduced[spin][innrg] += pvpR_reduced_thread[innrg];
        }
        delete[] pvpR_reduced_thread;
#endif

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
#ifdef _OPENMP
    } // end omp
#endif

#ifdef __MKL
    mkl_set_num_threads(mkl_threads);
#endif

	ModuleBase::timer::tick("Gint_k","vlocal");
	return;
}