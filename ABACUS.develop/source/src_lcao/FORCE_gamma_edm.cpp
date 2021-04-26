#include "FORCE_gamma.h"
#include "../src_pw/global.h"

double Force_LCAO_gamma::set_EDM_element(
    const int &ii, const int &jj, 
    const bool with_energy,
    double*** coef1, double*** coef2, const int &is)
{
    double ene = 0.0;
    for (int ib = 0; ib < NBANDS; ib++)
    {
        const double w1 = wf.wg(is,ib);
        if(w1>0.0)
        {
            if(with_energy)
            {
                ene += w1 * wf.ekb[is][ib] * coef1[is][ib][ii] * coef2[is][ib][jj];
            }
            else
            {
                ene += w1 * coef1[is][ib][ii] * coef2[is][ib][jj];
            }
        }
    }
    return ene;
}

//#include "../src_develop/src_siao/selinv.h"
void Force_LCAO_gamma::set_EDM_gamma(matrix& dm, bool with_energy)
{
    TITLE("Force_LCAO_gamma","set_EDM_gamma");
    timer::tick("Force_LCAO_gamma","set_EDM");

#ifdef __SELINV
//xiaohui modified 2013-03-23, adding "/*"
/*    if(DIAGO_TYPE=="selinv")
    {
//      cout << " set_EDM_gamma()" << endl;
//      cout << " fill dm with density matrix or energy density matrix here." << endl;

        //--------------------------------------------------
        // job=1: density matrix for charge density (grid).
        // job=2: density matrix for force (dm2d).
        // job=3: energy density matrix for force (dm2d).
        //--------------------------------------------------
        const int ik = 0;
        int job = 0;
        if(with_energy) 
        {
            job=3;
        }
        else 
        {
            job=2;
        }
        Selinv::dm_ptr = dm;
        Selinv::using_SELINV(ik, job, LM.Hloc, LM.Sloc);    

        return;
    }
*/
#endif
    //xiaohui add 2015-03-24, like density matrix calculation
#ifdef __MPI //2015-09-06, xiaohui
    int nprocs,myid;
    MPI_Status status;
    MPI_Comm_size(DIAG_HPSEPS_WORLD,&nprocs);
    MPI_Comm_rank(DIAG_HPSEPS_WORLD,&myid);

    int local_band=NBANDS/DSIZE;
	if (DRANK<NBANDS%DSIZE) local_band=local_band+1;

	int *local_bands;
	local_bands = new int[DSIZE];
	ZEROS(local_bands, DSIZE);

	int lastband_in_proc = 0;
	int lastband_number = 0;
	int count_bands = 0;
	for (int i=0; i<DSIZE; i++)
	{
		if (i<NBANDS%DSIZE)
		{
			local_bands[i]=NBANDS/DSIZE+1;
		}
		else
		{
			local_bands[i]=NBANDS/DSIZE;
		}
		count_bands += local_bands[i];
		if (count_bands >= NBANDS)
		{
			lastband_in_proc = i;
			lastband_number = NBANDS - (count_bands - local_bands[i]);
			break;
		}
	}

	for(int ispin=0;ispin<NSPIN;ispin++)//need to optimize later! noted by zhengdy20210303
	{
		double** w1;
		w1 = new double*[NSPIN];
		for(int is=0; is<NSPIN; is++)
		{
			w1[is] = new double[local_band];
			ZEROS(w1[is], local_band);
		}

		int total_bands = 0;

		time_t time_wg2w1_start = time(NULL);
		double Time_Wg2w1_Start = (double)clock();

		int Total_Bands = 0;
		int Nbands = -1;
		for(int i=0; i <= lastband_in_proc; i++)
		{
			Nbands = local_bands[i];
			if(myid == i)
			{
				for (int ib=0; ib< Nbands; ib++)
				{
					if(with_energy)
					{
						w1[ispin][ib] = wf.wg(ispin, ib + Total_Bands) * wf.ekb[ispin][ib + Total_Bands];
					}
					else
					{
						w1[ispin][ib] = wf.wg(ispin, ib + Total_Bands);
					}
				}
			}
			Total_Bands += Nbands;
		}

		double** Z_wg;
		Z_wg = new double*[NSPIN];
		for(int is=0; is<NSPIN; is++)
		{
			Z_wg[is] = new double[local_band * NLOCAL];
			ZEROS(Z_wg[is], local_band * NLOCAL);
		}

		time_t time_Z_wg_start = time(NULL);
		double Time_Z_Wg_Start = (double)clock();

		if(myid <= lastband_in_proc)
		{
			for(int i=0; i<local_bands[myid]; i++)
			{
				for(int n=0; n<NLOCAL; n++)
				{
					Z_wg[ispin][i + n*local_bands[myid]] = ParaO.Z_LOC[ispin][i + n*local_bands[myid]]* w1[ispin][i];
				}
			}
		}

		for(int is=0; is<NSPIN; is++)
		{
			delete[] w1[is];
		}
		delete[] w1;
		delete[] local_bands;

		time_t time_Z_wg_end = time(NULL);
		double Time_Z_Wg_End =(double)clock();
		double time_Z_wg = difftime(time_Z_wg_start,time_Z_wg_end);

		double coll_time = 0.0;
		double Time_Cal_Start = (double)clock();

		int row_col = NLOCAL/300;
		if(row_col >= 1)
		{
			double** ZW_300;
			ZW_300 = new double*[NSPIN];
			for(int is=0; is<NSPIN; is++)
			{
				ZW_300[is] = new double[300 * local_band];
				ZEROS(ZW_300[is], 300 * local_band);
			}

			double** Z_300;
			Z_300 = new double*[NSPIN];
			for(int is=0; is<NSPIN; is++)
			{
				Z_300[is] = new double[300 * local_band];
				ZEROS(Z_300[is], 300 * local_band);
			}
			if(NLOCAL%300 != 0) row_col += 1;
			int row_global = 0;
			int row_index = 0;
			int col_global = 0;
			int col_index = 0;
			for(int row_count=1; row_count<=row_col; row_count++)
			{
				row_global = row_count*300;
				if(row_global <= NLOCAL)
				{
					for(int i_row=0; i_row<300; i_row++)
					{
						row_index = i_row +(row_count-1)*300;
						for(int ib=0; ib<local_band; ib++)
						{
							ZW_300[ispin][ib + i_row*local_band] = Z_wg[ispin][ib + row_index*local_band] ;
						}
					}

					for(int col_count=1; col_count<=row_col; col_count++)
					{
						col_global = col_count*300;
						if(col_global <= NLOCAL)
						{
							double **rho_300;
							rho_300 = new double*[NSPIN];
							for(int is=0; is<NSPIN; is++)
							{
								rho_300[is] = new double[300*300];
								ZEROS(rho_300[is], 300*300);
							}

							for(int i_col=0; i_col<300; i_col++)
							{
								col_index = i_col +(col_count-1)*300;
								for(int ib=0; ib<local_band; ib++)
								{
									Z_300[ispin][ib + i_col*local_band] = ParaO.Z_LOC[ispin][ib + col_index*local_band] ;
								}
							}

							double alpha=1;
							double beta=0;
							char transa = 'T';
							char transb = 'N';
							int M = 300;
							int N = 300;
							int K = local_band;
							int lda = local_band;
							int ldb = local_band;
							int ldc = 300;
							dgemm_(&transa, &transb, &M, &N, &K,
								&alpha, &ZW_300[ispin][0], &lda, &Z_300[ispin][0], &ldb,
								&beta, &rho_300[ispin][0], &ldc);

							double Time_Coll_Start = (double)clock();
							MPI_Barrier(DIAG_HPSEPS_WORLD);
							Parallel_Reduce::reduce_double_all( rho_300[ispin], 300*300);
							double Time_Coll_End = (double)clock();

							coll_time += (Time_Coll_End-Time_Coll_Start);
							if(GAMMA_ONLY_LOCAL)
								{
								for(int i_row=0; i_row<300; i_row++)
								{
									row_index = i_row +(row_count-1)*300;
									int row_mu = ParaO.trace_loc_row[row_index];
									for(int i_col=0; i_col<300; i_col++)
									{
										col_index = i_col +(col_count-1)*300;
										int col_nu = ParaO.trace_loc_col[col_index];
										if(row_mu >= 0 && col_nu >= 0)
										{
											int index = row_mu * ParaO.ncol + col_nu;
											//dm[index] = rho_300[ispin][i_col + i_row*300];
											dm(ispin, index) = rho_300[ispin][i_row + i_col*300];
											}
									}
								}
							}
							for(int is=0; is<NSPIN; is++)
							{
								delete[] rho_300[is];
							}
							delete[] rho_300;
						}
						else
						{
							int col_remain = 0;
							int col_index2 = 0;
							col_remain = NLOCAL - (col_count - 1)*300;
							double** Z_remain;
							Z_remain = new double*[NSPIN];
							for(int is=0; is<NSPIN; is++)
							{
								Z_remain[is] = new double[col_remain * local_band];
								ZEROS(Z_remain[is], col_remain * local_band);
							}

							double** rho_300_remain;
							rho_300_remain = new double*[NSPIN];
							for(int is=0; is<NSPIN; is++)
							{
								rho_300_remain[is] = new double[300*col_remain];
								ZEROS(rho_300_remain[is], 300*col_remain);
							}
							for(int i_col=0; i_col<col_remain; i_col++)
							{
								col_index2 = i_col +(col_count-1)*300;
								for(int ib=0; ib<local_band; ib++)
								{
									Z_remain[ispin][ib + i_col*local_band] = ParaO.Z_LOC[ispin][ib + col_index2*local_band] ;
								}
							}

							double alpha=1;
							double beta=0;
							char transa = 'T';
							char transb = 'N';
							int M = 300;
							int N = col_remain;
							int K = local_band;
							int lda = local_band;
							int ldb = local_band;
							int ldc = 300;
							dgemm_(&transa, &transb, &M, &N, &K,
								&alpha, &ZW_300[ispin][0], &lda, &Z_remain[ispin][0], &ldb,
									&beta, &rho_300_remain[ispin][0], &ldc);

							double Time_Coll_Start = (double)clock();
							MPI_Barrier(DIAG_HPSEPS_WORLD);
							Parallel_Reduce::reduce_double_all( rho_300_remain[ispin], 300*col_remain);
							double Time_Coll_End = (double)clock();

							coll_time += (Time_Coll_End-Time_Coll_Start);
							if(GAMMA_ONLY_LOCAL)
								{
								for(int i_row=0; i_row<300; i_row++)
								{
									row_index = i_row +(row_count-1)*300;
									int row_mu = ParaO.trace_loc_row[row_index];
									for(int i_col=0; i_col<col_remain; i_col++)
									{
										col_index = i_col +(col_count-1)*300;
										int col_nu = ParaO.trace_loc_col[col_index];
										if(row_mu >= 0 && col_nu >= 0)
										{
											int index = row_mu * ParaO.ncol + col_nu;
											//dm[index] = rho_300_remain[ispin][i_col + i_row*col_remain];
											dm(ispin, index) = rho_300_remain[ispin][i_row + i_col*300];
										}
									}
								}
							}
							for(int is=0; is<NSPIN; is++)
							{
								delete[] Z_remain[is];
								delete[] rho_300_remain[is];
							}
							delete[] Z_remain;
							delete[] rho_300_remain;
						}
					}
				}
				else
				{
					int row_remain = 0;
					int row_index2 = 0;
					row_remain = NLOCAL - (row_count - 1)*300;
					double** Z_row_remain;
					Z_row_remain = new double*[NSPIN];
					for(int is=0; is<NSPIN; is++)
					{
						Z_row_remain[is] = new double[row_remain * local_band];
						ZEROS(Z_row_remain[is], row_remain * local_band);
					}

					for(int i_row=0; i_row<row_remain; i_row++)
					{
						row_index2 = i_row +(row_count-1)*300;
						for(int ib=0; ib<local_band; ib++)
						{
							Z_row_remain[ispin][ib + i_row*local_band] = Z_wg[ispin][ib + row_index2*local_band] ;
						}
					}

					for(int col_count=1; col_count<=row_col; col_count++)
					{
						col_global = col_count*300;
						if(col_global <= NLOCAL)
						{
							double** rho_remain_300;
							rho_remain_300 = new double*[NSPIN];
							for(int is=0; is<NSPIN; is++)
							{
								rho_remain_300[is] = new double[row_remain * 300];
								ZEROS(rho_remain_300[is], row_remain * 300);
							}

							for(int i_col=0; i_col<300; i_col++)
							{
								col_index = i_col +(col_count-1)*300;
								for(int ib=0; ib<local_band; ib++)
								{
									Z_300[ispin][ib + i_col*local_band] = ParaO.Z_LOC[ispin][ib + col_index*local_band] ;
								}
							}

							double alpha=1;
							double beta=0;
							char transa = 'T';
							char transb = 'N';
							int M = row_remain;
							int N = 300;
							int K = local_band;
							int lda = local_band;
							int ldb = local_band;
							int ldc = row_remain;
							dgemm_(&transa, &transb, &M, &N, &K,
								&alpha, &Z_row_remain[ispin][0], &lda, &Z_300[ispin][0], &ldb,
								&beta, &rho_remain_300[ispin][0], &ldc);

							double Time_Coll_Start = (double)clock();
							MPI_Barrier(DIAG_HPSEPS_WORLD);
							Parallel_Reduce::reduce_double_all( rho_remain_300[ispin], row_remain*300);
							double Time_Coll_End = (double)clock();

							coll_time += (Time_Coll_End-Time_Coll_Start);
							if(GAMMA_ONLY_LOCAL)
								{
								for(int i_row=0; i_row<row_remain; i_row++)
								{
									row_index = i_row +(row_count-1)*300;
									int row_mu = ParaO.trace_loc_row[row_index];
									for(int i_col=0; i_col<300; i_col++)
									{
										col_index = i_col +(col_count-1)*300;
										int col_nu = ParaO.trace_loc_col[col_index];
										if(row_mu >= 0 && col_nu >= 0)
										{
											int index = row_mu * ParaO.ncol + col_nu;
											dm(ispin, index) = rho_remain_300[ispin][i_row + i_col*row_remain];
										}
									}
								}
							}
							for(int is=0; is<NSPIN; is++)
							{
								delete[] rho_remain_300[is];
							}
							delete[] rho_remain_300;
						}
						else
						{
							int col_remain = 0;
							int col_index2 = 0;
							col_remain = NLOCAL - (col_count - 1)*300;
							double** Z_col_remain;
							Z_col_remain = new double*[NSPIN];
							for(int is=0; is<NSPIN; is++)
							{
								Z_col_remain[is] = new double[col_remain * local_band];
								ZEROS(Z_col_remain[is], col_remain * local_band);
							}
							double** rho_remain_remain;
							rho_remain_remain = new double*[NSPIN];
							for(int is=0; is<NSPIN; is++)
							{
								rho_remain_remain[is] = new double[row_remain * col_remain];
								ZEROS(rho_remain_remain[is], row_remain * col_remain);
							}
							for(int i_col=0; i_col<col_remain; i_col++)
							{
								col_index2 = i_col +(col_count-1)*300;
								for(int ib=0; ib<local_band; ib++)
								{
									Z_col_remain[ispin][ib + i_col*local_band] = ParaO.Z_LOC[ispin][ib + col_index2*local_band] ;
								}
							}

							double alpha=1;
							double beta=0;
							char transa = 'T';
							char transb = 'N';
							int M = row_remain;
							int N = col_remain;
							int K = local_band;
							int lda = local_band;
							int ldb = local_band;
							int ldc = row_remain;
							dgemm_(&transa, &transb, &M, &N, &K,
								&alpha, &Z_row_remain[ispin][0], &lda, &Z_col_remain[ispin][0], &ldb,
								&beta, &rho_remain_remain[ispin][0], &ldc);

							double Time_Coll_Start = (double)clock();
							MPI_Barrier(DIAG_HPSEPS_WORLD);
							Parallel_Reduce::reduce_double_all( rho_remain_remain[ispin], row_remain*col_remain);
							double Time_Coll_End = (double)clock();

							coll_time += (Time_Coll_End-Time_Coll_Start);
							if(GAMMA_ONLY_LOCAL)
								{
								for(int i_row=0; i_row<row_remain; i_row++)
								{
									row_index = i_row +(row_count-1)*300;
									int row_mu = ParaO.trace_loc_row[row_index];
									for(int i_col=0; i_col<col_remain; i_col++)
									{
										col_index = i_col +(col_count-1)*300;
										int col_nu = ParaO.trace_loc_col[col_index];
										if(row_mu >= 0 && col_nu >= 0)
										{
											int index = row_mu * ParaO.ncol + col_nu;
											dm(ispin, index) = rho_remain_remain[ispin][i_row + i_col*row_remain];
										}
									}
								}
							}
							for(int is=0; is<NSPIN; is++)
							{
								delete[] Z_col_remain[is];
								delete[] rho_remain_remain[is];
							}
							delete[] Z_col_remain;
							delete[] rho_remain_remain;
								
						}
					}
					for(int is=0; is<NSPIN; is++)
					{
						delete[] Z_row_remain[is];
					}
					delete[] Z_row_remain;
				}
			}
			for(int is=0; is<NSPIN; is++)
			{
				delete[] ZW_300[is];
				delete[] Z_300[is];
			}
			delete[] ZW_300;
			delete[] Z_300;
		}
		else
		{
			double** rho_NLOCAL_NLOCAL;
			rho_NLOCAL_NLOCAL = new double*[NSPIN];
			for(int is=0; is<NSPIN; is++)
			{
				rho_NLOCAL_NLOCAL[is] = new double[NLOCAL*NLOCAL];
				ZEROS(rho_NLOCAL_NLOCAL[is], NLOCAL*NLOCAL);
			}
			for(int i_row=0; i_row<NLOCAL; i_row++)
			{
				for(int i_col=0; i_col<NLOCAL; i_col++)
				{
					for(int ib=0; ib<local_band; ib++)
					{
						rho_NLOCAL_NLOCAL[ispin][i_col + i_row*NLOCAL] += Z_wg[ispin][ib + i_row*local_band]  * ParaO.Z_LOC[ispin][ib + i_col*local_band];
					}
				}
			}
			double Time_Coll_Start = (double)clock();
			MPI_Barrier(DIAG_HPSEPS_WORLD);
			Parallel_Reduce::reduce_double_all( rho_NLOCAL_NLOCAL[ispin], NLOCAL*NLOCAL);
			double Time_Coll_End = (double)clock();
			coll_time += (Time_Coll_End-Time_Coll_Start);
			if(GAMMA_ONLY_LOCAL)
			{
				for(int i_row=0; i_row<NLOCAL; i_row++)
				{
					int row_mu = ParaO.trace_loc_row[i_row];
					for(int i_col=0; i_col<NLOCAL; i_col++)
					{
						int col_nu = ParaO.trace_loc_col[i_col];
						if(row_mu >= 0 && col_nu >= 0)
						{
							int index = row_mu * ParaO.ncol + col_nu;
							dm(ispin, index) = rho_NLOCAL_NLOCAL[ispin][i_col + i_row*NLOCAL];
						}
					}
				}
			}
			for(int is=0; is<NSPIN; is++)
			{
				delete[] rho_NLOCAL_NLOCAL[is];
			}
			delete[] rho_NLOCAL_NLOCAL;
		}

		for(int is=0; is<NSPIN; is++)
		{
			delete[] Z_wg[is];
		}
		delete[] Z_wg;
	}//end ispin
#endif //2015-09-06, xiaohui
#ifndef __MPI //2015-09-06, xiaohui
    // orbital 1
	for(int ispin=0;ispin<NSPIN;ispin++)
	{
		for(int i=0; i<NLOCAL; i++)
		{
			const int mu = ParaO.trace_loc_row[i];
			if(mu<0) continue;
			const int ii = GridT.trace_lo[i];
			// orbital 2
			for(int j=0; j<NLOCAL; j++)
			{
				const int nu = ParaO.trace_loc_col[j];
				if(nu<0) continue;
				const int jj = GridT.trace_lo[j];
				// energy density matrix
				const int index = mu * ParaO.ncol + nu;
				double ene = 0.0;

				// (1) density matrix can be generated from LOWF.WFC_GAMMA directly.
				if(ii>=0 && jj>=0)
				{
					ene = this->set_EDM_element(ii, jj, with_energy, LOWF.WFC_GAMMA, LOWF.WFC_GAMMA, ispin);
				}

				// (2)
				else if(ii>=0 && jj<0)
				{
					const int a4 = LOWF.trace_aug[j];
					assert(a4>=0);
					ene = this->set_EDM_element(ii, a4, with_energy,  LOWF.WFC_GAMMA,  LOWF.WFC_GAMMA_aug, ispin);
				}
				else if(ii<0 && jj>=0)
				{
					const int a3 = LOWF.trace_aug[i];
					assert(a3>=0);
					// mohan fix serious bug 2011-07-01 (ii, a3) -> (a3, jj) !!!!!!!!!!!!
					ene = this->set_EDM_element(a3, jj, with_energy, LOWF.WFC_GAMMA_aug, LOWF.WFC_GAMMA, ispin);
				}
				else if(ii<0 && jj<0)
				{
					const int a3 = LOWF.trace_aug[i];
					const int a4 = LOWF.trace_aug[j];
					assert(a3>=0);
					assert(a4>=0);
					ene = this->set_EDM_element(a3, a4, with_energy, LOWF.WFC_GAMMA_aug, LOWF.WFC_GAMMA_aug, ispin);
				}

				dm(ispin, index) = ene;
				//dm[index] = 1.0;// mohan tmp
			}// end j
		}// end i
	}//end ispin
#endif //2015-09-06, xiaohui
    timer::tick("Force_LCAO_gamma","set_EDM");
    return;
}


// force due to the overlap matrix.
// need energy density matrix here.
void Force_LCAO_gamma::cal_foverlap(
	const bool isforce, 
	const bool isstress, 
	matrix& foverlap, 
	matrix& soverlap)
{
    TITLE("Force_LCAO_gamma","cal_foverlap");
    timer::tick("Force_LCAO_gamma","cal_foverlap",'G');

    // set energy density matrix.
    if(INPUT.new_dm>0)
    {
        timer::tick("Force_LCAO_gamma","cal_edm_2d",'H');

        matrix wgEkb;
        wgEkb.create(NSPIN, NBANDS);

        for(int is=0; is<NSPIN; is++)
        {
            for(int ib=0; ib<NBANDS; ib++)
            {
                wgEkb(is,ib) = wf.wg(is,ib) * wf.ekb[is][ib];
            }
        }

        Wfc_Dm_2d wfc_edm_2d;
        wfc_edm_2d.init();
        wfc_edm_2d.wfc_gamma=LOC.wfc_dm_2d.wfc_gamma;
        wfc_edm_2d.cal_dm(wgEkb);

        timer::tick("Force_LCAO_gamma","cal_edm_2d",'H');

        for(int i=0; i<NLOCAL; i++)
        {
            const int iat = ucell.iwt2iat[i];
            for(int j=0; j<NLOCAL; j++)
            {
                const int mu = ParaO.trace_loc_row[j];
                const int nu = ParaO.trace_loc_col[i];
                if(mu>=0 && nu>=0)
                {
                    const int index = mu * ParaO.ncol + nu;
                    double sum = 0.0;
                    for(int is=0; is<NSPIN; ++is)
                    {
                        sum += wfc_edm_2d.dm_gamma[is](nu, mu);
                    }
                    sum *= 2.0;

					if(isforce)
					{
						foverlap(iat,0) += sum * LM.DSloc_x[index];
						foverlap(iat,1) += sum * LM.DSloc_y[index];
						foverlap(iat,2) += sum * LM.DSloc_z[index];
					}

                    if(isstress)
                    {
                        soverlap(0,0) += sum/2.0 * LM.DSloc_11[index];
                        soverlap(0,1) += sum/2.0 * LM.DSloc_12[index];
                        soverlap(0,2) += sum/2.0 * LM.DSloc_13[index];
                        soverlap(1,1) += sum/2.0 * LM.DSloc_22[index];
                        soverlap(1,2) += sum/2.0 * LM.DSloc_23[index];
                        soverlap(2,2) += sum/2.0 * LM.DSloc_33[index];
                    }
                }

            }
        }
    }
    else
    {
        timer::tick("Force_LCAO_gamma","cal_edm_grid",'H');
        matrix edm2d;
		edm2d.create(NSPIN, ParaO.nloc);

        bool with_energy = true;

        // for different spins. mohan fix bugs 2012-07-19

		this->set_EDM_gamma(edm2d, with_energy);

        timer::tick("Force_LCAO_gamma","cal_edm_grid",'H');

        //summation \sum_{i,j} E(i,j)*dS(i,j)
        //BEGIN CALCULATION OF FORCE OF EACH ATOM

        for(int i=0; i<NLOCAL; i++)
        {
            const int iat = ucell.iwt2iat[i];
            for(int j=0; j<NLOCAL; j++)
            {
                const int mu = ParaO.trace_loc_row[j];
                const int nu = ParaO.trace_loc_col[i];

                if (mu >= 0 && nu >= 0 )
                {
                    const int index = mu * ParaO.ncol + nu;

                    //================================================================
                    // here is the normal order, the force of each atom is calculated
                    // according to each 'column' in DSloc_x, y, z
                    // because the DSloc_x,y,z are anti-symmetry matrix.
                    //================================================================

                    double sum = 0.0;
                    for(int is=0; is<NSPIN; ++is)
                    {
                        sum += edm2d(is,index);
                    }
                    sum *= 2.0;

                    if(isforce)
					{
						foverlap(iat,0) += sum * LM.DSloc_x[index];
						foverlap(iat,1) += sum * LM.DSloc_y[index];
						foverlap(iat,2) += sum * LM.DSloc_z[index];
					}

                    if(isstress)
                    {
                        soverlap(0,0) += sum/2.0 * LM.DSloc_11[index];
                        soverlap(0,1) += sum/2.0 * LM.DSloc_12[index];
                        soverlap(0,2) += sum/2.0 * LM.DSloc_13[index];
                        soverlap(1,1) += sum/2.0 * LM.DSloc_22[index];
                        soverlap(1,2) += sum/2.0 * LM.DSloc_23[index];
                        soverlap(2,2) += sum/2.0 * LM.DSloc_33[index];
					}
                }
            }
        }
    }

    if(isstress)
    {
        for(int i=0;i<3;i++)
        {
            for(int j=0;j<3;j++)
            {
                if(i<j) soverlap(j,i) = soverlap(i,j);
				soverlap(i,j) *=  ucell.lat0 / ucell.omega;
            }
        }
    }
    timer::tick("Force_LCAO_gamma","cal_foverlap",'G');
    return;
}
