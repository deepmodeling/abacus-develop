#include "wf_local.h"
#include "../src_pw/global.h"

inline int globalIndex(int localIndex, int nblk, int nprocs, int myproc)
{
    int iblock, gIndex;
    iblock=localIndex/nblk;
    gIndex=(iblock*nprocs+myproc)*nblk+localIndex%nblk;
    return gIndex;
    //return (localIndex/nblk*nprocs+myproc)*nblk+localIndex%nblk;
}

inline int CTOT2q(
	int myid,
	int naroc[2],
	int nb,
	int dim0,
	int dim1,
	int iprow,
	int ipcol,
	double* work,
	double** CTOT)
{
    for(int j=0; j<naroc[1]; ++j)
    {
        int igcol=globalIndex(j, nb, dim1, ipcol);
        if(igcol>=NBANDS) continue;
        for(int i=0; i<naroc[0]; ++i)
        {
            int igrow=globalIndex(i, nb, dim0, iprow);
			//ofs_running << "i,j,igcol,igrow" << i<<" "<<j<<" "<<igcol<<" "<<igrow<<endl;
            if(myid==0) work[j*naroc[0]+i]=CTOT[igcol][igrow];
        }
    }
    return 0;
}

inline int CTOT2q_c(
	int myid,
	int naroc[2],
	int nb,
	int dim0,
	int dim1,
	int iprow,
	int ipcol,
	complex<double>* work,
	complex<double>** CTOT)
{
    for(int j=0; j<naroc[1]; ++j)
    {
        int igcol=globalIndex(j, nb, dim1, ipcol);
        if(igcol>=NBANDS) continue;
        for(int i=0; i<naroc[0]; ++i)
        {
            int igrow=globalIndex(i, nb, dim0, iprow);
			//ofs_running << "i,j,igcol,igrow" << i<<" "<<j<<" "<<igcol<<" "<<igrow<<endl;
            if(myid==0) work[j*naroc[0]+i]=CTOT[igcol][igrow];
        }
    }
    return 0;
}

// be called in local_orbital_wfc::allocate_k
int WF_Local::read_lowf_complex(complex<double> **c, const int &ik, const bool &newdm)
{
    TITLE("WF_Local","read_lowf_complex");
    timer::tick("WF_Local","read_lowf_complex");

    complex<double> **ctot;

    stringstream ss;
	// read wave functions
	// write is in ../src_pdiag/pdiag_basic.cpp
    ss << global_out_dir << "LOWF_K_" << ik+1 <<".dat";
//	cout << " name is = " << ss.str() << endl;

    ifstream ifs;

    int error = 0;

    if (DRANK==0)
    {
        ifs.open(ss.str().c_str());
        if (!ifs)
        {
            ofs_warning << " Can't open file:" << ss.str() << endl;
            error = 1;
        }
    }

#ifdef __MPI
    Parallel_Common::bcast_int(error);
#endif

    if (error==1) return 1;

    // otherwise, find the file.

    if (MY_RANK==0)
    {
		int ikr;
		double kx,ky,kz;
        int nbands, nlocal;
		READ_VALUE(ifs, ikr);
		ifs >> kx >> ky >> kz;
        READ_VALUE(ifs, nbands);
        READ_VALUE(ifs, nlocal);

		if(ikr!=ik+1)
		{
			ofs_warning << " ikr=" << ikr << " ik=" << ik << endl;
			ofs_warning << " k index is not correct" << endl;
			error = 4;
		}
		else if ( 
			abs(kx-kv.kvec_c[ik].x)>1.0e-5 ||
			abs(ky-kv.kvec_c[ik].y)>1.0e-5 ||
			abs(kz-kv.kvec_c[ik].z)>1.0e-5 )
		{	
			ofs_warning << " k vector is not correct" << endl;
			ofs_warning << " Read in kx=" << kx << " ky = " << ky << " kz = " << kz << endl;
			ofs_warning << " In fact, kx=" << kv.kvec_c[ik].x 
			 << " ky=" << kv.kvec_c[ik].y
			 << " kz=" << kv.kvec_c[ik].z << endl;
			 error = 4; 
		}
        else if (nbands!=NBANDS)
        {
            ofs_warning << " read in nbands=" << nbands;
            ofs_warning << " NBANDS=" << NBANDS << endl;
            error = 2;
        }
        else if (nlocal != NLOCAL)
        {
            ofs_warning << " read in nlocal=" << nlocal;
            ofs_warning << " NLOCAL=" << NLOCAL << endl;
            error = 3;
        }

        ctot = new complex<double>*[NBANDS];
        for (int i=0; i<NBANDS; i++)
        {
            ctot[i] = new complex<double>[NLOCAL];
        }

        for (int i=0; i<NBANDS; ++i)
        {
            int ib;
            READ_VALUE(ifs, ib);
			ib -= 1; // because in C++, ib should start from 0
			//------------------------------------------------
			// read the eigenvalues!
			// very important to determine the occupations.
			//------------------------------------------------
			READ_VALUE(ifs, wf.ekb[ik][ib]);
			READ_VALUE(ifs, wf.wg(ik,ib));
            assert( i==ib );
			double a, b;
            for (int j=0; j<NLOCAL; ++j)
            {
                ifs >> a >> b;
				ctot[i][j]=complex<double>(a,b);
				//cout << ctot[i][j] << " " << endl;
            }
        }
    }

#ifdef __MPI
    Parallel_Common::bcast_int(error);
#endif
	if(error==2) return 2;
	if(error==3) return 3;
	if(error==4) return 4;

	// mohan add 2012-02-15,
	// SGO.cal_totwfc();

	// distri_lowf need all processors.
	// otherwise, read in sucessfully.
    // if DRANK!=0, ctot is not used,
    // so it's save.
	
    //WF_Local::distri_lowf(ctot, SGO.totwfc[0]);
    if(newdm==0)
	{
		WF_Local::distri_lowf_complex(ctot, c); 
	}
	else
	{	
		WF_Local::distri_lowf_complex_new(ctot, ik);
	}
	
	// mohan add 2012-02-15,
	// still have bugs, but can solve it later.
	// distribute the wave functions again.
	// SGO.dis_subwfc();

    if (DRANK==0)
    {
        // delte the ctot
        for (int i=0; i<NBANDS; i++)
        {
            delete[] ctot[i];
        }
        delete[] ctot;
    }

	//---------
	// TEST
	//---------
	/*
	for(int i=0; i<NBANDS; ++i)
	{
		cout << " c band i=" << i+1 << endl;
		for(int j=0; j<NLOCAL; ++j)
		{
			cout << " " << c[i][j];
		}
		cout << endl;
	}
	*/


    timer::tick("WF_Local","read_lowf_complex");
	return 0;
}

int WF_Local::read_lowf(double **c, const int &is)
{
    TITLE("WF_Local","read_lowf");
    timer::tick("WF_Local","read_lowf");

    double **ctot;

    stringstream ss;
	if(GAMMA_ONLY_LOCAL)
	{
		// read wave functions
		// write is in ../src_pdiag/pdiag_basic.cpp
    	ss << global_out_dir << "LOWF_GAMMA_S" << is+1 <<".dat";
		cout << " name is = " << ss.str() << endl;
	}
	else
	{
		ss << global_out_dir << "LOWF_K.dat";
	}

    ifstream ifs;

    int error = 0;

    if (DRANK==0)
    {
        ifs.open(ss.str().c_str());
        if (!ifs)
        {
            ofs_warning << " Can't open file:" << ss.str() << endl;
            error = 1;
        }
    }

#ifdef __MPI
    Parallel_Common::bcast_int(error);
#endif

    if (error==1) return 1;

    // otherwise, find the file.

    if (MY_RANK==0)
    {
        int nbands, nlocal;
        READ_VALUE(ifs, nbands);
        READ_VALUE(ifs, nlocal);

        if (nbands!=NBANDS)
        {
            ofs_warning << " read in nbands=" << nbands;
            ofs_warning << " NBANDS=" << NBANDS << endl;
            error = 2;
        }
        else if (nlocal != NLOCAL)
        {
            ofs_warning << " read in nlocal=" << nlocal;
            ofs_warning << " NLOCAL=" << NLOCAL << endl;
            error = 3;
        }

        ctot = new double*[NBANDS];
        for (int i=0; i<NBANDS; i++)
        {
            ctot[i] = new double[NLOCAL];
        }

		//cout << "nbands" << NBANDS << endl;
        for (int i=0; i<NBANDS; i++)
        {
            int ib;
            READ_VALUE(ifs, ib);
			READ_VALUE(ifs, wf.ekb[CURRENT_SPIN][i]);
			READ_VALUE(ifs, wf.wg(CURRENT_SPIN,i));
            assert( (i+1)==ib);
			//cout << " ib=" << ib << endl;
            for (int j=0; j<NLOCAL; j++)
            {
                ifs >> ctot[i][j];
				//cout << ctot[i][j] << " ";
            }
			//cout << endl;
        }
    }


#ifdef __MPI
    Parallel_Common::bcast_int(error);
	Parallel_Common::bcast_double( wf.ekb[is], NBANDS);
	Parallel_Common::bcast_double( wf.wg.c, NSPIN*NBANDS);
#endif
	if(error==2) return 2;
	if(error==3) return 3;

	// mohan add 2012-02-15,
	// mohan comment out 2021-02-09
	//SGO.cal_totwfc();

	// distri_lowf need all processors.
	// otherwise, read in sucessfully.
    // if DRANK!=0, ctot is not used,
    // so it's save.

	if(INPUT.new_dm==0)
	{	
		WF_Local::distri_lowf(ctot, SGO.totwfc[0]);
	}
	else
	{
		WF_Local::distri_lowf_new(ctot, is);
	}
	
	// mohan add 2012-02-15,
	// still have bugs, but can solve it later.
	// distribute the wave functions again.
	// mohan comment out 2021-02-09
	// SGO.dis_subwfc();

    if (MY_RANK==0)
    {
        // delte the ctot
        for (int i=0; i<NBANDS; i++)
        {
            delete[] ctot[i];
        }
        delete[] ctot;
    }

    timer::tick("WF_Local","read_lowf");
    return 0;
}

void WF_Local::write_lowf(const string &name, double **ctot)
{
    TITLE("WF_Local","write_lowf");
    timer::tick("WF_Local","write_lowf");

    ofstream ofs;
    if (DRANK==0)
    {
        ofs.open(name.c_str());
        if (!ofs)
        {
            WARNING("Pdiag_Basic::write_lowf","Can't write local orbital wave functions.");
        }
        ofs << NBANDS << " (number of bands)" << endl;
        ofs << NLOCAL << " (number of orbitals)";
        ofs << setprecision(8);
        ofs << scientific;

        for (int i=0; i<NBANDS; i++)
        {
            // +1 to mean more clearly.
            // band index start from 1.
            ofs << "\n" << i+1 << " (band)";
			ofs << "\n" << wf.ekb[CURRENT_SPIN][i] << " (Ry)"; //mohan add 2012-03-26
			ofs << "\n" << wf.wg(CURRENT_SPIN,i) << " (Occupations)";
            for (int j=0; j<NLOCAL; j++)
            {
                if (j % 5 == 0) ofs << "\n";
                ofs << ctot[i][j] << " ";
            }
        }
        ofs.close();
    }

    timer::tick("WF_Local","write_lowf");
    return;
}

void WF_Local::write_lowf_complex(const string &name, complex<double> **ctot, const int &ik)
{
    TITLE("WF_Local","write_lowf_complex");
    timer::tick("WF_Local","write_lowf_complex");

    ofstream ofs;
    if (DRANK==0)
    {
        ofs.open(name.c_str());
        if (!ofs)
        {
            WARNING("Pdiag_Basic::write_lowf","Can't write local orbital wave functions.");
        }
        ofs << setprecision(25);
		ofs << ik+1 << " (index of k points)" << endl;
		ofs << kv.kvec_c[ik].x << " " << kv.kvec_c[ik].y << " " << kv.kvec_c[ik].z << endl;
        ofs << NBANDS << " (number of bands)" << endl;
        ofs << NLOCAL << " (number of orbitals)";
        ofs << scientific;

        for (int i=0; i<NBANDS; i++)
        {
            // +1 to mean more clearly.
            // band index start from 1.
            ofs << "\n" << i+1 << " (band)";
			ofs << "\n" << wf.ekb[ik][i] << " (Ry)";
			ofs << "\n" << wf.wg(ik,i) << " (Occupations)";
            for (int j=0; j<NLOCAL; j++)
            {
                if (j % 5 == 0) ofs << "\n";
                ofs << ctot[i][j].real() << " " << ctot[i][j].imag() << " ";
            }
        }
        ofs.close();
    }

    timer::tick("WF_Local","write_lowf_complex");
    return;
}

void WF_Local::distri_lowf_new(double **ctot, const int &is)
{
    TITLE("WF_Local","distri_lowf_new");
#ifdef __MPI

//1. alloc work array; set some parameters

	long maxnloc; // maximum number of elements in local matrix
	MPI_Reduce(&ParaO.nloc, &maxnloc, 1, MPI_LONG, MPI_MAX, 0, ParaO.comm_2D);
	MPI_Bcast(&maxnloc, 1, MPI_LONG, 0, ParaO.comm_2D);
	//reduce and bcast could be replaced by allreduce
	
    int nprocs, myid;
    MPI_Comm_size(ParaO.comm_2D, &nprocs);
    MPI_Comm_rank(ParaO.comm_2D, &myid);

	double *work=new double[maxnloc]; // work/buffer matrix
	int nb = 0;
	if(NB2D==0)
	{
		if(NLOCAL>0) nb = 1;
		if(NLOCAL>500) nb = 32;
		if(NLOCAL>1000) nb = 64;
	}
	else if(NB2D>0)
	{
		nb = NB2D; // mohan add 2010-06-28
	}
	int info;
	int naroc[2]; // maximum number of row or column
	
//2. copy from ctot to wfc_gamma
	for(int iprow=0; iprow<ParaO.dim0; ++iprow)
	{
		for(int ipcol=0; ipcol<ParaO.dim1; ++ipcol)
		{
//2.1 get and bcast local 2d matrix info
			const int coord[2]={iprow, ipcol};
			int src_rank;
			MPI_Cart_rank(ParaO.comm_2D, coord, &src_rank);
			if(myid==src_rank)
			{
				naroc[0]=ParaO.nrow;
				naroc[1]=ParaO.ncol;
			}
			info=MPI_Bcast(naroc, 2, MPI_INT, src_rank, ParaO.comm_2D);

//2.2 copy from ctot to work, then bcast work
			info=CTOT2q(myid, naroc, nb, ParaO.dim0, ParaO.dim1, iprow, ipcol, work, ctot);
			info=MPI_Bcast(work, maxnloc, MPI_DOUBLE, 0, ParaO.comm_2D);
			//ofs_running << "iprow, ipcow : " << iprow << ipcol << endl;
			//for (int i=0; i<maxnloc; ++i)
			//{
				//ofs_running << *(work+i)<<" ";
			//}
			//ofs_running << endl;
//2.3 copy from work to wfc_gamma
			const int inc=1;
			if(myid==src_rank)
			{
				LapackConnector::copy(ParaO.nloc, work, inc, LOC.wfc_dm_2d.wfc_gamma[is].c, inc);
			}
		}//loop ipcol
	}//loop	iprow

	delete[] work;
#else
	WARNING_QUIT("WF_Local::distri_lowf_new","check the code without MPI.");
#endif
    return;
}

void WF_Local::distri_lowf_complex_new(complex<double> **ctot, const int &ik)
{
    TITLE("WF_Local","distri_lowf_complex_new");
#ifdef __MPI

//1. alloc work array; set some parameters

	long maxnloc; // maximum number of elements in local matrix
	MPI_Reduce(&ParaO.nloc, &maxnloc, 1, MPI_LONG, MPI_MAX, 0, ParaO.comm_2D);
	MPI_Bcast(&maxnloc, 1, MPI_LONG, 0, ParaO.comm_2D);
	//reduce and bcast could be replaced by allreduce
	
    int nprocs, myid;
    MPI_Comm_size(ParaO.comm_2D, &nprocs);
    MPI_Comm_rank(ParaO.comm_2D, &myid);

	complex<double> *work=new complex<double>[maxnloc]; // work/buffer matrix
	int nb = 0;
	if(NB2D==0)
	{
		if(NLOCAL>0) nb = 1;
		if(NLOCAL>500) nb = 32;
		if(NLOCAL>1000) nb = 64;
	}
	else if(NB2D>0)
	{
		nb = NB2D; // mohan add 2010-06-28
	}
	int info;
	int naroc[2]; // maximum number of row or column
	
//2. copy from ctot to wfc_gamma
	for(int iprow=0; iprow<ParaO.dim0; ++iprow)
	{
		for(int ipcol=0; ipcol<ParaO.dim1; ++ipcol)
		{
//2.1 get and bcast local 2d matrix info
			const int coord[2]={iprow, ipcol};
			int src_rank;
			MPI_Cart_rank(ParaO.comm_2D, coord, &src_rank);
			if(myid==src_rank)
			{
				naroc[0]=ParaO.nrow;
				naroc[1]=ParaO.ncol;
			}
			info=MPI_Bcast(naroc, 2, MPI_INT, src_rank, ParaO.comm_2D);

//2.2 copy from ctot to work, then bcast work
			info=CTOT2q_c(myid, naroc, nb, ParaO.dim0, ParaO.dim1, iprow, ipcol, work, ctot);
			info=MPI_Bcast(work, maxnloc, MPI_DOUBLE_COMPLEX, 0, ParaO.comm_2D);
			//ofs_running << "iprow, ipcow : " << iprow << ipcol << endl;
			//for (int i=0; i<maxnloc; ++i)
			//{
				//ofs_running << *(work+i)<<" ";
			//}
			//ofs_running << endl;
//2.3 copy from work to wfc_k
			const int inc=1;
			if(myid==src_rank)
			{
				LapackConnector::copy(ParaO.nloc, work, inc, LOC.wfc_dm_2d.wfc_k[ik].c, inc);
			}
		}//loop ipcol
	}//loop	iprow

	delete[] work;
#else
	WARNING_QUIT("WF_Local::distri_lowf_new","check the code without MPI.");
#endif
    return;
}

void WF_Local::distri_lowf(double **ctot, double **c)
{
    TITLE("WF_Local","distri_lowf");
#ifdef __MPI

    MPI_Status status;
    for (int i=0; i<DSIZE; i++)
    {
        if (DRANK==0)
        {
            if (i==0)
            {
                // get the wave functions from 'ctot',
                // save them in the matrix 'c'.
                for (int iw=0; iw<NLOCAL; iw++)
                {
					// mohan update 2012-01-12
//                  const int mu_local = GridT.trace_lo[iw]; 
                    const int mu_local = SGO.trace_lo_tot[iw];

                    if (mu_local >= 0)
                    {
                        for (int ib=0; ib<NBANDS; ib++)
                        {
                            c[ib][mu_local] = ctot[ib][iw];
                        }
                    }
                }
            }
            else
            {
                int tag;
                // receive trace_lo2
                tag = i * 3;
                int* trace_lo2 = new int[NLOCAL];
                MPI_Recv(trace_lo2, NLOCAL, MPI_INT, i, tag, DIAG_WORLD, &status);

                // receive lgd2
                int lgd2 = 0;
                tag = i * 3 + 1;
                MPI_Recv(&lgd2, 1, MPI_INT, i, tag, DIAG_WORLD, &status);

                // send csend
                double* csend = new double[NBANDS*lgd2];
                ZEROS(csend, NBANDS*lgd2);

                for (int ib=0; ib<NBANDS; ib++)
                {
                    for (int iw=0; iw<NLOCAL; iw++)
                    {
                        const int mu_local = trace_lo2[iw];
                        if (mu_local>=0)
                        {
                            csend[mu_local*NBANDS+ib] = ctot[ib][iw];
                        }
                    }
                }

                tag = i * 3 + 2;
                MPI_Send(csend,NBANDS*lgd2,MPI_DOUBLE,i,tag,DIAG_WORLD);

                delete[] trace_lo2;
                delete[] csend;
            }
        }// end DRANK=0
        else if ( i == DRANK)
        {
            int tag;

            // send trace_lo
            tag = DRANK * 3;
			// mohan update 2012-01-12
            //MPI_Send(GridT.trace_lo, NLOCAL, MPI_INT, 0, tag, DIAG_WORLD);
            MPI_Send(SGO.trace_lo_tot, NLOCAL, MPI_INT, 0, tag, DIAG_WORLD);

            // send lgd
            tag = DRANK * 3 + 1;

			// mohan update 2012-01-12
			int lgdnow = SGO.lgd;
            MPI_Send(&lgdnow, 1, MPI_INT, 0, tag, DIAG_WORLD);

            // receive c
			ofs_running << " lgdnow=" << lgdnow << endl;
            double* crecv = new double[NBANDS*lgdnow];
            ZEROS(crecv, NBANDS*lgdnow);
            tag = DRANK * 3 + 2;
            MPI_Recv(crecv, NBANDS*lgdnow, MPI_DOUBLE, 0, tag, DIAG_WORLD, &status);

            for (int ib=0; ib<NBANDS; ib++)
            {
                for (int mu=0; mu<lgdnow; mu++)
                {
                    c[ib][mu] = crecv[mu*NBANDS+ib];
                }
            }

            delete[] crecv;
        }// end i==DRANK
        MPI_Barrier(DIAG_WORLD);
    }// end i

    //-----------
    // for test,
    //-----------
    /*
    ofs_running << " Wave Functions in local basis: " << endl;
    for(int i=0; i<NBANDS; i++)
    {
        for(int j=0; j<GridT.lgd; j++)
        {
            if(j%8==0) ofs_running << endl;
            if( abs(c[i][j]) > 1.0e-5  )
            {
                ofs_running << setw(15) << c[i][j];
            }
            else
            {
                ofs_running << setw(15) << "0";
            }
        }
    }
    ofs_running << endl;
    */
#else
	WARNING_QUIT("WF_Local::distri_lowf","check the code without MPI.");
#endif
    return;
}


void WF_Local::distri_lowf_complex(complex<double> **ctot, complex<double> **cc)
{
    TITLE("WF_Local","distri_lowf_complex");
#ifdef __MPI

    MPI_Status status;

    for (int i=0; i<DSIZE; i++)
    {
        if (DRANK==0)
        {
            if (i==0)
            {
                // get the wave functions from 'ctot',
                // save them in the matrix 'c'.
                for (int iw=0; iw<NLOCAL; iw++)
                {
                    const int mu_local = GridT.trace_lo[iw];
                    if (mu_local >= 0)
                    {
                        for (int ib=0; ib<NBANDS; ib++)
                        {
                            cc[ib][mu_local] = ctot[ib][iw];
                        }
                    }
                }
            }
            else
            {
				int tag;
                // receive lgd2
                int lgd2 = 0;
                tag = i * 3;
                MPI_Recv(&lgd2, 1, MPI_INT, i, tag, DIAG_WORLD, &status);
				if(lgd2==0)
				{

				}
				else
				{
					// receive trace_lo2
					tag = i * 3 + 1;
					int* trace_lo2 = new int[NLOCAL];
					MPI_Recv(trace_lo2, NLOCAL, MPI_INT, i, tag, DIAG_WORLD, &status);

//					ofs_running << " lgd2=" << lgd2 << " proc=" << i+1 << endl;
					// send csend
					complex<double>* csend = new complex<double>[NBANDS*lgd2];
					ZEROS(csend, NBANDS*lgd2);
					for (int ib=0; ib<NBANDS; ib++)
					{
						for (int iw=0; iw<NLOCAL; iw++)
						{
							const int mu_local = trace_lo2[iw];
							if (mu_local>=0)
							{
								csend[mu_local*NBANDS+ib] = ctot[ib][iw];
							}
						}
					}
					tag = i * 3 + 2;
					MPI_Send(csend,NBANDS*lgd2,mpicomplex,i,tag,DIAG_WORLD);
                	delete[] csend;
                	delete[] trace_lo2;
				}
            }
        }// end DRANK=0
        else if ( i == DRANK)
		{
			int tag;

			// send GridT.lgd
			tag = DRANK * 3;
			MPI_Send(&GridT.lgd, 1, MPI_INT, 0, tag, DIAG_WORLD);

			if(GridT.lgd != 0)
			{
				// send trace_lo
				tag = DRANK * 3 + 1;
				MPI_Send(GridT.trace_lo, NLOCAL, MPI_INT, 0, tag, DIAG_WORLD);

				// receive cc
				complex<double>* crecv = new complex<double>[NBANDS*GridT.lgd];
				ZEROS(crecv, NBANDS*GridT.lgd);

				tag = DRANK * 3 + 2;
				MPI_Recv(crecv, NBANDS*GridT.lgd, mpicomplex, 0, tag, DIAG_WORLD, &status);

				for (int ib=0; ib<NBANDS; ib++)
				{
					for (int mu=0; mu<GridT.lgd; mu++)
					{
						cc[ib][mu] = crecv[mu*NBANDS+ib];
					}
				}

				delete[] crecv;

			}
        }// end i==DRANK
        MPI_Barrier(DIAG_WORLD);
    }// end i

    //-----------
    // for test,
    //-----------
    /*
    ofs_running << " Wave Functions in local basis: " << endl;
    for(int i=0; i<NBANDS; i++)
    {
        for(int j=0; j<GridT.lgd; j++)
        {
            if(j%8==0) ofs_running << endl;
            if( abs(c[i][j]) > 1.0e-5  )
            {
                ofs_running << setw(15) << c[i][j];
            }
            else
            {
                ofs_running << setw(15) << "0";
            }
        }
    }
    ofs_running << endl;
    */
#else
	WARNING_QUIT("WF_Local::distri_lowf_complex","check the code without MPI.");
#endif
    return;
}


void WF_Local::distri_lowf_aug(double **ctot, double **c_aug)
{
    TITLE("WF_Local","distri_lowf_aug");

#ifdef __MPI
    MPI_Status status;

	//-------------------------------------------------
	// Do the distribution of augmented wave functions
	// for each of the processors
	//-------------------------------------------------
    for (int i=0; i<DSIZE; i++)
    {
		//----------------------------------------
		// Rank 0 processor, very special.
		//----------------------------------------
        if (DRANK==0)
        {
            if (i==0)
            {
				//-------------------------------------
                // get the wave functions from 'ctot',
                // save them in the matrix 'c_aug'.
				//-------------------------------------
                for (int iw=0; iw<NLOCAL; iw++)
                {
					const int mu = LOWF.trace_aug[iw];
					if( mu < 0 ) continue;
					for (int ib=0; ib<NBANDS; ib++)
					{
						c_aug[ib][mu] = ctot[ib][iw];
					}
				}
			}
            else
            {
				// send the needed data to processor 'i'
                int tag;
                // (1) receive 'trace_lo2'
                tag = i * 3;
                int* tmp_trace = new int[NLOCAL];
                MPI_Recv(tmp_trace, NLOCAL, MPI_INT, i, tag, DIAG_WORLD, &status);

//				ofs_running << " Recv trace from pro " << i << endl; 

                // (2) receive daug
                tag = i * 3 + 1;
                int daug = 0;
                MPI_Recv(&daug, 1, MPI_INT, i, tag, DIAG_WORLD, &status);

//				ofs_running << " Recv daug from pro " << i << endl;

                // ready to send csend, data number is NBANDS*daug
                double* csend = new double[NBANDS*daug];
                ZEROS(csend, NBANDS*daug);

                for (int ib=0; ib<NBANDS; ib++)
                {
                    for (int iw=0; iw<NLOCAL; iw++)
                    {
                        const int mu = tmp_trace[iw];
                        if (mu>=0)
                        {
                            csend[mu*NBANDS+ib] = ctot[ib][iw];
                        }
                    }
                }

				// (3) send the data to processor i.
                tag = i * 3 + 2;
                MPI_Send(csend,NBANDS*daug,MPI_DOUBLE,i,tag,DIAG_WORLD);

//				ofs_running << " Send ctot to pro " << i << endl;

                delete[] tmp_trace;
                delete[] csend;
            }
        }// end DRANK=0
		//------------------------------
		// other processors operations
		//------------------------------
        else if ( i == DRANK)
        {

            int tag;

			// mohan add 2010-09-26
			//--------------------------------------------
            // (1) LOWF.trace_aug : trace c_aug
			// tell processor 0 that which wave functions
			// they need
			//--------------------------------------------
            tag = DRANK * 3;
            MPI_Send(LOWF.trace_aug, NLOCAL, MPI_INT, 0, tag, DIAG_WORLD);

			//--------------------------------------------
			// (2) LOWF.daug: dimension of c_aug
			// and RANK=0 will receive the information.
			//--------------------------------------------
            tag = DRANK * 3 + 1;
            MPI_Send(&LOWF.daug, 1, MPI_INT, 0, tag, DIAG_WORLD);

			//--------------------------------------------
            // (3) receive the augmented wave functions
			//--------------------------------------------
            double* crecv = new double[NBANDS*LOWF.daug];
            ZEROS(crecv, NBANDS*LOWF.daug);

            tag = DRANK * 3 + 2;
            MPI_Recv(crecv, NBANDS*LOWF.daug, MPI_DOUBLE, 0, tag, DIAG_WORLD, &status);

			//--------------------------------------------
            // (4) copy the augmented wave functions
			//--------------------------------------------
            for (int ib=0; ib<NBANDS; ib++)
            {
                for (int mu=0; mu<LOWF.daug; mu++)
                {
                    c_aug[ib][mu] = crecv[mu*NBANDS+ib];
                }
            }

			//=================================================
			// mohan fix bug 2011-04-07 perfectly!!!!!!!!!!
			// (1) first notice that memory increasing for large system,
			// but not Si2 dimer;
			// (2) second calculting the leap memory for each iteration,
			// about 0.1*16G ~ 16 MB for Graphene QD 612 atoms,
			// so i choose the one closely ~~~~~~,
			// and found WFC_GAMMA_AUG might be the target.
			// (3) third I notice when I calculate 1600 atoms,
			// the DRANK==0 processor only 7% memory
			// while others as large as 25% memory.
			// it's really hard to find here.
			//=================================================
			delete[] crecv;
        }// end i==DRANK
        MPI_Barrier(DIAG_WORLD);
    }// end i

    //-----------
    // for test,
    //-----------
    /*
    ofs_running << " Wave Functions in local basis: " << endl;
    for(int i=0; i<NBANDS; i++)
    {
        for(int j=0; j<LOWF.daug; j++)
        {
            if(j%8==0) ofs_running << endl;
            if( abs(c[i][j]) > 1.0e-5  )
            {
                ofs_running << setw(15) << c[i][j];
            }
            else
            {
                ofs_running << setw(15) << "0";
            }
        }
    }
    ofs_running << endl;
    */
#else
	WARNING_QUIT("WF_Local::distri_lowf_aug","check code without MPI.");
#endif
    return;
}


void WF_Local::distri_lowf_aug_complex(complex<double> **ctot, complex<double> **c_aug)
{
    TITLE("WF_Local","distri_lowf_aug_complex");

#ifdef __MPI
    MPI_Status status;

	//-------------------------------------------------
	// Do the distribution of augmented wave functions
	// for each of the processors in Diag group
	//-------------------------------------------------
    for (int i=0; i<DSIZE; i++)
    {
		//----------------------------------------
		// Rank 0 processor, very special.
		//----------------------------------------
        if (DRANK==0)
        {
            if (i==0)
            {
				//-------------------------------------
                // get the wave functions from 'ctot',
                // save them in the matrix 'c_aug'.
				//-------------------------------------
                for (int iw=0; iw<NLOCAL; iw++)
                {
					const int mu = LOWF.trace_aug[iw];
					if( mu < 0 ) continue;
					for (int ib=0; ib<NBANDS; ib++)
					{
						c_aug[ib][mu] = ctot[ib][iw];
					}
				}
			}
            else
            {
				// send the needed data to processor 'i'
                int tag;
                // (1) receive 'trace_lo2'
                tag = i * 3;
                int* tmp_trace = new int[NLOCAL];
                MPI_Recv(tmp_trace, NLOCAL, MPI_INT, i, tag, DIAG_WORLD, &status);

                // (2) receive daug
                tag = i * 3 + 1;
                int daug = 0;
                MPI_Recv(&daug, 1, MPI_INT, i, tag, DIAG_WORLD, &status);

				if(daug!=0)
				{
					// ready to send csend, data number is NBANDS*daug
					complex<double>* csend = new complex<double>[NBANDS*daug];
					ZEROS(csend, NBANDS*daug);

					for (int ib=0; ib<NBANDS; ib++)
					{
						for (int iw=0; iw<NLOCAL; iw++)
						{
							const int mu = tmp_trace[iw];
							if (mu>=0)
							{
								csend[mu*NBANDS+ib] = ctot[ib][iw];
							}
						}
					}

					// (3) send the data to processor i.
					tag = i * 3 + 2;
					MPI_Send(csend,NBANDS*daug,mpicomplex,i,tag,DIAG_WORLD);
					//				ofs_running << " Send ctot to pro " << i << endl;
					delete[] csend;
				}
                delete[] tmp_trace;
            }
        }// end DRANK=0
		//------------------------------
		// other processors operations
		//------------------------------
        else if ( i == DRANK)
        {
            int tag;

			// mohan add 2010-09-26
			//--------------------------------------------
            // (1) LOWF.trace_aug : trace c_aug
			// tell processor 0 that which wave functions
			// they need
			//--------------------------------------------
            tag = DRANK * 3;
            MPI_Send(LOWF.trace_aug, NLOCAL, MPI_INT, 0, tag, DIAG_WORLD);
//			for(int i=0; i<NLOCAL; ++i)
//			{
//				ofs_running << " trace_aug = " << LOWF.trace_aug[i] << endl; 
//			}

			//--------------------------------------------
			// (2) LOWF.daug: dimension of c_aug
			// and RANK=0 will receive the information.
			//--------------------------------------------
            tag = DRANK * 3 + 1;
            MPI_Send(&LOWF.daug, 1, MPI_INT, 0, tag, DIAG_WORLD);

			if(LOWF.daug!=0)
			{
				//--------------------------------------------
				// (3) receive the augmented wave functions
				//--------------------------------------------
				complex<double>* crecv = new complex<double>[NBANDS*LOWF.daug];
				ZEROS(crecv, NBANDS*LOWF.daug);

				tag = DRANK * 3 + 2;
				MPI_Recv(crecv, NBANDS*LOWF.daug, mpicomplex, 0, tag, DIAG_WORLD, &status);

				//--------------------------------------------
				// (4) copy the augmented wave functions
				//--------------------------------------------
				for (int ib=0; ib<NBANDS; ib++)
				{
					for (int mu=0; mu<LOWF.daug; mu++)
					{
						c_aug[ib][mu] = crecv[mu*NBANDS+ib];
					}
				}

				//=================================================
				// mohan fix bug 2011-04-07 perfectly!!!!!!!!!!
				// (1) first notice that memory increasing for large system,
				// but not Si2 dimer;
				// (2) second calculting the leap memory for each iteration,
				// about 0.1*16G ~ 16 MB for Graphene QD 612 atoms,
				// so i choose the one closely ~~~~~~,
				// and found WFC_GAMMA_AUG might be the target.
				// (3) third I notice when I calculate 1600 atoms,
				// the DRANK==0 processor only 7% memory
				// while others as large as 25% memory.
				// it's really hard to find here.
				//=================================================
				delete[] crecv;
			}
        }// end i==DRANK
        MPI_Barrier(DIAG_WORLD);
    }// end i

    //-----------
    // for test,
    //-----------
	/*
    ofs_running << " Wave Functions in local basis: " << endl;
    for(int i=0; i<NBANDS; i++)
    {
        for(int j=0; j<LOWF.daug; j++)
        {
            if(j%8==0) ofs_running << endl;
            if( abs(c_aug[i][j]) > 1.0e-5  )
            {
                ofs_running << setw(15) << c_aug[i][j];
            }
            else
            {
                ofs_running << setw(15) << "0";
            }
        }
    }
    ofs_running << endl;
	*/
#else
	WARNING_QUIT("WF_Local::distri_lowf_aug_complex","check code without MPI.");
#endif
    return;
}
