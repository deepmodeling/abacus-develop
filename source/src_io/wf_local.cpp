#include "wf_local.h"
#include "../src_pw/global.h"
#include "../src_parallel/parallel_common.h"
#include "../module_base/timer.h"

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
        int igcol=Local_Orbital_wfc::globalIndex(j, nb, dim1, ipcol);
        if(igcol>=GlobalV::NBANDS) continue;
        for(int i=0; i<naroc[0]; ++i)
        {
            int igrow=Local_Orbital_wfc::globalIndex(i, nb, dim0, iprow);
			//GlobalV::ofs_running << "i,j,igcol,igrow" << i<<" "<<j<<" "<<igcol<<" "<<igrow<<std::endl;
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
	std::complex<double>* work,
	std::complex<double>** CTOT)
{
    for(int j=0; j<naroc[1]; ++j)
    {
        int igcol=Local_Orbital_wfc::globalIndex(j, nb, dim1, ipcol);
        if(igcol>=GlobalV::NBANDS) continue;
        for(int i=0; i<naroc[0]; ++i)
        {
            int igrow=Local_Orbital_wfc::globalIndex(i, nb, dim0, iprow);
			//ofs_running << "i,j,igcol,igrow" << i<<" "<<j<<" "<<igcol<<" "<<igrow<<std::endl;
            if(myid==0) work[j*naroc[0]+i]=CTOT[igcol][igrow];
        }
    }
    return 0;
}

// be called in local_orbital_wfc::allocate_k
int WF_Local::read_lowf_complex(std::complex<double>** ctot, const int& ik, 
    Local_Orbital_wfc &lowf)
{
    ModuleBase::TITLE("WF_Local","read_lowf_complex");
    ModuleBase::timer::tick("WF_Local","read_lowf_complex");

    lowf.wfc_k[ik].create(lowf.ParaV->ncol_bands, lowf.ParaV->nrow);
    std::stringstream ss;
	// read wave functions
	// write is in ../src_pdiag/pdiag_basic.cpp
    ss << GlobalV::global_out_dir << "LOWF_K_" << ik+1 <<".dat";
//	std::cout << " name is = " << ss.str() << std::endl;

    std::ifstream ifs;

    int error = 0;

    if (GlobalV::DRANK==0)
    {
        ifs.open(ss.str().c_str());
        if (!ifs)
        {
            GlobalV::ofs_warning << " Can't open file:" << ss.str() << std::endl;
            error = 1;
        }
    }

#ifdef __MPI
    Parallel_Common::bcast_int(error);
#endif

    if (error==1) return 1;

    // otherwise, find the file.

    if (GlobalV::MY_RANK==0)
    {
		int ikr;
		double kx,ky,kz;
        int nbands, nlocal;
		ModuleBase::GlobalFunc::READ_VALUE(ifs, ikr);
		ifs >> kx >> ky >> kz;
        ModuleBase::GlobalFunc::READ_VALUE(ifs, nbands);
        ModuleBase::GlobalFunc::READ_VALUE(ifs, nlocal);

		if(ikr!=ik+1)
		{
			GlobalV::ofs_warning << " ikr=" << ikr << " ik=" << ik << std::endl;
			GlobalV::ofs_warning << " k index is not correct" << std::endl;
			error = 4;
		}
		else if ( 
			abs(kx-GlobalC::kv.kvec_c[ik].x)>1.0e-5 ||
			abs(ky-GlobalC::kv.kvec_c[ik].y)>1.0e-5 ||
			abs(kz-GlobalC::kv.kvec_c[ik].z)>1.0e-5 )
		{	
			GlobalV::ofs_warning << " k std::vector is not correct" << std::endl;
			GlobalV::ofs_warning << " Read in kx=" << kx << " ky = " << ky << " kz = " << kz << std::endl;
			GlobalV::ofs_warning << " In fact, kx=" << GlobalC::kv.kvec_c[ik].x 
			 << " ky=" << GlobalC::kv.kvec_c[ik].y
			 << " kz=" << GlobalC::kv.kvec_c[ik].z << std::endl;
			 error = 4; 
		}
        else if (nbands!=GlobalV::NBANDS)
        {
            GlobalV::ofs_warning << " read in nbands=" << nbands;
            GlobalV::ofs_warning << " NBANDS=" << GlobalV::NBANDS << std::endl;
            error = 2;
        }
        else if (nlocal != GlobalV::NLOCAL)
        {
            GlobalV::ofs_warning << " read in nlocal=" << nlocal;
            GlobalV::ofs_warning << " NLOCAL=" << GlobalV::NLOCAL << std::endl;
            error = 3;
        }

        ctot = new std::complex<double>*[GlobalV::NBANDS];
        for (int i=0; i<GlobalV::NBANDS; i++)
        {
            ctot[i] = new std::complex<double>[GlobalV::NLOCAL];
        }

        for (int i=0; i<GlobalV::NBANDS; ++i)
        {
            int ib;
            ModuleBase::GlobalFunc::READ_VALUE(ifs, ib);
			ib -= 1; // because in C++, ib should start from 0
			//------------------------------------------------
			// read the eigenvalues!
			// very important to determine the occupations.
			//------------------------------------------------
			ModuleBase::GlobalFunc::READ_VALUE(ifs, GlobalC::wf.ekb[ik][ib]);
			ModuleBase::GlobalFunc::READ_VALUE(ifs, GlobalC::wf.wg(ik,ib));
            assert( i==ib );
			double a, b;
            for (int j=0; j<GlobalV::NLOCAL; ++j)
            {
                ifs >> a >> b;
				ctot[i][j]=std::complex<double>(a,b);
				//std::cout << ctot[i][j] << " " << std::endl;
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
	// GlobalC::SGO.cal_totwfc();

	// distri_lowf need all processors.
	// otherwise, read in sucessfully.
    // if GlobalV::DRANK!=0, ctot is not used,
    // so it's save.
	
    //WF_Local::distri_lowf(ctot, GlobalC::SGO.totwfc[0]);
	WF_Local::distri_lowf_complex_new(ctot, ik, lowf);
	
	// mohan add 2012-02-15,
	// still have bugs, but can solve it later.
	// distribute the wave functions again.
	// GlobalC::SGO.dis_subwfc();

    if (GlobalV::DRANK==0)
    {
        // delte the ctot
        for (int i=0; i<GlobalV::NBANDS; i++)
        {
            delete[] ctot[i];
        }
        delete[] ctot;
    }

	//---------
	// TEST
	//---------
	/*
	for(int i=0; i<GlobalV::NBANDS; ++i)
	{
		std::cout << " c band i=" << i+1 << std::endl;
		for(int j=0; j<GlobalV::NLOCAL; ++j)
		{
			std::cout << " " << c[i][j];
		}
		std::cout << std::endl;
	}
	*/


    ModuleBase::timer::tick("WF_Local","read_lowf_complex");
	return 0;
}

int WF_Local::read_lowf(double** ctot, const int& is,
    Local_Orbital_wfc &lowf)
{
    ModuleBase::TITLE("WF_Local","read_lowf");
    ModuleBase::timer::tick("WF_Local", "read_lowf");
    
    std::stringstream ss;
	if(GlobalV::GAMMA_ONLY_LOCAL)
	{
		// read wave functions
		// write is in ../src_pdiag/pdiag_basic.cpp
    	ss << GlobalV::global_out_dir << "LOWF_GAMMA_S" << is+1 <<".dat";
		std::cout << " name is = " << ss.str() << std::endl;
	}
	else
	{
		ss << GlobalV::global_out_dir << "LOWF_K.dat";
	}

    std::ifstream ifs;

    int error = 0;

    if (GlobalV::DRANK==0)
    {
        ifs.open(ss.str().c_str());
        if (!ifs)
        {
            GlobalV::ofs_warning << " Can't open file:" << ss.str() << std::endl;
            error = 1;
        }
    }

#ifdef __MPI
    Parallel_Common::bcast_int(error);
#endif

    if (error==1) return 1;

    // otherwise, find the file.

    if (GlobalV::MY_RANK==0)
    {
        int nbands, nlocal;
        ModuleBase::GlobalFunc::READ_VALUE(ifs, nbands);
        ModuleBase::GlobalFunc::READ_VALUE(ifs, nlocal);

        if (nbands!=GlobalV::NBANDS)
        {
            GlobalV::ofs_warning << " read in nbands=" << nbands;
            GlobalV::ofs_warning << " NBANDS=" << GlobalV::NBANDS << std::endl;
            error = 2;
        }
        else if (nlocal != GlobalV::NLOCAL)
        {
            GlobalV::ofs_warning << " read in nlocal=" << nlocal;
            GlobalV::ofs_warning << " NLOCAL=" << GlobalV::NLOCAL << std::endl;
            error = 3;
        }

        ctot = new double*[GlobalV::NBANDS];
        for (int i=0; i<GlobalV::NBANDS; i++)
        {
            ctot[i] = new double[GlobalV::NLOCAL];
        }

		//std::cout << "nbands" << GlobalV::NBANDS << std::endl;
        for (int i=0; i<GlobalV::NBANDS; i++)
        {
            int ib;
            ModuleBase::GlobalFunc::READ_VALUE(ifs, ib);
			ModuleBase::GlobalFunc::READ_VALUE(ifs, GlobalC::wf.ekb[GlobalV::CURRENT_SPIN][i]);
			ModuleBase::GlobalFunc::READ_VALUE(ifs, GlobalC::wf.wg(GlobalV::CURRENT_SPIN,i));
            assert( (i+1)==ib);
			//std::cout << " ib=" << ib << std::endl;
            for (int j=0; j<GlobalV::NLOCAL; j++)
            {
                ifs >> ctot[i][j];
				//std::cout << ctot[i][j] << " ";
            }
			//std::cout << std::endl;
        }
    }


#ifdef __MPI
    Parallel_Common::bcast_int(error);
	Parallel_Common::bcast_double( GlobalC::wf.ekb[is], GlobalV::NBANDS);
	Parallel_Common::bcast_double( GlobalC::wf.wg.c, GlobalV::NSPIN*GlobalV::NBANDS);
#endif
	if(error==2) return 2;
	if(error==3) return 3;

	// mohan add 2012-02-15,
	// mohan comment out 2021-02-09
	//GlobalC::SGO.cal_totwfc();

	// distri_lowf need all processors.
	// otherwise, read in sucessfully.
    // if GlobalV::DRANK!=0, ctot is not used,
    // so it's save.

	WF_Local::distri_lowf_new(ctot, is, lowf);
	
	// mohan add 2012-02-15,
	// still have bugs, but can solve it later.
	// distribute the wave functions again.
	// mohan comment out 2021-02-09
	// GlobalC::SGO.dis_subwfc();

    if (GlobalV::MY_RANK==0)
    {
        // delte the ctot
        for (int i=0; i<GlobalV::NBANDS; i++)
        {
            delete[] ctot[i];
        }
        delete[] ctot;
    }

    ModuleBase::timer::tick("WF_Local","read_lowf");
    return 0;
}

void WF_Local::write_lowf(const std::string &name, double **ctot)
{
    ModuleBase::TITLE("WF_Local","write_lowf");
    ModuleBase::timer::tick("WF_Local","write_lowf");

    std::ofstream ofs;
    if (GlobalV::DRANK==0)
    {
        ofs.open(name.c_str());
        if (!ofs)
        {
            ModuleBase::WARNING("Pdiag_Basic::write_lowf","Can't write local orbital wave functions.");
        }
        ofs << GlobalV::NBANDS << " (number of bands)" << std::endl;
        ofs << GlobalV::NLOCAL << " (number of orbitals)";
        ofs << std::setprecision(8);
        ofs << scientific;

        for (int i=0; i<GlobalV::NBANDS; i++)
        {
            // +1 to mean more clearly.
            // band index start from 1.
            ofs << "\n" << i+1 << " (band)";
			ofs << "\n" << GlobalC::wf.ekb[GlobalV::CURRENT_SPIN][i] << " (Ry)"; //mohan add 2012-03-26
			ofs << "\n" << GlobalC::wf.wg(GlobalV::CURRENT_SPIN,i) << " (Occupations)";
            for (int j=0; j<GlobalV::NLOCAL; j++)
            {
                if (j % 5 == 0) ofs << "\n";
                ofs << ctot[i][j] << " ";
            }
        }
        ofs.close();
    }

    ModuleBase::timer::tick("WF_Local","write_lowf");
    return;
}

void WF_Local::write_lowf_complex(const std::string &name, std::complex<double> **ctot, const int &ik)
{
    ModuleBase::TITLE("WF_Local","write_lowf_complex");
    ModuleBase::timer::tick("WF_Local","write_lowf_complex");

    std::ofstream ofs;
    if (GlobalV::DRANK==0)
    {
        ofs.open(name.c_str());
        if (!ofs)
        {
            ModuleBase::WARNING("Pdiag_Basic::write_lowf","Can't write local orbital wave functions.");
        }
        ofs << std::setprecision(25);
		ofs << ik+1 << " (index of k points)" << std::endl;
		ofs << GlobalC::kv.kvec_c[ik].x << " " << GlobalC::kv.kvec_c[ik].y << " " << GlobalC::kv.kvec_c[ik].z << std::endl;
        ofs << GlobalV::NBANDS << " (number of bands)" << std::endl;
        ofs << GlobalV::NLOCAL << " (number of orbitals)";
        ofs << scientific;

        for (int i=0; i<GlobalV::NBANDS; i++)
        {
            // +1 to mean more clearly.
            // band index start from 1.
            ofs << "\n" << i+1 << " (band)";
			ofs << "\n" << GlobalC::wf.ekb[ik][i] << " (Ry)";
			ofs << "\n" << GlobalC::wf.wg(ik,i) << " (Occupations)";
            for (int j=0; j<GlobalV::NLOCAL; j++)
            {
                if (j % 5 == 0) ofs << "\n";
                ofs << ctot[i][j].real() << " " << ctot[i][j].imag() << " ";
            }
        }
        ofs.close();
    }

    ModuleBase::timer::tick("WF_Local","write_lowf_complex");
    return;
}

void WF_Local::distri_lowf_new(double** ctot, const int& is,
    Local_Orbital_wfc &lowf)
{
    ModuleBase::TITLE("WF_Local","distri_lowf_new");
#ifdef __MPI

//1. alloc work array; set some parameters

	long maxnloc; // maximum number of elements in local matrix
	MPI_Reduce(&lowf.ParaV->nloc, &maxnloc, 1, MPI_LONG, MPI_MAX, 0, lowf.ParaV->comm_2D);
	MPI_Bcast(&maxnloc, 1, MPI_LONG, 0, lowf.ParaV->comm_2D);
	//reduce and bcast could be replaced by allreduce
	
    int nprocs, myid;
    MPI_Comm_size(lowf.ParaV->comm_2D, &nprocs);
    MPI_Comm_rank(lowf.ParaV->comm_2D, &myid);

	double *work=new double[maxnloc]; // work/buffer matrix
	int nb = 0;
	if(GlobalV::NB2D==0)
	{
		if(GlobalV::NLOCAL>0) nb = 1;
		if(GlobalV::NLOCAL>500) nb = 32;
		if(GlobalV::NLOCAL>1000) nb = 64;
	}
	else if(GlobalV::NB2D>0)
	{
		nb = GlobalV::NB2D; // mohan add 2010-06-28
	}
	int info;
	int naroc[2]; // maximum number of row or column
	
//2. copy from ctot to wfc_gamma
	for(int iprow=0; iprow<lowf.ParaV->dim0; ++iprow)
	{
		for(int ipcol=0; ipcol<lowf.ParaV->dim1; ++ipcol)
		{
//2.1 get and bcast local 2d matrix info
			const int coord[2]={iprow, ipcol};
			int src_rank;
			MPI_Cart_rank(lowf.ParaV->comm_2D, coord, &src_rank);
			if(myid==src_rank)
			{
				naroc[0]=lowf.ParaV->nrow;
				naroc[1]=lowf.ParaV->ncol;
			}
			info=MPI_Bcast(naroc, 2, MPI_INT, src_rank, lowf.ParaV->comm_2D);

//2.2 copy from ctot to work, then bcast work
			info=CTOT2q(myid, naroc, nb, lowf.ParaV->dim0, lowf.ParaV->dim1, iprow, ipcol, work, ctot);
			info=MPI_Bcast(work, maxnloc, MPI_DOUBLE, 0, lowf.ParaV->comm_2D);
			//GlobalV::ofs_running << "iprow, ipcow : " << iprow << ipcol << std::endl;
			//for (int i=0; i<maxnloc; ++i)
			//{
				//GlobalV::ofs_running << *(work+i)<<" ";
			//}
			//GlobalV::ofs_running << std::endl;
//2.3 copy from work to wfc_gamma
			const int inc=1;
			if(myid==src_rank)
			{
				BlasConnector::copy(lowf.ParaV->nloc, work, inc, lowf.wfc_gamma.at(is).c, inc);
			}
		}//loop ipcol
	}//loop	iprow

	delete[] work;
#else
	ModuleBase::WARNING_QUIT("WF_Local::distri_lowf_new","check the code without MPI.");
#endif
    return;
}

void WF_Local::distri_lowf_complex_new(std::complex<double>** ctot, const int& ik,
    Local_Orbital_wfc &lowf)
{
    ModuleBase::TITLE("WF_Local","distri_lowf_complex_new");
#ifdef __MPI

//1. alloc work array; set some parameters

	long maxnloc; // maximum number of elements in local matrix
	MPI_Reduce(&lowf.ParaV->nloc_wfc, &maxnloc, 1, MPI_LONG, MPI_MAX, 0, lowf.ParaV->comm_2D);
	MPI_Bcast(&maxnloc, 1, MPI_LONG, 0, lowf.ParaV->comm_2D);
	//reduce and bcast could be replaced by allreduce
	
    int nprocs, myid;
    MPI_Comm_size(lowf.ParaV->comm_2D, &nprocs);
    MPI_Comm_rank(lowf.ParaV->comm_2D, &myid);

	std::complex<double> *work=new std::complex<double>[maxnloc]; // work/buffer matrix
	int nb = 0;
	if(GlobalV::NB2D==0)
	{
		if(GlobalV::NLOCAL>0) nb = 1;
		if(GlobalV::NLOCAL>500) nb = 32;
		if(GlobalV::NLOCAL>1000) nb = 64;
	}
	else if(GlobalV::NB2D>0)
	{
		nb = GlobalV::NB2D; // mohan add 2010-06-28
	}
	int info;
	int naroc[2]; // maximum number of row or column
	
//2. copy from ctot to wfc_gamma
	for(int iprow=0; iprow<lowf.ParaV->dim0; ++iprow)
	{
		for(int ipcol=0; ipcol<lowf.ParaV->dim1; ++ipcol)
		{
//2.1 get and bcast local 2d matrix info
			const int coord[2]={iprow, ipcol};
			int src_rank;
			MPI_Cart_rank(lowf.ParaV->comm_2D, coord, &src_rank);
			if(myid==src_rank)
			{
				naroc[0]=lowf.ParaV->nrow;
				naroc[1]=lowf.ParaV->ncol_bands;
			}
			info=MPI_Bcast(naroc, 2, MPI_INT, src_rank, lowf.ParaV->comm_2D);

//2.2 copy from ctot to work, then bcast work
			info=CTOT2q_c(myid, naroc, nb, lowf.ParaV->dim0, lowf.ParaV->dim1, iprow, ipcol, work, ctot);
			info=MPI_Bcast(work, maxnloc, MPI_DOUBLE_COMPLEX, 0, lowf.ParaV->comm_2D);
			//ofs_running << "iprow, ipcow : " << iprow << ipcol << std::endl;
			//for (int i=0; i<maxnloc; ++i)
			//{
				//ofs_running << *(work+i)<<" ";
			//}
			//ofs_running << std::endl;
//2.3 copy from work to wfc_k
            const int inc = 1;
			if(myid==src_rank)
			{
				BlasConnector::copy(lowf.ParaV->nloc_wfc, work, inc, lowf.wfc_k.at(ik).c, inc);
			}
		}//loop ipcol
	}//loop	iprow

	delete[] work;
#else
	ModuleBase::WARNING_QUIT("WF_Local::distri_lowf_new","check the code without MPI.");
#endif
    return;
}

void WF_Local::distri_lowf(double **ctot, double **c)
{
    ModuleBase::TITLE("WF_Local","distri_lowf");
#ifdef __MPI

    MPI_Status status;
    for (int i=0; i<GlobalV::DSIZE; i++)
    {
        if (GlobalV::DRANK==0)
        {
            if (i==0)
            {
                // get the wave functions from 'ctot',
                // save them in the matrix 'c'.
                for (int iw=0; iw<GlobalV::NLOCAL; iw++)
                {
					// mohan update 2012-01-12
//                  const int mu_local = GlobalC::GridT.trace_lo[iw]; 
                    const int mu_local = GlobalC::SGO.trace_lo_tot[iw];

                    if (mu_local >= 0)
                    {
                        for (int ib=0; ib<GlobalV::NBANDS; ib++)
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
                int* trace_lo2 = new int[GlobalV::NLOCAL];
                MPI_Recv(trace_lo2, GlobalV::NLOCAL, MPI_INT, i, tag, DIAG_WORLD, &status);

                // receive lgd2
                int lgd2 = 0;
                tag = i * 3 + 1;
                MPI_Recv(&lgd2, 1, MPI_INT, i, tag, DIAG_WORLD, &status);

                // send csend
                double* csend = new double[GlobalV::NBANDS*lgd2];
                ModuleBase::GlobalFunc::ZEROS(csend, GlobalV::NBANDS*lgd2);

                for (int ib=0; ib<GlobalV::NBANDS; ib++)
                {
                    for (int iw=0; iw<GlobalV::NLOCAL; iw++)
                    {
                        const int mu_local = trace_lo2[iw];
                        if (mu_local>=0)
                        {
                            csend[mu_local*GlobalV::NBANDS+ib] = ctot[ib][iw];
                        }
                    }
                }

                tag = i * 3 + 2;
                MPI_Send(csend,GlobalV::NBANDS*lgd2,MPI_DOUBLE,i,tag,DIAG_WORLD);

                delete[] trace_lo2;
                delete[] csend;
            }
        }// end GlobalV::DRANK=0
        else if ( i == GlobalV::DRANK)
        {
            int tag;

            // send trace_lo
            tag = GlobalV::DRANK * 3;
			// mohan update 2012-01-12
            //MPI_Send(GlobalC::GridT.trace_lo, GlobalV::NLOCAL, MPI_INT, 0, tag, DIAG_WORLD);
            MPI_Send(GlobalC::SGO.trace_lo_tot, GlobalV::NLOCAL, MPI_INT, 0, tag, DIAG_WORLD);

            // send lgd
            tag = GlobalV::DRANK * 3 + 1;

			// mohan update 2012-01-12
			int lgdnow = GlobalC::SGO.lgd;
            MPI_Send(&lgdnow, 1, MPI_INT, 0, tag, DIAG_WORLD);

            // receive c
			GlobalV::ofs_running << " lgdnow=" << lgdnow << std::endl;
            double* crecv = new double[GlobalV::NBANDS*lgdnow];
            ModuleBase::GlobalFunc::ZEROS(crecv, GlobalV::NBANDS*lgdnow);
            tag = GlobalV::DRANK * 3 + 2;
            MPI_Recv(crecv, GlobalV::NBANDS*lgdnow, MPI_DOUBLE, 0, tag, DIAG_WORLD, &status);

            for (int ib=0; ib<GlobalV::NBANDS; ib++)
            {
                for (int mu=0; mu<lgdnow; mu++)
                {
                    c[ib][mu] = crecv[mu*GlobalV::NBANDS+ib];
                }
            }

            delete[] crecv;
        }// end i==GlobalV::DRANK
        MPI_Barrier(DIAG_WORLD);
    }// end i

    //-----------
    // for test,
    //-----------
    /*
    GlobalV::ofs_running << " Wave Functions in local basis: " << std::endl;
    for(int i=0; i<GlobalV::NBANDS; i++)
    {
        for(int j=0; j<GlobalC::GridT.lgd; j++)
        {
            if(j%8==0) GlobalV::ofs_running << std::endl;
            if( abs(c[i][j]) > 1.0e-5  )
            {
                GlobalV::ofs_running << std::setw(15) << c[i][j];
            }
            else
            {
                GlobalV::ofs_running << std::setw(15) << "0";
            }
        }
    }
    GlobalV::ofs_running << std::endl;
    */
#else
	ModuleBase::WARNING_QUIT("WF_Local::distri_lowf","check the code without MPI.");
#endif
    return;
}


void WF_Local::distri_lowf_complex(std::complex<double> **ctot, std::complex<double> **cc)
{
    ModuleBase::TITLE("WF_Local","distri_lowf_complex");
#ifdef __MPI

    MPI_Status status;

    for (int i=0; i<GlobalV::DSIZE; i++)
    {
        if (GlobalV::DRANK==0)
        {
            if (i==0)
            {
                // get the wave functions from 'ctot',
                // save them in the matrix 'c'.
                for (int iw=0; iw<GlobalV::NLOCAL; iw++)
                {
                    const int mu_local = GlobalC::GridT.trace_lo[iw];
                    if (mu_local >= 0)
                    {
                        for (int ib=0; ib<GlobalV::NBANDS; ib++)
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
					int* trace_lo2 = new int[GlobalV::NLOCAL];
					MPI_Recv(trace_lo2, GlobalV::NLOCAL, MPI_INT, i, tag, DIAG_WORLD, &status);

//					GlobalV::ofs_running << " lgd2=" << lgd2 << " proc=" << i+1 << std::endl;
					// send csend
					std::complex<double>* csend = new std::complex<double>[GlobalV::NBANDS*lgd2];
					ModuleBase::GlobalFunc::ZEROS(csend, GlobalV::NBANDS*lgd2);
					for (int ib=0; ib<GlobalV::NBANDS; ib++)
					{
						for (int iw=0; iw<GlobalV::NLOCAL; iw++)
						{
							const int mu_local = trace_lo2[iw];
							if (mu_local>=0)
							{
								csend[mu_local*GlobalV::NBANDS+ib] = ctot[ib][iw];
							}
						}
					}
					tag = i * 3 + 2;
					MPI_Send(csend,GlobalV::NBANDS*lgd2,mpicomplex,i,tag,DIAG_WORLD);
                	delete[] csend;
                	delete[] trace_lo2;
				}
            }
        }// end GlobalV::DRANK=0
        else if ( i == GlobalV::DRANK)
		{
			int tag;

			// send GlobalC::GridT.lgd
			tag = GlobalV::DRANK * 3;
			MPI_Send(&GlobalC::GridT.lgd, 1, MPI_INT, 0, tag, DIAG_WORLD);

			if(GlobalC::GridT.lgd != 0)
			{
				// send trace_lo
				tag = GlobalV::DRANK * 3 + 1;
				MPI_Send(GlobalC::GridT.trace_lo, GlobalV::NLOCAL, MPI_INT, 0, tag, DIAG_WORLD);

				// receive cc
				std::complex<double>* crecv = new std::complex<double>[GlobalV::NBANDS*GlobalC::GridT.lgd];
				ModuleBase::GlobalFunc::ZEROS(crecv, GlobalV::NBANDS*GlobalC::GridT.lgd);

				tag = GlobalV::DRANK * 3 + 2;
				MPI_Recv(crecv, GlobalV::NBANDS*GlobalC::GridT.lgd, mpicomplex, 0, tag, DIAG_WORLD, &status);

				for (int ib=0; ib<GlobalV::NBANDS; ib++)
				{
					for (int mu=0; mu<GlobalC::GridT.lgd; mu++)
					{
						cc[ib][mu] = crecv[mu*GlobalV::NBANDS+ib];
					}
				}

				delete[] crecv;

			}
        }// end i==GlobalV::DRANK
        MPI_Barrier(DIAG_WORLD);
    }// end i

    //-----------
    // for test,
    //-----------
    /*
    GlobalV::ofs_running << " Wave Functions in local basis: " << std::endl;
    for(int i=0; i<GlobalV::NBANDS; i++)
    {
        for(int j=0; j<GlobalC::GridT.lgd; j++)
        {
            if(j%8==0) GlobalV::ofs_running << std::endl;
            if( abs(c[i][j]) > 1.0e-5  )
            {
                GlobalV::ofs_running << std::setw(15) << c[i][j];
            }
            else
            {
                GlobalV::ofs_running << std::setw(15) << "0";
            }
        }
    }
    GlobalV::ofs_running << std::endl;
    */
#else
	ModuleBase::WARNING_QUIT("WF_Local::distri_lowf_complex","check the code without MPI.");
#endif
    return;
}
