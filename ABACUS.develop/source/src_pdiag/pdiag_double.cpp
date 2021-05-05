/*test generalized Stardand double precision symmetric eigenproblem*/
#include "pdiag_double.h"
#include "src_global/lapack_connector.h"
#include "src_pw/occupy.h"
#include "src_pw/global.h"
#include "src_lcao/local_orbital_charge.h"
#include "src_io/wf_local.h"


#ifdef __MPI
extern "C"
{
    #include "Cblacs.h"
//    #include "pblas.h"
//    #include "scalapack.h"
    #include "my_elpa.h"
	#include "src_global/scalapack_connector.h"
}
#include "GenELPA.h"
#include "pdgseps.h"
#include "pzgseps.h"
#include "src_global/lapack_connector.h"
#endif

#include "src_external/src_test/test_function.h"

inline int cart2blacs(MPI_Comm comm_2D, int nprows, int npcols, int N, int nblk, int lld, int *desc, int &mpi_comm_rows, int &mpi_comm_cols)
{
#ifdef __MPI
    int my_blacs_ctxt;
    int myprow, mypcol;
    int *usermap=new int[nprows*npcols];
    int info=0;
    for(int i=0; i<nprows; ++i)
    {
        for(int j=0; j<npcols; ++j)
        {
            int pcoord[2]={i, j};
            MPI_Cart_rank(comm_2D, pcoord, &usermap[i+j*nprows]);
        }
    }
    Cblacs_get(comm_2D, 0, &my_blacs_ctxt);
    Cblacs_gridmap(&my_blacs_ctxt, usermap, nprows, nprows, npcols);
    Cblacs_gridinfo(my_blacs_ctxt, &nprows, &npcols, &myprow, &mypcol);
    info=elpa_get_communicators(comm_2D, myprow, mypcol, &mpi_comm_rows, &mpi_comm_cols);
    delete[] usermap;
    int ISRC=0;
    descinit_(desc, &N, &N, &nblk, &nblk, &ISRC, &ISRC, &my_blacs_ctxt, &lld, &info);

    return my_blacs_ctxt;
#else
    return 0;
#endif
}

inline int q2ZLOC_WFC(int pos, int naroc[2], int nb,
                      int dim0, int dim1, int iprow, int ipcol,
                      int loc_size,
                      double* work, double* ZLOC, double** WFC)
{
    //OUT(ofs_running,"start q2ZLOC_WFC");
    for(int j=0; j<naroc[1]; ++j)
    {
        int igcol=globalIndex(j, nb, dim1, ipcol);
        if(igcol>=NBANDS) continue;
        int zcol=igcol-pos;
        for(int i=0; i<naroc[0]; ++i)
        {
            int igrow=globalIndex(i, nb, dim0, iprow);
            if(0<=zcol && zcol<loc_size)
            {
                ZLOC[igrow*loc_size+zcol]=work[j*naroc[0]+i];
            }
	        int mu_local=SGO.trace_lo_tot[igrow];
            if(mu_local>=0 && igrow<NBANDS)
            {
                WFC[igcol][mu_local]=work[j*naroc[0]+i];
            }
        }
    }
    //OUT(ofs_running,"WFC was done in q2ZLOC_WFC");
    return 0;
}

inline int q2ZLOC_WFC_WFCAUG(int pos, int naroc[2], int nb,
                             int dim0, int dim1, int iprow, int ipcol,
                             int loc_size,
                             double* work, double* ZLOC, double** WFC, double** WFCAUG)
{
	
    stringstream ss;
    for(int j=0; j<naroc[1]; ++j)
    {
        int igcol=globalIndex(j, nb, dim1, ipcol);
/*        ss<<"local column "<<j<<" nb "<<nb<<" dim1 "<<dim1<<" mypcol "<<ipcol<<" global column (NBANDS) "<<igcol;
		OUT(ofs_running,ss.str());
		ss.str("");*/
        if(igcol>=NBANDS) continue;
        int zcol=igcol-pos;
        for(int i=0; i<naroc[0]; ++i)
        {
            int igrow=globalIndex(i, nb, dim0, iprow);
/*	        ss<<"    local row "<<i<<" nb "<<nb<<" dim0 "<<dim0<<" myprow "<<iprow<<" global row (NLOCAL)"<<igrow;
			OUT(ofs_running,ss.str());
			ss.str("");*/
            if(0<=zcol && zcol<loc_size)
            {
                ZLOC[igrow*loc_size+zcol]=work[j*naroc[0]+i];
            }
	        int mu_local=SGO.trace_lo_tot[igrow];
            if(mu_local>=0)
            {
                WFC[igcol][mu_local]=work[j*naroc[0]+i];
            }
	        int mu_aug=LOWF.trace_aug[igrow];
            if(mu_aug>=0)
            {
                WFCAUG[igcol][mu_aug]=work[j*naroc[0]+i];
            }
        }
    }

    //OUT(ofs_running,"WFCAUG was done in q2ZLOC_WFC_WFCAUG");  
    return 0;
}

inline int q2ZLOC_WFC_CTOT(int myid, int pos, int naroc[2], int nb,
                           int dim0, int dim1, int iprow, int ipcol,
                           int loc_size,
                           double* work, double* ZLOC, double** WFC, double** CTOT)
{
    for(int j=0; j<naroc[1]; ++j)
    {
        int igcol=globalIndex(j, nb, dim1, ipcol);
        if(igcol>=NBANDS) continue;
        int zcol=igcol-pos;
        for(int i=0; i<naroc[0]; ++i)
        {
            int igrow=globalIndex(i, nb, dim0, iprow);
            if(0<=zcol && zcol<loc_size)
            {
                ZLOC[igrow*loc_size+zcol]=work[j*naroc[0]+i];
            }
	        int mu_local=SGO.trace_lo_tot[igrow];
            if(mu_local>=0)
            {
                WFC[igcol][mu_local]=work[j*naroc[0]+i];
            }
            if(myid==0) CTOT[igcol][igrow]=work[j*naroc[0]+i];
        }
    }
    return 0;
}

inline int q2ZLOC_WFC_WFCAUG_CTOT(int myid, int pos, int naroc[2], int nb,
                                  int dim0, int dim1, int iprow, int ipcol,
                                  int loc_size,
                                  double* work, double* ZLOC, double** WFC, double** WFCAUG, double** CTOT)
{
    for(int j=0; j<naroc[1]; ++j)
    {
        int igcol=globalIndex(j, nb, dim1, ipcol);
        if(igcol>=NBANDS) continue;
        int zcol=igcol-pos;
        for(int i=0; i<naroc[0]; ++i)
        {
            int igrow=globalIndex(i, nb, dim0, iprow);
            if(0<=zcol && zcol<loc_size)
            {
                ZLOC[igrow*loc_size+zcol]=work[j*naroc[0]+i];
            }
	        int mu_local=SGO.trace_lo_tot[igrow];
            if(mu_local>=0)
            {
                WFC[igcol][mu_local]=work[j*naroc[0]+i];
            }
	        int mu_aug=LOWF.trace_aug[igrow];
            if(mu_aug>=0)
            {
                WFCAUG[igcol][mu_aug]=work[j*naroc[0]+i];
            }
            if(myid==0) CTOT[igcol][igrow]=work[j*naroc[0]+i];
        }
    }
    return 0;
}

inline int q2WFC_complex(int naroc[2], int nb,
                      	 int dim0, int dim1, int iprow, int ipcol,
                         complex<double>* work, complex<double>** WFC)
{
    for(int j=0; j<naroc[1]; ++j)
    {
        int igcol=globalIndex(j, nb, dim1, ipcol);
        if(igcol>=NBANDS) continue;
        for(int i=0; i<naroc[0]; ++i)
        {
            int igrow=globalIndex(i, nb, dim0, iprow);
	        int mu_local=GridT.trace_lo[igrow];
            if(mu_local>=0)
            {
                WFC[igcol][mu_local]=work[j*naroc[0]+i];
            }
        }
    }
    return 0;
}

inline int q2WFC_WFCAUG_complex(int naroc[2], int nb,
                                int dim0, int dim1, int iprow, int ipcol,
                                complex<double>* work, complex<double>** WFC, complex<double>** WFCAUG)
{
    for(int j=0; j<naroc[1]; ++j)
    {
        int igcol=globalIndex(j, nb, dim1, ipcol);
        if(igcol>=NBANDS) continue;
        for(int i=0; i<naroc[0]; ++i)
        {
            int igrow=globalIndex(i, nb, dim0, iprow);
	        int mu_local=GridT.trace_lo[igrow];
            if(mu_local>=0)
            {
                WFC[igcol][mu_local]=work[j*naroc[0]+i];
            }
	        int mu_aug=LOWF.trace_aug[igrow];
            if(mu_aug>=0)
            {
                WFCAUG[igcol][mu_aug]=work[j*naroc[0]+i];
            }
        }
    }
    return 0;
}

inline int q2WFC_CTOT_complex(int myid, int naroc[2], int nb,
                              int dim0, int dim1, int iprow, int ipcol,
                           	  complex<double>* work, complex<double>** WFC, complex<double>** CTOT)
{
    for(int j=0; j<naroc[1]; ++j)
    {
        int igcol=globalIndex(j, nb, dim1, ipcol);
        if(igcol>=NBANDS) continue;
        for(int i=0; i<naroc[0]; ++i)
        {
            int igrow=globalIndex(i, nb, dim0, iprow);
	        int mu_local=GridT.trace_lo[igrow];
            if(mu_local>=0)
            {
                WFC[igcol][mu_local]=work[j*naroc[0]+i];
            }
            if(myid==0) CTOT[igcol][igrow]=work[j*naroc[0]+i];
        }
    }
    return 0;
}

inline int q2WFC_WFCAUG_CTOT_complex(int myid, int naroc[2], int nb,
                                     int dim0, int dim1, int iprow, int ipcol,
                                     complex<double>* work, complex<double>** WFC, complex<double>** WFCAUG, complex<double>** CTOT)
{
    for(int j=0; j<naroc[1]; ++j)
    {
        int igcol=globalIndex(j, nb, dim1, ipcol);
        if(igcol>=NBANDS) continue;
        for(int i=0; i<naroc[0]; ++i)
        {
            int igrow=globalIndex(i, nb, dim0, iprow);
	        int mu_local=GridT.trace_lo[igrow];
            if(mu_local>=0)
            {
                WFC[igcol][mu_local]=work[j*naroc[0]+i];
            }
	        int mu_aug=LOWF.trace_aug[igrow];
            if(mu_aug>=0)
            {
                WFCAUG[igcol][mu_aug]=work[j*naroc[0]+i];
            }
            if(myid==0) CTOT[igcol][igrow]=work[j*naroc[0]+i];
        }
    }
    return 0;
}

Pdiag_Double::Pdiag_Double()
{
	// default value of nb is 1,
	// but can change to larger value from input.
    nb = 1;
	MatrixInfo.row_set = new int[1];
	MatrixInfo.col_set = new int[1];
}

Pdiag_Double::~Pdiag_Double()
{
	delete[] MatrixInfo.row_set;
	delete[] MatrixInfo.col_set;
}

void Pdiag_Double::divide_HS_2d
(
#ifdef __MPI
	MPI_Comm DIAG_WORLD
#endif
)
{
	TITLE("Pdiag_Double","divide_HS_2d");
	assert(NLOCAL>0);
	assert(DSIZE>0);

#ifdef __MPI
	DIAG_HPSEPS_WORLD=DIAG_WORLD;
#endif

	if(DCOLOR!=0) return; // mohan add 2012-01-13

	// get the 2D index of computer.
	this->dim0 = (int)sqrt((double)DSIZE); //mohan update 2012/01/13
	//while (NPROC_IN_POOL%dim0!=0)
	while (DSIZE%dim0!=0)
	{
		this->dim0 = dim0 - 1;
	}
	assert(dim0 > 0);
	this->dim1=DSIZE/dim0;

	if(testpb)OUT(ofs_running,"dim0",dim0);
	if(testpb)OUT(ofs_running,"dim1",dim1);

#ifdef __MPI
	// mohan add 2011-04-16
	if(NB2D==0)
	{
		if(NLOCAL>0) this->nb = 1;
		if(NLOCAL>500) this->nb = 32;
		if(NLOCAL>1000) this->nb = 64;
	}
	else if(NB2D>0)
	{
		this->nb = NB2D; // mohan add 2010-06-28
	}
	OUT(ofs_running,"nb2d",nb);

	this->set_parameters();

	// call mpi_creat_cart
	this->mpi_creat_cart(&this->comm_2D,this->dim0,this->dim1);

	// call mat_2d
	this->mat_2d(this->comm_2D, NLOCAL, NLOCAL, this->nb, this->MatrixInfo);

	// mohan add 2010-06-29
	this->nrow = this->MatrixInfo.row_num;
	this->ncol = this->MatrixInfo.col_num;
	this->nloc = MatrixInfo.col_num * MatrixInfo.row_num;

	// init blacs context for genelpa
    if(KS_SOLVER=="genelpa" || KS_SOLVER=="scalapack_gvx")
    {
        blacs_ctxt=cart2blacs(comm_2D, dim0, dim1, NLOCAL, nb, nrow, desc, mpi_comm_rows, mpi_comm_cols);
    }
#else // single processor used.
	this->nb = NLOCAL;
	this->nrow = NLOCAL;
	this->ncol = NLOCAL;
	this->nloc = NLOCAL * NLOCAL;
	this->set_parameters();
	MatrixInfo.row_b = 1;
	MatrixInfo.row_num = NLOCAL;
	delete[] MatrixInfo.row_set;
	MatrixInfo.row_set = new int[NLOCAL];
	for(int i=0; i<NLOCAL; i++)
	{
		MatrixInfo.row_set[i]=i;
	}
	MatrixInfo.row_pos=0;

	MatrixInfo.col_b = 1;
	MatrixInfo.col_num = NLOCAL;
	delete[] MatrixInfo.col_set;
	MatrixInfo.col_set = new int[NLOCAL];
	for(int i=0; i<NLOCAL; i++)
	{
		MatrixInfo.col_set[i]=i;
	}
	MatrixInfo.col_pos=0;
#endif

	assert(nloc>0);
	if(testpb)OUT(ofs_running,"MatrixInfo.row_num",MatrixInfo.row_num);
	if(testpb)OUT(ofs_running,"MatrixInfo.col_num",MatrixInfo.col_num);
	if(testpb)OUT(ofs_running,"nloc",nloc);
	return;
}

void Pdiag_Double::diago_double_begin(
	const int &ik, // k-point index
	double **wfc, // wave functions
	matrix &wfc_2d, // wave functions in 2d
	double* h_mat, // hamiltonian matrix
	double* s_mat, // overlap matrix
	double* ekb) // eigenvalues for each k-point and band
{
	#ifdef TEST_DIAG
	{
		static int istep = 0;
		auto print_matrix_C = [&](const string &file_name, double*m)
		{
			ofstream ofs(file_name+"-C_"+TO_STRING(istep)+"_"+TO_STRING(MY_RANK));
			for(int ic=0; ic<ParaO.ncol; ++ic)
			{
				for(int ir=0; ir<ParaO.nrow; ++ir)
				{
					const int index=ic*ParaO.nrow+ir;
					if(abs(m[index])>1E-10)
						ofs<<m[index]<<"\t";
					else
						ofs<<0<<"\t";
				}
				ofs<<endl;
			}
		};
		auto print_matrix_F = [&](const string &file_name, double*m)
		{
			ofstream ofs(file_name+"-F_"+TO_STRING(istep)+"_"+TO_STRING(MY_RANK));
			for(int ir=0; ir<ParaO.nrow; ++ir)
			{
				for(int ic=0; ic<ParaO.ncol; ++ic)
				{
					const int index=ic*ParaO.nrow+ir;
					if(abs(m[index])>1E-10)
						ofs<<m[index]<<"\t";
					else
						ofs<<0<<"\t";
				}
				ofs<<endl;
			}
		};
		print_matrix_F("H_gamma", h_mat);
		print_matrix_F("S_gamma", s_mat);
		print_matrix_C("H_gamma", h_mat);
		print_matrix_C("S_gamma", s_mat);
		++istep;
	}
	#endif
	
#ifdef __MPI
	TITLE("Pdiag_Double","diago_begin");
	assert(this->loc_size > 0);
	assert(NLOCAL > 0);

	char uplo='U';
	const int inc=1;

    int nprocs, myid;
    MPI_Status status;
    MPI_Comm_size(comm_2D, &nprocs);
    MPI_Comm_rank(comm_2D, &myid);	

	// parallel diagonalize the
	// H | psi > = S | psi >
	// problem.
	int loc_pos;

	double* Stmp = LM.Sdiag;

    OUT(ofs_running,"start solver, KS_SOLVER",KS_SOLVER);
    if(KS_SOLVER=="hpseps")
    {
        double *eigen = new double[NLOCAL];
        ZEROS(eigen, NLOCAL);

        double* Z = new double[this->loc_size * NLOCAL];
        ZEROS(Z, this->loc_size * NLOCAL);

        Memory::record("Pdiag_Double","Z",loc_size * NLOCAL,"double");
        timer::tick("Diago_LCAO_Matrix","pdgseps",'G');
		LapackConnector::copy(nloc, s_mat, inc, Stmp, inc);
		pdgseps(comm_2D, NLOCAL, nb, h_mat, Stmp, Z, eigen, this->MatrixInfo, uplo, this->loc_size, loc_pos);
        timer::tick("Diago_LCAO_Matrix","pdgseps",'G');

        if(myid <= lastband_in_proc)
        {
            for(int i=0; i<loc_sizes[myid]; i++)
            {
                for(int n=0; n<NLOCAL; n++)
                {
                    Z_LOC[ik][n*loc_sizes[myid] + i] = Z[n*loc_sizes[myid] + i];
                }
            }
        }

        // the eigenvalues.
        //xiaohui modify 2014-06-15, move to the top
        LapackConnector::copy(NBANDS, eigen, inc, ekb, inc);
        delete[] eigen;
        //=====================================
        // gather the eigenvectors and
        // distribute them to each processor
        // Z is delete in gath_eig
        //=====================================

        //xiaohui modify 2014-06-18
        timer::tick("Diago_LCAO_Matrix","gath_eig",'G');
        this->gath_eig(DIAG_HPSEPS_WORLD, NLOCAL, wfc, Z);
        timer::tick("Diago_LCAO_Matrix","gath_eig",'G');
        //delete[] Z; //LiuXh 20171109
	}// HPSEPS method
    else if(KS_SOLVER=="genelpa")
    {
        double *eigen = new double[NLOCAL];
        ZEROS(eigen, NLOCAL);

        long maxnloc; // maximum number of elements in local matrix
        MPI_Reduce(&nloc, &maxnloc, 1, MPI_LONG, MPI_MAX, 0, comm_2D);
        MPI_Bcast(&maxnloc, 1, MPI_LONG, 0, comm_2D);
		wfc_2d.create(this->ncol,this->nrow);			// Fortran order

        double *work=new double[maxnloc]; // work/buffer matrix
        static int method;
        bool wantEigenVector=true;
        bool wantDebug=true;
        int info;
        int comm_2D_f=MPI_Comm_c2f(comm_2D);

        int THIS_REAL_ELPA_KERNEL_API=12;
        int useQR=0;						// may be changed to input parameter sometime

        if(CHR.get_new_e_iteration())
        {
            timer::tick("Diago_LCAO_Matrix","genelpa1",'G');
            method=0;
        	LapackConnector::copy(nloc, s_mat, inc, Stmp, inc);
            info=pdDecomposeRightMatrix2(NLOCAL, nrow, ncol, desc,
                                        Stmp, eigen, wfc_2d.c, work,
                                        comm_2D_f, mpi_comm_rows, mpi_comm_cols,
                                        method, THIS_REAL_ELPA_KERNEL_API, useQR);
            timer::tick("Diago_LCAO_Matrix","genelpa1",'G');
        }
        timer::tick("Diago_LCAO_Matrix","genelpa2",'G');
        info=pdSolveEigen2(NBANDS, NLOCAL, nrow, ncol, desc,
                          h_mat, Stmp, eigen, wfc_2d.c, work,
                          comm_2D_f, mpi_comm_rows, mpi_comm_cols, method,
                          THIS_REAL_ELPA_KERNEL_API, useQR,
                          wantEigenVector, wantDebug);
        timer::tick("Diago_LCAO_Matrix","genelpa2",'G');

    	OUT(ofs_running,"K-S equation was solved by genelpa2");
        LapackConnector::copy(NBANDS, eigen, inc, ekb, inc);
        delete[] eigen;
	    OUT(ofs_running,"eigenvalues were copied to ekb");

        if(NEW_DM==0)
        {
            // convert wave function to band distribution 
			// and calculate the density matrix in the tranditional way
            // redistribute eigenvectors to wfc / wfc_aug

            timer::tick("Diago_LCAO_Matrix","gath_eig",'G');
            int pos=0;
            for(int i=0; i<myid; ++i)
            {
                pos+=loc_sizes[i];
            }
            int naroc[2]; // maximum number of row or column
            for(int iprow=0; iprow<dim0; ++iprow)
            {
                for(int ipcol=0; ipcol<dim1; ++ipcol)
                {
                    const int coord[2]={iprow, ipcol};
                    int src_rank;
                    MPI_Cart_rank(comm_2D, coord, &src_rank);
                    if(myid==src_rank)
                    {
                        LapackConnector::copy(nloc, wfc_2d.c, inc, work, inc);
                        naroc[0]=nrow;
                        naroc[1]=ncol;
                    }
                    info=MPI_Bcast(naroc, 2, MPI_INT, src_rank, comm_2D);
                    info=MPI_Bcast(work, maxnloc, MPI_DOUBLE, src_rank, comm_2D);

                    if(this->out_lowf)
                    {
                        double **ctot;
                        if(myid==0)
                        {
                            ctot = new double*[NBANDS];
                            for (int i=0; i<NBANDS; i++)
                            {
                                ctot[i] = new double[NLOCAL];
                                ZEROS(ctot[i], NLOCAL);
                            }
                            Memory::record("Pdiag_Basic","ctot",NBANDS*NLOCAL,"double");
                        }
						// mohan delete Bfield option 2021-02-12
						info=q2ZLOC_WFC_WFCAUG_CTOT(myid, pos, naroc, nb,
								dim0, dim1, iprow, ipcol, this->loc_size,
								work, Z_LOC[ik], wfc, LOWF.WFC_GAMMA_aug[CURRENT_SPIN], ctot);
						stringstream ss;
                        ss << global_out_dir << "LOWF_GAMMA_S" << CURRENT_SPIN+1 << ".dat";
                        // mohan add 2012-04-03, because we need the occupations for the
                        // first iteration.
                        WF_Local::write_lowf( ss.str(), ctot );//mohan add 2010-09-09
                        if(myid==0)
                        {
                            for (int i=0; i<NBANDS; i++)
                            {
                                delete[] ctot[i];
                            }
                            delete[] ctot;
                        }
                    }
                    else
                    {
						// mohan update 2021-02-12, delete Bfield option
						info=q2ZLOC_WFC_WFCAUG(pos, naroc, nb,
								dim0, dim1, iprow, ipcol, this->loc_size,
								work, Z_LOC[ik], wfc, LOWF.WFC_GAMMA_aug[CURRENT_SPIN]);
					}
                }
            }
        }

        delete[] work;
        timer::tick("Diago_LCAO_Matrix","gath_eig",'G');
    } // GenELPA method
	else if(KS_SOLVER=="lapack_gv")
	{
		wfc_2d.create(this->ncol, this->nrow, false);
		memcpy( wfc_2d.c, h_mat, sizeof(double)*this->ncol*this->nrow );
		matrix s_tmp(this->ncol, this->nrow, false);
		memcpy( s_tmp.c, s_mat, sizeof(double)*this->ncol*this->nrow );
		vector<double> ekb_tmp(NLOCAL,0);
		
		const char jobz='V', uplo='U';
		const int itype=1;
		int lwork=-1, info=0;
		vector<double> work(1,0);
		dsygv_(&itype, &jobz, &uplo, &NLOCAL, wfc_2d.c, &NLOCAL, 
			s_tmp.c, &NLOCAL, ekb_tmp.data(), work.data(), &lwork, &info);

		if(info)
		{
			throw runtime_error("info="+TO_STRING(info)+". "+TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));
		}

		lwork = work[0];
		work.resize(lwork);

		dsygv_(&itype, &jobz, &uplo, &NLOCAL, wfc_2d.c, &NLOCAL, 
			s_tmp.c, &NLOCAL, ekb_tmp.data(), work.data(), &lwork, &info);

		if(info)
		{
			throw runtime_error("info="+TO_STRING(info)+". "+TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));
		}
		memcpy( ekb, ekb_tmp.data(), sizeof(double)*NBANDS ); 
		
		if(NEW_DM==0)
		{
			throw domain_error("NEW_DM must be 1. "+TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));
		}
	}
	else if(KS_SOLVER=="lapack_gvx")
	{
		matrix h_tmp(this->ncol, this->nrow, false);
		memcpy( h_tmp.c, h_mat, sizeof(double)*this->ncol*this->nrow );
		matrix s_tmp(this->ncol, this->nrow, false);
		memcpy( s_tmp.c, s_mat, sizeof(double)*this->ncol*this->nrow );
		wfc_2d.create(this->ncol, this->nrow, false);
		
		const char jobz='V', range='I', uplo='U';
		const int itype=1, il=1, iu=NBANDS;
		int M=0, lwork=-1, info=0;
		const double abstol=0;
		vector<double> work(1,0);
		vector<int> iwork(5*NLOCAL,0);
		vector<int> ifail(NLOCAL,0);

		dsygvx_(&itype, &jobz, &range, &uplo,
			&NLOCAL, h_tmp.c, &NLOCAL, s_tmp.c, &NLOCAL, NULL, NULL, &il, &iu, &abstol,
			&M, ekb, wfc_2d.c, &NLOCAL, work.data(), &lwork, iwork.data(), ifail.data(), &info);

		if(info)
		{
			throw runtime_error("info="+TO_STRING(info)+". "+TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));
		}

		lwork = work[0];
		work.resize(lwork);
		dsygvx_(&itype, &jobz, &range, &uplo,
			&NLOCAL, h_tmp.c, &NLOCAL, s_tmp.c, &NLOCAL, NULL, NULL, &il, &iu, &abstol,
			&M, ekb, wfc_2d.c, &NLOCAL, work.data(), &lwork, iwork.data(), ifail.data(), &info);

		if(info)
		{
			throw runtime_error("info="+TO_STRING(info)+". "+TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));
		}
		if(M!=NBANDS)
		{
			throw runtime_error("M="+TO_STRING(M)+". NBANDS="+TO_STRING(NBANDS)+". "+TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));
		}
		
		if(NEW_DM==0)
		{
			throw domain_error("NEW_DM must be 1. "+TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));
		}
	}
	else if(KS_SOLVER=="scalapack_gvx")
	{
		matrix h_tmp(this->ncol, this->nrow, false);
		memcpy( h_tmp.c, h_mat, sizeof(double)*this->ncol*this->nrow );
		matrix s_tmp(this->ncol, this->nrow, false);
		memcpy( s_tmp.c, s_mat, sizeof(double)*this->ncol*this->nrow );
		wfc_2d.create(this->ncol, this->nrow, false);
		
		const char jobz='V', range='I', uplo='U';
		const int itype=1, il=1, iu=NBANDS, one=1;
		int M=0, NZ=0, lwork=-1, liwork=-1, info=0;
		const double abstol=0, orfac=-1;
		vector<double> work(1,0);
		vector<int> iwork(1,0);
		vector<int> ifail(NLOCAL,0);
		vector<int> iclustr(2*DSIZE);
		vector<double> gap(DSIZE);
		
		pdsygvx_(&itype, &jobz, &range, &uplo,
			&NLOCAL, h_tmp.c, &one, &one, desc, s_tmp.c, &one, &one, desc,
			NULL, NULL, &il, &iu, &abstol,
			&M, &NZ, ekb, &orfac, wfc_2d.c, &one, &one, desc,
			work.data(), &lwork, iwork.data(), &liwork, ifail.data(), iclustr.data(), gap.data(), &info);

		ofs_running<<"lwork="<<work[0]<<"\t"<<"liwork="<<iwork[0]<<endl;
		lwork = work[0];
		work.resize(lwork,0);
		liwork = iwork[0];
		iwork.resize(liwork,0);

		pdsygvx_(&itype, &jobz, &range, &uplo,
			&NLOCAL, h_tmp.c, &one, &one, desc, s_tmp.c, &one, &one, desc,
			NULL, NULL, &il, &iu, &abstol,
			&M, &NZ, ekb, &orfac, wfc_2d.c, &one, &one, desc,
			work.data(), &lwork, iwork.data(), &liwork, ifail.data(), iclustr.data(), gap.data(), &info);

		ofs_running<<"M="<<M<<"\t"<<"NZ="<<NZ<<endl;

		if(info)
		{
			throw runtime_error("info="+TO_STRING(info)+". "+TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));
		}
		if(M!=NBANDS)
		{
			throw runtime_error("M="+TO_STRING(M)+". NBANDS="+TO_STRING(NBANDS)+". "+TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));
		}
		if(M!=NZ)
		{
			throw runtime_error("M="+TO_STRING(M)+". NZ="+TO_STRING(NZ)+". "+TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));
		}
		if(NEW_DM==0)
		{
			throw domain_error("NEW_DM must be 1. "+TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));
		}
	}	
    //delete[] Stmp; //LiuXh 20171109
#endif

#ifdef TEST_DIAG
	{
		static int istep = 0;
		{
			ofstream ofs("ekb_"+TO_STRING(istep)+"_"+TO_STRING(MY_RANK));
			for(int ib=0; ib<NBANDS; ++ib)
				ofs<<ekb[ib]<<endl;
		}
		{
			ofstream ofs("wfc-C_"+TO_STRING(istep)+"_"+TO_STRING(MY_RANK));
			ofs<<wfc_2d<<endl;
		}
		{
			ofstream ofs("wfc-F_"+TO_STRING(istep)+"_"+TO_STRING(MY_RANK));
			ofs<<transpose(wfc_2d)<<endl;
		}
		++istep;
	}
#endif

	return;
}


void Pdiag_Double::diago_complex_begin(const int &ik, complex<double> **wfc, ComplexMatrix &wfc_2d,
	complex<double>* ch_mat, complex<double>* cs_mat, double *ekb)
{
    #ifdef TEST_DIAG
   	{
		static int istep = 0;
		auto print_matrix_C = [&](const string &file_name, complex<double>*m)
		{
			ofstream ofs(file_name+"-C_"+TO_STRING(istep)+"_"+TO_STRING(MY_RANK));
			for(int ic=0; ic<ParaO.ncol; ++ic)
			{
				for(int ir=0; ir<ParaO.nrow; ++ir)
				{
					const int index=ic*ParaO.nrow+ir;
					if(std::norm(m[index])>1E-10)
                    {
                        if(std::imag(m[index])>1E-10)
                            ofs<<m[index]<<"\t";
                        else
                            ofs<<std::real(m[index])<<"\t";
                    }
					else
						ofs<<0<<"\t";
				}
				ofs<<endl;
			}
		};
		auto print_matrix_F = [&](const string &file_name, complex<double>*m)
		{
			ofstream ofs(file_name+"-F_"+TO_STRING(istep)+"_"+TO_STRING(MY_RANK));
			for(int ir=0; ir<ParaO.nrow; ++ir)
			{
				for(int ic=0; ic<ParaO.ncol; ++ic)
				{
					const int index=ic*ParaO.nrow+ir;
					if(std::norm(m[index])>1E-10)
                    {
                        if(std::imag(m[index])>1E-10)
                            ofs<<m[index]<<"\t";
                        else
                            ofs<<std::real(m[index])<<"\t";
                    }
					else
						ofs<<0<<"\t";
				}
				ofs<<endl;
			}
		};
		print_matrix_F("H_gamma", ch_mat);
		print_matrix_F("S_gamma", cs_mat);
		print_matrix_C("H_gamma", ch_mat);
		print_matrix_C("S_gamma", cs_mat);
		++istep;
	}
    #endif

#ifdef __MPI
	TITLE("Pdiag_Double","diago_complex_begin");

	char uplo='U';
	const int inc=1;

    int nprocs, myid;
    MPI_Status status;
    MPI_Comm_size(comm_2D, &nprocs);
    MPI_Comm_rank(comm_2D, &myid);

	// parallel diagonalize the
	// H | psi > = S | psi >
	// problem.
	int loc_pos;

	// because the output Stmp will be different from Sloc2, so we need to copy that.
	complex<double>* Stmp = LM.Sdiag2;

	if(KS_SOLVER=="hpseps")
	{
        double *eigen = new double[NLOCAL];
        ZEROS(eigen, NLOCAL);

        assert(loc_size > 0);
        complex<double>* Z = new complex<double>[this->loc_size * NLOCAL];
        ZEROS(Z, this->loc_size * NLOCAL);

        Memory::record("Pdiag_Double","Z",loc_size * NLOCAL,"cdouble");
		int nbands_tmp = NBANDS;
        timer::tick("Diago_LCAO_Matrix","pzgseps",'G');
		LapackConnector::copy(nloc, cs_mat, inc, Stmp, inc);
    	pzgseps(comm_2D, NLOCAL, nb, nbands_tmp, ch_mat, Stmp, Z, eigen, this->MatrixInfo, uplo, this->loc_size, loc_pos);
        timer::tick("Diago_LCAO_Matrix","pzgseps",'G');
        // the eigenvalues.
        LapackConnector::copy(NBANDS, eigen, inc, ekb, inc);
        delete[] eigen;

        // Z is delete in gath_eig
        timer::tick("Diago_LCAO_Matrix","gath_eig_complex",'G');
        this->gath_eig_complex(DIAG_HPSEPS_WORLD, NLOCAL, wfc, Z, ik);
        timer::tick("Diago_LCAO_Matrix","gath_eig_complex",'G');
        //delete[] Z; //LiuXh 20180329, fix bug of 'double free()'
        //this->gath_full_eig_complex(DIAG_WORLD, NLOCAL, c, Z);
	} // HPSEPS method
    else if(KS_SOLVER=="genelpa")
    {
        double *eigen = new double[NLOCAL];
        ZEROS(eigen, NLOCAL);
        long maxnloc; // maximum number of elements in local matrix
        MPI_Reduce(&nloc, &maxnloc, 1, MPI_LONG, MPI_MAX, 0, comm_2D);
        MPI_Bcast(&maxnloc, 1, MPI_LONG, 0, comm_2D);
		wfc_2d.create(this->ncol,this->nrow);			// Fortran order
        complex<double> *work=new complex<double>[maxnloc]; // work/buffer matrix
        bool wantEigenVector=true;
        bool wantDebug=true;
        int info;
        int comm_2D_f=MPI_Comm_c2f(comm_2D);

        int THIS_REAL_ELPA_KERNEL_API=9;
        timer::tick("Diago_LCAO_Matrix","genelpa",'G');
        LapackConnector::copy(nloc, cs_mat, inc, Stmp, inc);
        int method=0;
        //info=pzSolveGenEigen2(NBANDS, NLOCAL, nrow, ncol, desc,
        //                      ch_mat, Stmp, eigen, wfc_2d.c, work,
        //                      comm_2D_f, blacs_ctxt,
        //                      method, THIS_REAL_ELPA_KERNEL_API,
        //                      wantEigenVector, wantDebug);
        //timer::tick("Diago_LCAO_Matrix","genelpa",'G');
        timer::tick("Diago_LCAO_Matrix","genelpa1",'G');
        info=pzDecomposeRightMatrix2(NLOCAL, nrow, ncol, desc,
                                    Stmp, eigen, wfc_2d.c, work,
                                    comm_2D_f, mpi_comm_rows, mpi_comm_cols,
                                    method, THIS_REAL_ELPA_KERNEL_API);
        timer::tick("Diago_LCAO_Matrix","genelpa1",'G');
        timer::tick("Diago_LCAO_Matrix","genelpa2",'G');
        info=pzSolveEigen2(NBANDS, NLOCAL, nrow, ncol, desc,
                    ch_mat, Stmp, eigen, wfc_2d.c, work,
                    comm_2D_f, mpi_comm_rows, mpi_comm_cols, method,
                    THIS_REAL_ELPA_KERNEL_API,
                    wantEigenVector, wantDebug);
        timer::tick("Diago_LCAO_Matrix","genelpa2",'G');

        // the eigenvalues.
        LapackConnector::copy(NBANDS, eigen, inc, ekb, inc);
        delete[] eigen;

        //change eigenvector matrix from block-cycle distribute matrix to column-divided distribute matrix
        timer::tick("Diago_LCAO_Matrix","gath_eig_complex",'G');
        int naroc[2]; // maximum number of row or column
        for(int iprow=0; iprow<dim0; ++iprow)
        {
            for(int ipcol=0; ipcol<dim1; ++ipcol)
            {
                const int coord[2]={iprow, ipcol};
                int src_rank;
                MPI_Cart_rank(comm_2D, coord, &src_rank);
                if(myid==src_rank)
                {
                    LapackConnector::copy(nloc, wfc_2d.c, inc, work, inc);
                    naroc[0]=nrow;
                    naroc[1]=ncol;
                }
                info=MPI_Bcast(naroc, 2, MPI_INT, src_rank, comm_2D);
                info=MPI_Bcast(work, maxnloc, MPI_DOUBLE_COMPLEX, src_rank, comm_2D);

                if(this->out_lowf)
                {
                    complex<double> **ctot;
                    if(myid==0)
                    {
                        ctot = new complex<double>*[NBANDS];
                        for (int i=0; i<NBANDS; i++)
                        {
                            ctot[i] = new complex<double>[NLOCAL];
                            ZEROS(ctot[i], NLOCAL);
                        }
                        Memory::record("Pdiag_Basic","ctot",NBANDS*NLOCAL,"cdouble");
                    }
					// mohan update 2021-02-12, delete BFIELD option
					info=q2WFC_WFCAUG_CTOT_complex(myid, naroc, nb,
							dim0, dim1, iprow, ipcol,
							work, wfc, LOWF.WFC_K_aug[ik], ctot);
                    stringstream ss;
	                ss << global_out_dir << "LOWF_K_" << ik+1 << ".dat";
                    // mohan add 2012-04-03, because we need the occupations for the
                    // first iteration.
                    WF_Local::write_lowf_complex( ss.str(), ctot, ik );//mohan add 2010-09-09
                    if(myid==0)
                    {
                        for (int i=0; i<NBANDS; i++)
                        {
                            delete[] ctot[i];
                        }
                        delete[] ctot;
                    }
                }
                else
                {
					// mohan update 2021-02-12, delte BFIELD option
					info=q2WFC_WFCAUG_complex(naroc, nb,
							dim0, dim1, iprow, ipcol,
							work, wfc, LOWF.WFC_K_aug[ik]);
				}
            }
        }
        delete[] work;
        timer::tick("Diago_LCAO_Matrix","gath_eig_complex",'G');
    } // GenELPA method
	else if(KS_SOLVER=="scalapack_gvx")
	{
		ComplexMatrix h_tmp(this->ncol, this->nrow, false);
		memcpy( h_tmp.c, ch_mat, sizeof(complex<double>)*this->ncol*this->nrow );
		ComplexMatrix s_tmp(this->ncol, this->nrow, false);
		memcpy( s_tmp.c, cs_mat, sizeof(complex<double>)*this->ncol*this->nrow );
		wfc_2d.create(this->ncol, this->nrow, false);
		
		const char jobz='V', range='I', uplo='U';
		const int itype=1, il=1, iu=NBANDS, one=1;
		int M=0, NZ=0, lwork=-1, lrwork=-1, liwork=-1, info=0;
		const double abstol=0, orfac=-1;
		vector<complex<double>> work(1,0);
		vector<double> rwork(1,0);
		vector<int> iwork(1,0);
		vector<int> ifail(NLOCAL,0);
		vector<int> iclustr(2*DSIZE);
		vector<double> gap(DSIZE);
		
		pzhegvx_(&itype, &jobz, &range, &uplo,
			&NLOCAL, h_tmp.c, &one, &one, desc, s_tmp.c, &one, &one, desc,
			NULL, NULL, &il, &iu, &abstol,
			&M, &NZ, ekb, &orfac, wfc_2d.c, &one, &one, desc,
			work.data(), &lwork, rwork.data(), &lrwork, iwork.data(), &liwork, ifail.data(), iclustr.data(), gap.data(), &info);
		ofs_running<<"lwork="<<work[0]<<"\t"<<"liwork="<<iwork[0]<<endl;
		lwork = work[0].real();
		work.resize(lwork,0);
		lrwork = rwork[0];
		rwork.resize(lrwork,0);
		liwork = iwork[0];
		iwork.resize(liwork,0);
		pzhegvx_(&itype, &jobz, &range, &uplo,
			&NLOCAL, h_tmp.c, &one, &one, desc, s_tmp.c, &one, &one, desc,
			NULL, NULL, &il, &iu, &abstol,
			&M, &NZ, ekb, &orfac, wfc_2d.c, &one, &one, desc,
			work.data(), &lwork, rwork.data(), &lrwork, iwork.data(), &liwork, ifail.data(), iclustr.data(), gap.data(), &info);
		ofs_running<<"M="<<M<<"\t"<<"NZ="<<NZ<<endl;

		if(info)
			throw runtime_error("info="+TO_STRING(info)+". "+TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));
		if(M!=NBANDS)
			throw runtime_error("M="+TO_STRING(M)+". NBANDS="+TO_STRING(NBANDS)+". "+TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));
		if(M!=NZ)
			throw runtime_error("M="+TO_STRING(M)+". NZ="+TO_STRING(NZ)+". "+TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));
		
//		if(NEW_DM==0)
//			throw domain_error("NEW_DM must be 1. "+TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));
		// the follow will be deleted after finish newdm
		{
			//change eigenvector matrix from block-cycle distribute matrix to column-divided distribute matrix
			timer::tick("Diago_LCAO_Matrix","gath_eig_complex",'G');

			long maxnloc; // maximum number of elements in local matrix
			MPI_Reduce(&nloc, &maxnloc, 1, MPI_LONG, MPI_MAX, 0, comm_2D);
			MPI_Bcast(&maxnloc, 1, MPI_LONG, 0, comm_2D);
			complex<double> *work=new complex<double>[maxnloc]; // work/buffer matrix

			int naroc[2]; // maximum number of row or column
			for(int iprow=0; iprow<dim0; ++iprow)
			{
				for(int ipcol=0; ipcol<dim1; ++ipcol)
				{
					const int coord[2]={iprow, ipcol};
					int src_rank;
					MPI_Cart_rank(comm_2D, coord, &src_rank);
					if(myid==src_rank)
					{
						LapackConnector::copy(nloc, wfc_2d.c, inc, work, inc);
						naroc[0]=nrow;
						naroc[1]=ncol;
					}
					info=MPI_Bcast(naroc, 2, MPI_INT, src_rank, comm_2D);
					info=MPI_Bcast(work, maxnloc, MPI_DOUBLE_COMPLEX, src_rank, comm_2D);

					if(this->out_lowf)
					{
						complex<double> **ctot;
						if(myid==0)
						{
							ctot = new complex<double>*[NBANDS];
							for (int i=0; i<NBANDS; i++)
							{
								ctot[i] = new complex<double>[NLOCAL];
								ZEROS(ctot[i], NLOCAL);
							}
							Memory::record("Pdiag_Basic","ctot",NBANDS*NLOCAL,"cdouble");
						}
						// mohan update 2021-02-12, delete BFIELD option
						info=q2WFC_WFCAUG_CTOT_complex(myid, naroc, nb,
								dim0, dim1, iprow, ipcol,
								work, wfc, LOWF.WFC_K_aug[ik], ctot);
						stringstream ss;
						ss << global_out_dir << "LOWF_K_" << ik+1 << ".dat";
						// mohan add 2012-04-03, because we need the occupations for the
						// first iteration.
						WF_Local::write_lowf_complex( ss.str(), ctot, ik );//mohan add 2010-09-09
						if(myid==0)
						{
							for (int i=0; i<NBANDS; i++)
							{
								delete[] ctot[i];
							}
							delete[] ctot;
						}
					}
					else
					{
						// mohan update 2021-02-12, delte BFIELD option
						info=q2WFC_WFCAUG_complex(naroc, nb,
								dim0, dim1, iprow, ipcol,
								work, wfc, LOWF.WFC_K_aug[ik]);
					}
				}
			}
			delete[] work;
			timer::tick("Diago_LCAO_Matrix","gath_eig_complex",'G');
		}
	}		

#endif
	return;
}



#ifdef __MPI
void Pdiag_Double::readin(const string &fa, const string &fb, const int &nlocal_tot, double *eigen, double *eigvr)
{
    TITLE("Pdiag_Double","readin");

    int coord[2];
    int dim[2];
    int period[2];
    int i,j,tmp1,tmp2;
    int k,loc_size,loc_pos;
    double time1,time2;

    MPI_Comm comm=DIAG_HPSEPS_WORLD,comm_2D,comm_col,comm_row,newcomm;

    dim[0]=(int)sqrt((double)DSIZE);

    while (DSIZE%dim[0]!=0)
    {
        dim[0]=dim[0]-1;
    }
    dim[1]=DSIZE/dim[0];

    // call mpi_creat_cart
    this->mpi_creat_cart(&comm_2D,dim[0],dim[1]);

    // call mat_2d
    this->mat_2d(comm_2D,nlocal_tot,nlocal_tot,nb,MatrixInfo);

    loc_size=nlocal_tot/DSIZE;
    if (DRANK<nlocal_tot%DSIZE) loc_size=loc_size+1;

    ofs_running << " loc_size = " << loc_size;

    /*Distribute the matrix*/
    const long nloc = MatrixInfo.col_num * MatrixInfo.row_num;

    double *A = new double[nloc];
    double *B = new double[nloc];
    double *Z = new double[loc_size*nlocal_tot];
    ZEROS(A, nloc);
    ZEROS(B, nloc);
    ZEROS(Z, loc_size * nlocal_tot);

    ofs_running << "\n Data distribution of H." << endl;
    this->data_distribution(comm_2D,fa,nlocal_tot,nb,A,MatrixInfo);
    ofs_running << "\n Data distribution of S." << endl;
    this->data_distribution(comm_2D,fb,nlocal_tot,nb,B,MatrixInfo);

    time1=MPI_Wtime();
    // call pdgseps
    char uplo = 'U';
    pdgseps(comm_2D,nlocal_tot,nb,A,B,Z,eigen,MatrixInfo,uplo,loc_size,loc_pos);
    time2=MPI_Wtime();
    OUT(ofs_running,"time1",time1);
    OUT(ofs_running,"time2",time2);

    //this->gath_eig(comm,n,eigvr,Z);

    ofs_running << "\n " << setw(6) << "Band" << setw(25) << "Ry" << setw(25) << " eV" << endl;
    for(int i=0; i<nlocal_tot; i++)
    {
        ofs_running << " " << setw(6) << i << setw(25) << eigen[i] << setw(25)<< eigen[i] * 13.6058 << endl;
    }

    delete[] A;
    delete[] B;
    delete[] Z;
}
#endif

