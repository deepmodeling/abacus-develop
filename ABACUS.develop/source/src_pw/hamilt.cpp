#include "global.h"
#include "hamilt.h"
#include "diago_cg.h"
#include "diago_david.h"

Hamilt::Hamilt() {}
Hamilt::~Hamilt() {}


void Hamilt::diago(
    const int &istep,
    const int &iter,
    const int &ik,
    const double *precondition,
    double &avg_iter)
{
	TITLE("Hamilt","diago");
	timer::tick("Hamilt","diago",'F');
    double avg = 0.0;

	// set ik0 because of mem_saver.
	// if mem_saver is not used, ik0=ik, as usual.
	// but if mem_saver is used, ik0=0.
	int ik0 = ik;

	if(CALCULATION=="nscf" && wf.mem_saver==1)
	{
		if(BASIS_TYPE=="pw")
		{
			// generate PAOs first, then diagonalize to get
			// inital wavefunctions.
			wf.diago_PAO_in_pw_k2(ik, wf.evc[0]);	
		}
		else if(BASIS_TYPE=="lcao_in_pw")
		{
			wf.LCAO_in_pw_k(ik, wf.wanf2[0]);
		}
		ik0 = 0;
	}

    if(BASIS_TYPE=="lcao_in_pw")
    {
		if(KS_SOLVER=="lapack")
		{
			assert(NLOCAL >= NBANDS);
        	this->cinitcgg(ik, NLOCAL, NBANDS, wf.wanf2[ik0], wf.evc[ik0], wf.ekb[ik]);
		}
		else
		{
			ofs_warning << " The diago_type " << KS_SOLVER << " not implemented yet." << endl; //xiaohui add 2013-09-02
			WARNING_QUIT("Hamilt::diago","no implemt yet.");
		}
    }
    else
    {
        int ntry = 0;
        int notconv = 0;
        do
        {	
	   		if(KS_SOLVER=="cg")
            {			
                if ( iter > 0 || istep > 0 ||  ntry > 0)
                {
                    this->cinitcgg( ik ,NBANDS, NBANDS, wf.evc[ik0], wf.evc[ik0], wf.ekb[ik]);
                    avg_iter += 1.0;
                }
                Diago_CG cg;
    	
				bool reorder = true;

				if(NPOL==1) cg.diag(wf.evc[ik0], wf.ekb[ik], kv.ngk[ik], wf.npwx,
						NBANDS, precondition, ETHR,
						DIAGO_CG_MAXITER, reorder, notconv, avg);
				else{
					cg.diag(wf.evc[ik0], wf.ekb[ik], wf.npwx*NPOL, wf.npwx*NPOL,
						NBANDS, precondition, ETHR,
						DIAGO_CG_MAXITER, reorder, notconv, avg);
				}
				// P.S. : nscf is the flag about reorder.
				// if cinitcgg is done once,
				// we don't need to reorder the eigenvectors order.
				// if cinitcgg has not been called,
				// we need to reorder the eigenvectors.
            }
	   		else if(KS_SOLVER=="dav")
        	{
				Diago_David david;
				if(NPOL==1) 
				{
					david.diag(wf.evc[ik0], wf.ekb[ik], kv.ngk[ik],
						NBANDS, precondition, DIAGO_DAVID_NDIM,
				 		ETHR, DIAGO_CG_MAXITER, notconv, avg);
				}
				else
				{
					david.diag(wf.evc[ik0], wf.ekb[ik], wf.npwx*NPOL,
						NBANDS, precondition, DIAGO_DAVID_NDIM,
						ETHR, DIAGO_CG_MAXITER, notconv, avg);
				}
        	}
        	else
        	{
				WARNING_QUIT("calculate_bands","Check KS_SOLVER !");
        	}
            avg_iter += avg;
            ++ntry;
        }
        while ( this->test_exit_cond(ntry, notconv) );

        if ( notconv > max(5, NBANDS/4) )
        {
            cout << "\n notconv = " << notconv;
            cout << "\n Hamilt::diago', too many bands are not converged! \n";
        }
    }

	timer::tick("Hamilt","diago",'F');
    return;
}


bool Hamilt::test_exit_cond(const int &ntry, const int &notconv)
{
    //================================================================
    // If this logical function is true, need to do cinitcgg and cg
    // again.
    //================================================================

	bool scf = true;
	if(CALCULATION=="nscf") scf=false;

    // If ntry <=5, try to do it better, if ntry > 5, exit.
    const bool f1 = (ntry <= 5);

    // In non-self consistent calculation, do until totally converged.
    const bool f2 = ( (!scf && (notconv > 0)) );

    // if self consistent calculation, if not converged > 5,
    // using cinitcgg and cg method again. ntry++
    const bool f3 = ( ( scf && (notconv > 5)) );
    return  ( f1 && ( f2 || f3 ) );
}

void Hamilt::cinitcgg(
    const int ik,
    const int nstart,
    const int n_band,
    const ComplexMatrix &psi,
    ComplexMatrix &evc,
    double *en)
{
	if(nstart < n_band)
	{
		WARNING_QUIT("cinitcgg","nstart < n_band!");
	}

    if(BASIS_TYPE=="pw" || BASIS_TYPE=="lcao_in_pw")
    {
        this->hpw.cinitcgg(ik, nstart, n_band, psi, evc, en);
    }
    else
    {
		WARNING_QUIT("cinitcgg","Check parameters: BASIS_TYPE. ");
    }
    return;
}



void Hamilt::cdiaghg(
	const int nstart,
	const int nbands,
	const ComplexMatrix &hc,
	const ComplexMatrix &sc,
	const int ldh, // nstart
	double *e,
	ComplexMatrix &hvec)
{
    TITLE("Hamilt","cdiaghg");
	timer::tick("Hamilt","cdiaghg");
    //====================================================================
    // calculates eigenvalues and eigenvectors of the generalized problem
    // Hv=eSv, with H hermitean matrix, S overlap matrix .
    // On output both matrix are unchanged
    // LAPACK version - uses both ZHEGV and ZHEGVX
    //=====================================================================

    int lwork;
    //========================================
    // int ILAENV();
    // ILAENV returns optimal block size "nb"
    //========================================

    ComplexMatrix sdum(nstart, ldh);
    ComplexMatrix hdum;

    sdum = sc;

    const bool all_eigenvalues = (nstart == nbands);

    int nb = LapackConnector::ilaenv(1, "ZHETRD", "U", nstart, -1, -1, -1);
//  int nb = ILAENV(1,  "ZHETRD", "U", n, -1, -1, -1);
//  int nb = 32;

    if (nb < 1)
    {
        nb = std::max(1, nstart);
    }
    
	if (nb == 1 || nb >= nstart)
    {
        lwork = 2 * nstart; // mohan modify 2009-08-02
    }
    else
    {
        lwork = (nb + 1) * nstart;
    }

    complex<double> *work = new complex<double>[lwork];
	ZEROS(work, lwork);
	
    //=====================================================================
    // input s and (see below) h are copied so that they are not destroyed
    //=====================================================================

    int info = 0;
    int rwork_dim;
    if (all_eigenvalues)
    {
        rwork_dim = 3*nstart-2;
    }
    else
    {
        rwork_dim = 7*nstart;
    }
	
    double *rwork = new double[rwork_dim];
    ZEROS( rwork, rwork_dim );

    if (all_eigenvalues)
    {
        //===========================
        // calculate all eigenvalues
        //===========================
        hvec = hc;
        LapackConnector::zhegv(1, 'V', 'U', nstart, hvec , ldh, sdum, ldh, e, work , lwork , rwork, info);
    }
    else
    {
        //=====================================
        // calculate only m lowest eigenvalues
        //=====================================
        double *work2 = new double[2*lwork];
        int *iwork = new int [5*nstart];
        int *ifail = new int[nstart];

        ZEROS(work2,2*lwork);
        ZEROS(rwork,7*nstart);
        ZEROS(iwork,5*nstart);
        ZEROS(ifail,nstart);

        hdum.create(nstart, ldh);
        hdum = hc;

    	//=============================
    	// Number of calculated bands
    	//=============================
    	int mm = nbands;

        LapackConnector::zhegvx
        (
                1,      //INTEGER
                'V',    //CHARACTER*1
                'I',    //CHARACTER*1
                'U',    //CHARACTER*1
                nstart, //INTEGER
                hdum,   //COMPLEX*16 array
                ldh,    //INTEGER
                sdum,   //COMPLEX*16 array
                ldh,    //INTEGER
           		0.0,    //DOUBLE PRECISION
                0.0,    //DOUBLE PRECISION
                1,      //INTEGER
                nbands, //INTEGER
                0.0,    //DOUBLE PRECISION
                mm,     //INTEGER
                e,      //DOUBLE PRECISION array
                hvec,   //COMPLEX*16 array
                ldh,    //INTEGER
                work,   //DOUBLE array, dimension (MAX(1,LWORK))
                lwork,  //INTEGER
                rwork , //DOUBLE PRECISION array, dimension (7*N)
                iwork,  //INTEGER array, dimension (5*N)
                ifail,  //INTEGER array, dimension (N)
                info    //INTEGER
        );

        delete[] work2;
        delete[] iwork;
        delete[] ifail;
    }
    delete[] rwork;
    delete[] work;

	timer::tick("Hamilt","cdiaghg");
    return;
}
