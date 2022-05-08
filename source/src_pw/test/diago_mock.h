#include "../hamilt_pw.h"
#include "../hamilt.h"
#include<random>
#include "../../module_base/lapack_connector.h"
#include "../../module_base/blas_connector.h"

Hamilt_PW::Hamilt_PW() {};
Hamilt_PW::~Hamilt_PW() {};
Hamilt::Hamilt() {};
Hamilt::~Hamilt() {};

namespace DIAGOTEST
{
    ModuleBase::ComplexMatrix hmatrix;
    int npw;

    void readh(std::ifstream &inf, ModuleBase::ComplexMatrix &hm)
    {
        int npw;
        inf >> npw;
        hm.create(npw,npw);
        for(int i=0;i<npw;i++)
        {
            for(int j=i;j<npw;j++)
            {
                inf >> hm(i,j);
                if (i != j) {hm(j,i) = conj(hm(i,j));}                
            }
            double real = hm(i,i).real();
            hm(i,i) = std::complex<double> {real,0.0};
        }
    }

	//totaly same as the original code.
	void diagH_LAPACK(
		const int nstart,
		const int nbands,
		const ModuleBase::ComplexMatrix &hc,
		const ModuleBase::ComplexMatrix &sc,
		const int ldh, // nstart
		double *e,
		ModuleBase::ComplexMatrix &hvec)
	{
	    int lwork=0;
	
	    ModuleBase::ComplexMatrix sdum(nstart, ldh);
	    ModuleBase::ComplexMatrix hdum;
	
	    sdum = sc;
	
	    const bool all_eigenvalues = (nstart == nbands);
	
	    //workspace query
	    int nb = LapackConnector::ilaenv(1, "ZHETRD", "U", nstart, -1, -1, -1);
	
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
	
	    std::complex<double> *work = new std::complex<double>[lwork];
		ModuleBase::GlobalFunc::ZEROS(work, lwork);
		
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
	    ModuleBase::GlobalFunc::ZEROS( rwork, rwork_dim );
	
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
	        int *iwork = new int [5*nstart];
	        int *ifail = new int[nstart];
	
	        ModuleBase::GlobalFunc::ZEROS(rwork,7*nstart);
	        ModuleBase::GlobalFunc::ZEROS(iwork,5*nstart);
	        ModuleBase::GlobalFunc::ZEROS(ifail,nstart);
	
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
	
	        delete[] iwork;
	        delete[] ifail;
	    }
	    delete[] rwork;
	    delete[] work;
	
		return;
	}
}


class HPsi
{
    /**
     * This calss used to produce the Hermite matrix, the initial 
     * guess wave function, and the precondition by the random 
     * number. The elements of Hermite matrix and wave function are
     * between -1.0 to 1.0, and the preconddition is between 1.0 to 2.0.
     * 
     * The parameters in construct function or function create()
     * are same:
     *  - int nband/nbd: number of calculated bands
     *  - int npw: number of plane wave
     *  - int sparsity: the sparsity of Halmit matrix, between 0 and 10. 
     *                  (0 means no sparsity, 10 means a diagonal matrix)
     * 
     * After instantiation a HPsi, one can use below functions:
     *  - hamilt(): return the Hermite matrix (type: ModuleBase::ComplexMatrix)
     *  - psi(): return the wavefunction (type: ModuleBase::ComplexMatrix) 
     *  - precond(): return the precondition (type: double Pointer) 
     * 
     */

    public:
    HPsi(int nband,int npw, int sparsity=7):
        nband(nband),npw(npw),sparsity(sparsity) {genhmatrix();genpsi();genprecondition();}
    HPsi(){};
    ~HPsi() {delete [] precondition;}

    void create(int nbd,int npw, int sparsity=7)
    {
        this->nband = nbd;
        this->npw = npw;
        this->sparsity = sparsity;
        genhmatrix();genpsi();genprecondition();
    }

    //return the matrix
    ModuleBase::ComplexMatrix hamilt() {return hmatrix;}
    ModuleBase::ComplexMatrix psi() {return psimatrix;}
    double* precond() {return precondition;}

    //generate the Hermite matrix
    void genhmatrix()
    {
        hmatrix.create(npw,npw);
        std::default_random_engine e(100);
        std::uniform_int_distribution<unsigned> u(min,max);
        if (sparsity < 0) sparsity = 0;
        if (sparsity > 10) sparsity = 10;
        for(int i=0;i<npw;i++)
        {
            for(int j=0;j<=i;j++)
            {
                double mincoef = 0.0;
                double realp= pow(-1.0,u(e)%2) * static_cast<double>(u(e))/max;
                //double imagp= pow(-1.0,u(e)%2) * static_cast<double>(u(e))/max;
                if (u(e) % 10 > (sparsity-1)) mincoef = 1.0;
                if(i==j)
                {
                    hmatrix(i,j) = std::complex<double>{realp,0.0};
                }
                else
                {
                    //hmatrix(i,j) = mincoef*std::complex<double>{realp,imagp};
                    hmatrix(i,j) = mincoef*std::complex<double>{realp,0.0};
                    hmatrix(j,i) = conj(hmatrix(i,j));
                }
            }
        }
    }

    //generate the psi matrix
    void genpsi()
    {
        psimatrix.create(nband,npw);
        std::default_random_engine e(10);
        std::uniform_int_distribution<unsigned> u(min,max);
        for(int i=0;i<nband;i++)
        {
            for(int j=0;j<npw;j++)
            {
                double realp=pow(-1.0,u(e)%2) * static_cast<double>(u(e))/max;
                double imagp=pow(-1.0,u(e)%2) * static_cast<double>(u(e))/max;
                psimatrix(i,j) = std::complex<double>{realp,imagp};
            }
        }
    }

    //generate precondition
    void genprecondition()
    {
        precondition = new double [npw];
        std::default_random_engine e(1000);
        std::uniform_int_distribution<unsigned> u(min,max);
        for(int i=0;i<npw;i++) 
        {
            precondition[i] = 1.0 + static_cast<double>(u(e))/max;
        }
    }

    private:
    int npw;
    int nband;
    ModuleBase::ComplexMatrix hmatrix;
    ModuleBase::ComplexMatrix psimatrix;
    double* precondition;
    int min=0;
    int max=9999;
    int sparsity;
};

//totally same as the original function
void Hamilt_PW::h_1psi( const int npw_in, const std::complex < double> *psi, 
		       std::complex<double> *hpsi, std::complex < double> *spsi)
{
    this->h_psi(psi, hpsi);

    for (int i=0;i<npw_in;i++)
    {
        spsi[i] = psi[i];
    }
    return;
}

//totally same as the original function
void Hamilt_PW::s_1psi
(
    const int dim,
    const std::complex<double> *psi,
    std::complex<double> *spsi
)
{
    for (int i=0; i<dim; i++)
    {
        spsi[i] = psi[i];
    }
    return;
}

//Mock function h_psi
void Hamilt_PW::h_psi(const std::complex<double> *psi_in, std::complex<double> *hpsi, const int m)
{
    for(int i=0;i<DIAGOTEST::npw;i++){
        hpsi[i] = 0.0;
        for(int j=0;j<DIAGOTEST::npw;j++)
        {
            hpsi[i] += DIAGOTEST::hmatrix(i,j) * psi_in[j];
        }
    }
}

//totaly same as the original code.
void Hamilt_PW::diagH_subspace(
		const int ik,
		const int n_start,
		const int n_band,
                const ModuleBase::ComplexMatrix &psi,
		ModuleBase::ComplexMatrix &evc,
                double *en)
{
	std::complex<double> one;
	std::complex<double> zero;
        one = std::complex<double>{1.0, 0.0};
        zero = std::complex<double>{0.0, 0.0};
    	ModuleBase::ComplexMatrix hc(n_start,n_band);
    	ModuleBase::ComplexMatrix sc(n_start,n_band);
    	ModuleBase::ComplexMatrix hvec(n_start,n_band);

	std::complex<double> *aux=new std::complex<double> [DIAGOTEST::npw*n_band];

	this->h_psi(psi.c, aux, n_band);
	char trans1 = 'C';
	char trans2 = 'N';
	zgemm_(&trans1,&trans2,&n_band,&n_band,&(DIAGOTEST::npw),&one,psi.c,&(DIAGOTEST::npw),aux,&(DIAGOTEST::npw),&zero,hc.c,&n_band);
	hc=transpose(hc,false);

	zgemm_(&trans1,&trans2,&n_band,&n_band,&(DIAGOTEST::npw),&one,psi.c,&(DIAGOTEST::npw),psi.c,&(DIAGOTEST::npw),&zero,sc.c,&n_band);
	sc=transpose(sc,false);

	delete []aux;

	DIAGOTEST::diagH_LAPACK(n_band, n_band, hc, sc, n_band, en, hvec);

	char transa = 'N';
	char transb = 'T';
	ModuleBase::ComplexMatrix evctmp(n_band,DIAGOTEST::npw,false);
	zgemm_(&transa,&transb,&(DIAGOTEST::npw),&n_band,&n_band,&one,psi.c,&(DIAGOTEST::npw),hvec.c,&n_band,&zero,evctmp.c,&(DIAGOTEST::npw));
	for(int ib=0; ib<n_band; ib++)
	{
		for(int ig=0; ig<DIAGOTEST::npw; ig++)
		{
			evc(ib,ig) = evctmp(ib,ig);
		}
	}
}