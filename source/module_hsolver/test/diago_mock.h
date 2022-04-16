#include "../../src_pw/hamilt_pw.h"
#include "../../src_pw/hamilt.h"
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
    double* precond() {return precondition;}
    ModuleBase::ComplexMatrix hamilt() {return hmatrix;}
    //ModuleBase::ComplexMatrix psi() {return psimatrix;}
    psi::Psi<std::complex<double>> psi()
    {
        PW_Basis* pbas;
        int* ngk = nullptr;
        psi::Psi<std::complex<double>> psitmp(ngk,1,nband,npw);
        for(int i=0;i<nband;i++)
	    {
		    for(int j=0;j<npw;j++) psitmp(0,i,j) = psimatrix(i,j);
	    }   
        return psitmp;
    };

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
                if (int (u(e) % 10) > int (sparsity-1)) mincoef = 1.0;
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

