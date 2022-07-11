#include<random>
#include "../../module_base/lapack_connector.h"
#include "../../module_base/blas_connector.h"
#include "mpi.h"
#include "src_parallel/parallel_reduce.h"
#include "module_hamilt/hamilt_pw.h"

namespace DIAGOTEST
{
    ModuleBase::ComplexMatrix hmatrix;
    ModuleBase::ComplexMatrix hmatrix_local;
    int npw;
    int* npw_local; //number of plane wave distributed to each process

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

    template<class T>
    void divide_psi(T *psi, T *psi_local)
    {
        int nprocs=1, mypnum=0;
#ifdef __MPI        
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &mypnum);
#endif                

        if(mypnum == 0)
        {                
            for(int j=0;j<npw_local[0];j++) psi_local[j] = psi[j];
#ifdef __MPI                                        
            int start_point = npw_local[0];
            for(int j=1;j<nprocs;j++)
            {
                if (std::is_same<T, double>::value) 
                    MPI_Send(&(psi[start_point]),npw_local[j],MPI_DOUBLE,j,0,MPI_COMM_WORLD);
                else if(std::is_same<T, std::complex<double>>::value) 
                    MPI_Send(&(psi[start_point]),npw_local[j],MPI_DOUBLE_COMPLEX,j,0,MPI_COMM_WORLD);
                start_point += npw_local[j];
            }  
        }
        else
        {
            int recv_len = mypnum < (npw%nprocs) ? npw/nprocs + 1  : npw/nprocs;
            if (std::is_same<T, double>::value) 
                MPI_Recv(psi_local, npw_local[mypnum],MPI_DOUBLE,0,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            else if(std::is_same<T, std::complex<double>>::value)
                MPI_Recv(psi_local, npw_local[mypnum],MPI_DOUBLE_COMPLEX,0,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#endif                 
        }     
    }
 

#ifdef __MPI
    void cal_division(int &npw)
    {
        int nprocs, mypnum;
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &mypnum);
        for(int i=0;i<nprocs;i++)
        {
            if(i<npw%nprocs) npw_local[i] = npw/nprocs + 1;
            else npw_local[i] = npw/nprocs;
        }
    }

    void divide_hpsi(psi::Psi<std::complex<double>> &psi, psi::Psi<std::complex<double>> &psi_local)
    {
        int nprocs, mypnum;
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &mypnum);

        int npw = psi.get_nbasis();
        int nbands = psi.get_nbands();
        int nk = psi.get_nk();

        hmatrix_local.create(npw,npw_local[mypnum]);
        psi_local.resize(nk,nbands,npw_local[mypnum]);
        for(int i=0;i<npw;i++)
        {
            divide_psi<std::complex<double>>(&(hmatrix.c[i*npw]),&(hmatrix_local.c[i*npw_local[mypnum]]));
            if(i<nbands) 
            {
                for(int k=0;k<nk;k++)
                {
                    psi.fix_k(k);
                    psi_local.fix_k(k);
                    divide_psi<std::complex<double>>(psi.get_pointer(i),psi_local.get_pointer(i));
                }
            }
        }
    } 
#endif

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
        Structure_Factor* sf;
        int* ngk = nullptr;
        psi::Psi<std::complex<double>> psitmp(1,nband,npw,ngk);
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
                //double imagp=pow(-1.0,u(e)%2) * static_cast<double>(u(e))/max;
                double imagp = 0.0;
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
void hamilt::HamiltPW::sPsi
(
    const std::complex<double> *psi, 
    std::complex<double> *spsi, 
    const size_t size
)const
{
    for (size_t i=0; i<size; i++)
    {
        spsi[i] = psi[i];
    }
    return;
}

//Mock function h_psi
void hamilt::HamiltPW::hPsi(const std::complex<double> *psi_in, std::complex<double> *hpsi_local, const size_t size)const
{
    int m = size / DIAGOTEST::npw;
    int nprocs=1, mypnum=0;
#ifdef __MPI    
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &mypnum);
#endif        

    std::complex<double> *hpsi = new std::complex<double>[DIAGOTEST::npw];
    for(int i=0;i<DIAGOTEST::npw;i++)
    {
        hpsi[i] = 0.0;
        for(int j=0;j<(DIAGOTEST::npw_local[mypnum]);j++)
        {
            hpsi[i] += DIAGOTEST::hmatrix_local(i,j) * psi_in[j];
        }
    }
    Parallel_Reduce::reduce_complex_double_pool(hpsi, DIAGOTEST::npw);
    DIAGOTEST::divide_psi<std::complex<double>>(hpsi,hpsi_local);
    delete [] hpsi;
}

void hamilt::HamiltPW::updateHk(const int ik)
{
    return;
}

hamilt::HamiltPW::HamiltPW()
{
}

hamilt::HamiltPW::~HamiltPW()
{
}
