#include<random>
#include "../../module_base/lapack_connector.h"
#include "../../module_base/blas_connector.h"
#include "mpi.h"
#include "module_base/parallel_reduce.h"
#include "module_hamilt_pw/hamilt_pwdft/hamilt_pw.h"

namespace DIAGOTEST
{
    std::vector<double> hmatrix_d;
    std::vector<double> hmatrix_local_d;
    std::vector<std::complex<double>> hmatrix;
    std::vector<std::complex<double>> hmatrix_local;
    std::vector<std::complex<float>> hmatrix_f;
    std::vector<std::complex<float>> hmatrix_local_f;
    int h_nr;
    int h_nc;
    int npw;
    int* npw_local; //number of plane wave distributed to each process

    void readh(std::ifstream& inf, std::vector<double>& hm)
    {
        int npw_;
        inf >> npw_;
        DIAGOTEST::npw = npw_;
        hm.resize(npw * npw);
        h_nr = npw;
        h_nc = npw;
        for (int i = 0;i < npw;i++)
        {
            for (int j = i;j < npw;j++)
            {
                inf >> hm[i * npw + j];
                if (i != j) { hm[j * npw + i] = hm[i * npw + j]; }
            }
        }
    }

    void readh(std::ifstream& inf, std::vector<std::complex<double>>& hm)
    {
        int npw_;
        inf >> npw_;
        DIAGOTEST::npw = npw_;
        hm.resize(npw * npw);
        h_nr = npw;
        h_nc = npw;
        for(int i=0;i<npw;i++)
        {
            for(int j=i;j<npw;j++)
            {
                inf >> hm[i * npw + j];
                if (i != j) {hm[j * npw + i] = conj(hm[i * npw + j]);}                
            }
            double real = hm[i * npw + i].real();
            hm[i * npw + i] = std::complex<double> {real,0.0};
        }
    }

    void readh(std::ifstream &inf, std::vector<std::complex<float>> &hm)
    {
        int npw_;
        inf >> npw_;
        DIAGOTEST::npw = npw_;
        hm.resize(npw * npw);
        h_nr = npw;
        h_nc = npw;
        for(int i=0;i<npw;i++)
        {
            for(int j=i;j<npw;j++)
            {
                inf >> hm[i * npw + j];
                if (i != j) {hm[j * npw + i] = conj(hm[i * npw + j]);}                
            }
            float real = hm[i * npw + i].real();
            hm[i * npw + i] = std::complex<float> {real,0.0};
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
                else if (std::is_same<T, float>::value)
                    MPI_Send(&(psi[start_point]), npw_local[j], MPI_FLOAT, j, 0, MPI_COMM_WORLD);
                else if (std::is_same<T, std::complex<float>>::value)
                    MPI_Send(&(psi[start_point]), npw_local[j], MPI_C_FLOAT_COMPLEX, j, 0, MPI_COMM_WORLD);
                start_point += npw_local[j];
            }
        }
        else
        {
            int recv_len = mypnum < (npw%nprocs) ? npw/nprocs + 1  : npw/nprocs;
            if (std::is_same<T, double>::value) 
                MPI_Recv(psi_local, npw_local[mypnum],MPI_DOUBLE,0,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            else if(std::is_same<T, std::complex<double>>::value)
                MPI_Recv(psi_local, npw_local[mypnum], MPI_DOUBLE_COMPLEX, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            else if (std::is_same<T, float>::value)
                MPI_Recv(psi_local, npw_local[mypnum], MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            else if (std::is_same<T, std::complex<float>>::value)
                MPI_Recv(psi_local, npw_local[mypnum], MPI_C_FLOAT_COMPLEX, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#endif                 
        }     
    }
    template void divide_psi<double>(double* psi, double* psi_local);
    template void divide_psi<std::complex<double>>(std::complex<double>* psi, std::complex<double>* psi_local);
    template void divide_psi<std::complex<float>>(std::complex<float>* psi, std::complex<float>* psi_local);
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

    template<typename T>
    void divide_hpsi(psi::Psi<T>& psi, psi::Psi<T>& psi_local, std::vector<T>& hmatrix, std::vector<T>& hmatrix_local)
    {
        int nprocs, mypnum;
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &mypnum);

        int npw = psi.get_nbasis();
        int nbands = psi.get_nbands();
        int nk = psi.get_nk();

        hmatrix_local.resize(npw * npw_local[mypnum]);
        h_nc = npw_local[mypnum];
        h_nr = npw;
        psi_local.resize(nk,nbands,npw_local[mypnum]);
        for(int i=0;i<npw;i++)
        {
            divide_psi<T>(&(hmatrix[i * npw]), &(hmatrix_local[i * npw_local[mypnum]]));
            if(i<nbands) 
            {
                for(int k=0;k<nk;k++)
                {
                    psi.fix_k(k);
                    psi_local.fix_k(k);
                    divide_psi<T>(psi.get_pointer(i), psi_local.get_pointer(i));
                }
            }
        }
    }
#endif
}

#include "module_base/macros.h"
template<typename T>
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
     *  - hamilt(): return the Hermite matrix (type: std::vector<T>)
     *  - precond(): return the precondition (type: Real Pointer)
     * 
     */
    using Real = typename GetTypeReal<T>::type;
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
    Real* precond() { return precondition; }
    std::vector<T> hamilt() { return hmatrix; }
    //ModuleBase::ComplexMatrix psi() {return psimatrix;}
    psi::Psi<T> psi()
    {
        Structure_Factor* sf;
        int* ngk = nullptr;
        psi::Psi<T> psitmp(1, nband, npw, ngk);
        for(int i=0;i<nband;i++)
	    {
		    for(int j=0;j<npw;j++) psitmp(0,i,j) = psimatrix[i * npw + j];
	    }   
        return psitmp;
    };

    //generate the Hermite matrix
    void genhmatrix();

    //generate the psi matrix
    void genpsi();

    //generate precondition
    void genprecondition()
    {
        precondition = new Real[npw];
        std::default_random_engine e(1000);
        std::uniform_int_distribution<unsigned> u(min, max);
        for (int i = 0;i < npw;i++)
        {
            precondition[i] = 1.0 + static_cast<Real>(u(e)) / max;
        }
    }

private:
    int npw;
    int nband;
    std::vector<T> hmatrix;
    std::vector<T> psimatrix;
    Real* precondition;
    int min = 0;
    int max = 9999;
    int sparsity;
};

template<>
void HPsi<double>::genhmatrix()
{
    hmatrix.resize(npw * npw);
    DIAGOTEST::h_nr = npw;
    DIAGOTEST::h_nc = npw;
    std::default_random_engine e(100);
    std::uniform_int_distribution<unsigned> u(min, max);
    if (sparsity < 0) sparsity = 0;
    if (sparsity > 10) sparsity = 10;
    for (int i = 0;i < npw;i++)
    {
        for (int j = 0;j <= i;j++)
        {
            double mincoef = 0.0;
            double realp = pow(-1.0, u(e) % 2) * static_cast<double>(u(e)) / max;
            if (int(u(e) % 10) > int(sparsity - 1)) mincoef = 1.0;
            if (i == j)
                hmatrix[i * npw + j] = realp;
            else
                hmatrix[j * npw + i] = hmatrix[i * npw + j] = mincoef * realp;
        }
    }
}
template<>
void HPsi<std::complex<double>>::genhmatrix()
{
    hmatrix.resize(npw * npw);
    DIAGOTEST::h_nr = npw;
    DIAGOTEST::h_nc = npw;
    std::default_random_engine e(100);
    std::uniform_int_distribution<unsigned> u(min, max);
    if (sparsity < 0) sparsity = 0;
    if (sparsity > 10) sparsity = 10;
    for (int i = 0;i < npw;i++)
    {
        for (int j = 0;j <= i;j++)
        {
            double mincoef = 0.0;
            double realp = pow(-1.0, u(e) % 2) * static_cast<double>(u(e)) / max;
            //double imagp= pow(-1.0,u(e)%2) * static_cast<double>(u(e))/max;
            if (int(u(e) % 10) > int(sparsity - 1)) mincoef = 1.0;
            if (i == j)
            {
                hmatrix[i * npw + j] = std::complex<double>{ realp,0.0 };
            }
            else
            {
                //hmatrix(i,j) = mincoef*std::complex<double>{realp,imagp};
                hmatrix[i * npw + j] = mincoef * std::complex<double>{ realp, 0.0 };
                hmatrix[j * npw + i] = conj(hmatrix[i * npw + j]);
            }
        }
    }
}
template<>
void HPsi<std::complex<float>>::genhmatrix()
{
    hmatrix.resize(npw * npw);
    DIAGOTEST::h_nr = npw;
    DIAGOTEST::h_nc = npw;
    std::default_random_engine e(100);
    std::uniform_int_distribution<unsigned> u(min, max);
    if (sparsity < 0) sparsity = 0;
    if (sparsity > 10) sparsity = 10;
    for (int i = 0;i < npw;i++)
    {
        for (int j = 0;j <= i;j++)
        {
            float mincoef = 0.0;
            float realp = pow(-1.0, u(e) % 2) * static_cast<float>(u(e)) / max;
            if (int(u(e) % 10) > int(sparsity - 1)) mincoef = 1.0;
            if (i == j)
            {
                hmatrix[i * npw + j] = std::complex<float>{ realp,0.0 };
            }
            else
            {
                hmatrix[i * npw + j] = mincoef * std::complex<float>{ realp, 0.0 };
                hmatrix[j * npw + i] = conj(hmatrix[i * npw + j]);
            }
        }
    }
}

template<>
void HPsi<double>::genpsi()
{
    {
        psimatrix.resize(nband * npw);
        std::default_random_engine e(10);
        std::uniform_int_distribution<unsigned> u(min, max);
        for (int i = 0;i < nband;i++)
            for (int j = 0;j < npw;j++)
                psimatrix[i * npw + j] = pow(-1.0, u(e) % 2) * static_cast<double>(u(e)) / max;
    }
}
template<>
void HPsi<std::complex<double>>::genpsi()
{
    {
        psimatrix.resize(nband * npw);
        std::default_random_engine e(10);
        std::uniform_int_distribution<unsigned> u(min, max);
        for (int i = 0;i < nband;i++)
        {
            for (int j = 0;j < npw;j++)
            {
                double realp = pow(-1.0, u(e) % 2) * static_cast<double>(u(e)) / max;
                //double imagp=pow(-1.0,u(e)%2) * static_cast<double>(u(e))/max;
                double imagp = 0.0;
                psimatrix[i * npw + j] = std::complex<double>{ realp,imagp };
            }
        }
    }
}
template<>
void HPsi<std::complex<float>>::genpsi()
{
    {
        psimatrix.resize(nband * npw);
        std::default_random_engine e(10);
        std::uniform_int_distribution<unsigned> u(min, max);
        for (int i = 0;i < nband;i++)
        {
            for (int j = 0;j < npw;j++)
            {
                float realp = pow(-1.0, u(e) % 2) * static_cast<float>(u(e)) / max;
                float imagp = 0.0;
                psimatrix[i * npw + j] = std::complex<float>{ realp,imagp };
            }
        }
    }
}

template class HPsi<double>;
template class HPsi<std::complex<double>>;
template class HPsi<std::complex<float>>;

//totally same as the original function
template<> void hamilt::Hamilt<double>::sPsi
(
    const double* psi,
    double* spsi,
    const size_t size
)const
{
    for (size_t i = 0; i < size; i++)
    {
        spsi[i] = psi[i];
    }
    return;
}
template<> void hamilt::Hamilt<std::complex<double>>::sPsi
(
    const std::complex<double>* psi,
    std::complex<double>* spsi,
    const size_t size
)const
{
    for (size_t i = 0; i < size; i++)
    {
        spsi[i] = psi[i];
    }
    return;
}
template<> void hamilt::Hamilt<std::complex<float>>::sPsi
(
    const std::complex<float>* psi,
    std::complex<float>* spsi,
    const size_t size
)const
{
    for (size_t i = 0; i < size; i++)
    {
        spsi[i] = psi[i];
    }
    return;
}

//Mock function h_psi
#include "module_hamilt_pw/hamilt_pwdft/operator_pw/operator_pw.h"
template<typename T>
class OperatorMock : public hamilt::Operator<T>
{
    ~OperatorMock()
    {
        if(this->hpsi != nullptr) 
        {
            delete this->hpsi;
            this->hpsi = nullptr;
        }
    }
    virtual void act
    (
        const int nbands,
        const int nbasis,
        const int npol,
        const T* tmpsi_in,
        T* tmhpsi,
        const int ngk_ik = 0)const;
};
template<>
void OperatorMock<double>::act(
    const int nbands,
    const int nbasis,
    const int npol,
    const double* tmpsi_in,
    double* tmhpsi,
    const int ngk_ik)const
{
    int nprocs = 1, mypnum = 0;
#ifdef __MPI    
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &mypnum);
#endif        

    double* hpsi0 = new double[DIAGOTEST::npw];
    for (int m = 0; m < nbands; m++)
    {
        for (int i = 0;i < DIAGOTEST::npw;i++)
        {
            hpsi0[i] = 0.0;
            for (int j = 0;j < (DIAGOTEST::npw_local[mypnum]);j++)
            {
                hpsi0[i] += DIAGOTEST::hmatrix_local_d[i * DIAGOTEST::h_nc + j] * tmpsi_in[j];
            }
        }
        Parallel_Reduce::reduce_pool(hpsi0, DIAGOTEST::npw);
        DIAGOTEST::divide_psi<double>(hpsi0, tmhpsi);
        tmhpsi += nbasis;
        tmpsi_in += nbasis;
    }
    delete[] hpsi0;
}
template<>
void OperatorMock<std::complex<double>>::act(
    const int nbands,
    const int nbasis,
    const int npol,
    const std::complex<double>* tmpsi_in,
    std::complex<double>* tmhpsi,
    const int ngk_ik)const
{
    int nprocs = 1, mypnum = 0;
#ifdef __MPI    
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &mypnum);
#endif        

    std::complex<double>* hpsi0 = new std::complex<double>[DIAGOTEST::npw];
    for (int m = 0; m < nbands; m++)
    {
        for (int i = 0;i < DIAGOTEST::npw;i++)
        {
            hpsi0[i] = 0.0;
            for (int j = 0;j < (DIAGOTEST::npw_local[mypnum]);j++)
            {
                hpsi0[i] += DIAGOTEST::hmatrix_local[i * DIAGOTEST::h_nc + j] * tmpsi_in[j];
            }
        }
        Parallel_Reduce::reduce_pool(hpsi0, DIAGOTEST::npw);
        DIAGOTEST::divide_psi<std::complex<double>>(hpsi0, tmhpsi);
        tmhpsi += nbasis;
        tmpsi_in += nbasis;
    }
    delete[] hpsi0;
}
template<>
void OperatorMock<std::complex<float>>::act(
    const int nbands,
    const int nbasis,
    const int npol,
    const std::complex<float>* tmpsi_in,
    std::complex<float>* tmhpsi,
    const int ngk_ik)const
{
    int nprocs = 1, mypnum = 0;
#ifdef __MPI    
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &mypnum);
#endif        

    std::complex<float>* hpsi0 = new std::complex<float>[DIAGOTEST::npw];
    for (int m = 0; m < nbands; m++)
    {
        for (int i = 0;i < DIAGOTEST::npw;i++)
        {
            hpsi0[i] = 0.0;
            for (int j = 0;j < (DIAGOTEST::npw_local[mypnum]);j++)
            {
                hpsi0[i] += DIAGOTEST::hmatrix_local_f[i * DIAGOTEST::h_nc + j] * tmpsi_in[j];
            }
        }
        Parallel_Reduce::reduce_pool(hpsi0, DIAGOTEST::npw);
        DIAGOTEST::divide_psi<std::complex<float>>(hpsi0, tmhpsi);
        tmhpsi += nbasis;
        tmpsi_in += nbasis;
    }
    delete[] hpsi0;
}
template<> void hamilt::HamiltPW<double>::updateHk(const int ik)
{
    return;
}

template<> hamilt::HamiltPW<double>::HamiltPW(elecstate::Potential* pot_in, ModulePW::PW_Basis_K* wfc_basis, K_Vectors* pkv)
{
    this->ops = new OperatorMock<double>;
}

template<> hamilt::HamiltPW<double>::~HamiltPW()
{
    delete this->ops;
}

template<> void hamilt::HamiltPW<std::complex<double>>::updateHk(const int ik)
{
    return;
}

template<> hamilt::HamiltPW<std::complex<double>>::HamiltPW(elecstate::Potential* pot_in, ModulePW::PW_Basis_K* wfc_basis, K_Vectors* pkv)
{
    this->ops = new OperatorMock<std::complex<double>>;
}

template<> hamilt::HamiltPW<std::complex<double>>::~HamiltPW()
{
    delete this->ops;
}

template<> void hamilt::HamiltPW<std::complex<float>>::updateHk(const int ik)
{
    return;
}

template<> hamilt::HamiltPW<std::complex<float>>::HamiltPW(elecstate::Potential* pot_in, ModulePW::PW_Basis_K* wfc_basis, K_Vectors* pkv)
{
    this->ops = new OperatorMock<std::complex<float>>;
}

template<> hamilt::HamiltPW<std::complex<float>>::~HamiltPW()
{
    delete this->ops;
}
