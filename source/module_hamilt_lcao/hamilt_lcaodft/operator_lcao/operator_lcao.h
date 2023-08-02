#ifndef OPERATORLCAO_H
#define OPERATORLCAO_H
#include "module_base/vector3.h"
#include "module_hamilt_general/matrixblock.h"
#include "module_hamilt_general/operator.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_matrix.h"

namespace hamilt
{

template <typename T>
class OperatorLCAO : public Operator<T>
{
  public:
    OperatorLCAO(LCAO_Matrix* LM_in, const std::vector<ModuleBase::Vector3<double>>& kvec_d_in)
        : LM(LM_in), kvec_d(kvec_d_in){};
    virtual ~OperatorLCAO()
    {
        if (this->allocated_smatrix)
            delete[] this->smatrix_k;
    };

    /* Function init(k) is used for update HR and HK ,
    data pointers of HR and HK are not passed by this function, but passed by constructors of every derived classes.
    No need to override init() in base class, but must override in derived class */
    virtual void init(const int ik_in) override;

    void refresh_h();

    /* Function getHR() is designed to update HR matrix only, it will loop all contributeHR() functions in chain table.
    Detail of this function is still in developed, HR matrix has two form: HR_all_spin and HR_one_spin.
    For NSPIN=2 case, HR_one_spin for spin-up and spin-down is not constructed at same time.
    */
    // void getHR(T* hr_pointer);

    /* Function contributeHR() is defined in derived class, for constructing <phi_{\mu, R}|H|phi_{\nu, 0}>
     */
    virtual void contributeHR()
    {
        return;
    }

    /* Function matrixHk() is used for get information of HK matrix and SK matrix for diagolization.
    Matrixes HK and SK come from LCAO_Matrix class.
    Gamma_only case (T = double), SK would not changed during one SCF loop, a template triangle matrix SK_temp is used
    for accelerating. General case (T = std::complex<double>), only pointers of HK and SK saved in OperatorLCAO
    */
    void matrixHk(MatrixBlock<T>& hk_in, MatrixBlock<T>& sk_in)
    {
        this->get_hs_pointers();
#ifdef __MPI
        hk_in = MatrixBlock<T>{hmatrix_k,
                               (size_t)this->LM->ParaV->nrow,
                               (size_t)this->LM->ParaV->ncol,
                               this->LM->ParaV->desc};
        sk_in = MatrixBlock<T>{smatrix_k,
                               (size_t)this->LM->ParaV->nrow,
                               (size_t)this->LM->ParaV->ncol,
                               this->LM->ParaV->desc};
#else
        hk_in = MatrixBlock<T>{hmatrix_k, (size_t)this->LM->ParaV->nrow, (size_t)this->LM->ParaV->ncol, nullptr};
        sk_in = MatrixBlock<T>{smatrix_k, (size_t)this->LM->ParaV->nrow, (size_t)this->LM->ParaV->ncol, nullptr};
#endif
    }

    /* Function contributeHk() is defined in derived class, for constructing <phi_{\mu}|H|phi_{\nu}>(K)
     */
    virtual void contributeHk(int ik)
    {
        return;
    }

    // protected:
    //  Hamiltonian matrix which are stored in LCAO_Matrix and calculated in OperatorLCAO
    LCAO_Matrix* LM = nullptr;
    const std::vector<ModuleBase::Vector3<double>>& kvec_d;

  protected:
    bool new_e_iteration = true;

  private:
    void get_hs_pointers();

    // there are H and S matrix for each k point in reciprocal space
    // type double for gamma_only case, type complex<double> for multi-k-points case
    T* hmatrix_k = nullptr;
    T* smatrix_k = nullptr;

    // only used for Gamma_only case
    bool allocated_smatrix = false;

    // fixed HR matrix folding to HK
    void folding_fixed(const int ik, const std::vector<ModuleBase::Vector3<double>>& kvec_d);

    // const std::vector<ModuleBase::Vector3<double>>& kvec_d;
};

} // end namespace hamilt

#endif