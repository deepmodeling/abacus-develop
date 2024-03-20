#ifndef OPERATORLCAO_H
#define OPERATORLCAO_H
#include "module_base/vector3.h"
#include "module_hamilt_general/matrixblock.h"
#include "module_hamilt_general/operator.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_matrix.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"

namespace hamilt
{

template <typename TK, typename TR>
class OperatorLCAO : public Operator<TK>
{
  public:
    OperatorLCAO(
        LCAO_Matrix* LM_in, 
        const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
        HContainer<TR>* hR_in,
        std::vector<TK>* hK_in)
        : LM(LM_in), kvec_d(kvec_d_in), hR(hR_in), hK(hK_in){};
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
    // void getHR(TR* hr_pointer);

    /* Function contributeHR() is defined in derived class, for constructing <phi_{\mu, R}|H|phi_{\nu, 0}>
     */
    virtual void contributeHR()
    {
        return;
    }

    /* Function matrixHk() is used for get information of HK matrix and SK matrix for diagolization.
    Matrixes HK and SK come from LCAO_Matrix class.
    Gamma_only case (TK = double), SK would not changed during one SCF loop, a template triangle matrix SK_temp is used
    for accelerating. General case (TK = std::complex<double>), only pointers of HK and SK saved in OperatorLCAO
    */
    void matrixHk(MatrixBlock<TK>& hk_in, MatrixBlock<TK>& sk_in)
    {
        this->get_hs_pointers();
#ifdef __MPI
        hk_in = MatrixBlock<TK>{hmatrix_k,
                               (size_t)this->LM->ParaV->nrow,
                               (size_t)this->LM->ParaV->ncol,
                               this->LM->ParaV->desc};
        sk_in = MatrixBlock<TK>{smatrix_k,
                               (size_t)this->LM->ParaV->nrow,
                               (size_t)this->LM->ParaV->ncol,
                               this->LM->ParaV->desc};
#else
        hk_in = MatrixBlock<TK>{hmatrix_k, (size_t)this->LM->ParaV->nrow, (size_t)this->LM->ParaV->ncol, nullptr};
        sk_in = MatrixBlock<TK>{smatrix_k, (size_t)this->LM->ParaV->nrow, (size_t)this->LM->ParaV->ncol, nullptr};
#endif
    }

    /* Function contributeHk() is defined in derived class, for constructing <phi_{\mu}|H|phi_{\nu}>(K)
     */
    virtual void contributeHk(int ik);

    /**
     * @brief set_HR_fixed() is used for pass HR_fixed matrix to the next node in sub-chain table
     * not used in base class, only be override in fixed Hamiltonian Operators (e.g. Ekinetic and Nonlocal)
    */
    virtual void set_HR_fixed(void*)
    {
        return;
    }

    /**
     * @brief reset hr_done status
    */
    void set_hr_done(bool hr_done_in);

    // protected:
    //  Hamiltonian matrix which are stored in LCAO_Matrix and calculated in OperatorLCAO
    LCAO_Matrix* LM = nullptr;
    const std::vector<ModuleBase::Vector3<double>>& kvec_d;

  protected:
    bool new_e_iteration = true;

    // Real space Hamiltonian pointer
    hamilt::HContainer<TR>* hR = nullptr;

    // vector of HK matrix for current k point in reciprocal space
    std::vector<TK>* hK = nullptr;

  private:
    void get_hs_pointers();

    // there are H and S matrix for each k point in reciprocal space
    // type double for gamma_only case, type complex<double> for multi-k-points case
    TK* hmatrix_k = nullptr;
    TK* smatrix_k = nullptr;

    // only used for Gamma_only case
    bool allocated_smatrix = false;
    
    // if HR is calculated
    bool hr_done = false;
};

} // end namespace hamilt

#endif