#ifndef Diago_David_CUH
#define Diago_David_CUH

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/complexmatrix.h"

class Diago_David_CUDA
{
public:

    Diago_David_CUDA();
    ~Diago_David_CUDA();
    
    // SchmitOrth in Davison Method.
    static void SchmitOrth(
        const int& npw,
        const int n_band,
        const int m,
        const ModuleBase::ComplexMatrix &psi,
        double2* psi_m,
        double2* spsi);
    
    // Use LapackConnector::zhegvx operations.
    static void diag_zhegvx(
        const int& n,
        const int& m,
        const ModuleBase::ComplexMatrix &hc,
        const ModuleBase::ComplexMatrix &sc,
        const int& ldh,
        double* e,
        ModuleBase::ComplexMatrix &vc);

    // Solving for eigenvalues using davison's method 
    void diag(
        ModuleBase::ComplexMatrix &psi,
        double *en,
        const int &npw,
        const int &nband,
        const double *precondition,
        const int order,
        const double &eps,
        const int &maxiter,
        int &notconv,
        double &avg_iter);

private:

    int test_david;

    // calculate gradient
    void cal_grad(
        const int& npw,
        const int& nbase,
        const int& notconv,
        ModuleBase::ComplexMatrix &basis,
        ModuleBase::ComplexMatrix &hp,
        ModuleBase::ComplexMatrix &sp,
        const ModuleBase::ComplexMatrix &vc,
        const int* unconv,
        const double* precondition,
        const double* en,
        double2* hpsi,
        double2* spsi,
        double2* ppsi,
        double2* respsi);

    void cal_elem(
        const int& npw,
        int& nbase,
        int& nbase_x,
        const int& notconv,
        const ModuleBase::ComplexMatrix &basis,
        const ModuleBase::ComplexMatrix &hp,
        const ModuleBase::ComplexMatrix &sp,
        ModuleBase::ComplexMatrix &hc,
        ModuleBase::ComplexMatrix &sc);

    // Update hp, sp, basis and reduce Hamiltonian.
    void refresh(
        const int& npw,
        const int& nband,
        int& nbase,
        const double* en,
        const ModuleBase::ComplexMatrix &psi,
        ModuleBase::ComplexMatrix &basis,
        ModuleBase::ComplexMatrix &hp,
        ModuleBase::ComplexMatrix &sp,
        ModuleBase::ComplexMatrix &hc,
        ModuleBase::ComplexMatrix &sc,
        ModuleBase::ComplexMatrix &vc);

    // Calculate residuals.
    void cal_err(
        const int& npw,
        const int& nband,
        const int& nbase,
        const ModuleBase::ComplexMatrix &vc,
        const ModuleBase::ComplexMatrix &hp,
        const ModuleBase::ComplexMatrix &basis,
        const double* en,
        double2* respsi);

};

#endif
