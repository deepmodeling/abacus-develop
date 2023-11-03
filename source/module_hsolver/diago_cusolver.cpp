#include "diago_cusolver.h"

#include "module_base/global_variable.h"
#include "module_base/lapack_connector.h"
#include "module_base/timer.h"
#include "module_base/tool_quit.h"
extern "C"
{
#include "module_base/blacs_connector.h"
#include "module_base/scalapack_connector.h"
}
#include "module_hsolver/kernels/cuda/diag_cusolver.cuh"
typedef hamilt::MatrixBlock<double> matd;
typedef hamilt::MatrixBlock<std::complex<double>> matcd;
using namespace std;
namespace hsolver
{
    template<>
    int DiagoCusolver<double>::DecomposedState = 0;
    template<>
    int DiagoCusolver<std::complex<double>>::DecomposedState = 0;
    template<>
    void DiagoCusolver<std::complex<double>>::diag(hamilt::Hamilt<std::complex<double>>* phm_in, psi::Psi<std::complex<double>>& psi, Real* eigenvalue_in)
    {
        ModuleBase::TITLE("DiagoCusolver", "diag");
        matcd h_mat,s_mat;
        phm_in->matrix(h_mat, s_mat);

        std::vector<double> eigen(GlobalV::NLOCAL, 0.0);
        bool isReal=false;
        Diag_Cusolver_gvd dc;
        this->DecomposedState=0; 
        ModuleBase::timer::tick("DiagoCusolver", "cusolver");
        //es.generalized_eigenvector(h_mat.p, s_mat.p, this->DecomposedState, eigen.data(), psi.get_pointer());
        //   subroutines that are related to calculating generalized eigenvalues and eigenvectors for dense matrix pairs:
//  - Dngvd_double : dense double type matrix
//  - Dngvd_complex : dense complex type matrix
//      Input Parameters
//          N: the number of rows of the matrix 
//          M: the number of cols of the matrix  
//          A: the hermitian matrix A in A x=lambda B (column major) 
//          B: the SPD matrix B in A x=lambda B (column major) 
//      Output Parameter
//          W: generalized eigenvalues
//          V: generalized eigenvectors (column major)
        dc.Dngvd_complex(h_mat.col, h_mat.row, h_mat.p, s_mat.p, eigen.data(), psi.get_pointer());
        ModuleBase::timer::tick("DiagoElpa", "elpa_solve");
        const int inc = 1;
        BlasConnector::copy(GlobalV::NBANDS, eigen.data(), inc, eigenvalue_in, inc);

    }

    void DiagoCusolver<double>::diag(hamilt::Hamilt<double>* phm_in, psi::Psi<double>& psi, Real* eigenvalue_in)
    {
        ModuleBase::TITLE("DiagoCusolver", "diag");
        printf("which have diagosolver\n");
        matd h_mat,s_mat;
        phm_in->matrix(h_mat, s_mat);

        std::vector<double> eigen(GlobalV::NLOCAL, 0.0);
        bool isReal=false;
        Diag_Cusolver_gvd dc;
        this->DecomposedState=0; 
        ModuleBase::timer::tick("DiagoCusolver", "cusolver");
        //es.generalized_eigenvector(h_mat.p, s_mat.p, this->DecomposedState, eigen.data(), psi.get_pointer());
        //   subroutines that are related to calculating generalized eigenvalues and eigenvectors for dense matrix pairs:
        //  - Dngvd_double : dense double type matrix
        //  - Dngvd_complex : dense complex type matrix
        //      Input Parameters
        //          N: the number of rows of the matrix 
        //          M: the number of cols of the matrix  
        //          A: the hermitian matrix A in A x=lambda B (column major) 
        //          B: the SPD matrix B in A x=lambda B (column major) 
        //      Output Parameter
        //          W: generalized eigenvalues
        //          V: generalized eigenvectors (column major)
        dc.Dngvd_double(h_mat.col, h_mat.row, h_mat.p, s_mat.p, eigen.data(), psi.get_pointer());
        ModuleBase::timer::tick("DiagoElpa", "elpa_solve");
        const int inc = 1;
        BlasConnector::copy(GlobalV::NBANDS, eigen.data(), inc, eigenvalue_in, inc);
    }
}