#include "diago_elpa.h"

#include "module_base/global_variable.h"
#include "module_base/lapack_connector.h"
#include "module_base/timer.h"
#include "module_base/tool_quit.h"
extern "C"
{
#include "module_base/blacs_connector.h"
#include "module_base/scalapack_connector.h"
}
#include "genelpa/elpa_solver.h"

typedef hamilt::MatrixBlock<double> matd;
typedef hamilt::MatrixBlock<std::complex<double>> matcd;

namespace hsolver
{
    template<>
    int DiagoElpa<double>::DecomposedState = 0;
    template<>
    int DiagoElpa<std::complex<double>>::DecomposedState = 0;
    template<>
    void DiagoElpa<std::complex<double>>::diag(hamilt::Hamilt<std::complex<double>>* phm_in, psi::Psi<std::complex<double>>& psi, Real* eigenvalue_in)
{
    ModuleBase::TITLE("DiagoElpa", "diag");
#ifdef __MPI
    matcd h_mat, s_mat;
    phm_in->matrix(h_mat, s_mat);

    std::vector<double> eigen(GlobalV::NLOCAL, 0.0);

    bool isReal=false;
    const MPI_Comm COMM_DIAG=MPI_COMM_WORLD; // use all processes
    ELPA_Solver es((const bool)isReal, COMM_DIAG, (const int)GlobalV::NBANDS, (const int)h_mat.row, (const int)h_mat.col, (const int*)h_mat.desc);
    this->DecomposedState=0; // for k pointer, the decomposed s_mat can not be reused
    ModuleBase::timer::tick("DiagoElpa", "elpa_solve");
    es.generalized_eigenvector(h_mat.p, s_mat.p, this->DecomposedState, eigen.data(), psi.get_pointer());
    ModuleBase::timer::tick("DiagoElpa", "elpa_solve");
    es.exit();

    const int inc = 1;
    BlasConnector::copy(GlobalV::NBANDS, eigen.data(), inc, eigenvalue_in, inc);
#else
    ModuleBase::WARNING_QUIT("DiagoElpa", "DiagoElpa only can be used with macro __MPI");
#endif
}

    template<>
    void DiagoElpa<double>::diag(hamilt::Hamilt<double>* phm_in, psi::Psi<double>& psi, Real* eigenvalue_in)
{
    ModuleBase::TITLE("DiagoElpa", "diag");
#ifdef __MPI
    matd h_mat, s_mat;
    phm_in->matrix(h_mat, s_mat);

    std::vector<double> eigen(GlobalV::NLOCAL, 0.0);

    bool isReal=true;
    MPI_Comm COMM_DIAG=MPI_COMM_WORLD; // use all processes
    //ELPA_Solver es(isReal, COMM_DIAG, GlobalV::NBANDS, h_mat.row, h_mat.col, h_mat.desc);
    ELPA_Solver es((const bool)isReal, COMM_DIAG, (const int)GlobalV::NBANDS, (const int)h_mat.row, (const int)h_mat.col, (const int*)h_mat.desc);
    ModuleBase::timer::tick("DiagoElpa", "elpa_solve");
    es.generalized_eigenvector(h_mat.p, s_mat.p, this->DecomposedState, eigen.data(), psi.get_pointer());
    ModuleBase::timer::tick("DiagoElpa", "elpa_solve");
    es.exit();

    const int inc = 1;
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "K-S equation was solved by genelpa2");
    BlasConnector::copy(GlobalV::NBANDS, eigen.data(), inc, eigenvalue_in, inc);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "eigenvalues were copied to ekb");
#else
    ModuleBase::WARNING_QUIT("DiagoElpa", "DiagoElpa only can be used with macro __MPI");
#endif
}

#ifdef __MPI
    template<typename T>
    bool DiagoElpa<T>::ifElpaHandle(const bool& newIteration, const bool& ifNSCF)
{
    int doHandle = false;
    if (newIteration)
        doHandle = true;
    if (ifNSCF)
        doHandle = true;
    return doHandle;
}
#endif

} // namespace hsolver
