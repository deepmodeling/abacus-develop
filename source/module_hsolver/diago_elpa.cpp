#include "diago_elpa.h"

#include "module_parameter/parameter.h"
#include "genelpa/elpa_solver.h"
#include "module_base/blacs_connector.h"
#include "module_base/global_variable.h"
#include "module_base/lapack_connector.h"
#include "module_base/scalapack_connector.h"
#include "module_base/timer.h"
#include "module_base/tool_quit.h"

typedef hamilt::MatrixBlock<double> matd;
typedef hamilt::MatrixBlock<std::complex<double>> matcd;

namespace hsolver {
#ifdef __MPI
template <>
MPI_Comm DiagoElpa<double>::setmpicomm() {
    if (this->elpa_num_thread == -1) {
        return MPI_COMM_WORLD;
    } else {
        int _num = 0;
        MPI_Comm_size(MPI_COMM_WORLD, &_num);
        if (elpa_num_thread > _num || elpa_num_thread <= 0) {
            return MPI_COMM_WORLD;
        } else {
            lastmpinum++;
            int* _ranks = new int[elpa_num_thread];
            for (int i = 0; i < elpa_num_thread; i++) {
                _ranks[i] = (lastmpinum + i) % _num;
            }
            MPI_Group _tempgroup, _oldgroup;
            MPI_Comm_group(MPI_COMM_WORLD, &_oldgroup);
            MPI_Group_incl(_oldgroup, elpa_num_thread, _ranks, &_tempgroup);
            MPI_Comm _new_comm;
            MPI_Comm_create(MPI_COMM_WORLD, _tempgroup, &_new_comm);
            delete[] _ranks;
            return _new_comm;
        }
    }
}
template <>
MPI_Comm DiagoElpa<std::complex<double>>::setmpicomm() {
    if (this->elpa_num_thread == -1) {
        return MPI_COMM_WORLD;
    } else {
        int _num = 0;
        MPI_Comm_size(MPI_COMM_WORLD, &_num);
        if (elpa_num_thread > _num || elpa_num_thread <= 0) {
            return MPI_COMM_WORLD;
        } else {
            lastmpinum++;
            int* _ranks = new int[elpa_num_thread];
            for (int i = 0; i < elpa_num_thread; i++) {
                _ranks[i] = (lastmpinum + i) % _num;
            }
            MPI_Group _tempgroup, _oldgroup;
            MPI_Comm_group(MPI_COMM_WORLD, &_oldgroup);
            MPI_Group_incl(_oldgroup, elpa_num_thread, _ranks, &_tempgroup);
            MPI_Comm _new_comm;
            MPI_Comm_create(MPI_COMM_WORLD, _tempgroup, &_new_comm);
            delete[] _ranks;
            return _new_comm;
        }
    }
}
#endif
template <>
void DiagoElpa<std::complex<double>>::diag(
    hamilt::Hamilt<std::complex<double>>* phm_in,
    psi::Psi<std::complex<double>>& psi,
    Real* eigenvalue_in) {
    ModuleBase::TITLE("DiagoElpa", "diag");
#ifdef __MPI
    matcd h_mat, s_mat;
    phm_in->matrix(h_mat, s_mat);

    std::vector<double> eigen(PARAM.globalv.nlocal, 0.0);

    bool isReal = false;
    MPI_Comm COMM_DIAG = setmpicomm(); // set mpi_comm needed
    ELPA_Solver es((const bool)isReal,
                   COMM_DIAG,
                   (const int)PARAM.inp.nbands,
                   (const int)h_mat.row,
                   (const int)h_mat.col,
                   (const int*)h_mat.desc);
    this->DecomposedState
        = 0; // for k pointer, the decomposed s_mat can not be reused
    ModuleBase::timer::tick("DiagoElpa", "elpa_solve");
    es.generalized_eigenvector(h_mat.p,
                               s_mat.p,
                               this->DecomposedState,
                               eigen.data(),
                               psi.get_pointer());
    ModuleBase::timer::tick("DiagoElpa", "elpa_solve");
    es.exit();

    const int inc = 1;
    BlasConnector::copy(PARAM.inp.nbands, eigen.data(), inc, eigenvalue_in, inc);
#else
    ModuleBase::WARNING_QUIT("DiagoElpa",
                             "DiagoElpa only can be used with macro __MPI");
#endif
}

template <>
void DiagoElpa<double>::diag(hamilt::Hamilt<double>* phm_in,
                             psi::Psi<double>& psi,
                             Real* eigenvalue_in) {
    ModuleBase::TITLE("DiagoElpa", "diag");
#ifdef __MPI
    matd h_mat, s_mat;
    phm_in->matrix(h_mat, s_mat);

    std::vector<double> eigen(PARAM.globalv.nlocal, 0.0);

    bool isReal = true;
    MPI_Comm COMM_DIAG = setmpicomm(); // set mpi_comm needed
    // ELPA_Solver es(isReal, COMM_DIAG, PARAM.inp.nbands, h_mat.row, h_mat.col,
    // h_mat.desc);
    ELPA_Solver es((const bool)isReal,
                   COMM_DIAG,
                   (const int)PARAM.inp.nbands,
                   (const int)h_mat.row,
                   (const int)h_mat.col,
                   (const int*)h_mat.desc);
    ModuleBase::timer::tick("DiagoElpa", "elpa_solve");
    es.generalized_eigenvector(h_mat.p,
                               s_mat.p,
                               this->DecomposedState,
                               eigen.data(),
                               psi.get_pointer());
    ModuleBase::timer::tick("DiagoElpa", "elpa_solve");
    es.exit();

    const int inc = 1;
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,
                                "K-S equation was solved by genelpa2");
    BlasConnector::copy(PARAM.inp.nbands, eigen.data(), inc, eigenvalue_in, inc);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,
                                "eigenvalues were copied to ekb");
#else
    ModuleBase::WARNING_QUIT("DiagoElpa",
                             "DiagoElpa only can be used with macro __MPI");
#endif
}


#ifdef __MPI
template <>
void DiagoElpa<std::complex<double>>::diag_pool(hamilt::MatrixBlock<std::complex<double>>& h_mat,
    hamilt::MatrixBlock<std::complex<double>>& s_mat,
    psi::Psi<std::complex<double>>& psi,
    Real* eigenvalue_in,
    MPI_Comm& comm)
{
    std::vector<double> eigen(PARAM.globalv.nlocal, 0.0);
    bool isReal = false;
    ELPA_Solver es((const bool)isReal,
                   comm,
                   (const int)PARAM.inp.nbands,
                   (const int)h_mat.row,
                   (const int)h_mat.col,
                   (const int*)h_mat.desc);
    this->DecomposedState
        = 0; // for k pointer, the decomposed s_mat can not be reused
    ModuleBase::timer::tick("DiagoElpa", "elpa_solve");
    es.generalized_eigenvector(h_mat.p,
                               s_mat.p,
                               this->DecomposedState,
                               eigen.data(),
                               psi.get_pointer());
    ModuleBase::timer::tick("DiagoElpa", "elpa_solve");
    es.exit();
    const int inc = 1;
    BlasConnector::copy(PARAM.inp.nbands, eigen.data(), inc, eigenvalue_in, inc);
}

template <>
void DiagoElpa<double>::diag_pool(hamilt::MatrixBlock<double>& h_mat,
    hamilt::MatrixBlock<double>& s_mat,
    psi::Psi<double>& psi,
    Real* eigenvalue_in,
    MPI_Comm& comm)
{
    std::vector<double> eigen(PARAM.globalv.nlocal, 0.0);

    bool isReal = true;
    // ELPA_Solver es(isReal, COMM_DIAG, PARAM.inp.nbands, h_mat.row, h_mat.col,
    // h_mat.desc);
    ELPA_Solver es((const bool)isReal,
                   comm,
                   (const int)PARAM.inp.nbands,
                   (const int)h_mat.row,
                   (const int)h_mat.col,
                   (const int*)h_mat.desc);
    ModuleBase::timer::tick("DiagoElpa", "elpa_solve");
    es.generalized_eigenvector(h_mat.p,
                               s_mat.p,
                               this->DecomposedState,
                               eigen.data(),
                               psi.get_pointer());
    ModuleBase::timer::tick("DiagoElpa", "elpa_solve");
    es.exit();

    const int inc = 1;
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,
                                "K-S equation was solved by genelpa2");
    BlasConnector::copy(PARAM.inp.nbands, eigen.data(), inc, eigenvalue_in, inc);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,
                                "eigenvalues were copied to ekb");
}
#endif



#ifdef __MPI
template <typename T>
bool DiagoElpa<T>::ifElpaHandle(const bool& newIteration, const bool& ifNSCF) {
    int doHandle = false;
    if (newIteration) {
        doHandle = true;
    }
    if (ifNSCF) {
        doHandle = true;
    }
    return doHandle;
}
#endif

} // namespace hsolver
