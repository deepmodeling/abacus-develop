#include "diago_elpa.h"

#include "module_base/global_variable.h"
#include "module_base/lapack_connector.h"
#include "module_base/timer.h"
extern "C"
{
#include "module_base/blacs_connector.h"
#include "my_elpa.h"
#include "module_base/scalapack_connector.h"
}

typedef hamilt::MatrixBlock<double> matd;
typedef hamilt::MatrixBlock<std::complex<double>> matcd;

namespace hsolver
{
bool DiagoElpa::is_already_decomposed = false;
#ifdef __MPI
inline int set_elpahandle(elpa_t &handle,
                          const int *desc,
                          const int local_nrows,
                          const int local_ncols,
                          const int nbands)
{
    int error;
    int nprows, npcols, myprow, mypcol;
    Cblacs_gridinfo(desc[1], &nprows, &npcols, &myprow, &mypcol);
    elpa_init(20210430);
    handle = elpa_allocate(&error);
    elpa_set_integer(handle, "na", desc[2], &error);
    elpa_set_integer(handle, "nev", nbands, &error);

    elpa_set_integer(handle, "local_nrows", local_nrows, &error);

    elpa_set_integer(handle, "local_ncols", local_ncols, &error);

    elpa_set_integer(handle, "nblk", desc[4], &error);

    elpa_set_integer(handle, "mpi_comm_parent", MPI_Comm_c2f(MPI_COMM_WORLD), &error);

    elpa_set_integer(handle, "process_row", myprow, &error);

    elpa_set_integer(handle, "process_col", mypcol, &error);

    elpa_set_integer(handle, "blacs_context", desc[1], &error);

    elpa_set_integer(handle, "cannon_for_generalized", 0, &error);
    /* Setup */
    elpa_setup(handle); /* Set tunables */
    return 0;
}
#endif

void DiagoElpa::diag(hamilt::Hamilt *phm_in, psi::Psi<std::complex<double>> &psi, double *eigenvalue_in)
{
    matcd h_mat, s_mat;
    phm_in->matrix(h_mat, s_mat);

    std::vector<double> eigen(GlobalV::NLOCAL, 0.0);

    static elpa_t handle;
    static bool has_set_elpa_handle = false;
    if (!has_set_elpa_handle)
    {
        set_elpahandle(handle, h_mat.desc, h_mat.row, h_mat.col, GlobalV::NBANDS);
        has_set_elpa_handle = true;
    }

    // compare to old code from pplab, there is no need to copy Sloc2 to another memory,
    // just change Sloc2, which is a temporary matrix
    // size_t nloc = h_mat.col * h_mat.row,
    // BlasConnector::copy(nloc, s_mat, inc, Stmp, inc);

    ModuleBase::timer::tick("DiagoElpa", "elpa_solve");
    int elpa_derror;
    elpa_generalized_eigenvectors_dc(handle,
                                     reinterpret_cast<double _Complex *>(h_mat.p),
                                     reinterpret_cast<double _Complex *>(s_mat.p),
                                     eigen.data(),
                                     reinterpret_cast<double _Complex *>(psi.get_pointer()),
                                     0,
                                     &elpa_derror);
    ModuleBase::timer::tick("DiagoElpa", "elpa_solve");

    // the eigenvalues.
    const int inc = 1;
    BlasConnector::copy(GlobalV::NBANDS, eigen.data(), inc, eigenvalue_in, inc);
}

void DiagoElpa::diag(hamilt::Hamilt *phm_in, psi::Psi<double> &psi, double *eigenvalue_in)
{
    matd h_mat, s_mat;
    phm_in->matrix(h_mat, s_mat);

    std::vector<double> eigen(GlobalV::NLOCAL, 0.0);

    static elpa_t handle;
    static bool has_set_elpa_handle = false;
    if (!has_set_elpa_handle)
    {
        set_elpahandle(handle, h_mat.desc, h_mat.row, h_mat.col, GlobalV::NBANDS);
        has_set_elpa_handle = true;
    }

    // compare to old code from pplab, there is no need to copy Sloc2 to another memory,
    // just change Sloc2, which is a temporary matrix
    // change this judgement to HamiltLCAO
    /*int is_already_decomposed;
    if(ifElpaHandle(GlobalC::CHR.get_new_e_iteration(), (GlobalV::CALCULATION=="nscf")))
    {
        ModuleBase::timer::tick("DiagoElpa","decompose_S");
        BlasConnector::copy(pv->nloc, s_mat, inc, Stmp, inc);
        is_already_decomposed=0;
        ModuleBase::timer::tick("DiagoElpa","decompose_S");
    }
    else
    {
        is_already_decomposed=1;
    }*/

    ModuleBase::timer::tick("DiagoElpa", "elpa_solve");
    int elpa_error;
    elpa_generalized_eigenvectors_d(handle,
                                    h_mat.p,
                                    s_mat.p,
                                    eigen.data(),
                                    psi.get_pointer(),
                                    DiagoElpa::is_already_decomposed,
                                    &elpa_error);
    ModuleBase::timer::tick("DiagoElpa", "elpa_solve");

    const int inc = 1;
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "K-S equation was solved by genelpa2");
    BlasConnector::copy(GlobalV::NBANDS, eigen.data(), inc, eigenvalue_in, inc);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "eigenvalues were copied to ekb");
}

bool DiagoElpa::ifElpaHandle(const bool &newIteration, const bool &ifNSCF)
{
    int doHandle = false;
    if (newIteration)
        doHandle = true;
    if (ifNSCF)
        doHandle = true;
    return doHandle;
}

} // namespace hsolver