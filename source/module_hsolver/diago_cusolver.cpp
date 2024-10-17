#include "diago_cusolver.h"

#include "module_parameter/parameter.h"
#include "module_base/blacs_connector.h"
#include "module_base/global_variable.h"
#include "module_base/scalapack_connector.h"
#include "module_base/timer.h"

#include <memory>
#include <type_traits>
#include <vector>

using complex = std::complex<double>;

// Namespace for the diagonalization solver
namespace hsolver
{
// this struct is used for collecting matrices from all processes to root process
template <typename T>
struct Matrix_g
{
    std::shared_ptr<T> p;
    size_t row;
    size_t col;
    std::shared_ptr<int> desc;
};

// Initialize the DecomposedState variable for real and complex numbers
template <typename T>
int DiagoCusolver<T>::DecomposedState = 0;

template <typename T>
DiagoCusolver<T>::DiagoCusolver(const Parallel_Orbitals* ParaV)
{
    this->ParaV = ParaV;
}

template <typename T>
DiagoCusolver<T>::~DiagoCusolver()
{
}

// Wrapper for pdgemr2d and pzgemr2d
// static inline void Cpxgemr2d(
//    const int M, const int N,
//    double *a, const int ia, const int ja, const int *desca,
//    double *b, const int ib, const int jb, const int *descb,
//    const int blacs_ctxt)
//{
//    pdgemr2d_(&M, &N,
//              a, &ia, &ja, desca,
//              b, &ib, &jb, descb,
//              &blacs_ctxt);
//}
//
// static inline void Cpxgemr2d(
//    const int M, const int N,
//    complex *a, const int ia, const int ja, const int *desca,
//    complex *b, const int ib, const int jb, const int *descb,
//    const int blacs_ctxt)
//{
//    pzgemr2d_(&M, &N,
//              a, &ia, &ja, desca,
//              b, &ib, &jb, descb,
//              &blacs_ctxt);
//}

// Use Cpxgemr2d to collect matrices from all processes to root process
template <typename mat, typename matg>
static void gatherMatrix(const int myid, const int root_proc, const mat& mat_l, matg& mat_g)
{
    auto a = mat_l.p;
    const int* desca = mat_l.desc;
    int ctxt = desca[1];
    int nrows = desca[2];
    int ncols = desca[3];

    if (myid == root_proc)
    {
        mat_g.p.reset(new typename std::remove_reference<decltype(*a)>::type[nrows * ncols]);
    }
    else
    {
        mat_g.p.reset(new typename std::remove_reference<decltype(*a)>::type[1]);
    }

    // Set descb, which has all elements in the only block in the root process
    mat_g.desc.reset(new int[9]{1, ctxt, nrows, ncols, nrows, ncols, 0, 0, nrows});

    mat_g.row = nrows;
    mat_g.col = ncols;

    Cpxgemr2d(nrows, ncols, a, 1, 1, const_cast<int*>(desca), mat_g.p.get(), 1, 1, mat_g.desc.get(), ctxt);
}

// Convert the Psi to a 2D block storage format
template <typename T>
static void distributePsi(const int* desc_psi, T* psi, T* psi_g)
{
    int ctxt = desc_psi[1];
    int nrows = desc_psi[2];
    int ncols = desc_psi[3];
    int rsrc = desc_psi[6];
    int csrc = desc_psi[7];

    int descg[9] = {1, ctxt, nrows, ncols, nrows, ncols, rsrc, csrc, nrows};
    int descl[9];

    std::copy(desc_psi, desc_psi + 9, descl);

    Cpxgemr2d(nrows, ncols, psi_g, 1, 1, descg, psi, 1, 1, descl, ctxt);
}

// Diagonalization function
template <typename T>
void DiagoCusolver<T>::diag(hamilt::Hamilt<T>* phm_in, psi::Psi<T>& psi, Real* eigenvalue_in)
{
    // Output the title for the current operation
    ModuleBase::TITLE("DiagoCusolver", "diag");

    // Create matrices for the Hamiltonian and overlap
    hamilt::MatrixBlock<T> h_mat;
    hamilt::MatrixBlock<T> s_mat;
    phm_in->matrix(h_mat, s_mat);

#ifdef __MPI
    // global matrix
    Matrix_g<T> h_mat_g;
    Matrix_g<T> s_mat_g;

    // get the context and process information
    int ctxt = ParaV->blacs_ctxt;
    int nprows = 0;
    int npcols = 0;
    int myprow = 0;
    int mypcol = 0;
    Cblacs_gridinfo(ctxt, &nprows, &npcols, &myprow, &mypcol);
    const int num_procs = nprows * npcols;
    const int myid = Cblacs_pnum(ctxt, myprow, mypcol);
    const int root_proc = Cblacs_pnum(ctxt, ParaV->desc[6], ParaV->desc[7]);
#endif

    // Allocate memory for eigenvalues
    std::vector<double> eigen(PARAM.globalv.nlocal, 0.0);

    // Start the timer for the cusolver operation
    ModuleBase::timer::tick("DiagoCusolver", "cusolver");

#ifdef __MPI
    if (num_procs > 1)
    {
        // gather matrices from processes to root process
        gatherMatrix(myid, root_proc, h_mat, h_mat_g);
        gatherMatrix(myid, root_proc, s_mat, s_mat_g);
    }
#endif

    // Call the dense diagonalization routine
#ifdef __MPI
    if (num_procs > 1)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        // global psi for distribute
        int psi_len = myid == root_proc ? h_mat_g.row * h_mat_g.col : 1;
        std::vector<T> psi_g(psi_len);
        if (myid == root_proc)
        {
            this->dc.Dngvd(h_mat_g.col, h_mat_g.row, h_mat_g.p.get(), s_mat_g.p.get(), eigen.data(), psi_g.data());
        }

        MPI_Barrier(MPI_COMM_WORLD);

        // broadcast eigenvalues to all processes
        MPI_Bcast(eigen.data(), PARAM.inp.nbands, MPI_DOUBLE, root_proc, MPI_COMM_WORLD);

        // distribute psi to all processes
        distributePsi(this->ParaV->desc_wfc, psi.get_pointer(), psi_g.data());
    }
    else
    {
        // Be careful that h_mat.row * h_mat.col != psi.get_nbands() * psi.get_nbasis() under multi-k situation
        std::vector<T> eigenvectors(h_mat.row * h_mat.col);
        this->dc.Dngvd(h_mat.row, h_mat.col, h_mat.p, s_mat.p, eigen.data(), eigenvectors.data());
        const int size = psi.get_nbands() * psi.get_nbasis();
        BlasConnector::copy(size, eigenvectors.data(), 1, psi.get_pointer(), 1);
    }
#else
    std::vector<T> eigenvectors(h_mat.row * h_mat.col);
    this->dc.Dngvd(h_mat.row, h_mat.col, h_mat.p, s_mat.p, eigen.data(), eigenvectors.data());
    const int size = psi.get_nbands() * psi.get_nbasis();
    BlasConnector::copy(size, eigenvectors.data(), 1, psi.get_pointer(), 1);
#endif
    // Stop the timer for the cusolver operation
    ModuleBase::timer::tick("DiagoCusolver", "cusolver");

    // Copy the eigenvalues to the output arrays
    const int inc = 1;
    BlasConnector::copy(PARAM.inp.nbands, eigen.data(), inc, eigenvalue_in, inc);
}

// Explicit instantiation of the DiagoCusolver class for real and complex numbers
template class DiagoCusolver<double>;
template class DiagoCusolver<complex>;

} // namespace hsolver
