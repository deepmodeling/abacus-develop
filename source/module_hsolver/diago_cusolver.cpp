#include "diago_cusolver.h"
#include "module_base/global_variable.h"
#include "module_base/timer.h"
#include "module_base/blacs_connector.h"
#include "module_base/scalapack_connector.h"

#include <type_traits>

using complex = std::complex<double>;

// Namespace for the diagonalization solver
namespace hsolver
{
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
    static inline void Cpxgemr2d(
        const int M, const int N,
        double *a, const int ia, const int ja, const int *desca,
        double *b, const int ib, const int jb, const int *descb,
        const int blacs_ctxt) 
    {
        pdgemr2d_(&M, &N,
                  a, &ia, &ja, desca,
                  b, &ib, &jb, descb,
                  &blacs_ctxt);
    }

    static inline void Cpxgemr2d(
        const int M, const int N,
        complex *a, const int ia, const int ja, const int *desca,
        complex *b, const int ib, const int jb, const int *descb,
        const int blacs_ctxt) 
    {
        pzgemr2d_(&M, &N,
                  a, &ia, &ja, desca,
                  b, &ib, &jb, descb,
                  &blacs_ctxt);
    }

    // Use Cpxgemr2d to collect matrices from all processes to root process
    template <typename mat>
    static void gatherMatrix(const int myid,
                             const int root_proc,
                             const mat& mat_l,
                             mat& mat_g)
    {
        auto a = mat_l.p;
        decltype(a) b;
        const int* desca = mat_l.desc;
        int ctxt = desca[1];
        int nrows = desca[2];
        int ncols = desca[3];

        if (myid == root_proc)
            b = new typename std::remove_reference<decltype(*a)>::type[nrows * ncols];
        else
            b = new typename std::remove_reference<decltype(*a)>::type[1];

        // Set descb, which has all elements in the only block in the root process
        int descb[9] = {1, ctxt, nrows, ncols, nrows, ncols, 0, 0, nrows};

        mat_g.desc = descb;
        mat_g.row = nrows;
        mat_g.col = ncols;
        mat_g.p = b;

        Cpxgemr2d(nrows, ncols, a, 1, 1, desca, b, 1, 1, descb, ctxt);
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
    void DiagoCusolver<T>::diag(hamilt::Hamilt<T>* phm_in,
                                psi::Psi<T>& psi,
                                Real* eigenvalue_in)
    {
        // Output the title for the current operation
        ModuleBase::TITLE("DiagoCusolver", "diag");

        // Create matrices for the Hamiltonian and overlap
        hamilt::MatrixBlock<T> h_mat, s_mat;
        phm_in->matrix(h_mat, s_mat);

#ifdef __MPI
        // global matrix
        hamilt::MatrixBlock<T> h_mat_g, s_mat_g;

        // global psi for distribute
        T* psi_g;

        // get the context and process information
        int ctxt = ParaV->blacs_ctxt;
        int nprows, npcols, myprow, mypcol;
        Cblacs_gridinfo(ctxt, &nprows, &npcols, &myprow, &mypcol);
        int myid = Cblacs_pnum(ctxt, myprow, mypcol);
        const int root_proc = Cblacs_pnum(ctxt, ParaV->desc[6], ParaV->desc[7]);

#endif

        // Allocate memory for eigenvalues
        std::vector<double> eigen(GlobalV::NLOCAL, 0.0);

        // Start the timer for the cusolver operation
        ModuleBase::timer::tick("DiagoCusolver", "cusolver");

#ifdef __MPI
        // gather matrices from processes to root process
        gatherMatrix(myid, root_proc, h_mat, h_mat_g);
        gatherMatrix(myid, root_proc, s_mat, s_mat_g);
#endif

        // Call the dense diagonalization routine
#ifdef __MPI
        MPI_Barrier(MPI_COMM_WORLD);
        if (myid == root_proc)
        {
            psi_g = new T[h_mat_g.row * h_mat_g.col];
            this->dc.Dngvd(h_mat_g.col, h_mat_g.row, h_mat_g.p, s_mat_g.p, eigen.data(), psi_g);
        }
        else
        {
            psi_g = new T[1];
        }
        MPI_Barrier(MPI_COMM_WORLD);
        // broadcast eigenvalues to all processes
        MPI_Bcast(eigen.data(), GlobalV::NBANDS, MPI_DOUBLE, root_proc, MPI_COMM_WORLD);

        // distribute psi to all processes
        distributePsi(this->ParaV->desc_wfc, psi.get_pointer(), psi_g);
#else
        // Call the dense diagonalization routine
        this->dc.Dngvd(h_mat.row, h_mat.col, h_mat.p, s_mat.p, eigen.data(), psi.get_pointer());
#endif
        // Stop the timer for the cusolver operation
        ModuleBase::timer::tick("DiagoCusolver", "cusolver");

        // Copy the eigenvalues and eigenvectors to the output arrays
        const int inc = 1;
        BlasConnector::copy(GlobalV::NBANDS, eigen.data(), inc, eigenvalue_in, inc);

        // Free allocated memory
#ifdef __MPI
        delete[] h_mat_g.p;
        delete[] s_mat_g.p;
        delete[] psi_g;
#endif
    }

    // Explicit instantiation of the DiagoCusolver class for real and complex numbers
    template class DiagoCusolver<double>;
    template class DiagoCusolver<complex>;

} // namespace hsolver
