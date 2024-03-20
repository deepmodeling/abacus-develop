#include "diago_cusolver.h"
#include "module_base/global_variable.h"
#include "module_base/lapack_connector.h"
#include "module_base/timer.h"

#include <type_traits>

extern "C"
{
#include "module_hsolver/genelpa/Cblacs.h"
#include "module_hsolver/genelpa/scalapack.h"
}

// Define matrix types for real and complex numbers
using complex = std::complex<double>;
using matd = hamilt::MatrixBlock<double>;
using matcd = hamilt::MatrixBlock<complex>;

// Namespace for the diagonalization solver
namespace hsolver
{
    // Initialize the DecomposedState variable for real and complex numbers
    template <>
    int DiagoCusolver<double>::DecomposedState = 0;
    template <>
    int DiagoCusolver<complex>::DecomposedState = 0;

    template <>
    DiagoCusolver<double>::DiagoCusolver(const Parallel_Orbitals* ParaV)
    {
        this->ParaV = ParaV;
    }

    template <>
    DiagoCusolver<complex>::DiagoCusolver(const Parallel_Orbitals* ParaV)
    {
        this->ParaV = ParaV;
    }

    template <>
    DiagoCusolver<double>::~DiagoCusolver()
    {
    }

    template <>
    DiagoCusolver<complex>::~DiagoCusolver()
    {
    }

    static inline void Cpxgemr2d(int M, int N, double *a, int ia, int ja,
                                 int *desca, double *b, int ib, int jb,
                                 int *descb, int blacs_ctxt) {
      pdgemr2d_(&M, &N, a, &ia, &ja, desca, b, &ib, &jb, descb, &blacs_ctxt);
    }

    static inline void Cpxgemr2d(int M, int N, complex *a, int ia, int ja,
                                 int *desca, complex *b, int ib, int jb,
                                 int *descb, int blacs_ctxt) {
      pzgemr2d_(&M, &N, reinterpret_cast<double __complex__ *>(a), &ia, &ja,
                desca, reinterpret_cast<double __complex__ *>(b), &ib, &jb,
                descb, &blacs_ctxt);
    }

    // Use pdgemr2d to collect matrices from all processes to root process
    template <typename mat>
    static void gatherMatrix(const mat& mat_l, mat& mat_g)
    {
        auto a = mat_l.p;
        decltype(a) b;
        int* desca = mat_l.desc;
        int ctxt = desca[1];
        int nrows = desca[2];
        int ncols = desca[3];
        int nlrows = desca[4];
        int nlcols = desca[5];
        int nprows, npcols, myprow, mypcol;

        Cblacs_gridinfo(ctxt, &nprows, &npcols, &myprow, &mypcol);
        int myid = Cblacs_pnum(ctxt, myprow, mypcol);

        const int ROOT_PROC = Cblacs_pnum(ctxt, desca[6], desca[7]);

        mat_g.row = nrows;
        mat_g.col = ncols;

        if (myid == ROOT_PROC)
            b = new typename std::remove_reference<decltype(*a)>::type[nrows * ncols];
        else
            b = new typename std::remove_reference<decltype(*a)>::type[1];

        // Set descb, which has all elements in the only block in the root process
        int descb[9] = {1, ctxt, nrows, ncols, nrows, ncols, 0, 0, nrows};
        mat_g.desc = descb;
        mat_g.p = b;

        Cpxgemr2d(nrows, ncols, a, 1, 1, desca, b, 1, 1, descb, ctxt);
    }

    // Convert the Psi to a 2D block storage format
    template <typename T>
    static void distributePsi(const int* desc_psi, T* psi, T* psi_g)
    {
        ModuleBase::timer::tick("DiagoCusolver", "distribute");
        int ctxt = desc_psi[1];
        int nrows = desc_psi[2];
        int ncols = desc_psi[3];
        int rsrc = desc_psi[6];
        int csrc = desc_psi[7];

        int descg[9] = {1, ctxt, nrows, ncols, nrows, ncols, rsrc, csrc, nrows};
        int descl[9];

        std::copy(desc_psi, desc_psi + 9, descl);

        Cpxgemr2d(nrows, ncols, psi_g, 1, 1, descg, psi, 1, 1, descl, ctxt);

        ModuleBase::timer::tick("DiagoCusolver", "distribute");
    }

    // Diagonalization function for complex numbers
    template <>
    void DiagoCusolver<complex>::diag(hamilt::Hamilt<complex>* phm_in, psi::Psi<complex>& psi, Real* eigenvalue_in)
    {
        // Output the title for the current operation
        ModuleBase::TITLE("DiagoCusolver", "diag");

        // Create matrices for the Hamiltonian and overlap
        matcd h_mat, s_mat;
        phm_in->matrix(h_mat, s_mat);

#ifdef __MPI
        // global matrix
        matcd h_mat_g, s_mat_g;

        // global psi for distribute
        complex* psi_g;

        // get the context and process information
        int ctxt = ParaV->blacs_ctxt;
        int nprows, npcols, myprow, mypcol;
        Cblacs_gridinfo(ctxt, &nprows, &npcols, &myprow, &mypcol);
        int myid = Cblacs_pnum(ctxt, myprow, mypcol);
        const int ROOT_PROC = Cblacs_pnum(ctxt, ParaV->desc[6], ParaV->desc[7]);

#endif

        // Allocate memory for eigenvalues
        std::vector<double> eigen(GlobalV::NLOCAL, 0.0);

        // Start the timer for the cusolver operation
        ModuleBase::timer::tick("DiagoCusolver", "cusolver");

#ifdef __MPI
        // gather matrices from processes to root process
        gatherMatrix(h_mat, h_mat_g);
        gatherMatrix(s_mat, s_mat_g);
#endif

        // Call the dense complex diagonalization routine
#ifdef __MPI
        MPI_Barrier(MPI_COMM_WORLD);
        if (myid == ROOT_PROC)
        {
            psi_g = new complex[h_mat_g.row * h_mat_g.col];
            this->dc.Dngvd_complex(h_mat_g.col, h_mat_g.row, h_mat_g.p, s_mat_g.p, eigen.data(), psi_g);
        }
        else
        {
            psi_g = new complex[1];
        }
        MPI_Barrier(MPI_COMM_WORLD);
        // broadcast eigenvalues to all processes
        MPI_Bcast(eigen.data(), GlobalV::NBANDS, MPI_DOUBLE, ROOT_PROC, MPI_COMM_WORLD);

        // distribute psi to all processes
        distributePsi(this->ParaV->desc_wfc, psi.get_pointer(), psi_g);
#else
        // Call the dense complex diagonalization routine
        this->dc.Dngvd_complex(h_mat.row, h_mat.col, h_mat.p, s_mat.p, eigen.data(), psi.get_pointer());
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

    // Diagonalization function for real numbers
    template <>
    void DiagoCusolver<double>::diag(hamilt::Hamilt<double>* phm_in, psi::Psi<double>& psi, Real* eigenvalue_in)
    {
        // Output the title for the current operation
        ModuleBase::TITLE("DiagoCusolver", "diag");

        // Create matrices for the Hamiltonian and overlap
        matd h_mat, s_mat;
        phm_in->matrix(h_mat, s_mat);

#ifdef __MPI
        // global matrix
        matd h_mat_g, s_mat_g;

        // global psi for distribute
        double *psi_g;

        // get the context and process information
        int ctxt = ParaV->blacs_ctxt;
        int nprows, npcols, myprow, mypcol;
        Cblacs_gridinfo(ctxt, &nprows, &npcols, &myprow, &mypcol);
        int myid = Cblacs_pnum(ctxt, myprow, mypcol);
        const int ROOT_PROC = Cblacs_pnum(ctxt, ParaV->desc[6], ParaV->desc[7]);

#endif

        // Allocate memory for eigenvalues
        std::vector<double> eigen(GlobalV::NLOCAL, 0.0);

        // Start the timer for the cusolver operation
        ModuleBase::timer::tick("DiagoCusolver", "cusolver");

#ifdef __MPI
        // gather matrices from processes to root process
        gatherMatrix(h_mat, h_mat_g);
        gatherMatrix(s_mat, s_mat_g);
#endif

        // Call the dense double diagonalization routine
#ifdef __MPI
        MPI_Barrier(MPI_COMM_WORLD);
        if (myid == ROOT_PROC)
        {
            psi_g = new double[h_mat_g.row * h_mat_g.col];
            this->dc.Dngvd_double(h_mat_g.col, h_mat_g.row, h_mat_g.p, s_mat_g.p, eigen.data(), psi_g);
        }
        else
        {
            psi_g = new double[1];
        }
        MPI_Barrier(MPI_COMM_WORLD);
        // broadcast eigenvalues to all processes
        MPI_Bcast(eigen.data(), GlobalV::NBANDS, MPI_DOUBLE, ROOT_PROC, MPI_COMM_WORLD);

        // distribute psi to all processes
        distributePsi(this->ParaV->desc_wfc, psi.get_pointer(), psi_g);
#else
        this->dc.Dngvd_double(h_mat.col, h_mat.row, h_mat.p, s_mat.p, eigen.data(), psi.get_pointer());
#endif

        // Stop the timer for the cusolver operation
        ModuleBase::timer::tick("DiagoCusolver", "cusolver");

        // Copy the eigenvalues to the output array
        const int inc = 1;
        BlasConnector::copy(GlobalV::NBANDS, eigen.data(), inc, eigenvalue_in, inc);

        // Free allocated memory
#ifdef __MPI
        delete[] h_mat_g.p;
        delete[] s_mat_g.p;
        delete[] psi_g;
#endif
    }

} // namespace hsolver
