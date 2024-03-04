#include "diago_cusolver.h"
#include "module_base/global_variable.h"
#include "module_base/lapack_connector.h"
#include "module_base/timer.h"
#include "module_base/tool_quit.h"
#include "mpi.h"

extern "C"
{
#include "module_hsolver/genelpa/Cblacs.h"
#include "module_hsolver/genelpa/scalapack.h"
}


// Define matrix types for real and complex numbers
typedef hamilt::MatrixBlock<double> matd;
typedef hamilt::MatrixBlock<std::complex<double>> matcd;

static inline void Cpdgemr2d(int M,
                            int N,
                            double* a,
                            int ia,
                            int ja,
                            int* desca,
                            double* b,
                            int ib,
                            int jb,
                            int* descb,
                            int blacs_ctxt)
{
    pdgemr2d_(&M, &N, a, &ia, &ja, desca, b, &ib, &jb, descb, &blacs_ctxt);
}

// use pdgemr2d to collect matrix from all processes to root process
void gatherMatrix(matd &mat_l, matd &mat_g)
{
    double *a;a = mat_l.p;
    int *desca;desca = mat_l.desc;
    double *b;
    // int *descb = mat_g.desc;
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
        b = new double[nrows * ncols];
    else
        b = new double[1];

    // set descb, which has all elements in the only block in the root process
    int descb[9] = {1, ctxt, nrows, ncols, nrows, ncols, 0, 0, nrows};
    mat_g.desc = descb;
    mat_g.p = b;

    Cpdgemr2d(nrows, ncols, a, 1, 1, desca, b, 1, 1, descb, ctxt);
}

// convert the Psi to a 2D block storage format
void distributePsi(int ctxt, int nrows, int ncols, int rsrc, int csrc, double *psi, double *psi_g)
{
    int descg[9] = {1, ctxt, nrows, ncols, nrows, ncols, rsrc, csrc, nrows};

    // not sure whether nprows and npcols equal to 1
    // nrows = nbasis
    int descl[9] = {1, ctxt, nrows, ncols, 1, 1, rsrc, csrc, nrows};
    Cpdgemr2d(nrows, ncols, psi_g, 1, 1, descg, psi, 1, 1, descl, ctxt);
}

// Namespace for the diagonalization solver
namespace hsolver
{
    // Initialize the DecomposedState variable for real and complex numbers
    template <>
    int DiagoCusolver<double>::DecomposedState = 0;
    template <>
    int DiagoCusolver<std::complex<double>>::DecomposedState = 0;

    // Diagonalization function for complex numbers
    template <>
    void DiagoCusolver<std::complex<double>>::diag(hamilt::Hamilt<std::complex<double>>* phm_in,
                                                   psi::Psi<std::complex<double>>& psi,
                                                   Real* eigenvalue_in)
    {
        // Output the title for the current operation
        ModuleBase::TITLE("DiagoCusolver", "diag");

        // Create matrices for the Hamiltonian and overlap
        matcd h_mat, s_mat;
        phm_in->matrix(h_mat, s_mat);

        // Calculate the size based on the number of bands and basis functions
        int size = psi.get_nbands() * psi.get_nbasis();

        // Allocate memory for eigenvalues and eigenvectors
        std::vector<double> eigen(GlobalV::NLOCAL, 0.0);
        std::complex<double>* eigenvectors = new std::complex<double>[h_mat.row * h_mat.col];

        // Start the timer for the cusolver operation
        ModuleBase::timer::tick("DiagoCusolver", "cusolver");

        // Call the dense complex diagonalization routine
        this->dc.Dngvd_complex(h_mat.row, h_mat.col, h_mat.p, s_mat.p, eigen.data(), eigenvectors);

        // Stop the timer for the cusolver operation
        ModuleBase::timer::tick("DiagoCusolver", "cusolver");

        // Copy the eigenvalues and eigenvectors to the output arrays
        const int inc = 1;
        BlasConnector::copy(GlobalV::NBANDS, eigen.data(), inc, eigenvalue_in, inc);
        BlasConnector::copy(size, eigenvectors, inc, psi.get_pointer(), inc);

        // Free allocated memory
        delete[] eigenvectors;
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
        double *psi_g;
        int ctxt = h_mat.desc[1];
        int nprows, npcols, myprow, mypcol;
        Cblacs_gridinfo(ctxt, &nprows, &npcols, &myprow, &mypcol);
        int myid = Cblacs_pnum(ctxt, myprow, mypcol);
        const int ROOT_PROC = Cblacs_pnum(ctxt, h_mat.desc[6], h_mat.desc[7]);
        #endif

        // Allocate memory for eigenvalues
        std::vector<double> eigen(GlobalV::NLOCAL, 0.0);

        // Start the timer for the cusolver operation
        ModuleBase::timer::tick("DiagoCusolver", "cusolver");

        #ifdef __MPI
        // gather matrixs from processes
        gatherMatrix(h_mat, h_mat_g);
        gatherMatrix(s_mat, s_mat_g);
        #endif

        // Call the dense double diagonalization routine
        #ifdef  __MPI
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
        MPI_Bcast(eigen.data(), GlobalV::NBANDS, MPI_DOUBLE, ROOT_PROC, MPI_COMM_WORLD);
        distributePsi(ctxt, h_mat_g.row, h_mat_g.col, 0, 0, psi.get_pointer(), psi_g);
        #else
        this->dc.Dngvd_double(h_mat.col, h_mat.row, h_mat.p, s_mat.p, eigen.data(), psi.get_pointer());
        #endif

        
        // Stop the timer for the cusolver operation
        ModuleBase::timer::tick("DiagoCusolver", "cusolver");

        // Copy the eigenvalues to the output array
        const int inc = 1;
        BlasConnector::copy(GlobalV::NBANDS, eigen.data(), inc, eigenvalue_in, inc);
        #ifdef __MPI
        delete[] h_mat_g.p;
        delete[] s_mat_g.p;
        delete[] psi_g;
        #endif 
    }

} // namespace hsolver
