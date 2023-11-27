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

// Define matrix types for real and complex numbers
typedef hamilt::MatrixBlock<double> matd;
typedef hamilt::MatrixBlock<std::complex<double>> matcd;

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
    void DiagoCusolver<double>::diag(hamilt::Hamilt<double>* phm_in, psi::Psi<double>& psi, Real* eigenvalue_in)
    {
        // Output the title for the current operation
        ModuleBase::TITLE("DiagoCusolver", "diag");

        // Create matrices for the Hamiltonian and overlap
        matd h_mat, s_mat;
        phm_in->matrix(h_mat, s_mat);

        // Allocate memory for eigenvalues
        std::vector<double> eigen(GlobalV::NLOCAL, 0.0);

        // Start the timer for the cusolver operation
        ModuleBase::timer::tick("DiagoCusolver", "cusolver");

        // Call the dense double diagonalization routine
        this->dc.Dngvd_double(h_mat.col, h_mat.row, h_mat.p, s_mat.p, eigen.data(), psi.get_pointer());

        // Stop the timer for the cusolver operation
        ModuleBase::timer::tick("DiagoCusolver", "cusolver");

        // Copy the eigenvalues to the output array
        const int inc = 1;
        BlasConnector::copy(GlobalV::NBANDS, eigen.data(), inc, eigenvalue_in, inc);
    }

} // namespace hsolver
