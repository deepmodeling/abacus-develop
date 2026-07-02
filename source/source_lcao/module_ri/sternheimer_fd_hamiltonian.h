#ifndef STERNHEIMER_FD_HAMILTONIAN_H
#define STERNHEIMER_FD_HAMILTONIAN_H

#include <complex>
#include <vector>

namespace ModuleRI
{

class SternheimerFDHamiltonian
{
  public:
    using Complex = std::complex<double>;
    using Vector = std::vector<Complex>;
    using Matrix = std::vector<Vector>;

    struct Eigenpairs
    {
        std::vector<double> eigenvalues;
        std::vector<Vector> eigenvectors;
    };

    struct Grid
    {
        int nx = 0;
        int ny = 0;
        int nz = 0;
        double hx = 0.0;
        double hy = 0.0;
        double hz = 0.0;
        bool periodic = true;

        int size() const;
    };

    // kinetic_prefactor = 0.5 gives the Hartree convention -1/2 laplacian.
    // ABACUS PW eigenvalues are in Ry, so ABACUS bindings pass 1.0.
    SternheimerFDHamiltonian(Grid grid, std::vector<double> local_potential, double kinetic_prefactor = 0.5);

    const Grid& grid() const;
    const std::vector<double>& local_potential() const;
    double kinetic_prefactor() const;

    void apply(const Vector& psi, Vector& hpsi) const;

    Matrix dense_matrix(int max_size = 4096) const;

    // Diagonalizes the zero-order Hermitian FD Hamiltonian on small debug grids.
    // The shifted Sternheimer operator H - eps + i omega is non-Hermitian and
    // must be handled as a linear operator, not by this routine.
    Eigenpairs diagonalize_dense(int max_size = 4096) const;

    int index(int ix, int iy, int iz) const;

  private:
    int shifted_index(int ix, int iy, int iz) const;

    Grid grid_;
    std::vector<double> local_potential_;
    double kinetic_prefactor_ = 0.5;
};

} // namespace ModuleRI

#endif
