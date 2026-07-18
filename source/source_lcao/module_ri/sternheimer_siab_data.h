#ifndef STERNHEIMER_SIAB_DATA_H
#define STERNHEIMER_SIAB_DATA_H

#include <complex>
#include <string>
#include <vector>

namespace module_ri
{
namespace sternheimer_siab
{

struct PrimitiveBlock
{
    std::string element;
    int atom_index;
    int l;
    int m;
    int n_primitive;
    int offset;
};

struct ReferenceRow
{
    int occupied_state;
    int auxiliary_channel;
    double frequency_ha;
    double occupation;
    double frequency_weight;
    double norm;
    std::vector<std::complex<double>> q;
    int frequency_index;
};

// The eight fields below are required-core v1; Task 4 may append deterministic optional provenance fields.
struct Provenance
{
    std::string abacus_commit;
    std::string auxiliary_basis_sha256;
    std::vector<double> cell_bohr; ///< complete 3x3 row-major lattice vectors in Bohr
    double ecut_ry;
    std::string kernel;
    std::string orbital_sha256;
    std::string pseudopotential_sha256;
    std::string spin_convention;

    // Task 4 production fields.  Keeping these unset preserves the exact Task 3 canonical fixture.
    std::string executable_sha256;
    double exx_pca_thr = -1.0;
    int sternheimer_nfreq = 0;
    std::vector<double> frequency_ha;
    std::vector<double> frequency_weights_ha;
    int mpi_ranks = 0;
    int omp_threads = 0;
};

} // namespace sternheimer_siab
} // namespace module_ri

#endif
