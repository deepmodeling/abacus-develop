#ifndef RPA_LRI_DETAIL_H
#define RPA_LRI_DETAIL_H

#include "source_psi/psi.h"

#include <complex>
#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <string>
#include <vector>
#include <mpi.h>

class Parallel_Orbitals;
class UnitCell;
class Numerical_Orbital_Lm;
namespace RI
{
template <typename T>
class Tensor;
}

namespace RpaLriDetail
{
constexpr int LIBRPA_COULOMB_V1_MARKER = -20129433;
constexpr int LIBRPA_LRICOEF_V1_MARKER = -10267453;
constexpr int LIBRPA_SHRINK_SINVS_V1_MARKER = -30241621;
constexpr int LIBRPA_KS_EIGENVECTOR_V1_MARKER = -12345679;
constexpr int LIBRPA_KS_EIGENVECTOR_V1_KIND_COMPLEX_DOUBLE = 28;
constexpr int LIBRPA_COULOMB_V1_COMPLEX_FLAG = 1;

static_assert(sizeof(std::complex<double>) == 2 * sizeof(double),
              "LibRPA v1 binary output expects complex<double> as two doubles.");

void trim_malloc_cache();

bool debug_dump_exx_ao_enabled();

std::size_t coulomb_atom_pair_index(const std::size_t I, const std::size_t J, const std::size_t natoms);

void checked_write(std::ofstream& ofs, const void* data, const std::size_t bytes, const std::string& filename);

unsigned long long checked_mul_u64(const unsigned long long lhs,
                                   const unsigned long long rhs,
                                   const std::string& context);

unsigned long long checked_add_u64(const unsigned long long lhs,
                                   const unsigned long long rhs,
                                   const std::string& context);

std::int64_t checked_i64_from_u64(const unsigned long long value, const std::string& context);

std::int32_t checked_i32_from_size(const std::size_t value, const std::string& context);

std::int32_t checked_i32_from_int(const int value, const std::string& context);

int sum_int_vector(const std::vector<int>& values);

double real_as_double(const std::complex<double>& value);

std::vector<std::vector<int>> collect_abfs_l_nchi(
    const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& abfs);

std::vector<std::vector<int>> collect_wfc_l_nchi(const UnitCell& ucell);

int basis_size_from_shell_counts(const std::vector<int>& shell_counts);

void write_librpa_split_basis_file(const UnitCell& ucell,
                                   const std::vector<int>& type_sizes,
                                   const std::vector<std::vector<int>>& type_shell_counts,
                                   const std::string& filename);

template <typename Value>
void write_scalar(std::ofstream& ofs, const Value& value, const std::string& filename);
double real_as_double(double value);
bool has_valid_matrix_shape(const RI::Tensor<double>& tensor);
void abort_output(MPI_Comm comm, const std::string& message);

#ifdef __MPI
struct KSEigenvectorMpiLayout
{
    MPI_Datatype filetype = MPI_C_DOUBLE_COMPLEX;
    bool free_filetype = false;
    unsigned long long local_count = 0;
    unsigned long long max_local_count = 0;
};

std::string mpi_error_string(const int error_code);

void collective_mpi_check(const MPI_Comm mpi_comm, const int local_error, const std::string& context);

void collective_require(const MPI_Comm mpi_comm, const bool local_condition, const std::string& context);

MPI_Aint checked_mpi_aint_from_u64(const unsigned long long value, const std::string& context);

MPI_Offset checked_mpi_offset_from_u64(const unsigned long long value, const std::string& context);

KSEigenvectorMpiLayout make_ks_eigenvector_mpi_layout(const MPI_Comm mpi_comm,
                                                      const Parallel_Orbitals& parav,
                                                      const int nbands,
                                                      const int nbasis_wfc,
                                                      const bool is_soc,
                                                      const int spinor_component);

int free_ks_eigenvector_mpi_layout(KSEigenvectorMpiLayout& layout);

template <typename T>
void write_ks_eigenvector_v1_mpi(const MPI_Comm mpi_comm,
                                 const Parallel_Orbitals& parav,
                                 const psi::Psi<T>& psi,
                                 const int nks_tot,
                                 const int nspins,
                                 const int nspin_abacus,
                                 const int nbands,
                                 const int nbasis_wfc,
                                 const std::string& final_name);

#endif
} // namespace RpaLriDetail
#endif
