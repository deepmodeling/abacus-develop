#include "exx_lri.h"
#include "rpa_lri.h"
#include "rpa_lri_detail.h"
#include "source_base/global_function.h"
#include "source_basis/module_ao/elem_basis_idx_orb.h"
#include "source_estate/elecstate_lcao.h"
#include "source_io/module_parameter/input_parameter.h"
#include "source_lcao/module_ri/module_exx_symmetry/symm_rotation.h"

#include <complex>
#include <stdexcept>

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::out_eigen_vector_v1(const Parallel_Orbitals& parav, const psi::Psi<T>& psi)
{
#ifdef __MPI
    const int nks_tot = this->runtime.input.nspin == 2 ? p_kv->get_nks() / 2 : p_kv->get_nks();
    const int npsin_tmp = this->runtime.input.nspin == 2 ? 2 : 1;
    ModuleBase::timer::start("RPA_LRI", "out_eigen_vector_v1_mpi_io");
    try
    {
        const MPI_Comm io_comm = parav.comm();
        if (io_comm == MPI_COMM_NULL)
        {
            throw std::runtime_error("KS eigenvector MPI-IO writer has no wavefunction communicator.");
        }
        int communicator_relation = MPI_UNEQUAL;
        RpaLriDetail::collective_mpi_check(io_comm,
                                           MPI_Comm_compare(io_comm, this->mpi_comm, &communicator_relation),
                                           "Failed to compare KS eigenvector MPI communicators");
        RpaLriDetail::collective_require(io_comm,
                                         communicator_relation == MPI_IDENT || communicator_relation == MPI_CONGRUENT,
                                         "KS eigenvector wavefunction and RPA communicators are inconsistent");
        RpaLriDetail::write_ks_eigenvector_v1_mpi(io_comm,
                                                  parav,
                                                  psi,
                                                  nks_tot,
                                                  npsin_tmp,
                                                  this->runtime.input.nspin,
                                                  this->runtime.input.nbands,
                                                  this->runtime.nlocal,
                                                  this->runtime.input.rpa_outdir + "KS_wfc_0.dat");
    }
    catch (...)
    {
        ModuleBase::timer::end("RPA_LRI", "out_eigen_vector_v1_mpi_io");
        throw;
    }
    ModuleBase::timer::end("RPA_LRI", "out_eigen_vector_v1_mpi_io");
#endif
}

template void RPA_LRI<double, double>::out_eigen_vector_v1(const Parallel_Orbitals& parav, const psi::Psi<double>& psi);
template void RPA_LRI<std::complex<double>, double>::out_eigen_vector_v1(const Parallel_Orbitals& parav,
                                                                         const psi::Psi<std::complex<double>>& psi);
