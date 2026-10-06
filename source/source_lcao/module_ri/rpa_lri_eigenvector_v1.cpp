#include "exx_lri.h"
#include "rpa_lri.h"
#include "rpa_lri_detail.h"
#include "source_base/global_function.h"
#include "source_basis/module_ao/elem_basis_idx_orb.h"
#include "source_estate/elecstate_lcao.h"
#include "source_io/module_parameter/input_parameter.h"
#include "source_lcao/module_ri/module_exx_symmetry/symm_rotation.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::out_eigen_vector_v1(const Parallel_Orbitals& parav, const psi::Psi<T>& psi)
{
    const int nks_tot = this->runtime.input.nspin == 2 ? p_kv->get_nks() / 2 : p_kv->get_nks();
    const int npsin_tmp = this->runtime.input.nspin == 2 ? 2 : 1;
    const int nbands = parav.get_wfc_global_nbands();
    const int nbasis = parav.get_wfc_global_nbasis();
    const std::complex<double> zero(0.0, 0.0);
#ifdef __MPI
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
#endif
#ifdef __MPI
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
#ifndef __MPI
    struct KSEigenRecord
    {
        std::int32_t ik = 0;
        std::int64_t payload_offset = 0;
        std::vector<std::complex<double>> payload;
    };

    std::vector<KSEigenRecord> records;
    records.reserve(static_cast<std::size_t>(nks_tot));

    for (int ik = 0; ik < nks_tot; ik++)
    {
        std::vector<ModuleBase::ComplexMatrix> is_wfc_ib_iw(npsin_tmp);
        for (int is = 0; is < npsin_tmp; is++)
        {
            is_wfc_ib_iw[is].create(this->runtime.input.nbands, this->runtime.nlocal);
            for (int ib_global = 0; ib_global < this->runtime.input.nbands; ++ib_global)
            {
                std::vector<std::complex<double>> wfc_iks(this->runtime.nlocal, zero);

                const int ib_local = parav.global2local_col(ib_global);

                if (ib_local >= 0)
                {
                    for (int ir = 0; ir < psi.get_nbasis(); ir++)
                    {
                        wfc_iks[parav.local2global_row(ir)] = psi(ik + nks_tot * is, ib_local, ir);
                    }
                }

                for (int iw = 0; iw < this->runtime.nlocal; iw++)
                {
                    is_wfc_ib_iw[is](ib_global, iw) = wfc_iks[iw];
                }
            }
        }

        if (this->runtime.rank == 0)
        {
            KSEigenRecord record;
            record.ik = static_cast<std::int32_t>(ik + 1);
            record.payload.reserve(static_cast<std::size_t>(npsin_tmp)
                                   * static_cast<std::size_t>(this->runtime.input.nbands)
                                   * static_cast<std::size_t>(this->runtime.nlocal));
            if (this->runtime.input.nspin == 4)
            {
                if (this->runtime.nlocal % 2 != 0)
                {
                    throw std::runtime_error("SOC KS eigenvector output expects an even basis size.");
                }
                const int nlocal_ao = this->runtime.nlocal / 2;
                for (int isoc = 0; isoc < 2; ++isoc)
                {
                    for (int ib = 0; ib < this->runtime.input.nbands; ++ib)
                    {
                        for (int iw = 0; iw < nlocal_ao; ++iw)
                        {
                            record.payload.push_back(is_wfc_ib_iw[0](ib, iw * 2 + isoc));
                        }
                    }
                }
            }
            else
            {
                for (int is = 0; is < npsin_tmp; ++is)
                {
                    for (int ib = 0; ib < this->runtime.input.nbands; ++ib)
                    {
                        for (int iw = 0; iw < this->runtime.nlocal; ++iw)
                        {
                            record.payload.push_back(is_wfc_ib_iw[is](ib, iw));
                        }
                    }
                }
            }
            records.push_back(std::move(record));
        }
    }

    if (this->runtime.rank == 0)
    {
        const std::string out_name = this->runtime.input.rpa_outdir + "KS_wfc_0.dat";
        const std::int64_t record_bytes
            = static_cast<std::int64_t>(sizeof(std::int32_t)) + static_cast<std::int64_t>(sizeof(std::int64_t));
        std::int64_t offset = 6 * static_cast<std::int64_t>(sizeof(std::int32_t))
                              + static_cast<std::int64_t>(records.size()) * record_bytes;
        for (auto& record: records)
        {
            record.payload_offset = offset;
            offset += static_cast<std::int64_t>(record.payload.size() * sizeof(std::complex<double>));
        }

        std::ofstream ofs(out_name.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
        if (!ofs.good())
        {
            throw std::runtime_error("Failed to open " + out_name);
        }

        const std::int32_t marker = RpaLriDetail::LIBRPA_KS_EIGENVECTOR_V1_MARKER;
        const std::int32_t kind = RpaLriDetail::LIBRPA_KS_EIGENVECTOR_V1_KIND_COMPLEX_DOUBLE;
        const std::int32_t nkpoints_local
            = RpaLriDetail::checked_i32_from_size(records.size(), "KS eigenvector k-point count");
        const std::int32_t nspins = RpaLriDetail::checked_i32_from_int(npsin_tmp, "KS eigenvector spin count");
        const std::int32_t nstates
            = RpaLriDetail::checked_i32_from_int(this->runtime.input.nbands, "KS eigenvector state count");
        const std::int32_t nbasis_wfc
            = RpaLriDetail::checked_i32_from_int(this->runtime.nlocal, "KS eigenvector basis count");
        RpaLriDetail::write_scalar(ofs, marker, out_name);
        RpaLriDetail::write_scalar(ofs, kind, out_name);
        RpaLriDetail::write_scalar(ofs, nkpoints_local, out_name);
        RpaLriDetail::write_scalar(ofs, nspins, out_name);
        RpaLriDetail::write_scalar(ofs, nstates, out_name);
        RpaLriDetail::write_scalar(ofs, nbasis_wfc, out_name);
        for (const auto& record: records)
        {
            RpaLriDetail::write_scalar(ofs, record.ik, out_name);
            RpaLriDetail::write_scalar(ofs, record.payload_offset, out_name);
        }
        for (const auto& record: records)
        {
            RpaLriDetail::checked_write(ofs,
                                        record.payload.data(),
                                        record.payload.size() * sizeof(std::complex<double>),
                                        out_name);
        }
        ofs.close();
    }
#endif
}

template void RPA_LRI<double, double>::out_eigen_vector_v1(const Parallel_Orbitals& parav, const psi::Psi<double>& psi);
template void RPA_LRI<std::complex<double>, double>::out_eigen_vector_v1(const Parallel_Orbitals& parav,
                                                                         const psi::Psi<std::complex<double>>& psi);
