#ifndef WRITE_HS_R_H
#define WRITE_HS_R_H

#include "source_base/matrix.h"
#include "source_basis/module_nao/two_center_bundle.h"
#include "source_cell/klist.h"
#include "source_hamilt/hamilt.h"
#include "source_lcao/LCAO_HS_arrays.hpp"
#include "source_lcao/module_dftu/dftu.h" // mohan add 20251107

#ifdef __EXX
#include "source_hamilt/module_xc/exx_info.h" // use GlobalC::exx_info
#include "RI/global/Tensor.h" // for RI::Tensor
#endif

namespace ModuleIO
{
void output_dHR(const int& istep,
                const ModuleBase::matrix& v_eff,
                const UnitCell& ucell,
                const Parallel_Orbitals& pv,
                LCAO_HS_Arrays& HS_Arrays,
                const Grid_Driver& grid, // mohan add 2024-04-06
                const TwoCenterBundle& two_center_bundle,
                const LCAO_Orbitals& orb,
                const K_Vectors& kv,
                const bool& binary = false,
                const double& sparse_threshold = 1e-10);

void output_dSR(const int& istep,
                const UnitCell& ucell,
                const Parallel_Orbitals& pv,
                LCAO_HS_Arrays& HS_Arrays,
                const Grid_Driver& grid, // mohan add 2024-04-06
                const TwoCenterBundle& two_center_bundle,
                const LCAO_Orbitals& orb,
                const K_Vectors& kv,
                const bool& binary = false,
                const double& sparse_thr = 1e-10);

void output_TR(const int istep,
               const UnitCell& ucell,
               const Parallel_Orbitals& pv,
               LCAO_HS_Arrays& HS_Arrays,
               const Grid_Driver& grid,
               const TwoCenterBundle& two_center_bundle,
               const LCAO_Orbitals& orb,
               const std::string& TR_filename = "trs1_nao.csr",
               const bool& binary = false,
               const double& sparse_threshold = 1e-10);

template <typename TK>
void output_SR(Parallel_Orbitals& pv,
               const Grid_Driver& grid,
               hamilt::Hamilt<TK>* p_ham,
               const std::string& SR_filename = "srs1_nao.csr",
               const bool& binary = false,
               const double& sparse_threshold = 1e-10);

/// Generate filename for HR/SR CSR output.
std::string hsr_gen_fname(const std::string& prefix,
                          const int ispin,
                          const bool append,
                          const int istep);

/// Generate filename for derivative matrices (dH/dR, dS/dR).
std::string dhr_gen_fname(const std::string& prefix,
                          const int ispin,
                          const bool append,
                          const int istep);

/// Write a single HContainer to CSR file with header.
template <typename TR>
void write_hcontainer_csr(const std::string& fname,
                          const UnitCell* ucell,
                          const int precision,
                          hamilt::HContainer<TR>* mat_serial,
                          const int istep,
                          const int ispin,
                          const int nspin,
                          const std::string& label);

/// Write H(R) and S(R) in CSR format, unified with write_dmr interface.
template <typename TR>
void write_hsr(const std::vector<hamilt::HContainer<TR>*>& hr_vec,
               const hamilt::HContainer<TR>* sr,
               const UnitCell* ucell,
               const int precision,
               const Parallel_2D& paraV,
               const bool append,
               const int* iat2iwt,
               const int nat,
               const int istep);

/// Write real-space matrix in CSR format (generic interface).
template <typename TR>
void write_matrix_r(const std::string& matrix_label,
                    const std::string& description,
                    const std::vector<hamilt::HContainer<TR>*>& matrices,
                    const UnitCell* ucell,
                    const int precision,
                    const Parallel_2D& paraV,
                    const bool append,
                    const int* iat2iwt,
                    const int nat,
                    const int istep);

} // namespace ModuleIO

#endif
