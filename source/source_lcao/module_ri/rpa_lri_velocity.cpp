#include "exx_lri.h"
#include "rpa_lri.h"
#include "rpa_lri_detail.h"
#include "source_base/global_function.h"
#include "source_basis/module_ao/elem_basis_idx_orb.h"
#include "source_estate/elecstate_lcao.h"
#include "source_io/module_parameter/input_parameter.h"
#include "source_lcao/module_lr/utils/spectrum_mo.hpp"
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
void RPA_LRI<T, Tdata>::out_velocity(const UnitCell& ucell,
                                     const Grid_Driver& gd,
                                     const TwoCenterBundle& two_center_bundle,
                                     const Parallel_Orbitals& parav,
                                     const psi::Psi<T>& psi,
                                     const elecstate::ElecState* pelec)
{
    ModuleBase::TITLE("DFT_RPA_interface", "out_velocity");
    ModuleBase::timer::start("RPA_LRI", "out_velocity");

    Parallel_2D parac;
    LR_Util::setup_2d_division(parac,
                               parav.get_block_size(),
                               this->runtime.nlocal,
                               this->runtime.input.nbands
#ifdef __MPI
                               ,
                               parav.blacs_ctxt
#endif
    );

    const int nk = this->runtime.input.nspin == 2 ? p_kv->get_nks() / 2 : p_kv->get_nks();
    const int nspin_tmp = this->runtime.input.nspin == 2 ? 2 : 1;
    const int nbands = parav.get_wfc_global_nbands();
    const int nbasis = parav.get_wfc_global_nbasis();

    std::vector<int> nocc(2, nbands);
    std::vector<int> nvirt(2, 0);
    const std::vector<std::complex<double>> velocity_mo = LR_Util::cal_velocity_mo(ucell,
                                                                                   gd,
                                                                                   two_center_bundle,
                                                                                   parav,
                                                                                   parac,
                                                                                   *this->p_kv,
                                                                                   psi,
                                                                                   nk,
                                                                                   nspin_tmp,
                                                                                   this->runtime.nlocal,
                                                                                   nocc,
                                                                                   nvirt);
    if (this->runtime.rank == 0)
    {
        LR_Util::output_spectrum_mo_librpa(velocity_mo,
                                           this->runtime.input.rpa_outdir + "velocity_matrix.txt",
                                           nk,
                                           nspin_tmp,
                                           nbands,
                                           nbasis,
                                           nbands,
                                           *this->p_kv);
    }
    ModuleBase::timer::end("RPA_LRI", "out_velocity");
}

template void RPA_LRI<double, double>::out_velocity(const UnitCell& ucell,
                                                    const Grid_Driver& gd,
                                                    const TwoCenterBundle& two_center_bundle,
                                                    const Parallel_Orbitals& parav,
                                                    const psi::Psi<double>& psi,
                                                    const elecstate::ElecState* pelec);
template void RPA_LRI<std::complex<double>, double>::out_velocity(const UnitCell& ucell,
                                                                  const Grid_Driver& gd,
                                                                  const TwoCenterBundle& two_center_bundle,
                                                                  const Parallel_Orbitals& parav,
                                                                  const psi::Psi<std::complex<double>>& psi,
                                                                  const elecstate::ElecState* pelec);
