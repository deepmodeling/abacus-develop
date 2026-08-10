#ifndef WRITE_H_TERMS_H
#define WRITE_H_TERMS_H

#include "source_basis/module_nao/two_center_bundle.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/klist.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_estate/module_charge/charge.h"
#include "source_estate/module_pot/potential_new.h"
#include "source_lcao/lcao_domain.h"
#include "source_hamilt/module_hcontainer/hcontainer.h"

#include <complex>
#include <vector>

template <typename T, typename Tdata>
class Exx_LRI_Interface;

namespace ModuleIO
{

struct WriteHParams
{
    const UnitCell* ucell = nullptr;
    const Grid_Driver* gd = nullptr;
    const Parallel_Orbitals* pv = nullptr;
    const TwoCenterBundle* two_center_bundle = nullptr;
    const LCAO_Orbitals* orb = nullptr;
    const K_Vectors* kv = nullptr;
    const elecstate::Potential* pot = nullptr;   // used by write_h_vl (local pp only)
    const Charge* chg = nullptr;                 // used by write_h_vh, write_h_vxc
    const ModulePW::PW_Basis* rho_basis = nullptr; // used by write_h_vh
    int nrxx = 0;                                // used by write_h_vxc
    int nspin = 1;
    int istep = 0;
    bool append = false;
    const int* iat2iwt = nullptr;
    int nat = 0;
    bool also_hR = false; // H(k) is always written; H(R) (CSR) only when this is true
#ifdef __EXX
    // The gamma-only (TK==double) exx interfaces used by write_h_exx. 
    // Deliberately NOT templated on TK, because it would force WriteHParams, WriteDHParams and
    //      every free function taking them to become templates as well -- a large, purely
    //      mechanical change for a case nobody needs.
    // Multi-k + EXX is therefore rejected up front (see write_h_exx)
    // instead of silently producing output with the EXX term missing.
    Exx_LRI_Interface<double, double>* exd = nullptr;
    Exx_LRI_Interface<double, std::complex<double>>* exc = nullptr;
#endif
};

void write_h_t(WriteHParams& params);

void write_h_vnl(WriteHParams& params);

void write_h_vl(WriteHParams& params);

void write_h_vh(WriteHParams& params);

void write_h_vxc(WriteHParams& params);

#ifdef __EXX
// Build V^EXX(R) into a real HContainer via add_HexxR (from exd/exc->get_Hexxs()) and write it.
// exd (real Hexx) and exc (complex Hexx) are mutually exclusive; picked by info_ri.real_number.
void write_h_exx(WriteHParams& params);
#endif

} // namespace ModuleIO

#endif
