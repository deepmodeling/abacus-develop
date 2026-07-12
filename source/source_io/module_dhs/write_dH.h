#ifndef WRITE_DH_H
#define WRITE_DH_H

#include "source_basis/module_nao/two_center_bundle.h"
#include "source_cell/klist.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_estate/module_pot/potential_new.h"
#include "source_lcao/LCAO_domain.h"
#include "source_lcao/module_hcontainer/hcontainer.h"

#include <array>
#include <complex>
#include <vector>

///for lack of make_unique in c++11 
template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args &&... args)
{
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

template <typename T, typename Tdata>
class Exx_LRI_Interface;

namespace ModuleIO
{

struct WriteDHParams
{
    const UnitCell* ucell = nullptr;
    const Grid_Driver* gd = nullptr;
    const Parallel_Orbitals* pv = nullptr;
    const TwoCenterBundle* two_center_bundle = nullptr;
    const LCAO_Orbitals* orb = nullptr;
    const K_Vectors* kv = nullptr;
    const ModuleBase::matrix* v_eff = nullptr;
    const int* iat2iwt = nullptr;
    elecstate::Potential* pot = nullptr;
    // Dedicated single-component potentials for the Veff-based dH terms. pelec->pot mixes
    // V^L + V^H + V^XC in get_eff_v(), so cal_dH would read the wrong potential for the
    // separate V^H / V^XC (and V^L) outputs. Each of these is built with exactly one
    // component registered ("local" / "hartree" / "xc"); see ctrl_scf_lcao.
    elecstate::Potential* pot_vl = nullptr;
    elecstate::Potential* pot_vh = nullptr;
    elecstate::Potential* pot_vxc = nullptr;
    int nat = 0;
    int nspin = 1;
    int istep = 0;
    bool gamma_only = false;
    bool append = false;
    bool also_dhR = false; // whether to write the real-space dH(R) in addition to the k-space dH(k)
    // per-spin real-space DM (size nspin for nspin=1/2). Used by the Veff Hellmann-Feynman
    // terms: V^H needs the total density (sum over spins), V^XC the spin-resolved densities.
    std::vector<const hamilt::HContainer<double>*> dmR;
    const Charge* chg = nullptr; // ground-state charge for XC Hellmann-Feynman (FDM)
#ifdef __EXX
    // gamma (TK==double) exx interfaces used by write_dH_exx; exactly one is set depending on
    // GlobalC::exx_info.info_ri.real_number (exd: real Hexx, exc: complex Hexx).
    Exx_LRI_Interface<double, double>* exd = nullptr;
    Exx_LRI_Interface<double, std::complex<double>>* exc = nullptr;
#endif
};

// Returns 0-based atom indices to output (converted from the 1-based user-facing values stored at param[2+]); 
// empty vector (param.size() <= 2) means all atoms.  Out-of-range checking is done in
// write_dh_perI where nat is available: indices >= nat are warned about and silently skipped.
inline std::vector<int> dh_atom_filter(const std::vector<int>& param)
{
    if (param.size() <= 2)  // param elements: [on/off][precition][iat1][iat2][...]
        return {};
    return std::vector<int>(param.begin() + 2, param.end());
}

// Shared writer for the per-atom-I dH terms. For every differentiated atom I it writes:
//   - dH(R) in CSR real-space format  ({rprefix}{x,y,z}_iat{I}...)
//   - dH(k) dense matrices            ({kprefix}{x,y,z}_iat{I}...) folded like H(k),
//     so they can be compared directly with the H(k) term matrices (*_nao.txt).
// g[d] are nat per-I HContainers for direction d=0..2 (already filled by an operator's cal_dH).
// atom_filter: if non-empty, only the listed 0-based atom indices are written; empty = all atoms.
void write_dh_perI(WriteDHParams& params,
    int ispin,
    const std::string& rprefix,
    const std::string& kprefix,
    const std::string& label,
    std::array<std::vector<hamilt::HContainer<double>*>, 3>& g,
    const std::vector<int>& atom_filter = {});

void write_dH_components(WriteDHParams& params);

bool write_dH_t(WriteDHParams& params);

bool write_dH_vnl(WriteDHParams& params);

bool write_dH_vl(WriteDHParams& params);

bool write_dH_vh(WriteDHParams& params);

bool write_dH_vh_pulay(WriteDHParams& params);

bool write_dH_vxc(WriteDHParams& params);

bool write_dH_vxc_pulay(WriteDHParams& params);

bool write_dH_sum(WriteDHParams& params);

#ifdef __EXX
bool write_dH_exx(WriteDHParams& params);
#endif

} // namespace ModuleIO

#endif
