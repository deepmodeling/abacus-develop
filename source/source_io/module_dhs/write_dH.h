#ifndef WRITE_DH_H
#define WRITE_DH_H

#include "source_basis/module_nao/two_center_bundle.h"
#include "source_cell/klist.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_estate/module_pot/potential_new.h"
#include "source_lcao/LCAO_domain.h"
#include "source_lcao/module_hcontainer/hcontainer.h"

#include <vector>

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
    int nat = 0;
    int nspin = 1;
    int istep = 0;
    bool gamma_only = false;
    bool append = false;
    const hamilt::HContainer<double>* dmR = nullptr;
};

bool any_dh_term_enabled();

// Shared writer for the per-atom-I dH terms. For every differentiated atom I it writes:
//   - dH(R) in CSR real-space format  ({rprefix}{x,y,z}_iat{I}...)
//   - dH(k) dense matrices            ({kprefix}{x,y,z}_iat{I}...) folded like H(k),
//     so they can be compared directly with the H(k) term matrices (*_nao.txt).
// gx/gy/gz are nat per-I HContainers (already filled by an operator's cal_dH).
void write_dh_perI(WriteDHParams& params,
                   int ispin,
                   const std::string& rprefix,
                   const std::string& kprefix,
                   const std::string& label,
                   std::vector<hamilt::HContainer<double>*>& gx,
                   std::vector<hamilt::HContainer<double>*>& gy,
                   std::vector<hamilt::HContainer<double>*>& gz);

void write_dH_components(WriteDHParams& params);

bool write_dH_t(WriteDHParams& params);

bool write_dH_vnl(WriteDHParams& params);

bool write_dH_vl(WriteDHParams& params);

bool write_dH_vh(WriteDHParams& params);

bool write_dH_vxc(WriteDHParams& params);

bool write_dH_sum(WriteDHParams& params);

} // namespace ModuleIO

#endif
