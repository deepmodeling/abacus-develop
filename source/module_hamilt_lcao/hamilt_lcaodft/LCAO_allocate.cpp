#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_domain.h"
#include "module_base/timer.h"
#include "module_parameter/parameter.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_hamilt_lcao/module_deepks/LCAO_deepks.h"

namespace LCAO_domain
{

void divide_HS_in_frag(const bool isGamma, const UnitCell& ucell, Parallel_Orbitals& pv,const int& nks, const LCAO_Orbitals& orb) {
    ModuleBase::TITLE("LCAO_domain", "divide_HS_in_frag");

    //(1), (2): set up matrix division have been moved into ESolver_KS_LCAO::init_basis_lcao
    // just pass `ParaV` as pointer is enough
#ifdef __DEEPKS
    // wenfei 2021-12-19
    // preparation for DeePKS

    if (PARAM.inp.deepks_out_labels || PARAM.inp.deepks_scf) {
        // allocate relevant data structures for calculating descriptors
        std::vector<int> na;
        na.resize(ucell.ntype);
        for (int it = 0; it < ucell.ntype; it++) {
            na[it] = ucell.atoms[it].na;
        }

        GlobalC::ld.init(orb,
                         ucell.nat,
                         ucell.ntype,
                         pv,
                         na);

        if (PARAM.inp.deepks_scf) {
            if (isGamma) {
                GlobalC::ld.allocate_V_delta(ucell.nat);
            } else {
                GlobalC::ld.allocate_V_delta(ucell.nat, nks);
            }
        }
    }
#endif
    return;
}

}
