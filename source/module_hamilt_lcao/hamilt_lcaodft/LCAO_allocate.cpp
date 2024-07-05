#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_domain.h"
#include "module_base/timer.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_hamilt_lcao/module_deepks/LCAO_deepks.h"

namespace LCAO_domain
{

void divide_HS_in_frag(LCAO_Matrix& lm, const bool isGamma, Parallel_Orbitals& pv,const int& nks) {
    ModuleBase::TITLE("LCAO_Matrix", "divide_HS_in_frag");

    //(1), (2): set up matrix division have been moved into ORB_control
    // just pass `ParaV` as pointer is enough
    lm.ParaV = &pv;
    // (3) allocate for S, H_fixed, H, and S_diag
    if (isGamma) {
        LCAO_domain::allocate_HS_gamma(lm, lm.ParaV->nloc);
    } else {
        LCAO_domain::allocate_HS_k(lm, lm.ParaV->nloc);
    }
#ifdef __DEEPKS
    // wenfei 2021-12-19
    // preparation for DeePKS

    if (GlobalV::deepks_out_labels || GlobalV::deepks_scf) {
        // allocate relevant data structures for calculating descriptors
        std::vector<int> na;
        na.resize(GlobalC::ucell.ntype);
        for (int it = 0; it < GlobalC::ucell.ntype; it++) {
            na[it] = GlobalC::ucell.atoms[it].na;
        }

        GlobalC::ld.init(GlobalC::ORB,
                         GlobalC::ucell.nat,
                         GlobalC::ucell.ntype,
                         pv,
                         na);

        if (GlobalV::deepks_scf) {
            if (isGamma) {
                GlobalC::ld.allocate_V_delta(GlobalC::ucell.nat);
            } else {
                GlobalC::ld.allocate_V_delta(GlobalC::ucell.nat, nks);
            }
        }
    }
#endif
    return;
}

void allocate_HS_gamma(LCAO_Matrix& lm, const long& nloc) {
    ModuleBase::TITLE("LCAO_Matrix", "allocate_HS_gamma");

    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "nloc", nloc);

    if (nloc == 0) {
        return;
    }

    // because we initilize in the constructor function
    // with dimension '1', so here we reconstruct these
    // matrices

    lm.Sloc.resize(nloc);
    lm.Hloc_fixed.resize(nloc);
    lm.Hloc.resize(nloc);

    ModuleBase::GlobalFunc::ZEROS(lm.Sloc.data(), nloc);
    ModuleBase::GlobalFunc::ZEROS(lm.Hloc_fixed.data(), nloc);
    ModuleBase::GlobalFunc::ZEROS(lm.Hloc.data(), nloc);

    return;
}

void allocate_HS_k(LCAO_Matrix& lm, const long& nloc) {
    ModuleBase::TITLE("LCAO_Matrix", "allocate_HS_k");

    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "nloc", nloc);

    if (nloc == 0) {
        return; // mohan fix bug 2012-05-25
    }

    // because we initilize in the constructor function
    // with dimension '1', so here we reconstruct these
    // matrices
    lm.Sloc2.resize(nloc);
    lm.Hloc_fixed2.resize(nloc);
    lm.Hloc2.resize(nloc);

    ModuleBase::GlobalFunc::ZEROS(lm.Sloc2.data(), nloc);
    ModuleBase::GlobalFunc::ZEROS(lm.Hloc_fixed2.data(), nloc);
    ModuleBase::GlobalFunc::ZEROS(lm.Hloc2.data(), nloc);

    return;
}

}
