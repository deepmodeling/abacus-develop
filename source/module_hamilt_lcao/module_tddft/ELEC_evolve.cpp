#include "ELEC_evolve.h"

#include "LCAO_evolve.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
#include "module_elecstate/module_charge/symmetry_rho.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"
#include "module_hamilt_lcao/module_dftu/dftu.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

ELEC_evolve::ELEC_evolve(){};
ELEC_evolve::~ELEC_evolve(){};

double ELEC_evolve::td_force_dt;
bool ELEC_evolve::td_vext;
std::vector<int> ELEC_evolve::td_vext_dire_case;
bool ELEC_evolve::out_dipole;
bool ELEC_evolve::out_efield;
double ELEC_evolve::td_print_eij; // the threshold to output Eij elements
int ELEC_evolve::td_edm; // 0: new edm method   1: old edm method

// this routine only serves for TDDFT using LCAO basis set
void ELEC_evolve::evolve_psi(const int& istep,
                             hamilt::Hamilt<double>* phm,
                             Local_Orbital_wfc& lowf,
                             psi::Psi<std::complex<double>>* psi,
                             psi::Psi<std::complex<double>>* psi_laststep,
                             std::complex<double>** Hk_laststep,
                             ModuleBase::matrix& ekb,
                             int htype,
                             int propagator)
{
    ModuleBase::TITLE("ELEC_evolve", "eveolve_psi");
    ModuleBase::timer::tick("ELEC_evolve", "evolve_psi");

    // pool parallization in future -- mohan note 2021-02-09
    for (int ik = 0; ik < GlobalC::kv.nks; ik++)
    {
        phm->updateHk(ik);

        ModuleBase::timer::tick("Efficience", "evolve_k");
        Evolve_LCAO_Matrix ELM(lowf.ParaV);
        psi->fix_k(ik);
        psi_laststep->fix_k(ik);
        if (Hk_laststep == nullptr)
        {
            ELM.evolve_complex_matrix(ik, phm, psi, psi_laststep, nullptr, &(ekb(ik, 0)), htype, propagator);
        }
        else
        {
            ELM.evolve_complex_matrix(ik, phm, psi, psi_laststep, Hk_laststep[ik], &(ekb(ik, 0)), htype, propagator);
        }
        ModuleBase::timer::tick("Efficience", "evolve_k");
    } // end k

    ModuleBase::timer::tick("ELEC_evolve", "evolve_psi");
    return;
}
