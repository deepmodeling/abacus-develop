#include "evolve_elec.h"

#include "evolve_psi.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
#include "module_elecstate/module_charge/symmetry_rho.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"
#include "module_hamilt_lcao/module_dftu/dftu.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

namespace module_tddft
{
Evolve_elec::Evolve_elec(){};
Evolve_elec::~Evolve_elec(){};

double Evolve_elec::td_force_dt;
bool Evolve_elec::td_vext;
std::vector<int> Evolve_elec::td_vext_dire_case;
bool Evolve_elec::out_dipole;
bool Evolve_elec::out_efield;
bool Evolve_elec::out_current;
double Evolve_elec::td_print_eij; // the threshold to output Eij elements
int Evolve_elec::td_edm;          // 0: new edm method   1: old edm method

// this routine only serves for TDDFT using LCAO basis set
void Evolve_elec::solve_psi(const int& istep,
                            const int nband,
                            const int nlocal,
                            hamilt::Hamilt<std::complex<double>>* phm,
                            Local_Orbital_wfc& lowf,
                            psi::Psi<std::complex<double>>* psi,
                            psi::Psi<std::complex<double>>* psi_laststep,
                            std::complex<double>** Hk_laststep,
                            std::complex<double>** Sk_laststep,
                            ModuleBase::matrix& ekb,
                            int htype,
                            int propagator,
                            const int& nks)
{
    ModuleBase::TITLE("Evolve_elec", "eveolve_psi");
    ModuleBase::timer::tick("Evolve_elec", "evolve_psi");

    for (int ik = 0; ik < nks; ik++)
    {
        phm->updateHk(ik);

        ModuleBase::timer::tick("Efficience", "evolve_k");
        psi->fix_k(ik);
        psi_laststep->fix_k(ik);
        if (htype == 0)
        {
            evolve_psi(nband,
                       nlocal,
                       lowf.ParaV,
                       phm,
                       psi[0].get_pointer(),
                       psi_laststep[0].get_pointer(),
                       nullptr,
                       nullptr,
                       &(ekb(ik, 0)),
                       htype,
                       propagator);
        }
        else if (htype == 1)
        {
            evolve_psi(nband,
                       nlocal,
                       lowf.ParaV,
                       phm,
                       psi[0].get_pointer(),
                       psi_laststep[0].get_pointer(),
                       Hk_laststep[ik],
                       Sk_laststep[ik],
                       &(ekb(ik, 0)),
                       htype,
                       propagator);
        }
        else
        {
            std::cout << "method of htype is wrong" << std::endl;
        }

        ModuleBase::timer::tick("Efficience", "evolve_k");
    } // end k

    ModuleBase::timer::tick("Evolve_elec", "evolve_psi");
    return;
}
} // namespace module_tddft