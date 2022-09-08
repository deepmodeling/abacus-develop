#include "ELEC_evolve.h"

#include "../module_base/timer.h"
#include "../src_parallel/parallel_reduce.h"
#include "../src_pw/global.h"
#include "../src_pw/symmetry_rho.h"
#include "LCAO_diago.h"
#include "LCAO_evolve.h"
#include "dftu.h"
#include "module_hamilt/hamilt_lcao.h"

ELEC_evolve::ELEC_evolve(){};
ELEC_evolve::~ELEC_evolve(){};

int ELEC_evolve::tddft;
double ELEC_evolve::td_scf_thr;
double ELEC_evolve::td_dt;
double ELEC_evolve::td_force_dt;
int ELEC_evolve::td_val_elec_01;
int ELEC_evolve::td_val_elec_02;
int ELEC_evolve::td_val_elec_03;
int ELEC_evolve::td_vext;
int ELEC_evolve::td_vext_dire;
double ELEC_evolve::td_timescale;
int ELEC_evolve::td_vexttype;
int ELEC_evolve::td_vextout;
int ELEC_evolve::td_dipoleout;

// this routine only serves for TDDFT using LCAO basis set
void ELEC_evolve::evolve_psi(const int& istep,
                             hamilt::Hamilt* phm,
                             Local_Orbital_wfc& lowf,
                             psi::Psi<std::complex<double>>* psi,
                             psi::Psi<std::complex<double>>* psi_laststep)
{
    ModuleBase::TITLE("ELEC_evolve", "eveolve_psi");
    ModuleBase::timer::tick("ELEC_evolve", "evolve_psi");

    phm->constructHamilt();

    // pool parallization in future -- mohan note 2021-02-09
    for (int ik = 0; ik < GlobalC::kv.nks; ik++)
    {
        phm->updateHk(ik);

        ModuleBase::timer::tick("Efficience", "evolve_k");
        Evolve_LCAO_Matrix ELM(lowf.ParaV);
        psi->fix_k(ik);
        ELM.evolve_complex_matrix(ik, phm, lowf, psi, psi_laststep, GlobalC::wf.ekb[ik]);
        ModuleBase::timer::tick("Efficience", "evolve_k");
    } // end k

    ModuleBase::timer::tick("ELEC_evolve", "evolve_psi");
    return;
}
