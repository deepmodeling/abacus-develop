#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "global.h"
#include "electrons.h"
#include "../src_pw/symmetry_rho.h"
#include "../src_io/print_info.h"
#include "../src_io/wf_io.h"
#include "../src_io/chi0_hilbert.h"
#include "../src_io/chi0_standard.h"
#include "../src_io/epsilon0_pwscf.h"
#include "../src_io/epsilon0_vasp.h"
#include "../src_io/to_wannier90.h"
#include "../src_io/berryphase.h"
// new
#include "H_Ewald_pw.h"

double Electrons::avg_iter = 0;

Electrons::Electrons()
{
    iter = 0;
    test = 0;
    delta_total_energy = 0.0;
}

Electrons::~Electrons()
{
}

bool Electrons::check_stop_now(void)
{
    bool check_stop_now = false;

    if (check_stop_now)
    {
        conv_elec = false;
    }

    return check_stop_now;
} // END FUNCTION check_stop_now


void Electrons::init_mixstep_final_scf(void)
{
    ModuleBase::TITLE("electrons","init_mixstep_final_scf");

    GlobalC::CHR.irstep=0;
    GlobalC::CHR.idstep=0;
    GlobalC::CHR.totstep=0;

    return;
}
