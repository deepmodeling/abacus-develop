#include "esolver_ks.h"
#include <iostream>

#include "module_io/print_info.h"
#include "module_parameter/parameter.h"
#include "module_base/global_variable.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

namespace Print_functions
{

void print_wfcfft(
		const Input_para& inp, 
        ModulePW::PW_Basis_K& pw_wfc,
		std::ofstream& ofs)
{
    ofs << "\n\n\n\n";
    ofs << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
           ">>>>"
        << std::endl;
    ofs << " |                                                                 "
           "   |"
        << std::endl;
    ofs << " | Setup plane waves of wave functions:                            "
           "   |"
        << std::endl;
    ofs << " | Use the energy cutoff and the lattice vectors to generate the   "
           "   |"
        << std::endl;
    ofs << " | dimensions of FFT grid. The number of FFT grid on each "
           "processor   |"
        << std::endl;
    ofs << " | is 'nrxx'. The number of plane wave basis in reciprocal space "
           "is   |"
        << std::endl;
    ofs << " | different for charege/potential and wave functions. We also set "
           "   |"
        << std::endl;
    ofs << " | the 'sticks' for the parallel of FFT. The number of plane wave "
           "of  |"
        << std::endl;
    ofs << " | each k-point is 'npwk[ik]' in each processor                    "
           "   |"
        << std::endl;
    ofs << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
           "<<<<"
        << std::endl;
    ofs << "\n\n\n\n";
    ofs << "\n SETUP PLANE WAVES FOR WAVE FUNCTIONS" << std::endl;

    double ecut = inp.ecutwfc;
    if (std::abs(ecut - pw_wfc.gk_ecut * pw_wfc.tpiba2) > 1e-6)
    {
        ecut = pw_wfc.gk_ecut * pw_wfc.tpiba2;
        ofs << "Energy cutoff for wavefunc is incompatible with nx, ny, nz and "
               "it will be reduced!"
            << std::endl;
    }
    ModuleBase::GlobalFunc::OUT(ofs, "energy cutoff for wavefunc (unit:Ry)", ecut);
    ModuleBase::GlobalFunc::OUT(ofs,
                                "fft grid for wave functions",
                                pw_wfc.nx,
                                pw_wfc.ny,
                                pw_wfc.nz);
    ModuleBase::GlobalFunc::OUT(ofs, "number of plane waves", pw_wfc.npwtot);
    ModuleBase::GlobalFunc::OUT(ofs, "number of sticks", pw_wfc.nstot);

    ofs << "\n PARALLEL PW FOR WAVE FUNCTIONS" << std::endl;
    ofs << " " << std::setw(8) << "PROC" << std::setw(15) << "COLUMNS(POT)" << std::setw(15) << "PW" << std::endl;

    for (int i = 0; i < GlobalV::NPROC_IN_POOL; ++i)
    {
        ofs << " " << std::setw(8) << i + 1 << std::setw(15) << pw_wfc.nst_per[i] << std::setw(15)
            << pw_wfc.npw_per[i] << std::endl;
    }

    ofs << " --------------- sum -------------------" << std::endl;
    ofs << " " << std::setw(8) << GlobalV::NPROC_IN_POOL << std::setw(15) << pw_wfc.nstot << std::setw(15)
        << pw_wfc.npwtot << std::endl;
    ModuleBase::GlobalFunc::DONE(ofs, "INIT PLANEWAVE");
}

}
