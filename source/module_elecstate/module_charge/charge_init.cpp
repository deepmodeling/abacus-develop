#include <vector>

#include "charge.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/libm/libm.h"
#include "module_base/math_integral.h"
#include "module_base/math_sphbes.h"
#include "module_base/memory.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
#include "module_base/tool_threading.h"
#include "module_elecstate/magnetism.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_hamilt_pw/hamilt_pwdft/parallel_grid.h"
#include "module_io/rho_io.h"

void Charge::init_rho(elecstate::efermi& eferm_iout, const ModuleBase::ComplexMatrix& strucFac)
{
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "init_chg", GlobalV::init_chg);

    std::cout << " START CHARGE      : " << GlobalV::init_chg << std::endl;
    if (GlobalV::init_chg == "atomic") // mohan add 2007-10-17
    {
        this->atomic_rho(GlobalV::NSPIN, GlobalC::ucell.omega, rho, strucFac);
    }
    else if (GlobalV::init_chg == "file")
    {
        GlobalV::ofs_running << " try to read charge from file : ";
        for (int is = 0; is < GlobalV::NSPIN; ++is)
        {
            std::stringstream ssc;
            ssc << GlobalV::global_readin_dir << "SPIN" << is + 1 << "_CHG.cube";
            GlobalV::ofs_running << ssc.str() << std::endl;
            double& ef_tmp = eferm_iout.get_ef(is);
            if (ModuleIO::read_rho(
#ifdef __MPI
                    &(GlobalC::Pgrid),
#endif
                    is,
                    GlobalV::NSPIN,
                    ssc.str(),
                    this->rho[is],
                    this->rhopw->nx,
                    this->rhopw->ny,
                    this->rhopw->nz,
                    ef_tmp,
                    &(GlobalC::ucell),
                    this->prenspin))
            {
                GlobalV::ofs_running << " Read in the charge density: " << ssc.str() << std::endl;
            }
            else if (is > 0)
            {
                if (prenspin == 1)
                {
                    GlobalV::ofs_running << " Didn't read in the charge density but autoset it for spin " << is + 1
                                         << std::endl;
                    for (int ir = 0; ir < this->rhopw->nrxx; ir++)
                    {
                        this->rho[is][ir] = 0.0;
                    }
                }
                //
                else if (prenspin == 2)
                { // read up and down , then rearrange them.
                    if (is == 1)
                    {
                        ModuleBase::WARNING_QUIT("Charge::init_rho", "Incomplete charge density file!");
                    }
                    else if (is == 2)
                    {
                        GlobalV::ofs_running << " Didn't read in the charge density but would rearrange it later. "
                                             << std::endl;
                    }
                    else if (is == 3)
                    {
                        GlobalV::ofs_running << " rearrange charge density " << std::endl;
                        for (int ir = 0; ir < this->rhopw->nrxx; ir++)
                        {
                            this->rho[3][ir] = this->rho[0][ir] - this->rho[1][ir];
                            this->rho[0][ir] = this->rho[0][ir] + this->rho[1][ir];
                            this->rho[1][ir] = 0.0;
                            this->rho[2][ir] = 0.0;
                        }
                    }
                }
            }
            else
            {
                ModuleBase::WARNING_QUIT(
                    "init_rho",
                    "!!! Couldn't find the charge file !!! The default directory \n of SPIN1_CHG.cube is OUT.suffix, "
                    "or you must set read_file_dir \n to a specific directory. ");
            }
        }

        if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
        {
            for (int is = 0; is < GlobalV::NSPIN; is++)
            {
                std::stringstream ssc;
                ssc << GlobalV::global_readin_dir << "SPIN" << is + 1 << "_TAU.cube";
                GlobalV::ofs_running << " try to read kinetic energy density from file : " << ssc.str() << std::endl;
                // mohan update 2012-02-10, sunliang update 2023-03-09
                if (ModuleIO::read_rho(
#ifdef __MPI
                        &(GlobalC::Pgrid),
#endif
                        is,
                        GlobalV::NSPIN,
                        ssc.str(),
                        this->kin_r[is],
                        this->rhopw->nx,
                        this->rhopw->ny,
                        this->rhopw->nz,
                        eferm_iout.ef,
                        &(GlobalC::ucell),
                        this->prenspin))
                {
                    GlobalV::ofs_running << " Read in the kinetic energy density: " << ssc.str() << std::endl;
                }
            }
        }
    }
    else
    {
        ModuleBase::WARNING_QUIT("Charge::init_rho", "init_chg is wrong!");
    }

    // Peize Lin add 2020.04.04
    if (GlobalC::restart.info_load.load_charge && !GlobalC::restart.info_load.load_charge_finish)
    {
        for (int is = 0; is < GlobalV::NSPIN; ++is)
        {
            GlobalC::restart.load_disk("charge", is, this->nrxx, rho);
        }
        GlobalC::restart.info_load.load_charge_finish = true;
    }
}
