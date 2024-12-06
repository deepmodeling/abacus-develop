#include "cal_nelec_nband.h"
#include "module_base/constants.h"
#include "module_parameter/parameter.h"
#ifdef USE_PAW
#include "module_cell/module_paw/paw_cell.h"
#endif

namespace elecstate {

void cal_nelec(const Atom* atoms, const int& ntype, double& nelec)
{
    ModuleBase::TITLE("UnitCell", "cal_nelec");
    GlobalV::ofs_running << "\n SETUP THE ELECTRONS NUMBER" << std::endl;

    if (nelec == 0)
    {
        if (PARAM.inp.use_paw)
        {
#ifdef USE_PAW
            for (int it = 0; it < ntype; it++)
            {
                std::stringstream ss1, ss2;
                ss1 << " electron number of element " << GlobalC::paw_cell.get_zat(it) << std::endl;
                const int nelec_it = GlobalC::paw_cell.get_val(it) * atoms[it].na;
                nelec += nelec_it;
                ss2 << "total electron number of element " << GlobalC::paw_cell.get_zat(it);

                ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, ss1.str(), GlobalC::paw_cell.get_val(it));
                ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, ss2.str(), nelec_it);
            }
#endif
        }
        else
        {
            for (int it = 0; it < ntype; it++)
            {
                std::stringstream ss1, ss2;
                ss1 << "electron number of element " << atoms[it].label;
                const double nelec_it = atoms[it].ncpp.zv * atoms[it].na;
                nelec += nelec_it;
                ss2 << "total electron number of element " << atoms[it].label;

                ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, ss1.str(), atoms[it].ncpp.zv);
                ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, ss2.str(), nelec_it);
            }
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "AUTOSET number of electrons: ", nelec);
        }
    }
    if (PARAM.inp.nelec_delta != 0)
    {
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,
                                    "nelec_delta is NOT zero, please make sure you know what you are "
                                    "doing! nelec_delta: ",
                                    PARAM.inp.nelec_delta);
        nelec += PARAM.inp.nelec_delta;
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "nelec now: ", nelec);
    }
    return;
}

void cal_nbands(const int& nelec, const int& nlocal, const std::vector<double>& nelec_spin, int& nbands)
{
    if (PARAM.inp.esolver_type == "sdft") // qianrui 2021-2-20
    {
        return;
    }
    //=======================================
    // calculate number of bands (setup.f90)
    //=======================================
    double occupied_bands = static_cast<double>(nelec / ModuleBase::DEGSPIN);
    if (PARAM.inp.lspinorb == 1) {
        occupied_bands = static_cast<double>(nelec);
    }

    if ((occupied_bands - std::floor(occupied_bands)) > 0.0)
    {
        occupied_bands = std::floor(occupied_bands) + 1.0; // mohan fix 2012-04-16
    }

    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "occupied bands", occupied_bands);

    if (nbands == 0)
    {
        if (PARAM.inp.nspin == 1)
        {
            const int nbands1 = static_cast<int>(occupied_bands) + 10;
            const int nbands2 = static_cast<int>(1.2 * occupied_bands) + 1;
            nbands = std::max(nbands1, nbands2);
            if (PARAM.inp.basis_type != "pw") {
                nbands = std::min(nbands, nlocal);
            }
        }
        else if (PARAM.inp.nspin == 4)
        {
            const int nbands3 = nelec + 20;
            const int nbands4 = static_cast<int>(1.2 * nelec) + 1;
            nbands = std::max(nbands3, nbands4);
            if (PARAM.inp.basis_type != "pw") {
                nbands = std::min(nbands, nlocal);
            }
        }
        else if (PARAM.inp.nspin == 2)
        {
            const double max_occ = std::max(nelec_spin[0], nelec_spin[1]);
            const int nbands3 = static_cast<int>(max_occ) + 11;
            const int nbands4 = static_cast<int>(1.2 * max_occ) + 1;
            nbands = std::max(nbands3, nbands4);
            if (PARAM.inp.basis_type != "pw") {
                nbands = std::min(nbands, nlocal);
            }
        }
        ModuleBase::GlobalFunc::AUTO_SET("NBANDS", nbands);
    }
    // else if ( PARAM.inp.calculation=="scf" || PARAM.inp.calculation=="md" || PARAM.inp.calculation=="relax") //pengfei
    // 2014-10-13
    else
    {
        if (nbands < occupied_bands) {
            ModuleBase::WARNING_QUIT("unitcell", "Too few bands!");
        }
        if (PARAM.inp.nspin == 2)
        {
            if (nbands < nelec_spin[0])
            {
                ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "nelec_up", nelec_spin[0]);
                ModuleBase::WARNING_QUIT("ElecState::cal_nbands", "Too few spin up bands!");
            }
            if (nbands < nelec_spin[1])
            {
                ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "nelec_down", nelec_spin[1]);
                ModuleBase::WARNING_QUIT("ElecState::cal_nbands", "Too few spin down bands!");
            }
        }
    }

    // mohan add 2010-09-04
    // std::cout << "nbands(this-> = " <<nbands <<std::endl;
    if (nbands == occupied_bands)
    {
        if (PARAM.inp.smearing_method != "fixed")
        {
            ModuleBase::WARNING_QUIT("ElecState::cal_nbands", "for smearing, num. of bands > num. of occupied bands");
        }
    }

    // mohan update 2021-02-19
    // mohan add 2011-01-5
    if (PARAM.inp.basis_type == "lcao" || PARAM.inp.basis_type == "lcao_in_pw")
    {
        if (nbands > nlocal)
        {
            ModuleBase::WARNING_QUIT("ElecState::cal_nbandsc", "NLOCAL < NBANDS");
        }
        else
        {
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "NLOCAL", nlocal);
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "NBANDS", nbands);
        }
    }

    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "NBANDS", nbands);
}

}