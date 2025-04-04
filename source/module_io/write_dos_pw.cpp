#include "write_dos_pw.h"

#include "cal_dos.h"
#include "module_base/parallel_reduce.h"
#include "module_parameter/parameter.h"

void ModuleIO::write_dos_pw(const ModuleBase::matrix& ekb,
                            const ModuleBase::matrix& wg,
                            const K_Vectors& kv,
                            const int nbands,
                            const elecstate::efermi &energy_fermi,
                            const double& dos_edelta_ev,
                            const double& dos_scale,
                            const double& bcoeff,
                            std::ofstream& ofs_running)
{
    ModuleBase::TITLE("ModuleIO", "write_dos_pw");

    ofs_running << " DOS CALCULATIONS BEGINS" << std::endl; 

    ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
                   ">>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    ofs_running << " |                                            "
                   "                        |" << std::endl;
    ofs_running << " | DOS stands for Density of States. It represents the number of      |" << std::endl;
    ofs_running << " | available electronic states per unit energy range.                 |" << std::endl;
    ofs_running << " | By analyzing the DOS, we can gain insights into how electrons are  |" << std::endl;
    ofs_running << " | distributed among different energy levels within the material.     |" << std::endl;
    ofs_running << " |                                            "
                   "                        |" << std::endl;
    ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
                   ">>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;

    assert(nbands>0);

    ofs_running << std::setprecision(6);

    const int nspin0 = (PARAM.inp.nspin == 2) ? 2 : 1;

	if (PARAM.globalv.two_fermi == false)
	{
        ModuleBase::GlobalFunc::OUT(ofs_running, "Fermi energy (eV)", 
        energy_fermi.ef * ModuleBase::Ry_to_eV);
	}
	else
	{

        ModuleBase::GlobalFunc::OUT(ofs_running, "Spin up, Fermi energy (Ry)", 
        energy_fermi.ef_up * ModuleBase::Ry_to_eV);
        ModuleBase::GlobalFunc::OUT(ofs_running, "Spin dw, Fermi energy (Ry)", 
        energy_fermi.ef_dw * ModuleBase::Ry_to_eV);
	}

    // find energy range
    double emax = ekb(0, 0);
    double emin = ekb(0, 0);
    for (int ik = 0; ik < kv.get_nks(); ++ik)
    {
        for (int ib = 0; ib < nbands; ++ib)
        {
            emax = std::max(emax, ekb(ik, ib));
            emin = std::min(emin, ekb(ik, ib));
        }
    }

#ifdef __MPI
    Parallel_Reduce::gather_max_double_all(GlobalV::NPROC, emax);
    Parallel_Reduce::gather_min_double_all(GlobalV::NPROC, emin);
#endif

    emax *= ModuleBase::Ry_to_eV;
    emin *= ModuleBase::Ry_to_eV;

	if (PARAM.globalv.dos_setemax)
	{
		emax = PARAM.inp.dos_emax_ev;
	}
	if (PARAM.globalv.dos_setemin)
	{
		emin = PARAM.inp.dos_emin_ev;
	}

    if (!PARAM.globalv.dos_setemax && !PARAM.globalv.dos_setemin)
    {
        // scale up a little bit so the end peaks are displaced better
        double delta = (emax - emin) * dos_scale;
        emax = emax + delta / 2.0;
        emin = emin - delta / 2.0;
    }

    assert(dos_edelta_ev>0.0);

    ModuleBase::GlobalFunc::OUT(ofs_running, "Minimal energy is (eV)", emin);
    ModuleBase::GlobalFunc::OUT(ofs_running, "Maximal energy is (eV)", emax);
    ModuleBase::GlobalFunc::OUT(ofs_running, "Energy interval (eV)", dos_edelta_ev);


    for (int is = 0; is < nspin0; ++is)
    {
        // DOS_ispin contains not smoothed dos
        std::stringstream ss;
        ss << PARAM.globalv.global_out_dir << "DOS" << is + 1 << ".dat";

        std::stringstream ss1;
        ss1 << PARAM.globalv.global_out_dir << "DOS" << is + 1 << "_smear.dat";

        ModuleBase::GlobalFunc::OUT(ofs_running, "DOS file", ss.str());

		ModuleIO::cal_dos(is,
				ss.str(),
				ss1.str(),
				dos_edelta_ev,
				emax,
				emin,
				bcoeff,
				kv.get_nks(),
				kv.get_nkstot(),
				kv.wk,
				kv.isk,
				nbands,
				ekb,
				wg);
	}

    ofs_running << " DOS CALCULATIONS ENDS." << std::endl; 

}
