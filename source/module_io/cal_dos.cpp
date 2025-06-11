#include "cal_dos.h"

#include "module_base/constants.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/parallel_reduce.h"
#include "module_parameter/parameter.h"

void ModuleIO::prepare_dos(std::ofstream& ofs_running,
		const elecstate::efermi &energy_fermi,
        const ModuleBase::matrix& ekb,
        const int nks,
        const int nbands,
		const double& dos_edelta_ev,
		const double& dos_scale,
		double &emax,
		double &emin)
{
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

    ofs_running << std::setprecision(6);

    assert(nbands>0);

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
    emax = ekb(0, 0);
    emin = ekb(0, 0);
    for (int ik = 0; ik < nks; ++ik)
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

}

bool ModuleIO::cal_dos(const int& is,  // index for spin
		const std::string& fn,   // file name for DOS
		const double& de_ev,           // delta energy in ev
		const double& emax_ev, // maximal energy in eV
		const double& emin_ev, // minimal energy in ev.
		const double& bcoeff,
		const int& nks, // number of k points in this pool
		const int& nkstot, // number of total kpoints
		const std::vector<double>& wk, // weight of k points
		const std::vector<int>& isk,   // index of spin for each k-point
		const int& nbands,             // number of bands
		const ModuleBase::matrix& ekb, // energy for each k point and each band
		const ModuleBase::matrix& wg   // weight of k-points and bands 
		)
{
    ModuleBase::TITLE("ModuleIO", "cal_dos");

    std::ofstream ofs_dos;
    std::ofstream ofs_smear;

    if (GlobalV::MY_RANK == 0)
    {
        ofs_dos.open(fn.c_str());
    }

    std::vector<double> dos;
    std::vector<double> ene;
    std::vector<double> sum_elec;
    std::vector<double> dos_smear; // dos_smearing
    dos.clear();
    ene.clear();
    sum_elec.clear();
    dos_smear.clear();

#ifdef __MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    if (de_ev <= 0)
    {
        ModuleBase::WARNING("ModuleIO::cal_dos", "de <= 0 ");
        return false;
    }
    else if (emax_ev < emin_ev)
    {
        ModuleBase::WARNING("ModuleIO::cal_dos", "emax_ev < emin_ev");
        return false;
    }

    const int npoints = static_cast<int>(std::floor((emax_ev - emin_ev) / de_ev))+1;

    if (npoints <= 0)
    {
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "npoints", npoints);
        ModuleBase::WARNING("ModuleIO::cal_dos", "npoints <= 0");
        return false;
    }
    if (GlobalV::MY_RANK == 0)
    {
        ofs_dos << npoints << " # number of points" << std::endl;
        ofs_dos << "#" << std::setw(14) << "energy" 
                 << std::setw(15) << "elec_states" 
                 << std::setw(15) << "sum_states" 
                 << std::setw(15) << "states_smear" 
                 << std::setw(15) << "sum_states" << std::endl;
    }

    std::vector<double> e_mod(npoints, 0.0); 

    double sum = 0.0;
    double curr_energy = emin_ev;
    double e_old = 0.0;

    while (curr_energy < emax_ev)
    {
        double nstates = 0.0;
        e_old = curr_energy;
        curr_energy += de_ev;

        // nks is the number of k-points in the 'pool'
        for (int ik = 0; ik < nks; ik++)
        {
            // spin index
            if (is == isk[ik])
            {
                // band index
                for (int ib = 0; ib < nbands; ib++)
                {
                    //  compare et and e_old(curr_energy) in ev unit.
                    if (ekb(ik, ib) * ModuleBase::Ry_to_eV >= e_old 
                     && ekb(ik, ib) * ModuleBase::Ry_to_eV < curr_energy)
                    {
                        nstates += wk[ik] * nkstot; 
                    }
                }
            }
        }

#ifdef __MPI
        const int npool = GlobalV::KPAR * PARAM.inp.bndpar;
        Parallel_Reduce::reduce_double_allpool(npool, GlobalV::NPROC_IN_POOL, nstates);
#endif

        nstates = nstates / static_cast<double>(nkstot);
        sum += nstates;
        if (GlobalV::MY_RANK == 0)
        {
            dos.push_back(nstates);
            ene.push_back(curr_energy);
            sum_elec.push_back(sum);
        }
    }

    // Use Gaussian smearing to smooth the DOS
    if (GlobalV::MY_RANK == 0)
    {
        dos_smear.resize(dos.size());

        double b = sqrt(2.0) * bcoeff;
        for (int i = 0; i < dos.size() ; i++)
        {
            double Gauss = 0.0;

            for (int j = 0; j < dos.size(); j++)
            {
                double denergy = ene[j] - ene[i];
                double de2 = denergy * denergy;
                Gauss = exp(-de2 / b / b) / sqrt(ModuleBase::PI) / b;
                dos_smear[j] += dos[i] * Gauss;
            }
        }

        // mohan add 2025-06-08
        const double dos_thr = 1.0e-12; 
        double sum2 = 0.0;

        for (int i = 0; i < dos.size(); i++)
        {
            if(dos_smear[i]<dos_thr)
            {
                 dos_smear[i]=0.0;
            }
            sum2 += dos_smear[i] * de_ev;

            ofs_dos << std::setw(15) << ene[i] 
                 << std::setw(15) << dos[i]
                 << std::setw(15) << sum_elec[i]
                 << std::setw(15) << dos_smear[i] 
                 << std::setw(15) << sum2 << std::endl;
        }
    }

    if (GlobalV::MY_RANK == 0)
    {
        ofs_dos.close();
    }

    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Number of bands", nbands);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Total electronic states from DOS", sum);

    return true;
}
