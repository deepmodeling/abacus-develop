#include "cal_dos.h"

#include "source_base/constants.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/parallel_reduce.h"

#include <cmath>

void ModuleIO::prepare_dos(std::ofstream& ofs_running,
        const elecstate::Efermi &energy_fermi,
        const ModuleBase::matrix& ekb,
        const int nks,
        const int nbands,
        const double& dos_edelta_ev,
        const double& dos_scale,
        double &emax,
        double &emin,
        const bool dos_setemax,
        const double dos_emax_ev,
        const bool dos_setemin,
        const double dos_emin_ev,
        const bool two_fermi)
{
    ofs_running << "\n >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    ofs_running << " |                                                                    |" << std::endl;
    ofs_running << " |               #Calcualte Density of States (DOS)#                  |" << std::endl;                 
    ofs_running << " | DOS stands for Density of States. It represents the number of      |" << std::endl;
    ofs_running << " | available electronic states per unit energy range.                 |" << std::endl;
    ofs_running << " | By analyzing the DOS, we can gain insights into how electrons are  |" << std::endl;
    ofs_running << " | distributed among different energy levels within the material.     |" << std::endl;
    ofs_running << " |                                                                    |" << std::endl;
    ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;

    ofs_running << std::setprecision(6);

    assert(nbands>0);

    if (two_fermi == false)
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
    Parallel_Reduce::reduce_max(emax);
    Parallel_Reduce::reduce_min(emin);
#endif

    emax *= ModuleBase::Ry_to_eV;
    emin *= ModuleBase::Ry_to_eV;

    if (dos_setemax)
    {
        emax = dos_emax_ev;
    }
    if (dos_setemin)
    {
        emin = dos_emin_ev;
    }

    if (!dos_setemax && !dos_setemin)
    {
        // scale up a little bit so the end peaks are displaced better
        double delta = (emax - emin) * dos_scale;
        emax = emax + delta / 2.0;
        emin = emin - delta / 2.0;
    }

    assert(dos_edelta_ev>0.0);

    ModuleBase::GlobalFunc::OUT(ofs_running, "Minimal energy (eV)", emin);
    ModuleBase::GlobalFunc::OUT(ofs_running, "Maximal energy (eV)", emax);
    ModuleBase::GlobalFunc::OUT(ofs_running, "Energy interval (eV)", dos_edelta_ev);

}

bool ModuleIO::cal_dos(const int& is,  // index for spin
        const std::string& fn,   // file name for DOS
        const double& de_ev,   // delta energy in ev
        const double& emax_ev, // maximal energy in eV
        const double& emin_ev, // minimal energy in ev.
        const double& bcoeff,
        const int& nks, // number of k points in this pool
        const int& nkstot, // number of total kpoints
        const std::vector<double>& wk, // weight of k points
        const std::vector<int>& isk,   // index of spin for each k-point
        const int& nbands,             // number of bands
        const ModuleBase::matrix& ekb, // energy for each k point and each band
        const ModuleBase::matrix& wg,  // weight of k-points and bands
        const int istep, // ionic step
        const bool out_app_flag,
        const int bndpar)
{
    ModuleBase::TITLE("ModuleIO", "cal_dos");

    std::ofstream ofs_dos;

    if (GlobalV::MY_RANK == 0)
    {
        if(out_app_flag==true)
        {
            ofs_dos.open(fn.c_str(), std::ios::app);
        }
        else
        {
            ofs_dos.open(fn.c_str());
        }
        ofs_dos << istep+1 << "    # ionic step" << std::endl;
    }

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
        ofs_dos << std::fixed << std::setprecision(6);
        ofs_dos << "#" << std::setw(12) << "energy(eV)"
                 << std::setw(12) << "dos"
                 << std::setw(12) << "dos_int"
                 << std::setw(17) << "dos_smear(1/eV)"
                 << std::setw(15) << "dos_smear_int" << std::endl;
    }

    // Node-sampled grid shared with the PDOS text output: en = emin + n*de.
    // Each node carries the histogram count of bin [en, en+de) and the
    // Gaussian-smeared DOS evaluated directly at en from the exact
    // eigenvalue positions, using the same normalized-Gaussian convention
    // as cal_pdos (width parameter bcoeff), so that the sum of all PDOS
    // orbitals equals this smeared column point by point.
    std::vector<double> dos(npoints, 0.0);
    std::vector<double> ene(npoints, 0.0);
    std::vector<double> dos_smear(npoints, 0.0);
    std::vector<double> sum_elec(npoints, 0.0);

    const double b = std::sqrt(2.0) * bcoeff;
    const double gauss_norm = std::sqrt(ModuleBase::PI) * b;

    // en accumulates by de_ev so that adjacent bins [en, en+de) share the
    // exact same floating-point boundary and no eigenvalue is counted twice
    // or skipped; the printed labels still coincide with the PDOS nodes
    // emin + n*de at the 1e-6 eV output precision.
    double en = emin_ev;
    for (int n = 0; n < npoints; ++n)
    {
        ene[n] = en;
        const double bin_hi = en + de_ev;

        // nks is the number of k-points in the 'pool'
        for (int ik = 0; ik < nks; ik++)
        {
            // spin index
            if (is == isk[ik])
            {
                // band index
                for (int ib = 0; ib < nbands; ib++)
                {
                    const double e = ekb(ik, ib) * ModuleBase::Ry_to_eV;
                    if (e >= en && e < bin_hi)
                    {
                        dos[n] += wk[ik] * nkstot;
                    }
                    const double de = en - e;
                    dos_smear[n] += wk[ik] * std::exp(-de * de / b / b) / gauss_norm;
                }
            }
        }
        en = bin_hi;
    }

#ifdef __MPI
    const int npool = GlobalV::KPAR * bndpar;
    Parallel_Reduce::reduce_double_allpool(npool, GlobalV::NPROC_IN_POOL, dos.data(), npoints);
    Parallel_Reduce::reduce_double_allpool(npool, GlobalV::NPROC_IN_POOL, dos_smear.data(), npoints);
#endif

    double sum = 0.0;
    double sum2 = 0.0;
    const double dos_thr = 1.0e-12;

    if (GlobalV::MY_RANK == 0)
    {
        for (int n = 0; n < npoints; ++n)
        {
            dos[n] /= static_cast<double>(nkstot);
            sum += dos[n];
            sum_elec[n] = sum;

            if (dos_smear[n] < dos_thr)
            {
                dos_smear[n] = 0.0;
            }
            sum2 += dos_smear[n] * de_ev;

            ofs_dos << std::setw(13) << ene[n]
                 << " " << std::setw(12) << dos[n]
                 << " " << std::setw(12) << sum_elec[n]
                 << " " << std::setw(14) << dos_smear[n]
                 << " " << std::setw(14) << sum2 << std::endl;
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
