#include "cal_dos.h"

#include "module_base/constants.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/parallel_reduce.h"
#include "module_parameter/parameter.h"

bool ModuleIO::cal_dos(const int& is,  // index for spin
		const std::string& file_dos,   // file address for DOS
		const std::string& file_smear, // file address for DOS_smearing
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
        ofs_dos.open(file_dos.c_str());
        ofs_smear.open(file_smear.c_str());
    }

    std::vector<double> dos;
    std::vector<double> ene;
    std::vector<double> dos_smear; // dos_smearing
    dos.clear();
    ene.clear();
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

    const int npoints = static_cast<int>(std::floor((emax_ev - emin_ev) / de_ev));

    if (npoints <= 0)
    {
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "npoints", npoints);
        ModuleBase::WARNING("ModuleIO::cal_dos", "npoints <= 0");
        return false;
    }
    if (GlobalV::MY_RANK == 0)
    {
        ofs_dos << npoints << " # number of points" << std::endl;
        ofs_dos << nkstot << " # number of total k-points" << std::endl;
    }

    std::vector<double> e_mod(npoints, 0.0); 

    double sum = 0.0;
    double e_new = emin_ev;
    double e_old = 0.0;

    while (e_new < emax_ev)
    {
        double nstates = 0.0;
        e_old = e_new;
        e_new += de_ev;

        // nks is the number of k-points in the 'pool'
        for (int ik = 0; ik < nks; ik++)
        {
            // spin index
            if (is == isk[ik])
            {
                // band index
                for (int ib = 0; ib < nbands; ib++)
                {
                    //  compare et and e_old(e_new) in ev unit.
                    if (ekb(ik, ib) * ModuleBase::Ry_to_eV >= e_old 
                     && ekb(ik, ib) * ModuleBase::Ry_to_eV < e_new)
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
            ofs_dos << std::setw(15) << e_new 
                    << std::setw(15) << nstates << std::endl;
            dos.push_back(nstates);
            ene.push_back(e_new);
        }
    }

    // Use Gaussian smearing to smooth the DOS
    if (GlobalV::MY_RANK == 0)
    {
        dos_smear.resize(dos.size() - 1);

        double b = sqrt(2.0) * bcoeff;
        for (int i = 0; i < dos.size() - 1; i++)
        {
            double Gauss = 0.0;

            for (int j = 0; j < dos.size() - 1; j++)
            {
                double de = ene[j] - ene[i];
                double de2 = de * de;
                Gauss = exp(-de2 / b / b) / sqrt(ModuleBase::PI) / b;
                dos_smear[j] += dos[i] * Gauss;
            }
        }

        double sum2 = 0.0;
        ofs_smear << "#" << std::setw(14) << "energy" 
                 << std::setw(15) << "dos_smear" 
                 << std::setw(15) << "sum_elec" << std::endl;

        for (int i = 0; i < dos.size() - 1; i++)
        {
            sum2 += dos_smear[i];

            ofs_smear << std::setw(15) << ene[i] 
                 << std::setw(15) << dos_smear[i] 
                 << std::setw(15) << sum2 << std::endl;
        }
    }

    if (GlobalV::MY_RANK == 0)
    {
        ofs_dos.close();
        ofs_smear.close();
    }

    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "number of bands", nbands);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "sum up the electronic states", sum);

    return true;
}
