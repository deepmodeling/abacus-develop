#include "write_wfc_nao.h"
#include "module_base/timer.h"
#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"

void ModuleIO::write_wfc_nao(const std::string &name, double **ctot, const ModuleBase::matrix& ekb, const ModuleBase::matrix& wg)
{
    ModuleBase::TITLE("ModuleIO","write_wfc_nao");
    ModuleBase::timer::tick("ModuleIO","write_wfc_nao");

    std::ofstream ofs;
    if (GlobalV::DRANK==0)
    {
        if (GlobalV::out_app_flag)
        {
            ofs.open(name.c_str(), std::ofstream::app);
        }
        {
            ofs.open(name.c_str());
        }
        if (!ofs)
        {
            ModuleBase::WARNING("ModuleIO::write_wfc_nao","Can't write local orbital wave functions.");
        }
        ofs << GlobalV::NBANDS << " (number of bands)" << std::endl;
        ofs << GlobalV::NLOCAL << " (number of orbitals)";
        ofs << std::setprecision(8);
        ofs << std::scientific;

        for (int i=0; i<GlobalV::NBANDS; i++)
        {
            // +1 to mean more clearly.
            // band index start from 1.
            ofs << "\n" << i+1 << " (band)";
			ofs << "\n" << ekb(GlobalV::CURRENT_SPIN, i) << " (Ry)"; //mohan add 2012-03-26
			ofs << "\n" << wg(GlobalV::CURRENT_SPIN, i) << " (Occupations)";
            for (int j=0; j<GlobalV::NLOCAL; j++)
            {
                if (j % 5 == 0) ofs << "\n";
                ofs << ctot[i][j] << " ";
            }
        }
        ofs << std::endl;
        ofs.close();
    }

    ModuleBase::timer::tick("ModuleIO","write_wfc_nao");
    return;
}

void ModuleIO::write_wfc_nao_complex(const std::string &name, std::complex<double> **ctot, const int &ik, const ModuleBase::Vector3<double> &kvec_c, const ModuleBase::matrix& ekb, const ModuleBase::matrix& wg)
{
    ModuleBase::TITLE("ModuleIO","write_wfc_nao_complex");
    ModuleBase::timer::tick("ModuleIO","write_wfc_nao_complex");

    std::ofstream ofs;
    if (GlobalV::DRANK==0)
    {
        if (GlobalV::out_app_flag)
        {
            ofs.open(name.c_str(), std::ofstream::app);
        }
        {
            ofs.open(name.c_str());
        }
        if (!ofs)
        {
            ModuleBase::WARNING("ModuleIO::write_wfc_nao","Can't write local orbital wave functions.");
        }
        ofs << std::setprecision(25);
		ofs << ik+1 << " (index of k points)" << std::endl;
		ofs << kvec_c.x << " " << kvec_c.y << " " << kvec_c.z << std::endl;
        ofs << GlobalV::NBANDS << " (number of bands)" << std::endl;
        ofs << GlobalV::NLOCAL << " (number of orbitals)";
        ofs << std::scientific;

        for (int i=0; i<GlobalV::NBANDS; i++)
        {
            // +1 to mean more clearly.
            // band index start from 1.
            ofs << "\n" << i+1 << " (band)";
			ofs << "\n" << ekb(ik, i) << " (Ry)";
			ofs << "\n" << wg(ik,i) << " (Occupations)";
            for (int j=0; j<GlobalV::NLOCAL; j++)
            {
                if (j % 5 == 0) ofs << "\n";
                ofs << ctot[i][j].real() << " " << ctot[i][j].imag() << " ";
            }
        }
        ofs << std::endl;
        ofs.close();
    }

    ModuleBase::timer::tick("ModuleIO","write_wfc_nao_complex");
    return;
}
