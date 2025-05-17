#include "write_wfc_nao.h"

#include "module_parameter/parameter.h"
#include "module_base/memory.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"
#include "module_base/parallel_2d.h"
#include "module_base/scalapack_connector.h"
#include "module_base/global_variable.h"
#include "module_base/global_function.h"
#include "binstream.h"
#include "filename.h"

namespace ModuleIO
{

void wfc_nao_write2file(const std::string& name,
                        const double* ctot,
                        const int nlocal,
                        const int ik,
                        const ModuleBase::matrix& ekb,
                        const ModuleBase::matrix& wg,
                        const bool& writeBinary,
                        const bool& append_flag)
{
    ModuleBase::TITLE("ModuleIO", "wfc_nao_write2file");
    ModuleBase::timer::tick("ModuleIO", "wfc_nao_write2file");

    int nbands = ekb.nc;

    if (writeBinary)
    {
        Binstream ofs;
        if (append_flag)
        {
            ofs.open(name, "a");
        }
        else
        {
            ofs.open(name, "w");
        }
        if (!ofs)
        {
            ModuleBase::WARNING("ModuleIO::wfc_nao_write2file", "Can't write local orbital wave functions.");
        }

        ofs << nbands;
        ofs << nlocal;

        for (int i = 0; i < nbands; i++)
        {
            ofs << i + 1;
            ofs << ekb(ik, i);
            ofs << wg(ik, i);

            for (int j = 0; j < nlocal; j++)
            {
                ofs << ctot[i * nlocal + j];
            }
        }
        ofs.close();
    }
    else
    {
        std::ofstream ofs;
        if (append_flag)
        {
            ofs.open(name.c_str(), std::ofstream::app);
        }
        else
        {
            ofs.open(name.c_str());
        }
        if (!ofs)
        {
            ModuleBase::WARNING("ModuleIO::wfc_nao_write2file", "Can't write local orbital wave functions.");
        }
        ofs << nbands << " (number of bands)" << std::endl;
        ofs << nlocal << " (number of orbitals)";
        ofs << std::setprecision(8);
        ofs << std::scientific;

        for (int i = 0; i < nbands; i++)
        {
            // +1 to mean more clearly.
            // band index start from 1.
            ofs << "\n" << i + 1 << " (band)";
            ofs << "\n" << ekb(ik, i) << " (Ry)";
            ofs << "\n" << wg(ik, i) << " (Occupations)";
            for (int j = 0; j < nlocal; j++)
            {
                if (j % 5 == 0)
                {
                    ofs << "\n";
                }
                ofs << ctot[i * nlocal + j] << " ";
            }
        }
        ofs << std::endl;
        ofs.close();
    }

    ModuleBase::timer::tick("ModuleIO", "wfc_nao_write2file");
    return;
}

void wfc_nao_write2file_complex(const std::string& name,
                                const std::complex<double>* ctot,
                                const int nlocal,
                                const int& ik,
                                const ModuleBase::Vector3<double>& kvec_c,
                                const ModuleBase::matrix& ekb,
                                const ModuleBase::matrix& wg,
                                const bool& writeBinary,
                                const bool& append_flag)
{
    ModuleBase::TITLE("ModuleIO","wfc_nao_write2file_complex");
    ModuleBase::timer::tick("ModuleIO","wfc_nao_write2file_complex");

    int nbands = ekb.nc;

    if (writeBinary)
    {
        Binstream ofs;
        if (append_flag)
        {
            ofs.open(name, "a");
        }
        else
        {
            ofs.open(name, "w");
        }
        if (!ofs)
        {
            ModuleBase::WARNING("ModuleIO::wfc_nao_write2file_complex", "Can't write local orbital wave functions.");
        }
        ofs << ik + 1;
        ofs << kvec_c.x;
        ofs << kvec_c.y;
        ofs << kvec_c.z;
        ofs << nbands;
        ofs << nlocal;

        for (int i = 0; i < nbands; i++)
        {
            ofs << i + 1;
            ofs << ekb(ik, i);
            ofs << wg(ik, i);

            for (int j = 0; j < nlocal; j++)
            {
                ofs << ctot[i * nlocal + j].real() << ctot[i * nlocal + j].imag();
            }
        }
        ofs.close();
    }
    else
    {
        std::ofstream ofs;
        if (append_flag)
        {
            ofs.open(name.c_str(), std::ofstream::app);
        }
        else
        {
            ofs.open(name.c_str());
        }
        if (!ofs)
        {
            ModuleBase::WARNING("ModuleIO::wfc_nao_write2file_complex", "Can't write local orbital wave functions.");
        }
        ofs << std::setprecision(8);
        ofs << ik + 1 << " (index of k points)" << std::endl;
        ofs << kvec_c.x << " " << kvec_c.y << " " << kvec_c.z << std::endl;
        ofs << nbands << " (number of bands)" << std::endl;
        ofs << nlocal << " (number of orbitals)";
        ofs << std::scientific;

        for (int i = 0; i < nbands; i++)
        {
            // +1 to mean more clearly.
            // band index start from 1.
            ofs << "\n" << i + 1 << " (band)";
            ofs << "\n" << ekb(ik, i) << " (Ry)";
            ofs << "\n" << wg(ik, i) << " (Occupations)";
            for (int j = 0; j < nlocal; j++)
            {
                if (j % 5 == 0)
                {
                    ofs << "\n";
                }
                ofs << ctot[i * nlocal + j].real() << " " << ctot[i * nlocal + j].imag() << " ";
            }
        }
        ofs << std::endl;
        ofs.close();
    }

    ModuleBase::timer::tick("ModuleIO","wfc_nao_write2file_complex");
    return;
}

template <typename T>
void write_wfc_nao(const int out_type,
		const bool out_app_flag,
		const psi::Psi<T>& psi,
		const ModuleBase::matrix& ekb,
		const ModuleBase::matrix& wg,
		const std::vector<ModuleBase::Vector3<double>>& kvec_c,
		const std::vector<int> &ik2iktot,
		const int nkstot,
		const Parallel_Orbitals& pv,
		const int nspin,
		const int istep)
{
    if (!out_type)
    {
        return;
    }
    ModuleBase::TITLE("ModuleIO", "write_wfc_nao");
    ModuleBase::timer::tick("ModuleIO", "write_wfc_nao");
    int myid = 0;
    int nbands = 0;
    int nlocal = 0;

    // If using MPI, the nbasis and nbands in psi is the value on local rank, 
    // so get nlocal and nbands from pv->desc_wfc[2] and pv->desc_wfc[3]
#ifdef __MPI
    MPI_Comm_rank(pv.comm(), &myid);
    nlocal = pv.desc_wfc[2];
    nbands = pv.desc_wfc[3];
#else
    nlocal = psi.get_nbasis();
    nbands = psi.get_nbands();
#endif

    bool gamma_only = (std::is_same<T, double>::value);
    bool writeBinary = (out_type == 2);
    Parallel_2D pv_glb;
    int blk_glb = std::max(nlocal, nbands);
    std::vector<T> ctot(myid == 0 ? nbands * nlocal : 0);
    ModuleBase::Memory::record("ModuleIO::write_wfc_nao::glb", sizeof(T) * nlocal * nbands);

    for (int ik = 0; ik < psi.get_nk(); ik++)
    {
        psi.fix_k(ik);
#ifdef __MPI        
        pv_glb.set(nlocal, nbands, blk_glb, pv.blacs_ctxt);   
        Cpxgemr2d(nlocal,
                  nbands,
                  psi.get_pointer(),
                  1,
                  1,
                  const_cast<int*>(pv.desc_wfc),
                  ctot.data(),
                  1,
                  1,
                  pv_glb.desc,
                  pv_glb.blacs_ctxt);
#else
        for (int ib = 0; ib < nbands; ib++)
        {
            for (int i = 0; i < nlocal; i++)
            {
                ctot[ib * nlocal + i] = psi(ib,i);
            }
        }    
#endif

        if (myid == 0)
        {
            std::string fn = filename_output(PARAM.globalv.global_out_dir,"wf","nao",ik,ik2iktot,nspin,nkstot,
              out_type,out_app_flag,gamma_only,istep);

            bool append_flag = (istep > 0 && out_app_flag);
            if (std::is_same<double, T>::value)
            {
                wfc_nao_write2file(fn,
                                   reinterpret_cast<double*>(ctot.data()),
                                   nlocal,
                                   ik,
                                   ekb,
                                   wg,
                                   writeBinary,
                                   append_flag);
            }
            else
            {
                wfc_nao_write2file_complex(fn,
                                           reinterpret_cast<std::complex<double>*>(ctot.data()),
                                           nlocal,
                                           ik,
                                           kvec_c[ik],
                                           ekb,
                                           wg,
                                           writeBinary,
                                           append_flag);
            }
        }
    }
    ModuleBase::timer::tick("ModuleIO", "write_wfc_nao");
}

template void write_wfc_nao<double>(const int out_type,
		const bool out_app_flag,
		const psi::Psi<double>& psi,
		const ModuleBase::matrix& ekb,
		const ModuleBase::matrix& wg,
		const std::vector<ModuleBase::Vector3<double>>& kvec_c,
		const std::vector<int> &ik2iktot,
		const int nkstot,
		const Parallel_Orbitals& pv,
		const int nspin,
		const int istep);

template void write_wfc_nao<std::complex<double>>(const int out_type,
		const bool out_app_flag,
		const psi::Psi<std::complex<double>>& psi,
		const ModuleBase::matrix& ekb,
		const ModuleBase::matrix& wg,
		const std::vector<ModuleBase::Vector3<double>>& kvec_c,
		const std::vector<int> &ik2iktot,
		const int nkstot,
		const Parallel_Orbitals& pv,
		const int nspin,
		const int istep);

} // namespace ModuleIO
