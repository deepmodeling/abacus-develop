#include "write_wfc_lcao.h"

#include "module_base/memory.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"
#include "module_basis/module_ao/parallel_2d.h"
#include "module_base/scalapack_connector.h"
#include "module_base/global_variable.h"
#include "binstream.h"
#include "module_base/global_function.h"

namespace ModuleIO
{

std::string wfc_lcao_gen_fname(const int out_type,
                               const bool gamma_only,
                               const bool out_app_flag,
                               const int ik,
                               const int istep)
{
    // fn_out = "{GlobalV::global_out_dir}/WFC_LCAO_{K|GAMMA}{K index}{_ION} + {".txt"/".dat"}""
    std::string kgamma_block = (gamma_only) ? "_GAMMA" : "_K";
    std::string istep_block
        = (istep >= 0 && (!out_app_flag))
              ? "_ION" + std::to_string(istep + 1)
              : ""; // only when istep >= 0 and out_app_flag is true will write each wfc to a separate file
    std::string suffix_block = "";

    if (out_type == 1)
    {
        suffix_block = ".txt";
    }
    else if (out_type == 2)
    {
        suffix_block = ".dat";
    }
    else
    {
        std::cout << "WARNING: the out type of wave function is not 1 or 2. Write to a txt file." << std::endl;
        suffix_block = ".txt";
    }

    std::string fn_out
        = "WFC_LCAO" + kgamma_block + std::to_string(ik + 1) + istep_block + suffix_block;
    return fn_out;
}

void write_wfc_nao(const std::string &name, const double* ctot, const int nlocal, const ModuleBase::matrix& ekb, const ModuleBase::matrix& wg, bool writeBinary)
{
    ModuleBase::TITLE("ModuleIO", "write_wfc_nao");
    ModuleBase::timer::tick("ModuleIO", "write_wfc_nao");

    //if (GlobalV::DRANK == 0)
    {
        int nbands = ekb.nc;
        
        if (writeBinary)
        {
            Binstream ofs(name, "a");
            if (!ofs)
            {
                ModuleBase::WARNING("ModuleIO::write_wfc_nao", "Can't write local orbital wave functions.");
            }

            ofs << nbands;
            ofs << nlocal;

            for (int i = 0; i < nbands; i++)
            {
                ofs << i+1;
                ofs << ekb(GlobalV::CURRENT_SPIN, i);
                ofs << wg(GlobalV::CURRENT_SPIN, i);

                for (int j = 0; j < nlocal; j++)
                {
                    ofs << ctot[i*nlocal + j];
                }
            }
            ofs.close();
        }
        else
        {
            std::ofstream ofs;
            // if (GlobalV::out_app_flag)
            // {
            //     ofs.open(name.c_str(), std::ofstream::app);
            // }
            // else
            {   // the default value of `out_app_flag`is true, but usually there's no use to save each step's LCAO wave function.
                ofs.open(name.c_str());
            }
            if (!ofs)
            {
                ModuleBase::WARNING("ModuleIO::write_wfc_nao", "Can't write local orbital wave functions.");
            }
            ofs << nbands << " (number of bands)" << std::endl;
            ofs << nlocal << " (number of orbitals)";
            ofs << std::setprecision(8);
            ofs << std::scientific;

            for (int i=0; i<nbands; i++)
            {
                // +1 to mean more clearly.
                // band index start from 1.
                ofs << "\n" << i+1 << " (band)";
		    	ofs << "\n" << ekb(GlobalV::CURRENT_SPIN, i) << " (Ry)";
		    	ofs << "\n" << wg(GlobalV::CURRENT_SPIN,i) << " (Occupations)";
                for (int j=0; j<nlocal; j++)
                {
                    if (j % 5 == 0) ofs << "\n";
                    ofs << ctot[i*nlocal + j] << " ";
                }
            }
            ofs << std::endl;
            ofs.close();
        }
    }

    ModuleBase::timer::tick("ModuleIO", "write_wfc_nao");
    return;
}

void write_wfc_nao_complex(const std::string &name, const std::complex<double>* ctot, const int nlocal,const int &ik, const ModuleBase::Vector3<double> &kvec_c, const ModuleBase::matrix& ekb, const ModuleBase::matrix& wg, bool writeBinary)
{
    ModuleBase::TITLE("ModuleIO","write_wfc_nao_complex");
    ModuleBase::timer::tick("ModuleIO","write_wfc_nao_complex");

    
    //if (GlobalV::DRANK==0)
    {
        int nbands = ekb.nc;

        if (writeBinary)
        {
            Binstream ofs(name, "a");
            if (!ofs)
            {
                ModuleBase::WARNING("ModuleIO::write_wfc_nao", "Can't write local orbital wave functions.");
            }
            ofs << ik+1;
            ofs << kvec_c.x;
            ofs << kvec_c.y;
            ofs << kvec_c.z;
            ofs << nbands;
            ofs << nlocal;

            for (int i = 0; i < nbands; i++)
            {
                ofs << i+1;
                ofs << ekb(ik, i);
                ofs << wg(ik, i);

                for (int j = 0; j < nlocal; j++)
                {
                    ofs << ctot[i*nlocal + j].real() << ctot[i*nlocal + j].imag();
                }
            }
            ofs.close();
        }
        else
        {
            std::ofstream ofs;
            // if (GlobalV::out_app_flag)
            // {
            //     ofs.open(name.c_str(), std::ofstream::app);
            // }
            // else
            {   // the default value of `out_app_flag`is true, but usually there's no use to save each step's LCAO wave function.
                ofs.open(name.c_str());
            }
            if (!ofs)
            {
                ModuleBase::WARNING("ModuleIO::write_wfc_nao","Can't write local orbital wave functions.");
            }
            ofs << std::setprecision(25);
		    ofs << ik+1 << " (index of k points)" << std::endl;
		    ofs << kvec_c.x << " " << kvec_c.y << " " << kvec_c.z << std::endl;
            ofs << nbands << " (number of bands)" << std::endl;
            ofs << nlocal << " (number of orbitals)";
            ofs << std::scientific;

            for (int i=0; i<nbands; i++)
            {
                // +1 to mean more clearly.
                // band index start from 1.
                ofs << "\n" << i+1 << " (band)";
		    	ofs << "\n" << ekb(ik, i) << " (Ry)";
		    	ofs << "\n" << wg(ik,i) << " (Occupations)";
                for (int j=0; j<nlocal; j++)
                {
                    if (j % 5 == 0) ofs << "\n";
                    ofs << ctot[i*nlocal + j].real() << " " << ctot[i*nlocal + j].imag() << " ";
                }
            }
            ofs << std::endl;
            ofs.close();
        }
    }

    ModuleBase::timer::tick("ModuleIO","write_wfc_nao_complex");
    return;
}

template <typename T>
void write_wfc_lcao(const int out_type,
                    const psi::Psi<T>& psi,
                    const ModuleBase::matrix& ekb,
                    const ModuleBase::matrix& wg,
                    const std::vector<ModuleBase::Vector3<double>>& kvec_c,
                    const Parallel_Orbitals& pv,
                    const int istep)
{
    if (!out_type)
    {
        return;
    }
    ModuleBase::TITLE("ModuleIO", "write_wfc_lcao");
    ModuleBase::timer::tick("ModuleIO", "write_wfc_lcao");
    int myid = 0;
    int nbands;
    int nlocal;
    // If using MPI, the nbasis and nbands in psi is the value on local rank, 
    // so get nlocal and nbands from pv->desc_wfc[2] and pv->desc_wfc[3]
#ifdef __MPI
    MPI_Comm_rank(pv.comm_2D, &myid);
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
    ModuleBase::Memory::record("ModuleIO::write_wfc_lcao::glb", sizeof(T) * nlocal * nbands);

    for (int ik = 0; ik < psi.get_nk(); ik++)
    {
        psi.fix_k(ik);
#ifdef __MPI        
        pv_glb.set(nlocal, nbands, blk_glb, pv.comm_2D, pv.blacs_ctxt);   
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
            std::string fn = GlobalV::global_out_dir + wfc_lcao_gen_fname(out_type, gamma_only, GlobalV::out_app_flag, ik, istep);
            if (std::is_same<double, T>::value)
            {
                write_wfc_nao(fn, reinterpret_cast<double*>(ctot.data()), nlocal, ekb, wg, writeBinary);
            }
            else
            {
                write_wfc_nao_complex(fn,
                                      reinterpret_cast<std::complex<double>*>(ctot.data()),
                                      nlocal,
                                      ik,
                                      kvec_c[ik],
                                      ekb,
                                      wg,
                                      writeBinary);
            }
        }
    }
    ModuleBase::timer::tick("ModuleIO", "write_wfc_lcao");
}

template void write_wfc_lcao<double>(const int out_type,
                                     const psi::Psi<double>& psi,
                                     const ModuleBase::matrix& ekb,
                                     const ModuleBase::matrix& wg,
                                     const std::vector<ModuleBase::Vector3<double>>& kvec_c,
                                     const Parallel_Orbitals& pv,
                                     const int istep);

template void write_wfc_lcao<std::complex<double>>(const int out_type,
                                                   const psi::Psi<std::complex<double>>& psi,
                                                   const ModuleBase::matrix& ekb,
                                                   const ModuleBase::matrix& wg,
                                                   const std::vector<ModuleBase::Vector3<double>>& kvec_c,
                                                   const Parallel_Orbitals& pv,
                                                   const int istep);

} // namespace ModuleIO
