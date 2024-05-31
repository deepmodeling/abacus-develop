#include "write_wfc_lcao.h"

#include "module_base/memory.h"
#include "module_base/timer.h"
#include "module_basis/module_ao/parallel_2d.h"
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_wfc.h"
#include "write_wfc_nao.h"
#include "module_base/scalapack_connector.h"

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
        = GlobalV::global_out_dir + "WFC_LCAO" + kgamma_block + std::to_string(ik + 1) + istep_block + suffix_block;
    return fn_out;
}

template <typename T>
void write_wfc_lcao(const int out_type,
                    const psi::Psi<T>& psi,
                    const ModuleBase::matrix& ekb,
                    const ModuleBase::matrix& wg,
                    const std::vector<ModuleBase::Vector3<double>>& kvec_c,
                    const Parallel_Orbitals* pv,
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
#ifdef __MPI
    MPI_Comm_rank(pv->comm_2D, &myid);
    nbands = pv->desc_wfc[3];
    nlocal = pv->desc_wfc[2];
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
        pv_glb.set(nlocal, nbands, blk_glb, pv->comm_2D, pv->blacs_ctxt);
        Cpxgemr2d(nlocal,
                  nbands,
                  psi.get_pointer(),
                  1,
                  1,
                  const_cast<Parallel_Orbitals*>(pv)->desc_wfc,
                  ctot.data(),
                  1,
                  1,
                  pv_glb.desc,
                  pv_glb.blacs_ctxt);

        if (myid == 0)
        {
            std::string fn = wfc_lcao_gen_fname(out_type, gamma_only, GlobalV::out_app_flag, ik, istep);
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
                                     const Parallel_Orbitals* pv,
                                     const int istep);

template void write_wfc_lcao<std::complex<double>>(const int out_type,
                                                   const psi::Psi<std::complex<double>>& psi,
                                                   const ModuleBase::matrix& ekb,
                                                   const ModuleBase::matrix& wg,
                                                   const std::vector<ModuleBase::Vector3<double>>& kvec_c,
                                                   const Parallel_Orbitals* pv,
                                                   const int istep);

} // namespace ModuleIO
