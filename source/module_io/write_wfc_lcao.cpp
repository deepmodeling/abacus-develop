#include "write_wfc_lcao.h"

#include "module_base/memory.h"
#include "module_base/timer.h"
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_wfc.h"
#include "write_wfc_nao.h"


namespace ModuleIO
{

std::string wfc_lcao_gen_fname(const int& out_type,
                               const bool& gamma_only,
                               const bool& out_app_flag,
                               const int& ik,
                               const int& istep)
{
    // fn_out = "{GlobalV::global_out_dir}/WFC_LCAO_{K|GAMMA}{K index}{_ION} + {".txt"/".dat"}""
    std::string kgamma_block = (gamma_only) ? "_GAMMA" : "_K";
    std::string istep_block = (istep >=0 && (! out_app_flag)) ? "_ION" + std::to_string(istep+1): ""; // only when istep >= 0 and out_app_flag is true will write each wfc to a separate file
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

#ifdef __MPI
void wfc_lcao_bcast_work(const int& maxnloc, 
                         const int& src_rank, 
                         double* work, 
                         const Parallel_Orbitals* pv)
{

    int info = 0;
    info = MPI_Bcast(work, maxnloc, MPI_DOUBLE, src_rank, pv->comm_2D);
}

void wfc_lcao_bcast_work(const int& maxnloc,
                         const int& src_rank,
                         std::complex<double>* work,
                         const Parallel_Orbitals* pv)
{
    int info = 0;
    info = MPI_Bcast(work, maxnloc, MPI_DOUBLE_COMPLEX, src_rank, pv->comm_2D);
}
#endif

template <typename T>
void wfc_lcao_collect_wfc2d(const psi::Psi<T>& psi,
                   const int& nbands,
                   const int& nlocal,
                   const Parallel_Orbitals* pv,
                   const int& myid,
                   std::vector<std::vector<T>>& ctot)
{
    static_assert(std::is_same<T, double>::value || std::is_same<T, std::complex<double>>::value, "wfc_lcao_collect_wfc2d: type of psi must be double or std::complex<double>");
    const int inc = 1;
    int total_processes = 1;
#ifdef __MPI
    MPI_Comm_size(MPI_COMM_WORLD, &total_processes);
#endif
    if (myid == 0)
    {
        ctot.assign(nbands, std::vector<T>(nlocal));
    }
    ModuleBase::Memory::record("ModuleIO::write_wfc_lcao::ctot", sizeof(T) * nbands * nlocal);

    if (total_processes == 1 && myid == 0)
    {
        // only one process, and just transfer wfc_2d to ctot
        for (int i = 0; i < nbands; i++)
        {
            for (int j = 0; j < nlocal; j++)
            {
                ctot[i][j] = psi(i, j);
            }
        }
    }
    else
    {
#ifdef __MPI
        int info = 0;
        long maxnloc = pv->nloc_wfc; // maximum number of elements in local matrix
        // get the maximum number of elements in local matrix
        info = MPI_Reduce(&pv->nloc_wfc, &maxnloc, 1, MPI_LONG, MPI_MAX, 0, pv->comm_2D);
        info = MPI_Bcast(&maxnloc, 1, MPI_LONG, 0, pv->comm_2D);
        std::vector<T> work(maxnloc);
        ModuleBase::Memory::record("ModuleIO::write_wfc_lcao::work", sizeof(T) * maxnloc);

        int naroc[2] = {pv->nrow, pv->ncol_bands}; // maximum number of row or column
        for (int iprow = 0; iprow < pv->dim0; ++iprow)
        {
            for (int ipcol = 0; ipcol < pv->dim1; ++ipcol)
            {
                const int coord[2] = {iprow, ipcol};
                int src_rank = 0;
                info = MPI_Cart_rank(pv->comm_2D, coord, &src_rank);
                if (myid == src_rank)
                {
                    BlasConnector::copy(pv->nloc_wfc, psi.get_pointer(), inc, work.data(), inc);
                    naroc[0] = pv->nrow;
                    naroc[1] = pv->ncol_bands;
                }
                info = MPI_Bcast(naroc, 2, MPI_INT, src_rank, pv->comm_2D);
                wfc_lcao_bcast_work(maxnloc, src_rank, work.data(), pv);

                // info = MPI_Bcast(work, maxnloc, MPI_DOUBLE, src_rank, pv->comm_2D);

                if (myid == 0)
                {
                    for (int j = 0; j < naroc[1]; ++j)
                    {
                        int igcol = Local_Orbital_wfc::globalIndex(j, pv->nb, pv->dim1, ipcol);
                        if (igcol >= nbands)
                        {
                            continue;
                        }
                        for (int i = 0; i < naroc[0]; ++i)
                        {
                            int igrow = Local_Orbital_wfc::globalIndex(i, pv->nb, pv->dim0, iprow);
                            ctot[igcol][igrow] = work[j * naroc[0] + i];
                        }
                    }
                }
            } // loop ipcol
        } // loop iprow
#endif
    }
}

void wfc_lcao_write2file(const std::string &name, std::vector<std::vector<double>>& ctot, const int &ik, const ModuleBase::Vector3<double> &kvec_c, const ModuleBase::matrix& ekb, const ModuleBase::matrix& wg, bool writeBinary)
{
    write_wfc_nao(name, ctot, ekb, wg, writeBinary);
}

void wfc_lcao_write2file(const std::string &name, std::vector<std::vector<std::complex<double>>>& ctot, const int &ik, const ModuleBase::Vector3<double> &kvec_c, const ModuleBase::matrix& ekb, const ModuleBase::matrix& wg, bool writeBinary)
{
    write_wfc_nao_complex(name, ctot, ik, kvec_c, ekb, wg, writeBinary);
}

template <typename T>
void write_wfc_lcao(const int out_type,
                    const psi::Psi<T>& psi,
                    const ModuleBase::matrix& ekb,
                    const ModuleBase::matrix& wg,
                    const std::vector<ModuleBase::Vector3<double>>& kvec_c,
                    const Parallel_Orbitals* pv,
                    const int& istep)
{
    if (! out_type)
    {
        return;
    }
    ModuleBase::TITLE("ModuleIO", "write_wfc_lcao");
    ModuleBase::timer::tick("ModuleIO", "write_wfc_lcao");
    std::vector<std::vector<T>> ctot;
    int myid = 0;
    int nbands;
    int nlocal;
#ifdef __MPI
    MPI_Comm_rank(pv->comm_2D, &myid);
    nlocal = pv->desc[2];
#else
    nlocal = psi.get_nbasis();
#endif
    nbands = ekb.nc;

    bool gamma_only = false;
    if (std::is_same<T, double>::value)
    {
        gamma_only = true;
    }

    for (int ik=0;ik<psi.get_nk();ik++)
    {
        psi.fix_k(ik);
        wfc_lcao_collect_wfc2d(psi, nbands, nlocal, pv, myid, ctot);
        if (myid == 0)
        {
            std::string fn = wfc_lcao_gen_fname(out_type, gamma_only, GlobalV::out_app_flag, ik, istep);
            bool writeBinary = (out_type == 2);
            wfc_lcao_write2file(fn, ctot, ik, kvec_c[ik], ekb, wg, writeBinary);
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
                                       const int& istep);

template void write_wfc_lcao<std::complex<double>>(const int out_type,
                                       const psi::Psi<std::complex<double>>& psi,
                                       const ModuleBase::matrix& ekb,
                                       const ModuleBase::matrix& wg,
                                       const std::vector<ModuleBase::Vector3<double>>& kvec_c,
                                       const Parallel_Orbitals* pv,
                                       const int& istep);


} // namespace ModuleIO

