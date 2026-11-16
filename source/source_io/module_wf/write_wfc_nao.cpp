#include "write_wfc_nao.h"

#include "source_io/module_parameter/parameter.h"
#include "source_base/fs_compat.h"
#include "source_base/memory_recorder.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "source_base/parallel_blkcyc.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_base/global_variable.h"
#include "source_base/global_function.h"
#include "source_base/module_out/binstream.h"
#include "source_base/module_out/filename.h"
#include "source_base/formatter.h"
#include "source_io/module_output/provenance_header.h"
#include "source_io/module_output/spin_tag.h"

namespace
{

/**
 * @brief spin channel of a global k-point index, consistent with the s{N}
 *        block of filename_output.
 *
 * @param ik0 global k-point index (the ik2iktot value)
 * @param nspin number of spins
 * @param nkstot total number of k-points including the spin multiplicity
 * @return 0 for spin up, 1 for spin down, -1 when nspin != 2
 */
int wfc_spin_index(const int ik0, const int nspin, const int nkstot)
{
    if (nspin != 2)
    {
        return -1;
    }
    const int half_k = nkstot / 2;
    if (ik0 >= half_k)
    {
        return 1;
    }
    return 0;
}

/**
 * @brief build the provenance header that precedes the wavefunction data
 *        block, mirroring the bands1.txt/bands2.txt header style.
 *
 * In append mode (out_app_flag == 1) every accumulated frame carries its own
 * header, so readers can tell the ionic steps apart.
 *
 * @param istep zero-based ionic step index, printed one-based; the ionic-step
 *        line is omitted when istep < 0 (no step context in the run)
 * @param spin_tag label of the collinear spin channel, empty when nspin != 2
 * @return two or three comment lines terminated by newlines
 */
std::string build_wfc_nao_header(const int istep, const std::string& spin_tag)
{
    const std::string header = ModuleIO::provenance_header(ModuleIO::local_timestamp());
    if (istep < 0)
    {
        return header;
    }
    const int istep_one_based = istep + 1;
    const std::string step_line = FmtCore::format("# IONIC STEP %d%s\n", istep_one_based, spin_tag.c_str());
    return header + step_line;
}

} // namespace

namespace ModuleIO
{

void wfc_nao_write2file(const std::string& name,
                        const double* ctot,
                        const int nlocal,
                        const int ik,
                        const ModuleBase::matrix& ekb,
                        const ModuleBase::matrix& wg,
                        const bool& writeBinary,
                        const bool& append_flag,
                        const std::string& header)
{
    ModuleBase::TITLE("ModuleIO", "wfc_nao_write2file");
    ModuleBase::timer::start("ModuleIO", "wfc_nao_write2file");

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
            ModuleBase::WARNING_QUIT("ModuleIO::wfc_nao_write2file", "Can't write local orbital wave functions.");
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
            ModuleBase::WARNING_QUIT("ModuleIO::wfc_nao_write2file", "Can't write local orbital wave functions.");
        }
        // provenance header precedes the data block; the binary branch above
        // keeps its fixed record layout and ignores the header
        if (!header.empty())
        {
            ofs << header;
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

    ModuleBase::timer::end("ModuleIO", "wfc_nao_write2file");
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
                                const bool& append_flag,
                                const std::string& header)
{
    ModuleBase::TITLE("ModuleIO","wfc_nao_write2file_complex");
    ModuleBase::timer::start("ModuleIO","wfc_nao_write2file_complex");

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
            ModuleBase::WARNING_QUIT("ModuleIO::wfc_nao_write2file_complex", "Can't write local orbital wave functions.");
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
            ModuleBase::WARNING_QUIT("ModuleIO::wfc_nao_write2file_complex", "Can't write local orbital wave functions.");
        }
        // provenance header precedes the k-point record; the binary branch
        // above keeps its fixed record layout and ignores the header
        if (!header.empty())
        {
            ofs << header;
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

    ModuleBase::timer::end("ModuleIO","wfc_nao_write2file_complex");
    return;
}

template <typename T>
void write_wfc_nao(const int out_type,
        const int out_app_flag,
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
    ModuleBase::timer::start("ModuleIO", "write_wfc_nao");
    int myid = 0;
    int nbands = 0;
    int nlocal = 0;

    // If using MPI, the nbasis and nbands in psi is the value on local rank, 
    // so get nlocal and nbands from pv->desc_wfc[2] and pv->desc_wfc[3]
#ifdef __MPI
    MPI_Comm_rank(pv.para_2d().comm(), &myid);
    nlocal = pv.desc_wfc[2];
    nbands = pv.desc_wfc[3];
#else
    nlocal = psi.get_nbasis();
    nbands = psi.get_nbands();
#endif

    bool gamma_only = (std::is_same<T, double>::value);
    bool writeBinary = (out_type == 2);
    Parallel_Blkcyclic pv_glb;
    int blk_glb = std::max(nlocal, nbands);
    std::vector<T> ctot(myid == 0 ? nbands * nlocal : 0);
    ModuleBase::Memory::record("ModuleIO::write_wfc_nao::glb", sizeof(T) * nlocal * nbands);

    for (int ik = 0; ik < psi.get_nk(); ik++)
    {
        psi.fix_k(ik);
#ifdef __MPI
        pv_glb.set(nlocal, nbands, blk_glb, pv.para_2d().ictxt());
        Cpxgemr2d(nlocal,
                  nbands,
                  psi.get_pointer(),
                  1,
                  1,
                  const_cast<int*>(pv.desc_wfc),
                  ctot.data(),
                  1,
                  1,
                  const_cast<int*>(pv_glb.desc()),
                  pv_glb.ictxt());
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
            std::string wfc_dir;
            if(out_app_flag == 1)
            {
                wfc_dir = PARAM.globalv.global_out_dir;
            }
            else
            {
                wfc_dir = PARAM.globalv.global_wfc_dir;
            }
            // The mode-2 PSI/ subdirectory is created lazily at write time;
            // otherwise this is a no-op on the existing output directory.
            ModuleBase::ensure_directory(wfc_dir);
            std::string fn = filename_output(wfc_dir,"wf","nao",ik,ik2iktot,nspin,nkstot,
              out_type,out_app_flag,gamma_only,istep);

            bool append_flag = (istep > 0 && out_app_flag == 1);
            // spin channel of the current k-point and the per-frame
            // provenance header (version, timestamp, ionic step, spin tag)
            const int ik0 = ik2iktot[ik];
            const int is = wfc_spin_index(ik0, nspin, nkstot);
            const std::string spin_tag = ModuleIO::make_compact_spin_tag(is, nspin);
            const std::string header = build_wfc_nao_header(istep, spin_tag);
            if (std::is_same<double, T>::value)
            {
                wfc_nao_write2file(fn,
                                   reinterpret_cast<double*>(ctot.data()),
                                   nlocal,
                                   ik,
                                   ekb,
                                   wg,
                                   writeBinary,
                                   append_flag,
                                   header);
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
                                           append_flag,
                                           header);
            }
        }
    }
    ModuleBase::timer::end("ModuleIO", "write_wfc_nao");
}

template void write_wfc_nao<double>(const int out_type,
        const int out_app_flag,
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
        const int out_app_flag,
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
