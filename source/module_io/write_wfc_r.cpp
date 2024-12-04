//======================
// AUTHOR : Peize Lin
#include "module_parameter/parameter.h"
// DATE :   2021-11-21
//======================

//======================
// WARNING: These interfaces will be removed in the future! Do not use them!
// Taoni add 2024-10-08
//======================

#include "module_base/timer.h"
#include "module_base/tool_title.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "write_wfc_r.h"

#include <cstdlib>
#include <fstream>
#include <stdexcept>

namespace ModuleIO
{
// write ||wfc_r|| for all k-points and all bands
// Input: wfc_g(ik, ib, ig)
// loop order is for(z){for(y){for(x)}}
void write_psi_r_1(const UnitCell& ucell,
                   const psi::Psi<std::complex<double>>& wfc_g,
                   const ModulePW::PW_Basis_K* wfcpw,
                   const std::string& folder_name,
                   const bool& square,
                   const K_Vectors& kv)
{
    ModuleBase::TITLE("ModuleIO", "write_psi_r_1");
    ModuleBase::timer::tick("ModuleIO", "write_psi_r_1");

    const std::string outdir = PARAM.globalv.global_out_dir + folder_name + "/";
    ModuleBase::GlobalFunc::MAKE_DIR(outdir);
#ifdef __MPI
    std::vector<MPI_Request> mpi_requests;
#endif
    for (int ik = 0; ik < wfc_g.get_nk(); ++ik)
    {
        wfc_g.fix_k(ik);
        const int ik_out = (PARAM.inp.nspin != 2)
                               ? ik + GlobalC::Pkpoints.startk_pool[GlobalV::MY_POOL]
                               : ik - kv.get_nks() / 2 * kv.isk[ik] + kv.get_nkstot() / 2 * kv.isk[ik]
                                     + GlobalC::Pkpoints.startk_pool[GlobalV::MY_POOL];
        for (int ib = 0; ib < wfc_g.get_nbands(); ++ib)
        {
            const std::vector<std::complex<double>> wfc_r = cal_wfc_r(wfcpw, wfc_g, ik, ib);

            std::vector<double> wfc_real(wfc_r.size());
            std::vector<double> wfc_imag;
            if (square)
            {
                for (int ir = 0; ir < wfc_real.size(); ++ir)
                {
                    wfc_real[ir] = std::norm(wfc_r[ir]); // "std::norm(z)" returns |z|^2
                }
            }
            else
            {
                wfc_imag.resize(wfc_r.size());
                for (int ir = 0; ir < wfc_real.size(); ++ir)
                {
                    wfc_real[ir] = wfc_r[ir].real();
                    wfc_imag[ir] = wfc_r[ir].imag();
                }
            }
            const std::string file_name_base = outdir + "wfc_realspace_" + 
                                         ModuleBase::GlobalFunc::TO_STRING(ik_out) + "_" + 
                                         ModuleBase::GlobalFunc::TO_STRING(ib);
            const std::string file_name = square ?  file_name_base : file_name_base + "_imag";
#ifdef __MPI
            // Use write_chg_r_1 to output the real and imaginary parts of the wave function to file
            mpi_requests.push_back({});
            write_chg_r_1(ucell,wfcpw, wfc_real, file_name, mpi_requests.back());
#else
            write_chg_r_1(ucell,wfcpw, wfc_real, file_name);
            // if (!square)
            // write_chg_r_1(wfc_imag, file_name + "_imag", mpi_requests.back());
#endif
        }
    }
#ifdef __MPI
    MPI_Waitall(mpi_requests.size(), mpi_requests.data(), MPI_STATUSES_IGNORE);
#endif
    ModuleBase::timer::tick("ModuleIO", "write_psi_r_1");
}
// processes output pipeline:
//
//           t0  t1  t2  t3  t4  t5  t6  t7
//          -------------------------------->
//  rank0    k0  k1  k2  k3  k4  k5
//             \   \   \   \   \   \
    //  rank1        k0  k1  k2  k3  k4  k5
//                 \   \   \   \   \   \
    //  rank2            k0  k1  k2  k3  k4  k5

// Input: wfc_g(ib,ig)
// Output: wfc_r[ir]
std::vector<std::complex<double>> cal_wfc_r(const ModulePW::PW_Basis_K* wfcpw,
                                            const psi::Psi<std::complex<double>>& wfc_g,
                                            const int ik,
                                            const int ib)
{
    ModuleBase::timer::tick("ModuleIO", "cal_wfc_r");

    std::vector<std::complex<double>> wfc_r(wfcpw->nrxx);
    wfcpw->recip2real(&wfc_g(ib, 0), wfc_r.data(), ik);

    ModuleBase::timer::tick("ModuleIO", "cal_wfc_r");
    return wfc_r;
}

// Input: chg_r[ir]
void write_chg_r_1(const UnitCell& ucell,
                   const ModulePW::PW_Basis_K* wfcpw,
                   const std::vector<double>& chg_r,
                   const std::string& file_name
                   #ifdef __MPI
                   ,MPI_Request& mpi_request
                   #endif
                  )
{
    ModuleBase::timer::tick("ModuleIO", "write_chg_r_1");
    std::ofstream ofs;

#ifdef __MPI
    constexpr int mpi_tag = 100;
    if (GlobalV::RANK_IN_POOL == 0)
    {
#endif
        ofs.open(file_name);

        ofs << "calculated by ABACUS" << std::endl;
        ofs << ucell.lat0_angstrom << std::endl;
        ofs << ucell.latvec.e11 << " " << ucell.latvec.e12 << " " << ucell.latvec.e13
            << std::endl
            << ucell.latvec.e21 << " " << ucell.latvec.e22 << " " << ucell.latvec.e23
            << std::endl
            << ucell.latvec.e31 << " " << ucell.latvec.e32 << " " << ucell.latvec.e33
            << std::endl;

        for (int it = 0; it < ucell.ntype; ++it)
        {
            ofs << ucell.atoms[it].label << "\t";
        }
        ofs << std::endl;
        for (int it = 0; it < ucell.ntype; ++it)
        {
            ofs << ucell.atoms[it].na << "\t";
        }
        ofs << std::endl;

        ofs << "Direct" << std::endl;
        for (int it = 0; it < ucell.ntype; ++it)
        {
            for (int ia = 0; ia < ucell.atoms[it].na; ++ia)
            {
                ofs << ucell.atoms[it].taud[ia].x << " " << ucell.atoms[it].taud[ia].y << " "
                    << ucell.atoms[it].taud[ia].z << std::endl;
            }
        }
        ofs << std::endl;

        ofs << wfcpw->nx << " " << wfcpw->ny << " " << wfcpw->nz << std::endl;
#ifdef __MPI
    }
    else
    {
        char recv_tmp;
        MPI_Recv(&recv_tmp, 1, MPI_CHAR, GlobalV::RANK_IN_POOL - 1, mpi_tag, POOL_WORLD, MPI_STATUS_IGNORE);

        ofs.open(file_name, std::ofstream::app);
    }
#endif

    assert(wfcpw->nx * wfcpw->ny * wfcpw->nplane == chg_r.size());
    for (int iz = 0; iz < wfcpw->nplane; ++iz)
    {
        for (int iy = 0; iy < wfcpw->ny; ++iy)
        {
            for (int ix = 0; ix < wfcpw->nx; ++ix)
            {
                const int ir = (ix * wfcpw->ny + iy) * wfcpw->nplane + iz;
                ofs << chg_r[ir] << " ";
            }
            ofs << "\n";
        }
    }
    ofs.close();

#ifdef __MPI
    if (GlobalV::RANK_IN_POOL < GlobalV::NPROC_IN_POOL - 1)
    {
        const char send_tmp = 'c';
        MPI_Isend(&send_tmp, 1, MPI_CHAR, GlobalV::RANK_IN_POOL + 1, mpi_tag, POOL_WORLD, &mpi_request);
    }
    else
    {
        mpi_request = MPI_REQUEST_NULL;
    }
#endif
    ModuleBase::timer::tick("ModuleIO", "write_chg_r_1");
}
}; // namespace ModuleIO
