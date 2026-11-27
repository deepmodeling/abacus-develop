#include "to_w90_pw.h"
#include "source_base/parallel_comm.h" // use POOL_WORLD

#include "source_io/module_parameter/parameter.h"
#include "source_base/module_out/binstream.h"

void toW90_PW::cal_Mmn(
    const psi::Psi<std::complex<double>>& psi_pw,
    const ModulePW::PW_Basis_K* wfcpw
)
{
    std::ofstream mmn_file;

    if (GlobalV::MY_RANK == 0)
    {
        std::string fileaddress = PARAM.globalv.global_out_dir + wannier_file_name + ".mmn";
        mmn_file.open(fileaddress.c_str(), std::ios::out);

        time_t time_now = time(nullptr);
        mmn_file << " Created on " << ctime(&time_now);
        mmn_file << std::setw(12) << num_bands << std::setw(12) << cal_num_kpts << std::setw(12) << nntot << std::endl;
    }

    for (int ik = 0; ik < cal_num_kpts; ik++)
    {
        for (int ib = 0; ib < nntot; ib++)
        {
            int ikb = nnlist[ik][ib];
            ModuleBase::Vector3<double> phase_G = nncell[ik][ib];
            ModuleBase::ComplexMatrix Mmn;

            int cal_ik = ik + start_k_index;
            int cal_ikb = ikb + start_k_index;
            unkdotkb(psi_pw, wfcpw, cal_ik, cal_ikb, phase_G, Mmn);

            if (GlobalV::MY_RANK == 0)
            {
                mmn_file << std::setw(5) << ik + 1 << std::setw(5) << ikb + 1 << std::setw(5) << int(phase_G.x)
                         << std::setw(5) << int(phase_G.y) << std::setw(5) << int(phase_G.z) << std::endl;

                for (int n = 0; n < num_bands; n++)
                {
                    for (int m = 0; m < num_bands; m++)
                    {
                        mmn_file << std::setw(18) << std::setprecision(12) << std::showpoint << std::fixed
                                 << Mmn(m, n).real()
                                 << std::setw(18) << std::setprecision(12) << std::showpoint << std::fixed
                                 << Mmn(m, n).imag()
                                 // jingan test
                                 // << "    " << std::setw(12) << std::setprecision(9) << std::abs(Mmn(m, n))
                                 << std::endl;
                    } 
                }
            }

        }
    }

    if (GlobalV::MY_RANK == 0) mmn_file.close();

}


void toW90_PW::cal_Amn(
    const psi::Psi<std::complex<double>>& psi_pw, 
    const ModulePW::PW_Basis_K* wfcpw
)
{
    std::vector<ModuleBase::matrix> radial_in_q;
    gen_radial_function_in_q(radial_in_q);

    std::ofstream Amn_file;

    if (GlobalV::MY_RANK == 0)
    {
        time_t time_now = time(nullptr);
        std::string fileaddress = PARAM.globalv.global_out_dir + wannier_file_name + ".amn";
        Amn_file.open(fileaddress.c_str(), std::ios::out);
        Amn_file << " Created on " << ctime(&time_now);
        Amn_file << std::setw(12) << num_bands << std::setw(12) << cal_num_kpts << std::setw(12) << num_wannier
                 << std::endl;
    }

    for (int ik = start_k_index; ik < (cal_num_kpts + start_k_index); ik++)
    {
        ModuleBase::ComplexMatrix Amn;
        unkdotW_A(psi_pw, wfcpw, ik, radial_in_q, Amn);

        if (GlobalV::MY_RANK == 0)
        {
            for (int iw = 0; iw < num_wannier; iw++)
            {
                for (int ib_w = 0; ib_w < num_bands; ib_w++)
                {
                    Amn_file << std::setw(5) << ib_w + 1 << std::setw(5) << iw + 1 << std::setw(5)
                            << ik + 1 - start_k_index
                            << std::setw(18) << std::showpoint << std::fixed << std::setprecision(12)
                            << Amn(ib_w, iw).real()
                            << std::setw(18) << std::showpoint << std::fixed << std::setprecision(12)
                            << Amn(ib_w, iw).imag()
                            // jingan test
                            //<< "   " << std::setw(18) << std::setprecision(13) << std::abs(Amn(ib_w, iw))
                            << std::endl;
                }
            }
        }
    }

    if (GlobalV::MY_RANK == 0) Amn_file.close();

}

void toW90_PW::out_unk(
    const psi::Psi<std::complex<double>>& psi_pw,
    const ModulePW::PW_Basis_K* wfcpw,
    const ModulePW::PW_Basis_Big* bigpw
)
{
    const bool wvfn_formatted = out_wannier_wvfn_formatted;

#ifdef __MPI
    // which_ip: found iz belongs to which ip.
    std::vector<int> which_ip(wfcpw->nz);
    ModuleBase::GlobalFunc::ZEROS(which_ip.data(), wfcpw->nz);
    
    for (int ip = 0; ip < GlobalV::NPROC_IN_POOL; ip++)
    {
        int iz = wfcpw->startz[ip];
        for (int index = 0; index < wfcpw->numz[ip]; index++)
        {
            which_ip[iz+index] = ip;
        }
    }

    // only do in the first pool.
    std::vector<std::complex<double>> porter(wfcpw->nrxx);
    int nxy = wfcpw->nx * wfcpw->ny;
    std::vector<std::complex<double>> zpiece(nxy);

    if (GlobalV::MY_POOL == 0)
    {
        for (int ik = start_k_index; ik < (cal_num_kpts + start_k_index); ik++)
        {
            std::ofstream unkfile;
            Binstream unkfile_b;
            
            if (GlobalV::RANK_IN_POOL == 0)
            {
                
                std::stringstream name;
                if (nspin_ == 1 || nspin_ == 4)
                {
                    name << PARAM.globalv.global_out_dir << "UNK" << std::setw(5) << std::setfill('0')
                        << ik + 1 << ".1";
                }
                else if (nspin_ == 2)
                {
                    if (wannier_spin == "up")
                        name << PARAM.globalv.global_out_dir << "UNK" << std::setw(5) << std::setfill('0')
                            << ik + 1 - start_k_index << ".1";
                    else if (wannier_spin == "down")
                        name << PARAM.globalv.global_out_dir << "UNK" << std::setw(5) << std::setfill('0')
                            << ik + 1 - start_k_index << ".2";
                }
                if (wvfn_formatted)
                {
                    unkfile.open(name.str(), std::ios::out);
                    unkfile << std::setw(12) << wfcpw->nx << std::setw(12) << wfcpw->ny << std::setw(12) << wfcpw->nz
                            << std::setw(12) << ik + 1 << std::setw(12) << num_bands << std::endl;
                }
                else
                {
                    unkfile_b.open(name.str(), "w");
                    unkfile_b << int(20) << wfcpw->nx << wfcpw->ny << wfcpw->nz << ik + 1 << num_bands << 20;
                }
            }

            for (int ib_w = 0; ib_w < num_bands; ib_w++)
            {
                int ib = cal_band_index[ib_w];

                wfcpw->recip2real(&psi_pw(ik, ib, 0), porter.data(), ik);

                if (GlobalV::RANK_IN_POOL == 0)
                {
                    if (!wvfn_formatted)
                    {
                        unkfile_b << wfcpw->nz * wfcpw->ny * wfcpw->nx * 8 * 2; // sizeof(double) = 8
                    }
                }

                // save the rho one z by one z.
                for (int iz = 0; iz < wfcpw->nz; iz++)
                {
                    // tag must be different for different iz.
                    std::fill(zpiece.begin(), zpiece.end(), std::complex<double>(0.0, 0.0));
                    int tag = iz;
                    MPI_Status ierror;

                    // case 1: the first part of rho in processor 0.
                    if (which_ip[iz] == 0 && GlobalV::RANK_IN_POOL == 0)
                    {
                        for (int ir = 0; ir < nxy; ir++)
                        {
                            zpiece[ir] = porter[ir * wfcpw->nplane + iz - wfcpw->startz_current];
                        }
                    }
                    // case 2: > first part rho: send the rho to
                    // processor 0.
                    else if (which_ip[iz] == GlobalV::RANK_IN_POOL)
                    {
                        for (int ir = 0; ir < nxy; ir++)
                        {
                            zpiece[ir] = porter[ir * wfcpw->nplane + iz - wfcpw->startz_current];
                        }
                        MPI_Send(zpiece.data(), nxy, MPI_DOUBLE_COMPLEX, 0, tag, POOL_WORLD);
                    }

                    // case 2: > first part rho: processor 0 receive the rho
                    // from other processors
                    else if (GlobalV::RANK_IN_POOL == 0)
                    {
                        MPI_Recv(zpiece.data(), nxy, MPI_DOUBLE_COMPLEX, which_ip[iz], tag, POOL_WORLD, &ierror);
                    }

                    // write data
                    if (GlobalV::RANK_IN_POOL == 0)
                    {
                        if (wvfn_formatted)
                        {
                            for (int iy = 0; iy < wfcpw->ny; iy++)
                            {
                                for (int ix = 0; ix < wfcpw->nx; ix++)
                                {
                                    unkfile << std::setw(20) << std::setprecision(9)
                                            << std::setiosflags(std::ios::scientific)
                                            << zpiece[ix * wfcpw->ny + iy].real()
                                            << std::setw(20) << std::setprecision(9)
                                            << std::setiosflags(std::ios::scientific)
                                            << zpiece[ix * wfcpw->ny + iy].imag()
                                            << std::endl;
                                }
                            }
                        }
                        else
                        {
                            for (int iy = 0; iy < wfcpw->ny; iy++)
                            {
                                for (int ix = 0; ix < wfcpw->nx; ix++)
                                {
                                    unkfile_b << zpiece[ix * wfcpw->ny + iy].real()
                                              << zpiece[ix * wfcpw->ny + iy].imag();
                                }
                            }
                        }
                    }
                } // end iz
                if (GlobalV::RANK_IN_POOL == 0)
                {
                    if (!wvfn_formatted)
                    {
                        unkfile_b << wfcpw->nz * wfcpw->ny * wfcpw->nx * 8 * 2; // sizeof(double) = 8
                    }
                }
                MPI_Barrier(POOL_WORLD);
            } // ib_w

            if (GlobalV::RANK_IN_POOL == 0)
            {
                if (wvfn_formatted)
                    unkfile.close();
                else
                    unkfile_b.close();
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

#endif

}
