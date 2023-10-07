#include "to_wannier90_pw.h"

#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_base/math_integral.h"
#include "module_base/math_polyint.h"
#include "module_base/math_sphbes.h"
#include "module_base/math_ylmreal.h"
#include "module_base/parallel_reduce.h"
#include "binstream.h"

toWannier90_PW::toWannier90_PW(
    const bool &out_wannier_mmn, 
    const bool &out_wannier_amn, 
    const bool &out_wannier_unk, 
    const bool &out_wannier_eig,
    const bool &out_wannier_wvfn_formatted, 
    const std::string &nnkpfile,
    const std::string &wannier_spin
):toWannier90(out_wannier_mmn, out_wannier_amn, out_wannier_unk, out_wannier_eig, out_wannier_wvfn_formatted, nnkpfile, wannier_spin)
{

}

toWannier90_PW::~toWannier90_PW()
{
    
}

void toWannier90_PW::calculate(
    const ModuleBase::matrix& ekb,
    const ModulePW::PW_Basis_K* wfcpw,
    const ModulePW::PW_Basis_Big* bigpw,
    const K_Vectors& kv,
    const psi::Psi<std::complex<double>>* psi
)
{
    read_nnkp(kv);

    if (GlobalV::NSPIN == 2)
    {
        if (wannier_spin == "up")
        {
            start_k_index = 0;
        }
        else if (wannier_spin == "down")
        {
            start_k_index = num_kpts / 2;
        }
        else
        {
            ModuleBase::WARNING_QUIT("toWannier90::calculate", "Error wannier_spin set,is not \"up\" or \"down\" ");
        }
    }

    if (out_wannier_eig)
    {
        out_eig(ekb);
    }

    if (out_wannier_mmn)
    {
        cal_Mmn(*psi, wfcpw);
    }

    if (out_wannier_amn)
    {
        gen_radial_function();
        cal_Amn(*psi, wfcpw);
    }

    if (out_wannier_unk)
    {
        out_unk(*psi, wfcpw, bigpw);
    }

}

void toWannier90_PW::cal_Mmn(
    const psi::Psi<std::complex<double>>& psi_pw,
    const ModulePW::PW_Basis_K* wfcpw
)
{
    std::ofstream mmn_file;

    if (GlobalV::MY_RANK == 0)
    {
        std::string fileaddress = GlobalV::global_out_dir + wannier_file_name + ".mmn";
        mmn_file.open(fileaddress.c_str(), std::ios::out);

        time_t time_now = time(NULL);
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
                        mmn_file << std::setw(18) << std::setprecision(12) << std::showpoint << std::fixed << Mmn(m, n).real()
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


void toWannier90_PW::cal_Amn(
    const psi::Psi<std::complex<double>>& psi_pw, 
    const ModulePW::PW_Basis_K* wfcpw
)
{
    const int pwNumberMax = wfcpw->npwk_max;

    std::ofstream Amn_file;

    if (GlobalV::MY_RANK == 0)
    {
        time_t time_now = time(NULL);
        std::string fileaddress = GlobalV::global_out_dir + wannier_file_name + ".amn";
        Amn_file.open(fileaddress.c_str(), std::ios::out);
        Amn_file << " Created on " << ctime(&time_now);
        Amn_file << std::setw(12) << num_bands << std::setw(12) << cal_num_kpts << std::setw(12) << num_wannier
                 << std::endl;
    }

    ModuleBase::ComplexMatrix *trial_orbitals = new ModuleBase::ComplexMatrix[cal_num_kpts];
    for (int ik = 0; ik < cal_num_kpts; ik++)
    {
        trial_orbitals[ik].create(num_wannier, pwNumberMax);
        produce_trial_in_pw(psi_pw, ik, wfcpw, trial_orbitals[ik]);
    }

    for (int ik = start_k_index; ik < (cal_num_kpts + start_k_index); ik++)
    {
        for (int iw = 0; iw < num_wannier; iw++)
        {
            int index_band = 0;
            for (int ib = 0; ib < GlobalV::NBANDS; ib++)
            {
                if (!tag_cal_band[ib]) continue;
                index_band++;
                std::complex<double> amn(0.0, 0.0);
                std::complex<double> amn_tem(0.0, 0.0);

                if (GlobalV::NSPIN != 4)
                {
                    for (int ig = 0; ig < pwNumberMax; ig++)
                    {
                        int cal_ik = ik - start_k_index;
                        amn_tem = amn_tem + conj(psi_pw(ik, ib, ig)) * trial_orbitals[cal_ik](iw, ig);
                    }
                }
                else
                {
                    for (int ig = 0; ig < pwNumberMax; ig++)
                    {
                        int cal_ik = ik - start_k_index;
                        amn_tem += up_con[iw] * (psi_pw(ik, ib, ig)) * trial_orbitals[cal_ik](iw, ig)
                                 + dn_con[iw] * (psi_pw(ik, ib, ig+pwNumberMax)) * trial_orbitals[cal_ik](iw, ig);
                    }
                }
#ifdef __MPI
                MPI_Allreduce(&amn_tem, &amn, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
#else
                amn = amn_tem;
#endif
                if (GlobalV::MY_RANK == 0)
                {
                    Amn_file << std::setw(5) << index_band << std::setw(5) << iw + 1 << std::setw(5)
                             << ik + 1 - start_k_index << std::setw(18) << std::showpoint << std::fixed << std::setprecision(12)
                             << amn.real() << std::setw(18) << std::showpoint << std::fixed << std::setprecision(12)
                             << amn.imag()
                             // jingan test
                             //<< "   " << std::setw(18) << std::setprecision(13) << std::abs(amn)
                             << std::endl;
                }
            }
        }
    }

    if (GlobalV::MY_RANK == 0) Amn_file.close();

    delete[] trial_orbitals;

}

void toWannier90_PW::out_unk(
    const psi::Psi<std::complex<double>>& psi_pw,
    const ModulePW::PW_Basis_K* wfcpw,
    const ModulePW::PW_Basis_Big* bigpw
)
{
    const bool wvfn_formatted = out_wannier_wvfn_formatted;

#ifdef __MPI
    // num_z: how many planes on processor 'ip'
    int *num_z = new int[GlobalV::NPROC_IN_POOL];
    ModuleBase::GlobalFunc::ZEROS(num_z, GlobalV::NPROC_IN_POOL);
    for (int iz = 0; iz < bigpw->nbz; iz++)
    {
        int ip = iz % GlobalV::NPROC_IN_POOL;
        num_z[ip] += bigpw->bz;
    }

    // start_z: start position of z in
    // processor ip.
    int *start_z = new int[GlobalV::NPROC_IN_POOL];
    ModuleBase::GlobalFunc::ZEROS(start_z, GlobalV::NPROC_IN_POOL);
    for (int ip = 1; ip < GlobalV::NPROC_IN_POOL; ip++)
    {
        start_z[ip] = start_z[ip - 1] + num_z[ip - 1];
    }

    // which_ip: found iz belongs to which ip.
    int* which_ip = new int[wfcpw->nz];
    ModuleBase::GlobalFunc::ZEROS(which_ip, wfcpw->nz);
    for (int iz = 0; iz < wfcpw->nz; iz++)
    {
        for (int ip = 0; ip < GlobalV::NPROC_IN_POOL; ip++)
        {
            if (iz >= start_z[GlobalV::NPROC_IN_POOL - 1])
            {
                which_ip[iz] = GlobalV::NPROC_IN_POOL - 1;
                break;
            }
            else if (iz >= start_z[ip] && iz < start_z[ip + 1])
            {
                which_ip[iz] = ip;
                break;
            }
        }
    }

    // only do in the first pool.
    std::complex<double>* porter = new std::complex<double>[wfcpw->nrxx];
    int nxy = wfcpw->nx * wfcpw->ny;
    std::complex<double> *zpiece = new std::complex<double>[nxy];

    if (GlobalV::MY_POOL == 0)
    {
        for (int ik = start_k_index; ik < (cal_num_kpts + start_k_index); ik++)
        {
            std::ofstream unkfile;
            Binstream unkfile_b;
            
            if (GlobalV::MY_RANK == 0)
            {
                
                std::stringstream name;
                if (GlobalV::NSPIN == 1 || GlobalV::NSPIN == 4)
                {
                    name << GlobalV::global_out_dir << "UNK" << std::setw(5) << std::setfill('0') << ik + 1 << ".1";
                }
                else if (GlobalV::NSPIN == 2)
                {
                    if (wannier_spin == "up")
                        name << GlobalV::global_out_dir << "UNK" << std::setw(5) << std::setfill('0')
                            << ik + 1 - start_k_index << ".1";
                    else if (wannier_spin == "down")
                        name << GlobalV::global_out_dir << "UNK" << std::setw(5) << std::setfill('0')
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

            for (int ib = 0; ib < GlobalV::NBANDS; ib++)
            {
                if (!tag_cal_band[ib])
                    continue;

                wfcpw->recip2real(&psi_pw(ik, ib, 0), porter, ik);

                if (GlobalV::MY_RANK == 0)
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
                    ModuleBase::GlobalFunc::ZEROS(zpiece, nxy);
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
                        MPI_Send(zpiece, nxy, MPI_DOUBLE_COMPLEX, 0, tag, POOL_WORLD);
                    }

                    // case 2: > first part rho: processor 0 receive the rho
                    // from other processors
                    else if (GlobalV::RANK_IN_POOL == 0)
                    {
                        MPI_Recv(zpiece, nxy, MPI_DOUBLE_COMPLEX, which_ip[iz], tag, POOL_WORLD, &ierror);
                    }

                    // write data
                    if (GlobalV::MY_RANK == 0)
                    {
                        if (wvfn_formatted)
                        {
                            for (int iy = 0; iy < wfcpw->ny; iy++)
                            {
                                for (int ix = 0; ix < wfcpw->nx; ix++)
                                {
                                    unkfile << std::setw(20) << std::setprecision(9) << std::setiosflags(std::ios::scientific)
                                            << zpiece[ix * wfcpw->ny + iy].real() << std::setw(20) << std::setprecision(9)
                                            << std::setiosflags(std::ios::scientific) << zpiece[ix * wfcpw->ny + iy].imag()
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
                                    unkfile_b << zpiece[ix * wfcpw->ny + iy].real() << zpiece[ix * wfcpw->ny + iy].imag();
                                }
                            }
                        }
                    }
                } // end iz
                if (GlobalV::MY_RANK == 0)
                {
                    if (!wvfn_formatted)
                    {
                        unkfile_b << wfcpw->nz * wfcpw->ny * wfcpw->nx * 8 * 2; // sizeof(double) = 8
                    }
                }
                MPI_Barrier(POOL_WORLD);
            } // ib

            if (GlobalV::MY_RANK == 0)
            {
                if (wvfn_formatted)
                    unkfile.close();
                else
                    unkfile_b.close();
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    delete[] num_z;
    delete[] start_z;
    delete[] which_ip;
    delete[] porter;
    delete[] zpiece;

#endif

}


void toWannier90_PW::unkdotkb(
    const psi::Psi<std::complex<double>>& psi_pw, 
    const ModulePW::PW_Basis_K* wfcpw,
    const int& cal_ik,
    const int& cal_ikb,
    const ModuleBase::Vector3<double> G,
    ModuleBase::ComplexMatrix &Mmn
)
{
    Mmn.create(num_bands, num_bands);

    int count_m = -1;
    for (int m = 0; m < GlobalV::NBANDS; m++)
    {
        if (!tag_cal_band[m]) continue;
        count_m++;

        std::complex<double>* phase = new std::complex<double>[wfcpw->nmaxgr];
        ModuleBase::GlobalFunc::ZEROS(phase, wfcpw->nmaxgr);

        // get the phase value in realspace
        for (int ig = 0; ig < wfcpw->npwk[cal_ik]; ig++)
        {
            // if wfcpw->getgdirect(cal_ik, ig) == phase_G
            ModuleBase::Vector3<double> temp_G = wfcpw->getgdirect(cal_ik, ig) - G;
            if (temp_G.norm() < 1e-9)
            {
                phase[ig] = std::complex<double>(1.0, 0.0);
                break;
            }
        }

        wfcpw->recip2real(phase, phase, cal_ik);

        if (GlobalV::NSPIN == 4)
        {
            // (1) set value
            std::complex<double>* psir_up = new std::complex<double>[wfcpw->nmaxgr];
            std::complex<double>* psir_dn = new std::complex<double>[wfcpw->nmaxgr];
            ModuleBase::GlobalFunc::ZEROS(psir_up, wfcpw->nmaxgr);
            ModuleBase::GlobalFunc::ZEROS(psir_dn, wfcpw->nmaxgr);

            // (2) fft and get value
            int npw_ik = wfcpw->npwk[cal_ik];
            wfcpw->recip2real(&psi_pw(cal_ik, m, 0), psir_up, cal_ik);
            wfcpw->recip2real(&psi_pw(cal_ik, m, npw_ik), psir_dn, cal_ik);
            for (int ir = 0; ir < wfcpw->nrxx; ir++)
            {
                psir_up[ir] *= phase[ir];
                psir_dn[ir] *= phase[ir];
            }

            wfcpw->real2recip(psir_up, psir_up, cal_ikb);
            wfcpw->real2recip(psir_dn, psir_dn, cal_ikb);

            int count_n = -1;
            for (int n = 0; n < GlobalV::NBANDS; n++)
            {
                if (!tag_cal_band[n]) continue;
                count_n++;

                if (!gamma_only_wannier)
                {
                    std::complex<double> result_tem(0.0, 0.0);

                    // int npw_ikb = wfcpw->npwk[cal_ikb];
                    int pwNumberMax = wfcpw->npwk_max;

                    // Can be accelerated using lapack
                    for (int ig = 0; ig < pwNumberMax; ig++)
                    {
                        result_tem = result_tem + conj(psir_up[ig]) * psi_pw(cal_ikb, n, ig) + conj(psir_dn[ig]) * psi_pw(cal_ikb, n, ig+pwNumberMax);
                    }
#ifdef __MPI
                    Parallel_Reduce::reduce_complex_double_all(result_tem);
#endif
                    Mmn(count_m, count_n) = result_tem;
                }
                else
                {
                    // GlobalV::ofs_running << "gamma only test" << std::endl;
                }

            }

            delete[] psir_up;
            delete[] psir_dn;

        }
        else
        {
            // (1) set value
            std::complex<double>* psir = new std::complex<double>[wfcpw->nmaxgr];
            ModuleBase::GlobalFunc::ZEROS(psir, wfcpw->nmaxgr);

            // (2) fft and get value
            wfcpw->recip2real(&psi_pw(cal_ik, m, 0), psir, cal_ik);
            for (int ir = 0; ir < wfcpw->nrxx; ir++)
            {
                psir[ir] *= phase[ir];
            }

            wfcpw->real2recip(psir, psir, cal_ikb);

            int count_n = -1;
            for (int n = 0; n < GlobalV::NBANDS; n++)
            {
                if (!tag_cal_band[n]) continue;
                count_n++;

                if (!gamma_only_wannier)
                {
                    std::complex<double> result_tem(0.0, 0.0);

                    // Can be accelerated using lapack
                    for (int ig = 0; ig < wfcpw->npwk[cal_ikb]; ig++)
                    {
                        result_tem = result_tem + conj(psir[ig]) * psi_pw(cal_ikb, n, ig);
                    }

#ifdef __MPI
                    Parallel_Reduce::reduce_complex_double_all(result_tem);
#endif
                    Mmn(count_m, count_n) = result_tem;
                }
                else
                {
                    // GlobalV::ofs_running << "gamma only test" << std::endl;
                }

            }

            delete[] psir;
            
        }

        delete[] phase;

    }

}

void toWannier90_PW::gen_radial_function()
{
    r.create(num_wannier, mesh_r);
    dr.create(num_wannier, mesh_r);
    psi.create(num_wannier, mesh_r);
    psir.create(num_wannier, mesh_r);

    for (int i = 0; i < num_wannier; i++)
    {
        double x = 0;
        for (int ir = 0; ir < mesh_r; ir++)
        {
            x = x_min + ir * dx;
            r(i, ir) = exp(x) / alfa[i];
            dr(i, ir) = dx * r(i, ir);
        }
    }

    for (int i = 0; i < num_wannier; i++)
    {
        double alfa32 = pow(alfa[i], 3.0 / 2.0);
        double alfa_new = alfa[i];
        int wannier_index = i;

        if (rvalue[i] == 1)
        {
            for (int ir = 0; ir < mesh_r; ir++)
            {
                psi(wannier_index, ir) = 2.0 * alfa32 * exp(-alfa_new * r(wannier_index, ir));
            }
        }

        if (rvalue[i] == 2)
        {
            for (int ir = 0; ir < mesh_r; ir++)
            {
                psi(wannier_index, ir) = 1.0 / sqrt(8.0) * alfa32 * (2.0 - alfa_new * r(wannier_index, ir))
                                         * exp(-alfa_new * r(wannier_index, ir) * 0.5);
            }
        }

        if (rvalue[i] == 3)
        {
            for (int ir = 0; ir < mesh_r; ir++)
            {
                psi(wannier_index, ir)
                    = sqrt(4.0 / 27.0) * alfa32
                      * (1.0 - 2.0 / 3.0 * alfa_new * r(wannier_index, ir)
                         + 2.0 / 27.0 * pow(alfa_new, 2.0) * r(wannier_index, ir) * r(wannier_index, ir))
                      * exp(-alfa_new * r(wannier_index, ir) * 1.0 / 3.0);
            }
        }
    }

    for (int i = 0; i < num_wannier; i++)
    {
        for (int ir = 0; ir < mesh_r; ir++)
        {
            psir(i, ir) = psi(i, ir) * r(i, ir);
        }
    }

}

void toWannier90_PW::produce_trial_in_pw(
    const psi::Psi<std::complex<double>>& psi_pw,
    const int& ik,
    const ModulePW::PW_Basis_K* wfcpw,
    ModuleBase::ComplexMatrix& trial_orbitals_k
)
{
    const int npw = wfcpw->npwk[ik];
    const int npwx = wfcpw->npwk_max;
    const int total_lm = 16;
    ModuleBase::matrix ylm(total_lm, npw);

    double bs2, bs3, bs6, bs12;
    bs2 = 1.0 / sqrt(2.0);
    bs3 = 1.0 / sqrt(3.0);
    bs6 = 1.0 / sqrt(6.0);
    bs12 = 1.0 / sqrt(12.0);

    ModuleBase::Vector3<double> *gk = new ModuleBase::Vector3<double>[npw];
    for (int ig = 0; ig < npw; ig++)
    {
        gk[ig] = wfcpw->getgpluskcar(ik,ig);
    }

    ModuleBase::YlmReal::Ylm_Real(total_lm, npw, gk, ylm);

    for (int wannier_index = 0; wannier_index < num_wannier; wannier_index++)
    {
        if (L[wannier_index] >= 0)
        {
            get_trial_orbitals_lm_k(wannier_index,
                                    L[wannier_index],
                                    m[wannier_index],
                                    ylm,
                                    dr,
                                    r,
                                    psir,
                                    mesh_r,
                                    gk,
                                    npw,
                                    npwx,
                                    trial_orbitals_k);
        }
        else
        {
            if (L[wannier_index] == -1 && m[wannier_index] == 0)
            {
                get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 1, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                for (int ig = 0; ig < npwx; ig++)
                {
                    trial_orbitals_k(wannier_index, ig)
                        = bs2 * tem_array[ig] + bs2 * trial_orbitals_k(wannier_index, ig);
                }
                delete[] tem_array;
            }
            else if (L[wannier_index] == -1 && m[wannier_index] == 1)
            {
                get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 1, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                for (int ig = 0; ig < npwx; ig++)
                {
                    trial_orbitals_k(wannier_index, ig)
                        = bs2 * tem_array[ig] - bs2 * trial_orbitals_k(wannier_index, ig);
                }
                delete[] tem_array;
            }
            else if (L[wannier_index] == -2 && m[wannier_index] == 0)
            {
                get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_1 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_1[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 1, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_2 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_2[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 2, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                for (int ig = 0; ig < npwx; ig++)
                {
                    trial_orbitals_k(wannier_index, ig)
                        = bs3 * tem_array_1[ig] - bs6 * tem_array_2[ig] + bs2 * trial_orbitals_k(wannier_index, ig);
                }
                delete[] tem_array_1;
                delete[] tem_array_2;
            }
            else if (L[wannier_index] == -2 && m[wannier_index] == 1)
            {
                get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_1 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_1[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 1, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_2 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_2[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 2, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                for (int ig = 0; ig < npwx; ig++)
                {
                    trial_orbitals_k(wannier_index, ig)
                        = bs3 * tem_array_1[ig] - bs6 * tem_array_2[ig] - bs2 * trial_orbitals_k(wannier_index, ig);
                }
                delete[] tem_array_1;
                delete[] tem_array_2;
            }
            else if (L[wannier_index] == -2 && m[wannier_index] == 2)
            {
                get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 1, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                for (int ig = 0; ig < npwx; ig++)
                {
                    trial_orbitals_k(wannier_index, ig)
                        = bs3 * tem_array[ig] + 2 * bs6 * trial_orbitals_k(wannier_index, ig);
                }
                delete[] tem_array;
            }
            else if (L[wannier_index] == -3 && m[wannier_index] == 0)
            {
                get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_1 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_1[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 1, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_2 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_2[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 2, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_3 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_3[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                for (int ig = 0; ig < npwx; ig++)
                {
                    trial_orbitals_k(wannier_index, ig)
                        = 0.5
                          * (tem_array_1[ig] + tem_array_2[ig] + tem_array_3[ig] + trial_orbitals_k(wannier_index, ig));
                }
                delete[] tem_array_1;
                delete[] tem_array_2;
                delete[] tem_array_3;
            }
            else if (L[wannier_index] == -3 && m[wannier_index] == 1)
            {
                get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_1 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_1[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 1, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_2 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_2[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 2, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_3 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_3[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                for (int ig = 0; ig < npwx; ig++)
                {
                    trial_orbitals_k(wannier_index, ig)
                        = 0.5
                          * (tem_array_1[ig] + tem_array_2[ig] - tem_array_3[ig] - trial_orbitals_k(wannier_index, ig));
                }
                delete[] tem_array_1;
                delete[] tem_array_2;
                delete[] tem_array_3;
            }
            else if (L[wannier_index] == -3 && m[wannier_index] == 2)
            {
                get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_1 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_1[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 1, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_2 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_2[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 2, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_3 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_3[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                for (int ig = 0; ig < npwx; ig++)
                {
                    trial_orbitals_k(wannier_index, ig)
                        = 0.5
                          * (tem_array_1[ig] - tem_array_2[ig] + tem_array_3[ig] - trial_orbitals_k(wannier_index, ig));
                }
                delete[] tem_array_1;
                delete[] tem_array_2;
                delete[] tem_array_3;
            }
            else if (L[wannier_index] == -3 && m[wannier_index] == 3)
            {
                get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_1 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_1[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 1, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_2 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_2[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 2, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_3 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_3[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                for (int ig = 0; ig < npwx; ig++)
                {
                    trial_orbitals_k(wannier_index, ig)
                        = 0.5
                          * (tem_array_1[ig] - tem_array_2[ig] - tem_array_3[ig] + trial_orbitals_k(wannier_index, ig));
                }
                delete[] tem_array_1;
                delete[] tem_array_2;
                delete[] tem_array_3;
            }
            else if (L[wannier_index] == -4 && m[wannier_index] == 0)
            {
                get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_1 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_1[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 1, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_2 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_2[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 2, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                for (int ig = 0; ig < npwx; ig++)
                {
                    trial_orbitals_k(wannier_index, ig)
                        = bs3 * tem_array_1[ig] - bs6 * tem_array_2[ig] + bs2 * trial_orbitals_k(wannier_index, ig);
                }
                delete[] tem_array_1;
                delete[] tem_array_2;
            }
            else if (L[wannier_index] == -4 && m[wannier_index] == 1)
            {
                get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_1 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_1[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 1, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_2 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_2[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 2, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                for (int ig = 0; ig < npwx; ig++)
                {
                    trial_orbitals_k(wannier_index, ig)
                        = bs3 * tem_array_1[ig] - bs6 * tem_array_2[ig] - bs2 * trial_orbitals_k(wannier_index, ig);
                }
                delete[] tem_array_1;
                delete[] tem_array_2;
            }
            else if (L[wannier_index] == -4 && m[wannier_index] == 2)
            {
                get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_1 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_1[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 1, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                for (int ig = 0; ig < npwx; ig++)
                {
                    trial_orbitals_k(wannier_index, ig)
                        = bs3 * tem_array_1[ig] - 2 * bs6 * trial_orbitals_k(wannier_index, ig);
                }
                delete[] tem_array_1;
            }
            else if (L[wannier_index] == -4 && m[wannier_index] == 3)
            {
                get_trial_orbitals_lm_k(wannier_index, 1, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_1 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_1[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 2, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                for (int ig = 0; ig < npwx; ig++)
                {
                    trial_orbitals_k(wannier_index, ig)
                        = bs2 * tem_array_1[ig] + bs2 * trial_orbitals_k(wannier_index, ig);
                }
                delete[] tem_array_1;
            }
            else if (L[wannier_index] == -4 && m[wannier_index] == 4)
            {
                get_trial_orbitals_lm_k(wannier_index, 1, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_1 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_1[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 2, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                for (int ig = 0; ig < npwx; ig++)
                {
                    trial_orbitals_k(wannier_index, ig)
                        = -1.0 * bs2 * tem_array_1[ig] + bs2 * trial_orbitals_k(wannier_index, ig);
                }
                delete[] tem_array_1;
            }
            else if (L[wannier_index] == -5 && m[wannier_index] == 0)
            {
                get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_1 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_1[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 1, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_2 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_2[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 2, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_3 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_3[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 2, 3, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                for (int ig = 0; ig < npwx; ig++)
                {
                    trial_orbitals_k(wannier_index, ig) = bs6 * tem_array_1[ig] - bs2 * tem_array_2[ig]
                                                          - bs12 * tem_array_3[ig]
                                                          + 0.5 * trial_orbitals_k(wannier_index, ig);
                }
                delete[] tem_array_1;
                delete[] tem_array_2;
                delete[] tem_array_3;
            }
            else if (L[wannier_index] == -5 && m[wannier_index] == 1)
            {
                get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_1 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_1[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 1, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_2 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_2[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 2, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_3 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_3[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 2, 3, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                for (int ig = 0; ig < npwx; ig++)
                {
                    trial_orbitals_k(wannier_index, ig) = bs6 * tem_array_1[ig] + bs2 * tem_array_2[ig]
                                                          - bs12 * tem_array_3[ig]
                                                          + 0.5 * trial_orbitals_k(wannier_index, ig);
                }
                delete[] tem_array_1;
                delete[] tem_array_2;
                delete[] tem_array_3;
            }
            else if (L[wannier_index] == -5 && m[wannier_index] == 2)
            {
                get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_1 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_1[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 2, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_2 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_2[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 2, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_3 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_3[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 2, 3, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                for (int ig = 0; ig < npwx; ig++)
                {
                    trial_orbitals_k(wannier_index, ig) = bs6 * tem_array_1[ig] - bs2 * tem_array_2[ig]
                                                          - bs12 * tem_array_3[ig]
                                                          - 0.5 * trial_orbitals_k(wannier_index, ig);
                }
                delete[] tem_array_1;
                delete[] tem_array_2;
                delete[] tem_array_3;
            }
            else if (L[wannier_index] == -5 && m[wannier_index] == 3)
            {
                get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_1 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_1[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 2, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_2 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_2[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 2, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_3 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_3[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 2, 3, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                for (int ig = 0; ig < npwx; ig++)
                {
                    trial_orbitals_k(wannier_index, ig) = bs6 * tem_array_1[ig] + bs2 * tem_array_2[ig]
                                                          - bs12 * tem_array_3[ig]
                                                          - 0.5 * trial_orbitals_k(wannier_index, ig);
                }
                delete[] tem_array_1;
                delete[] tem_array_2;
                delete[] tem_array_3;
            }
            else if (L[wannier_index] == -5 && m[wannier_index] == 4)
            {
                get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_1 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_1[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_2 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_2[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 2, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                for (int ig = 0; ig < npwx; ig++)
                {
                    trial_orbitals_k(wannier_index, ig)
                        = bs6 * tem_array_1[ig] - bs2 * tem_array_2[ig] + bs3 * trial_orbitals_k(wannier_index, ig);
                }
                delete[] tem_array_1;
                delete[] tem_array_2;
            }
            else if (L[wannier_index] == -5 && m[wannier_index] == 5)
            {
                get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_1 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_1[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_2 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_2[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 2, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                for (int ig = 0; ig < npwx; ig++)
                {
                    trial_orbitals_k(wannier_index, ig)
                        = bs6 * tem_array_1[ig] + bs2 * tem_array_2[ig] + bs3 * trial_orbitals_k(wannier_index, ig);
                }
                delete[] tem_array_1;
                delete[] tem_array_2;
            }
        }
    }
    delete[] gk;
}

void toWannier90_PW::get_trial_orbitals_lm_k(
    const int wannier_index,
    const int orbital_L,
    const int orbital_m,
    ModuleBase::matrix &ylm,
    ModuleBase::matrix &dr,
    ModuleBase::matrix &r,
    ModuleBase::matrix &psir,
    const int mesh_r,
    ModuleBase::Vector3<double> *gk,
    const int npw,
    const int npwx,
    ModuleBase::ComplexMatrix &trial_orbitals_k
)
{

    double *psik = new double[npw];
    double *psir_tem = new double[mesh_r];
    double *r_tem = new double[mesh_r];
    double *dr_tem = new double[mesh_r];
    double *psik_tem = new double[GlobalV::NQX];
    ModuleBase::GlobalFunc::ZEROS(psir_tem, mesh_r);
    ModuleBase::GlobalFunc::ZEROS(r_tem, mesh_r);
    ModuleBase::GlobalFunc::ZEROS(dr_tem, mesh_r);

    for (int ir = 0; ir < mesh_r; ir++)
    {
        psir_tem[ir] = psir(wannier_index, ir);
        r_tem[ir] = r(wannier_index, ir);
        dr_tem[ir] = dr(wannier_index, ir);
    }

    this->integral(mesh_r, psir_tem, r_tem, dr_tem, orbital_L, psik_tem);

    for (int ig = 0; ig < npw; ig++)
    {
        psik[ig] = ModuleBase::PolyInt::Polynomial_Interpolation(psik_tem, GlobalV::NQX, GlobalV::DQ, gk[ig].norm() * GlobalC::ucell.tpiba);
    }

    std::complex<double> *sk = new std::complex<double>[npw];
    for (int ig = 0; ig < npw; ig++)
    {
        const double arg = (gk[ig] * R_centre[wannier_index]) * ModuleBase::TWO_PI;
        sk[ig] = std::complex<double>(cos(arg), -sin(arg));
    }

    double *wannier_ylm = new double[npw];
    for (int ig = 0; ig < npw; ig++)
    {
        int index = orbital_L * orbital_L + orbital_m;
        if (index == 2 || index == 3 || index == 5 || index == 6 || index == 14 || index == 15)
        {
            wannier_ylm[ig] = -1 * ylm(index, ig);
        }
        else
        {
            wannier_ylm[ig] = ylm(index, ig);
        }
    }

    std::complex<double> lphase = pow(ModuleBase::NEG_IMAG_UNIT, orbital_L);
    for (int ig = 0; ig < npwx; ig++)
    {
        if (ig < npw)
        {
            trial_orbitals_k(wannier_index, ig) = lphase * sk[ig] * wannier_ylm[ig] * psik[ig];
        }
        else
            trial_orbitals_k(wannier_index, ig) = std::complex<double>(0.0, 0.0);
    }

    std::complex<double> anorm(0.0, 0.0);
    for (int ig = 0; ig < npwx; ig++)
    {
        anorm = anorm + conj(trial_orbitals_k(wannier_index, ig)) * trial_orbitals_k(wannier_index, ig);
    }

    std::complex<double> anorm_tem(0.0, 0.0);
#ifdef __MPI
    MPI_Allreduce(&anorm, &anorm_tem, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, POOL_WORLD);
#else
    anorm_tem = anorm;
#endif

    for (int ig = 0; ig < npwx; ig++)
    {
        trial_orbitals_k(wannier_index, ig) = trial_orbitals_k(wannier_index, ig) / sqrt(anorm_tem);
    }

    delete[] psik;
    delete[] psir_tem;
    delete[] r_tem;
    delete[] dr_tem;
    delete[] psik_tem;
    delete[] sk;
    delete[] wannier_ylm;

    return;
}

void toWannier90_PW::integral(
    const int meshr,
    const double *psir,
    const double *r,
    const double *rab,
    const int &l,
    double *table
)
{
    const double pref = ModuleBase::FOUR_PI / sqrt(GlobalC::ucell.omega);

    double *inner_part = new double[meshr];
    for (int ir = 0; ir < meshr; ir++)
    {
        inner_part[ir] = psir[ir] * psir[ir];
    }

    double unit = 0.0;
    ModuleBase::Integral::Simpson_Integral(meshr, inner_part, rab, unit);
    delete[] inner_part;

    double *aux = new double[meshr];
    double *vchi = new double[meshr];
    for (int iq = 0; iq < GlobalV::NQX; iq++)
    {
        const double q = GlobalV::DQ * iq;
        ModuleBase::Sphbes::Spherical_Bessel(meshr, r, q, l, aux);
        for (int ir = 0; ir < meshr; ir++)
        {
            vchi[ir] = psir[ir] * aux[ir] * r[ir];
        }

        double vqint = 0.0;
        ModuleBase::Integral::Simpson_Integral(meshr, vchi, rab, vqint);

        table[iq] = vqint * pref;
    }
    delete[] aux;
    delete[] vchi;
    return;
}