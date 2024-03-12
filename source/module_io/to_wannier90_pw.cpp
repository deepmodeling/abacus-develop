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
    std::vector<ModuleBase::matrix> radial_in_q;
    gen_radial_function_in_q(radial_in_q);

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
                            << ik + 1 - start_k_index << std::setw(18) << std::showpoint << std::fixed << std::setprecision(12)
                            << Amn(ib_w, iw).real() << std::setw(18) << std::showpoint << std::fixed << std::setprecision(12)
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

void toWannier90_PW::out_unk(
    const psi::Psi<std::complex<double>>& psi_pw,
    const ModulePW::PW_Basis_K* wfcpw,
    const ModulePW::PW_Basis_Big* bigpw
)
{
    const bool wvfn_formatted = out_wannier_wvfn_formatted;

#ifdef __MPI
    // which_ip: found iz belongs to which ip.
    int* which_ip = new int[wfcpw->nz];
    ModuleBase::GlobalFunc::ZEROS(which_ip, wfcpw->nz);
    
    for (int ip = 0; ip < GlobalV::NPROC_IN_POOL; ip++)
    {
        int iz = wfcpw->startz[ip];
        for (int index = 0; index < wfcpw->numz[ip]; index++)
        {
            which_ip[iz+index] = ip;
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
            
            if (GlobalV::RANK_IN_POOL == 0)
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

            for (int ib_w = 0; ib_w < num_bands; ib_w++)
            {
                int ib = cal_band_index[ib_w];

                wfcpw->recip2real(&psi_pw(ik, ib, 0), porter, ik);

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
                    if (GlobalV::RANK_IN_POOL == 0)
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

    for (int m = 0; m < num_bands; m++)
    {
        int im = cal_band_index[m];
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
            // int npw_ik = wfcpw->npwk[cal_ik];
            int npwx = wfcpw->npwk_max;
            wfcpw->recip2real(&psi_pw(cal_ik, im, 0), psir_up, cal_ik);
            // wfcpw->recip2real(&psi_pw(cal_ik, im, npw_ik), psir_dn, cal_ik);
            wfcpw->recip2real(&psi_pw(cal_ik, im, npwx), psir_dn, cal_ik);
            for (int ir = 0; ir < wfcpw->nrxx; ir++)
            {
                psir_up[ir] *= phase[ir];
                psir_dn[ir] *= phase[ir];
            }

            wfcpw->real2recip(psir_up, psir_up, cal_ikb);
            wfcpw->real2recip(psir_dn, psir_dn, cal_ikb);

            for (int n = 0; n < num_bands; n++)
            {
                int in = cal_band_index[n];

                if (!gamma_only_wannier)
                {
                    std::complex<double> result_tem(0.0, 0.0);

                    // int npw_ikb = wfcpw->npwk[cal_ikb];
                    int pwNumberMax = wfcpw->npwk_max;

                    // Can be accelerated using lapack
                    for (int ig = 0; ig < pwNumberMax; ig++)
                    {
                        result_tem = result_tem + conj(psir_up[ig]) * psi_pw(cal_ikb, in, ig) + conj(psir_dn[ig]) * psi_pw(cal_ikb, in, ig+pwNumberMax);
                    }
#ifdef __MPI
                    Parallel_Reduce::reduce_all(result_tem);
#endif
                    Mmn(m, n) = result_tem;
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
            wfcpw->recip2real(&psi_pw(cal_ik, im, 0), psir, cal_ik);
            for (int ir = 0; ir < wfcpw->nrxx; ir++)
            {
                psir[ir] *= phase[ir];
            }

            wfcpw->real2recip(psir, psir, cal_ikb);

            for (int n = 0; n < num_bands; n++)
            {
                int in = cal_band_index[n];

                if (!gamma_only_wannier)
                {
                    std::complex<double> result_tem(0.0, 0.0);

                    // Can be accelerated using lapack
                    for (int ig = 0; ig < wfcpw->npwk[cal_ikb]; ig++)
                    {
                        result_tem = result_tem + conj(psir[ig]) * psi_pw(cal_ikb, in, ig);
                    }

#ifdef __MPI
                    Parallel_Reduce::reduce_all(result_tem);
#endif
                    Mmn(m, n) = result_tem;
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

void toWannier90_PW::gen_radial_function_in_q(std::vector<ModuleBase::matrix> &radial_in_q)
{
    // 径向函数傅里叶变换到q空间中
    radial_in_q.resize(num_wannier);

    double *r = new double[mesh_r];
    double *dr = new double[mesh_r];
    double *psi = new double[mesh_r];
    double *psir = new double[mesh_r];

    for (int wannier_index = 0; wannier_index < num_wannier; wannier_index++)
    {
        double x = 0;
        for (int ir = 0; ir < mesh_r; ir++)
        {
            x = x_min + ir * dx;
            r[ir] = exp(x) / alfa[wannier_index];
            dr[ir] = dx * r[ir];
        }

        double alfa32 = pow(alfa[wannier_index], 3.0 / 2.0);
        double alfa_new = alfa[wannier_index];

        if (rvalue[wannier_index] == 1)
        {
            for (int ir = 0; ir < mesh_r; ir++)
            {
                psi[ir] = 2.0 * alfa32 * exp(-alfa_new * r[ir]);
            }
        }

        if (rvalue[wannier_index] == 2)
        {
            for (int ir = 0; ir < mesh_r; ir++)
            {
                psi[ir] = 1.0 / sqrt(8.0) * alfa32 * (2.0 - alfa_new * r[ir]) * exp(-alfa_new * r[ir] * 0.5);
            }
        }

        if (rvalue[wannier_index] == 3)
        {
            for (int ir = 0; ir < mesh_r; ir++)
            {
                psi[ir] = sqrt(4.0 / 27.0) * alfa32
                        * (1.0 - 2.0 / 3.0 * alfa_new * r[ir] + 2.0 / 27.0 * pow(alfa_new, 2.0) * r[ir] * r[ir])
                        * exp(-alfa_new * r[ir] * 1.0 / 3.0);
            }
        }

        for (int ir = 0; ir < mesh_r; ir++)
        {
            psir[ir] = psi[ir] * r[ir];
        }

        auto &tmp_radial = radial_in_q[wannier_index];
        if (L[wannier_index] >= 0)
        {
            tmp_radial.create(1, GlobalV::NQX);
            integral(mesh_r, psir, r, dr, L[wannier_index], tmp_radial.c);
        }
        else
        {
            int tmp_size = 0;

            if (L[wannier_index] == -1 || L[wannier_index] == -2 || L[wannier_index] == -3) tmp_size = 2;

            if (L[wannier_index] == -4 || L[wannier_index] == -5) tmp_size = 3;

            tmp_radial.create(tmp_size, GlobalV::NQX);

            for (int tmp_L = 0; tmp_L < tmp_size; tmp_L++)
            {
                integral(mesh_r, psir, r, dr, tmp_L, tmp_radial.c+tmp_L*GlobalV::NQX);
            }
        }

    }

    delete[] r;
    delete[] dr;
    delete[] psi;
    delete[] psir;

}

void toWannier90_PW::produce_trial_in_pw(
    const psi::Psi<std::complex<double>>& psi_pw,
    const int& ik,
    const ModulePW::PW_Basis_K* wfcpw,
    const std::vector<ModuleBase::matrix> &radial_in_q,
    ModuleBase::ComplexMatrix& trial_orbitals_k
)
{
    const int npw = wfcpw->npwk[ik];
    const int npwx = wfcpw->npwk_max;

    trial_orbitals_k.create(num_wannier, npwx);

    const int total_lm = 16;
    ModuleBase::matrix ylm(total_lm, npw);

    ModuleBase::Vector3<double> *gk = new ModuleBase::Vector3<double>[npw];
    for (int ig = 0; ig < npw; ig++)
    {
        gk[ig] = wfcpw->getgpluskcar(ik, ig);
    }

    ModuleBase::YlmReal::Ylm_Real(total_lm, npw, gk, ylm);

    // 保持与Wannier90球谐函数定义一致
    std::vector<int> need_inv = {2, 3, 5, 6, 14, 15};
    for (auto index : need_inv)
    {
        for (int ig = 0; ig < npw; ig++)
        {
            ylm(index, ig) = -1.0 * ylm(index, ig);
        }
    }

    double bs2, bs3, bs6, bs12;
    bs2 = 1.0 / sqrt(2.0);
    bs3 = 1.0 / sqrt(3.0);
    bs6 = 1.0 / sqrt(6.0);
    bs12 = 1.0 / sqrt(12.0);

    std::complex<double> *sf = new std::complex<double>[npw];
    for (int wannier_index = 0; wannier_index < num_wannier; wannier_index++)
    {
        if (L[wannier_index] >= 0)
        {
            get_trial_orbitals_lm_k(L[wannier_index], m[wannier_index], ylm, gk, npw, radial_in_q[wannier_index].c, trial_orbitals_k.c + wannier_index*npwx);

            for (int ig = 0; ig < npw; ig++)
            {
                const double arg = (gk[ig] * R_centre[wannier_index]) * ModuleBase::TWO_PI;
                sf[ig] = std::complex<double>(cos(arg), -sin(arg));

                trial_orbitals_k(wannier_index, ig) = sf[ig] * trial_orbitals_k(wannier_index, ig);
            }

        }
        else
        {
            if (L[wannier_index] == -1)
            {
                double tmp_bs2 = 0;
                if (m[wannier_index] == 0) tmp_bs2 = bs2;
                if (m[wannier_index] == 1) tmp_bs2 = -bs2;

                std::complex<double> *orb_s = new std::complex<double>[npw];
                std::complex<double> *orb_px = new std::complex<double>[npw];

                get_trial_orbitals_lm_k(0, 0, ylm, gk, npw, radial_in_q[wannier_index].c, orb_s);
                get_trial_orbitals_lm_k(1, 1, ylm, gk, npw, radial_in_q[wannier_index].c+GlobalV::NQX, orb_px);

                for (int ig = 0; ig < npw; ig++)
                {
                    const double arg = (gk[ig] * R_centre[wannier_index]) * ModuleBase::TWO_PI;
                    sf[ig] = std::complex<double>(cos(arg), -sin(arg));

                    trial_orbitals_k(wannier_index, ig) = sf[ig] * (bs2 * orb_s[ig]  + tmp_bs2 * orb_px[ig]);
                }

                delete[] orb_s;
                delete[] orb_px;
            }

            if (L[wannier_index] == -2)
            {
                if (m[wannier_index] == 0 || m[wannier_index] == 1)
                {
                    double tmp_bs2 = bs2;
                    if (m[wannier_index] == 1) tmp_bs2 = -bs2;

                    std::complex<double> *orb_s = new std::complex<double>[npw];
                    std::complex<double> *orb_px = new std::complex<double>[npw];
                    std::complex<double> *orb_py = new std::complex<double>[npw];

                    get_trial_orbitals_lm_k(0, 0, ylm, gk, npw, radial_in_q[wannier_index].c, orb_s);
                    get_trial_orbitals_lm_k(1, 1, ylm, gk, npw, radial_in_q[wannier_index].c+GlobalV::NQX, orb_px);
                    get_trial_orbitals_lm_k(1, 2, ylm, gk, npw, radial_in_q[wannier_index].c+GlobalV::NQX, orb_py);

                    for (int ig = 0; ig < npw; ig++)
                    {
                        const double arg = (gk[ig] * R_centre[wannier_index]) * ModuleBase::TWO_PI;
                        sf[ig] = std::complex<double>(cos(arg), -sin(arg));

                        trial_orbitals_k(wannier_index, ig) = sf[ig] * (bs3 * orb_s[ig] - bs6 * orb_px[ig] + tmp_bs2 * orb_py[ig]);
                    }

                    delete[] orb_s;
                    delete[] orb_px;
                    delete[] orb_py;
                }
                else if (m[wannier_index] == 2)
                {
                    std::complex<double> *orb_s = new std::complex<double>[npw];
                    std::complex<double> *orb_px = new std::complex<double>[npw];
                    get_trial_orbitals_lm_k(0, 0, ylm, gk, npw, radial_in_q[wannier_index].c, orb_s);
                    get_trial_orbitals_lm_k(1, 1, ylm, gk, npw, radial_in_q[wannier_index].c+GlobalV::NQX, orb_px);

                    for (int ig = 0; ig < npw; ig++)
                    {
                        const double arg = (gk[ig] * R_centre[wannier_index]) * ModuleBase::TWO_PI;
                        sf[ig] = std::complex<double>(cos(arg), -sin(arg));

                        trial_orbitals_k(wannier_index, ig) = sf[ig] * (bs3 * orb_s[ig] + 2.0 * bs6 * orb_px[ig]);
                    }

                    delete[] orb_s;
                    delete[] orb_px;
                }
            }

            if (L[wannier_index] == -3)
            {
                double m_px = 1.0;
                double m_py = 1.0;
                double m_pz = 1.0;

                if (m[wannier_index] == 1)
                {
                    m_py = -1.0;
                    m_pz = -1.0;
                }
                else if (m[wannier_index] == 2)
                {
                    m_px = -1.0;
                    m_pz = -1.0;
                }
                else if (m[wannier_index] == 3)
                {
                    m_px = -1.0;
                    m_py = -1.0;
                }

                std::complex<double> *orb_s = new std::complex<double>[npw];
                std::complex<double> *orb_px = new std::complex<double>[npw];
                std::complex<double> *orb_py = new std::complex<double>[npw];
                std::complex<double> *orb_pz = new std::complex<double>[npw];

                get_trial_orbitals_lm_k(0, 0, ylm, gk, npw, radial_in_q[wannier_index].c, orb_s);
                get_trial_orbitals_lm_k(1, 1, ylm, gk, npw, radial_in_q[wannier_index].c+GlobalV::NQX, orb_px);
                get_trial_orbitals_lm_k(1, 2, ylm, gk, npw, radial_in_q[wannier_index].c+GlobalV::NQX, orb_py);
                get_trial_orbitals_lm_k(1, 0, ylm, gk, npw, radial_in_q[wannier_index].c+GlobalV::NQX, orb_pz);

                for (int ig = 0; ig < npw; ig++)
                {
                    const double arg = (gk[ig] * R_centre[wannier_index]) * ModuleBase::TWO_PI;
                    sf[ig] = std::complex<double>(cos(arg), -sin(arg));

                    trial_orbitals_k(wannier_index, ig) = sf[ig] * 0.5 * (orb_s[ig] + m_px * orb_px[ig] + m_py * orb_py[ig] + m_pz * orb_pz[ig]);
                }

                delete[] orb_s;
                delete[] orb_px;
                delete[] orb_py;
                delete[] orb_pz;
            }
            
            if (L[wannier_index] == -4)
            {
                if (m[wannier_index] == 0 || m[wannier_index] == 1)
                {
                    double tmp_bs2 = bs2;
                    if (m[wannier_index] == 1) tmp_bs2 = -bs2;

                    std::complex<double> *orb_s = new std::complex<double>[npw];
                    std::complex<double> *orb_px = new std::complex<double>[npw];
                    std::complex<double> *orb_py = new std::complex<double>[npw];

                    get_trial_orbitals_lm_k(0, 0, ylm, gk, npw, radial_in_q[wannier_index].c, orb_s);
                    get_trial_orbitals_lm_k(1, 1, ylm, gk, npw, radial_in_q[wannier_index].c+GlobalV::NQX, orb_px);
                    get_trial_orbitals_lm_k(1, 2, ylm, gk, npw, radial_in_q[wannier_index].c+GlobalV::NQX, orb_py);

                    for (int ig = 0; ig < npw; ig++)
                    {
                        const double arg = (gk[ig] * R_centre[wannier_index]) * ModuleBase::TWO_PI;
                        sf[ig] = std::complex<double>(cos(arg), -sin(arg));

                        trial_orbitals_k(wannier_index, ig) = sf[ig] * (bs3 * orb_s[ig] - bs6 * orb_px[ig] + tmp_bs2 * orb_py[ig]);
                    }

                    delete[] orb_s;
                    delete[] orb_px;
                    delete[] orb_py;
                }
                else if (m[wannier_index] == 2)
                {
                    std::complex<double> *orb_s = new std::complex<double>[npw];
                    std::complex<double> *orb_px = new std::complex<double>[npw];

                    get_trial_orbitals_lm_k(0, 0, ylm, gk, npw, radial_in_q[wannier_index].c, orb_s);
                    get_trial_orbitals_lm_k(1, 1, ylm, gk, npw, radial_in_q[wannier_index].c+GlobalV::NQX, orb_px);

                    for (int ig = 0; ig < npw; ig++)
                    {
                        const double arg = (gk[ig] * R_centre[wannier_index]) * ModuleBase::TWO_PI;
                        sf[ig] = std::complex<double>(cos(arg), -sin(arg));

                        trial_orbitals_k(wannier_index, ig) = sf[ig] * (bs3 * orb_s[ig] + 2.0 * bs6 * orb_px[ig]);
                    }

                    delete[] orb_s;
                    delete[] orb_px; 
                }
                else if (m[wannier_index] == 3 || m[wannier_index] == 4)
                {
                    double m_pz = 1.0;
                    if (m[wannier_index] == 4) m_pz = -1.0;

                    std::complex<double> *orb_pz = new std::complex<double>[npw];
                    std::complex<double> *orb_dz2 = new std::complex<double>[npw];

                    get_trial_orbitals_lm_k(1, 0, ylm, gk, npw, radial_in_q[wannier_index].c+GlobalV::NQX, orb_pz);
                    get_trial_orbitals_lm_k(2, 0, ylm, gk, npw, radial_in_q[wannier_index].c+2*GlobalV::NQX, orb_dz2);

                    for (int ig = 0; ig < npw; ig++)
                    {
                        const double arg = (gk[ig] * R_centre[wannier_index]) * ModuleBase::TWO_PI;
                        sf[ig] = std::complex<double>(cos(arg), -sin(arg));

                        trial_orbitals_k(wannier_index, ig) = sf[ig] * bs2 * (m_pz * orb_pz[ig] + orb_dz2[ig]);
                    }

                    delete[] orb_pz;
                    delete[] orb_dz2;
                }
            }

            if (L[wannier_index] == -5)
            {
                if (m[wannier_index] == 0 || m[wannier_index] == 1)
                {
                    double tmp_bs2 = -bs2;
                    double tmp_bs12 = -bs12;
                    double tmp_d = 0.5;

                    if (m[wannier_index] == 1)
                    {
                        tmp_bs2 = bs2;
                    }

                    std::complex<double> *orb_s = new std::complex<double>[npw];
                    std::complex<double> *orb_px = new std::complex<double>[npw];
                    std::complex<double> *orb_dz2 = new std::complex<double>[npw];
                    std::complex<double> *orb_dx2_y2 = new std::complex<double>[npw];

                    get_trial_orbitals_lm_k(0, 0, ylm, gk, npw, radial_in_q[wannier_index].c, orb_s);
                    get_trial_orbitals_lm_k(1, 1, ylm, gk, npw, radial_in_q[wannier_index].c+GlobalV::NQX, orb_px);
                    get_trial_orbitals_lm_k(2, 0, ylm, gk, npw, radial_in_q[wannier_index].c+2*GlobalV::NQX, orb_dz2);
                    get_trial_orbitals_lm_k(2, 3, ylm, gk, npw, radial_in_q[wannier_index].c+2*GlobalV::NQX, orb_dx2_y2);

                    for (int ig = 0; ig < npw; ig++)
                    {
                        const double arg = (gk[ig] * R_centre[wannier_index]) * ModuleBase::TWO_PI;
                        sf[ig] = std::complex<double>(cos(arg), -sin(arg));

                        trial_orbitals_k(wannier_index, ig) = sf[ig] * (bs6 * orb_s[ig] + tmp_bs2 * orb_px[ig] + tmp_bs12 * orb_dz2[ig] + tmp_d * orb_dx2_y2[ig]);
                    }

                    delete[] orb_s;
                    delete[] orb_px;
                    delete[] orb_dz2;
                    delete[] orb_dx2_y2;
                }
                else if (m[wannier_index] == 2 || m[wannier_index] == 3)
                {
                    double tmp_bs2 = -bs2;
                    double tmp_bs12 = -bs12;
                    double tmp_d = -0.5;

                    if (m[wannier_index] == 3)
                    {
                        tmp_bs2 = bs2;
                    }

                    std::complex<double> *orb_s = new std::complex<double>[npw];
                    std::complex<double> *orb_py = new std::complex<double>[npw];
                    std::complex<double> *orb_dz2 = new std::complex<double>[npw];
                    std::complex<double> *orb_dx2_y2 = new std::complex<double>[npw];

                    get_trial_orbitals_lm_k(0, 0, ylm, gk, npw, radial_in_q[wannier_index].c, orb_s);
                    get_trial_orbitals_lm_k(1, 2, ylm, gk, npw, radial_in_q[wannier_index].c+GlobalV::NQX, orb_py);
                    get_trial_orbitals_lm_k(2, 0, ylm, gk, npw, radial_in_q[wannier_index].c+2*GlobalV::NQX, orb_dz2);
                    get_trial_orbitals_lm_k(2, 3, ylm, gk, npw, radial_in_q[wannier_index].c+2*GlobalV::NQX, orb_dx2_y2);

                    for (int ig = 0; ig < npw; ig++)
                    {
                        const double arg = (gk[ig] * R_centre[wannier_index]) * ModuleBase::TWO_PI;
                        sf[ig] = std::complex<double>(cos(arg), -sin(arg));

                        trial_orbitals_k(wannier_index, ig) = sf[ig] * (bs6 * orb_s[ig] + tmp_bs2 * orb_py[ig] + tmp_bs12 * orb_dz2[ig] + tmp_d * orb_dx2_y2[ig]);
                    }

                    delete[] orb_s;
                    delete[] orb_py;
                    delete[] orb_dz2;
                    delete[] orb_dx2_y2;
                }
                else if (m[wannier_index] == 4 || m[wannier_index] == 5)
                {
                    double tmp_pz = -1.0;

                    if (m[wannier_index] == 5) tmp_pz = 1.0;

                    std::complex<double> *orb_s = new std::complex<double>[npw];
                    std::complex<double> *orb_pz = new std::complex<double>[npw];
                    std::complex<double> *orb_dz2 = new std::complex<double>[npw];

                    get_trial_orbitals_lm_k(0, 0, ylm, gk, npw, radial_in_q[wannier_index].c, orb_s);
                    get_trial_orbitals_lm_k(1, 0, ylm, gk, npw, radial_in_q[wannier_index].c+GlobalV::NQX, orb_pz);
                    get_trial_orbitals_lm_k(2, 0, ylm, gk, npw, radial_in_q[wannier_index].c+2*GlobalV::NQX, orb_dz2);

                    for (int ig = 0; ig < npw; ig++)
                    {
                        const double arg = (gk[ig] * R_centre[wannier_index]) * ModuleBase::TWO_PI;
                        sf[ig] = std::complex<double>(cos(arg), -sin(arg));

                        trial_orbitals_k(wannier_index, ig) = sf[ig] * (bs6 * orb_s[ig] + tmp_pz * bs2 * orb_pz[ig] + bs3 * orb_dz2[ig]);
                    }

                    delete[] orb_s;
                    delete[] orb_pz;
                    delete[] orb_dz2;
                }
            }

        }
    }

    for (int wannier_index = 0; wannier_index < num_wannier; wannier_index++)
    {
        std::complex<double> anorm(0.0, 0.0);
        for (int ig = 0; ig < npw; ig++)
        {
            anorm += conj(trial_orbitals_k(wannier_index, ig)) * trial_orbitals_k(wannier_index, ig);
        }

#ifdef __MPI
        Parallel_Reduce::reduce_all(anorm);
#endif
        for (int ig = 0; ig < npw; ig++)
        {
            trial_orbitals_k(wannier_index, ig) = trial_orbitals_k(wannier_index, ig) / sqrt(anorm);
        }
    }

    delete[] gk;
    delete[] sf;
}

void toWannier90_PW::get_trial_orbitals_lm_k(
    const int &orbital_L,
    const int &orbital_m,
    const ModuleBase::matrix &ylm,
    const ModuleBase::Vector3<double> *gk,
    const int &npw,
    double *radial_in_q_single,
    std::complex<double> *orbital_in_G_single
)
{
    for (int ig = 0; ig < npw; ig++)
    {
        orbital_in_G_single[ig] = ModuleBase::PolyInt::Polynomial_Interpolation(radial_in_q_single, GlobalV::NQX, GlobalV::DQ, gk[ig].norm() * GlobalC::ucell.tpiba);
    }

    std::complex<double> lphase = pow(ModuleBase::NEG_IMAG_UNIT, orbital_L);
    int index = orbital_L * orbital_L + orbital_m;
    for (int ig = 0; ig < npw; ig++)
    {
        orbital_in_G_single[ig] = lphase * ylm(index, ig) * orbital_in_G_single[ig];
    }

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

void toWannier90_PW::unkdotW_A(
    const psi::Psi<std::complex<double>>& psi_pw,
    const ModulePW::PW_Basis_K* wfcpw,
    const int& ik,
    const std::vector<ModuleBase::matrix> &radial_in_q,
    ModuleBase::ComplexMatrix &Amn
)
{
    Amn.create(num_bands, num_wannier);

    int npw = wfcpw->npwk[ik];
    int npwx = wfcpw->npwk_max;
    ModuleBase::ComplexMatrix trial_orbitals;
    produce_trial_in_pw(psi_pw, ik, wfcpw, radial_in_q, trial_orbitals);

    for (int iw = 0; iw < num_wannier; iw++)
    {
        for (int ib_w = 0; ib_w < num_bands; ib_w++)
        {
            int ib = cal_band_index[ib_w];

            if (GlobalV::NSPIN != 4)
            {
                for (int ig = 0; ig < npw; ig++)
                {
                    Amn(ib_w, iw) += conj(psi_pw(ik, ib, ig)) * trial_orbitals(iw, ig);
                }
            }
            else
            {
                for (int ig = 0; ig < npw; ig++)
                {
                    Amn(ib_w, iw) += up_con[iw] * conj(psi_pw(ik, ib, ig)) * trial_orbitals(iw, ig)
                                   + dn_con[iw] * conj(psi_pw(ik, ib, ig+npwx)) * trial_orbitals(iw, ig);
                }
            }
        }
    }

#ifdef __MPI
    Parallel_Reduce::reduce_all(Amn.c, Amn.size);
#endif

}