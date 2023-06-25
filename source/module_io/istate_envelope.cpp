#include "istate_envelope.h"

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/rho_io.h"
#include "module_io/write_wfc_pw.h"
#include "module_io/write_wfc_r.h"

IState_Envelope::IState_Envelope(const elecstate::ElecState* pes_in)
{pes = pes_in;}

IState_Envelope::~IState_Envelope()
{}

void IState_Envelope::begin(const psi::Psi<double>* psid,
                            const ModulePW::PW_Basis* rhopw,
                            const ModulePW::PW_Basis_K* wfcpw,
                            const ModulePW::PW_Basis_Big* bigpw,
                            Local_Orbital_wfc& lowf,
                            Gint_Gamma& gg,
                            int& out_wfc_pw,
                            int& out_wfc_r,
                            const K_Vectors& kv)
{
    ModuleBase::TITLE("IState_Envelope", "begin");

    std::cout << " perform |psi(band, r)| for selected bands." << std::endl;

    // (1) 
    // mohan update 2011-03-21
    // if ucell is odd, it's correct,
    // if ucell is even, it's also correct.
    // +1.0e-8 in case like (2.999999999+1)/2
    int fermi_band = static_cast<int>((GlobalV::nelec + 1) / 2 + 1.0e-8);
    int bands_below = GlobalV::NBANDS_ISTATE;
    int bands_above = GlobalV::NBANDS_ISTATE;

    std::cout << " number of electrons = " << GlobalV::nelec << std::endl;
    std::cout << " number of occupied bands = " << fermi_band << std::endl;
    std::cout << " plot band decomposed charge density below fermi surface with "
        << bands_below << " bands." << std::endl;

    std::cout << " plot band decomposed charge density above fermi surface with "
        << bands_above << " bands." << std::endl;

    // (2) cicle:

    // (2.1) calculate the selected density matrix
    // from wave functions.

    // (2.2) carry out the grid integration to
    // get the charge density.

    // (2.3) output the charge density in .cub format.
    this->bands_picked = new bool[GlobalV::NBANDS];
    ModuleBase::GlobalFunc::ZEROS(bands_picked, GlobalV::NBANDS);
    for (int ib = 0; ib < GlobalV::NBANDS; ib++)
    {
        if (ib >= fermi_band - bands_below)
        {
            if (ib < fermi_band + bands_above)
            {
                bands_picked[ib] = true;
            }
        }
    }

    //allocate grid wavefunction for gamma_only
    std::vector<double**> wfc_gamma_grid(GlobalV::NSPIN);
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        wfc_gamma_grid[is] = new double* [GlobalV::NBANDS];
        for (int ib = 0;ib < GlobalV::NBANDS; ++ib)
            wfc_gamma_grid[is][ib] = new double[gg.gridt->lgd];
    }

    //for pw-wfc in G space
    psi::Psi<std::complex<double>> pw_wfc_g;

    if (out_wfc_pw || out_wfc_r)
    {
        pw_wfc_g.resize(1, GlobalV::NBANDS, kv.ngk[0]);
    }


    for (int ib = 0; ib < GlobalV::NBANDS; ib++)
    {
        if (bands_picked[ib])
        {
            for (int is = 0; is < GlobalV::NSPIN; ++is)
            {
                std::cout << " Perform envelope function for band " << ib + 1 << std::endl;
                ModuleBase::GlobalFunc::ZEROS(pes->charge->rho[is], wfcpw->nrxx);

                psid->fix_k(is);
#ifdef __MPI
                lowf.wfc_2d_to_grid(-1, 0, psid->get_pointer(), wfc_gamma_grid[is], this->pes->ekb, this->pes->wg);
#else
                for (int i = 0;i < GlobalV::NBANDS;++i)
                {
                    for (int j = 0;j < GlobalV::NLOCAL;++j)
                        wfc_gamma_grid[is][i][j] = psid[0](i, j);
                }
#endif
                gg.cal_env(wfc_gamma_grid[is][ib], pes->charge->rho[is]);


                pes->charge->save_rho_before_sum_band(); //xiaohui add 2014-12-09
                std::stringstream ss;
                ss << GlobalV::global_out_dir << "BAND" << ib + 1 << "_s_" << is + 1 << "_ENV.cube";
                const double ef_tmp = this->pes->eferm.get_efval(is);
                ModuleIO::write_rho(
#ifdef __MPI
                    bigpw->bz,
                    bigpw->nbz,
                    rhopw->nplane,
                    rhopw->startz_current,
#endif
                    pes->charge->rho_save[is],
                    is,
                    GlobalV::NSPIN,
                    0,
                    ss.str(),
                    rhopw->nx,
                    rhopw->ny,
                    rhopw->nz,
                    ef_tmp,
                    &(GlobalC::ucell),
                    3);

                if (out_wfc_pw || out_wfc_r) //only for gamma_only now
                    this->set_pw_wfc(wfcpw, 0, ib, GlobalV::NSPIN,
                        pes->charge->rho_save, pw_wfc_g);
            }
        }
    }

    if (out_wfc_pw)
    {
        std::stringstream ssw;
        ssw << GlobalV::global_out_dir << "WAVEFUNC";
        std::cout << " write G-space wavefunction into \"" <<
            GlobalV::global_out_dir << "/" << ssw.str() << "\" files." << std::endl;
        ModuleIO::write_wfc_pw(ssw.str(), pw_wfc_g, kv, wfcpw);
    }
    if (out_wfc_r)
    {
        ModuleIO::write_psi_r_1(pw_wfc_g, wfcpw, "wfc_realspace", false, kv);
    }

    delete[] bands_picked;
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        for (int ib = 0;ib < GlobalV::NBANDS; ++ib)
            delete[] wfc_gamma_grid[is][ib];
        delete[] wfc_gamma_grid[is];
    }
    return;
}

void IState_Envelope::begin(const psi::Psi<std::complex<double>>* psi,
                            const ModulePW::PW_Basis* rhopw,
                            const ModulePW::PW_Basis_K* wfcpw,
                            const ModulePW::PW_Basis_Big* bigpw,
                            Local_Orbital_wfc& lowf,
                            Gint_k& gk,
                            int& out_wf,
                            int& out_wf_r,
                            const K_Vectors& kv)
{
    ModuleBase::TITLE("IState_Envelope", "begin");

    std::cout << " perform |psi(band, r)| for selected bands." << std::endl;

    // (1) 
    // mohan update 2011-03-21
    // if ucell is odd, it's correct,
    // if ucell is even, it's also correct.
    // +1.0e-8 in case like (2.999999999+1)/2
    int fermi_band = static_cast<int>((GlobalV::nelec + 1) / 2 + 1.0e-8);
    int bands_below = GlobalV::NBANDS_ISTATE;
    int bands_above = GlobalV::NBANDS_ISTATE;

    std::cout << " number of electrons = " << GlobalV::nelec << std::endl;
    std::cout << " number of occupied bands = " << fermi_band << std::endl;
    std::cout << " plot band decomposed charge density below fermi surface with "
        << bands_below << " bands." << std::endl;

    // (2) cicle:

    // (2.1) calculate the selected density matrix
    // from wave functions.

    // (2.2) carry out the grid integration to
    // get the charge density.

    // (2.3) output the charge density in .cub format.
    this->bands_picked = new bool[GlobalV::NBANDS];
    ModuleBase::GlobalFunc::ZEROS(bands_picked, GlobalV::NBANDS);
    for (int ib = 0; ib < GlobalV::NBANDS; ib++)
    {
        if (ib >= fermi_band - bands_below)
        {
            if (ib < fermi_band + bands_above)
            {
                bands_picked[ib] = true;
            }
        }
    }

    //for pw-wfc in G space
    psi::Psi<std::complex<double>> pw_wfc_g(kv.ngk.data());

    if (out_wf || out_wf_r)
    {
        pw_wfc_g.resize(kv.nks, GlobalV::NBANDS, wfcpw->npwk_max);
    }

    for (int ib = 0; ib < GlobalV::NBANDS; ib++)
    {
        if (bands_picked[ib])
        {
            const int nspin0 = (GlobalV::NSPIN == 2) ? 2 : 1;
            for (int ik = 0; ik < kv.nks; ++ik)    //the loop of nspin0 is included
            {
                const int ispin = kv.isk[ik];
                ModuleBase::GlobalFunc::ZEROS(pes->charge->rho[ispin], wfcpw->nrxx);
                std::cout << " Perform envelope function for kpoint " << ik << ",  band" << ib + 1 << std::endl;
                //  2d-to-grid conversion is unified into `wfc_2d_to_grid`.
                psi->fix_k(ik);
#ifdef __MPI
                // need to deal with NSPIN=4 !!!!
                lowf.wfc_2d_to_grid(-1,
                                    0,
                                    psi->get_pointer(),
                                    lowf.wfc_k_grid[ik],
                                    ik,
                                    this->pes->ekb,
                                    this->pes->wg,
                                    kv.kvec_c);
#else
                for (int i = 0;i < GlobalV::NBANDS;++i)
                {
                    for (int j = 0;j < GlobalV::NLOCAL;++j)
                        lowf.wfc_k_grid[ik][i][j] = psi[0](i, j);
                }
#endif
                //deal with NSPIN=4
                gk.cal_env_k(ik, lowf.wfc_k_grid[ik][ib], pes->charge->rho[ispin], kv.kvec_c, kv.kvec_d);

                std::stringstream ss;
                ss << GlobalV::global_out_dir << "BAND" << ib + 1 << "_k_" << ik / nspin0 + 1 << "_s_" << ispin + 1 << "_ENV.cube";
                const double ef_tmp = this->pes->eferm.get_efval(ispin);
                ModuleIO::write_rho(
#ifdef __MPI
                    bigpw->bz,
                    bigpw->nbz,
                    rhopw->nplane,
                    rhopw->startz_current,
#endif
                    pes->charge->rho[ispin],
                    ispin,
                    GlobalV::NSPIN,
                    0,
                    ss.str(),
                    rhopw->nx,
                    rhopw->ny,
                    rhopw->nz,
                    ef_tmp,
                    &(GlobalC::ucell),
                    3);

                if (out_wf || out_wf_r) //only for gamma_only now
                {
                    pw_wfc_g.fix_k(ik);
                    this->set_pw_wfc(wfcpw, ik, ib, GlobalV::NSPIN,
                        pes->charge->rho, pw_wfc_g);
                }
            }
        }
    }

    if (out_wf || out_wf_r)
    {
        if (out_wf)
        {
            std::stringstream ssw;
            ssw << GlobalV::global_out_dir << "WAVEFUNC";
            std::cout << " write G-space wavefunction into \"" <<
                GlobalV::global_out_dir << "/" << ssw.str() << "\" files." << std::endl;
            ModuleIO::write_wfc_pw(ssw.str(), pw_wfc_g, kv, wfcpw);
        }
        if (out_wf_r)
        {
            ModuleIO::write_psi_r_1(pw_wfc_g, wfcpw, "wfc_realspace", false, kv);
        }
    }

    delete[] bands_picked;
    return;
}

//for each band
void IState_Envelope::set_pw_wfc(const ModulePW::PW_Basis_K* wfcpw,
                                 const int& ik,
                                 const int& ib,
                                 const int& nspin,
                                 const double* const* const rho,
                                 psi::Psi<std::complex<double>>& wfc_g)
{
    if (ib == 0)//once is enough
        ModuleBase::TITLE("IState_Envelope", "set_pw_wfc");

    std::vector<std::complex<double>> Porter(wfcpw->nrxx);
    // here I refer to v_hartree, but I don't know how to deal with NSPIN=4
    const int nspin0 = (nspin == 2) ? 2 : 1;
    for (int is = 0; is < nspin0; is++)
        for (int ir = 0; ir < wfcpw->nrxx; ir++)
            Porter[ir] += std::complex<double>(rho[is][ir], 0.0);

    //call FFT
    wfcpw->real2recip(Porter.data(), &wfc_g(ib, 0), ik);
}
