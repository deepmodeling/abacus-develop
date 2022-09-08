#include "istate_envelope.h"
#include "../src_pw/global.h"
#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "src_io/wf_io.h"
#include "src_io/write_wfc_realspace.h"

IState_Envelope::IState_Envelope(const elecstate::ElecState* pes_in)
{pes = pes_in;}

IState_Envelope::~IState_Envelope()
{}


void IState_Envelope::begin(const psi::Psi<double>* psid, Local_Orbital_wfc& lowf, Gint_Gamma& gg, int& out_wfc_pw, int& out_wfc_r)
{
    ModuleBase::TITLE("IState_Envelope", "begin");

    std::cout << " perform |psi(band, r)| for selected bands." << std::endl;

    // (1) 
    // (1.1) allocate the space for GlobalC::LOWF.WFC_GAMMA

    // (1.2) read in LOWF_GAMMA.dat

    // mohan update 2011-03-21
    // if ucell is odd, it's correct,
    // if ucell is even, it's also correct.
    // +1.0e-8 in case like (2.999999999+1)/2
    int fermi_band = static_cast<int>((GlobalC::CHR.nelec + 1) / 2 + 1.0e-8);
    int bands_below = GlobalV::NBANDS_ISTATE;
    int bands_above = GlobalV::NBANDS_ISTATE;

    std::cout << " number of electrons = " << GlobalC::CHR.nelec << std::endl;
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
            wfc_gamma_grid[is][ib] = new double[GlobalC::GridT.lgd];
    }

    //for pw-wfc in G space
    psi::Psi<std::complex<double>> pw_wfc_g;

    if (out_wfc_pw || out_wfc_r)
    {
        pw_wfc_g.resize(1, GlobalV::NBANDS, GlobalC::kv.ngk[0]);
    }


    for (int ib = 0; ib < GlobalV::NBANDS; ib++)
    {
        if (bands_picked[ib])
        {
            for (int is = 0; is < GlobalV::NSPIN; ++is)
            {
                std::cout << " Perform envelope function for band " << ib + 1 << std::endl;
                ModuleBase::GlobalFunc::ZEROS(GlobalC::CHR.rho[is], GlobalC::wfcpw->nrxx);


                //---------------------------------------------------------
                // GlobalC::LOWF.WFC_GAMMA has been replaced by wfc_dm_2d.cpp 
                // and 2d-to-grid conversion is unified into `wfc_2d_to_grid`.
                //---------------------------------------------------------
                psid->fix_k(is);
#ifdef __MPI
                lowf.wfc_2d_to_grid(0, psid->get_pointer(), wfc_gamma_grid[is], this->pes->ekb, this->pes->wg);
#else
                for (int i = 0;i < GlobalV::NBANDS;++i)
                {
                    for (int j = 0;j < GlobalV::NLOCAL;++j)
                        wfc_gamma_grid[is][i][j] = psid[0](i, j);
                }
#endif
                gg.cal_env(wfc_gamma_grid[is][ib], GlobalC::CHR.rho[is]);


                GlobalC::CHR.save_rho_before_sum_band(); //xiaohui add 2014-12-09
                std::stringstream ss;
                ss << GlobalV::global_out_dir << "BAND" << ib + 1 << "_s_" << is + 1 << "_ENV";
                // 0 means definitely output charge density.
                bool for_plot = true;
                GlobalC::CHR.write_rho(GlobalC::CHR.rho_save[is], is, 0, ss.str(), 3, for_plot);

                if (out_wfc_pw || out_wfc_r) //only for gamma_only now
                    this->set_pw_wfc(GlobalC::wfcpw, 0, ib, GlobalV::NSPIN,
                        GlobalC::CHR.rho_save, pw_wfc_g);
            }
        }
    }

    if (out_wfc_pw || out_wfc_r)
    {
        if (out_wfc_pw)
        {
            std::stringstream ssw;
            ssw << GlobalV::global_out_dir << "WAVEFUNC";
            std::cout << " write G-space wavefunction into \"" <<
                GlobalV::global_out_dir << "/" << ssw.str() << "\" files." << std::endl;
            WF_io::write_wfc(ssw.str(), pw_wfc_g, &GlobalC::kv, GlobalC::wfcpw);
        }
        if (out_wfc_r)
        {
            Write_Wfc_Realspace::write_wfc_realspace_1(pw_wfc_g, "wfc_realspace", false);
        }
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

void IState_Envelope::begin(const psi::Psi<std::complex<double>>* psi, Local_Orbital_wfc& lowf, Gint_k& gk, int& out_wf, int& out_wf_r)
{
    ModuleBase::TITLE("IState_Envelope", "begin");

    std::cout << " perform |psi(band, r)| for selected bands." << std::endl;

    // (1) 
    // (1.1) allocate the space for GlobalC::LOWF.WFC_GAMMA

    // (1.2) read in LOWF_GAMMA.dat

    // mohan update 2011-03-21
    // if ucell is odd, it's correct,
    // if ucell is even, it's also correct.
    // +1.0e-8 in case like (2.999999999+1)/2
    int fermi_band = static_cast<int>((GlobalC::CHR.nelec + 1) / 2 + 1.0e-8);
    int bands_below = GlobalV::NBANDS_ISTATE;
    int bands_above = GlobalV::NBANDS_ISTATE;

    std::cout << " number of electrons = " << GlobalC::CHR.nelec << std::endl;
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
    psi::Psi<std::complex<double>> pw_wfc_g(GlobalC::kv.ngk.data());

    if (out_wf || out_wf_r)
    {
        pw_wfc_g.resize(GlobalC::kv.nks, GlobalV::NBANDS, GlobalC::wf.npwx);
    }

    for (int ib = 0; ib < GlobalV::NBANDS; ib++)
    {
        if (bands_picked[ib])
        {
            const int nspin0 = (GlobalV::NSPIN == 2) ? 2 : 1;
            for (int ik = 0; ik < GlobalC::kv.nks; ++ik)    //the loop of nspin0 is included
            {
                const int ispin = GlobalC::kv.isk[ik];
                ModuleBase::GlobalFunc::ZEROS(GlobalC::CHR.rho[ispin], GlobalC::wfcpw->nrxx);
                std::cout << " Perform envelope function for kpoint " << ik << ",  band" << ib + 1 << std::endl;
                //  2d-to-grid conversion is unified into `wfc_2d_to_grid`.
                psi->fix_k(ik);
#ifdef __MPI
                // need to deal with NSPIN=4 !!!!
                lowf.wfc_2d_to_grid(0, psi->get_pointer(), lowf.wfc_k_grid[ik], ik, this->pes->ekb, this->pes->wg);
#else
                for (int i = 0;i < GlobalV::NBANDS;++i)
                {
                    for (int j = 0;j < GlobalV::NLOCAL;++j)
                        lowf.wfc_k_grid[ik][i][j] = psi[0](i, j);
                }
#endif
                //deal with NSPIN=4
                gk.cal_env_k(ik, lowf.wfc_k_grid[ik][ib], GlobalC::CHR.rho[ispin]);

                std::stringstream ss;
                ss << GlobalV::global_out_dir << "BAND" << ib + 1 << "_k_" << ik / nspin0 + 1 << "_s_" << ispin + 1 << "_ENV";

                bool for_plot = true;   //if false, separate the output into spin up and spin down
                GlobalC::CHR.write_rho(GlobalC::CHR.rho[ispin], ispin, 0, ss.str(), 3, for_plot);

                if (out_wf || out_wf_r) //only for gamma_only now
                {
                    pw_wfc_g.fix_k(ik);
                    this->set_pw_wfc(GlobalC::wfcpw, ik, ib, GlobalV::NSPIN,
                        GlobalC::CHR.rho, pw_wfc_g);
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
            WF_io::write_wfc(ssw.str(), pw_wfc_g, &GlobalC::kv, GlobalC::wfcpw);
        }
        if (out_wf_r)
        {
            Write_Wfc_Realspace::write_wfc_realspace_1(pw_wfc_g, "wfc_realspace", false);
        }
    }

    delete[] bands_picked;
    return;
}

//for each band
void IState_Envelope::set_pw_wfc(ModulePW::PW_Basis_K* wfc_basis,
    const int& ik, const int& ib, const int& nspin,
    const double* const* const rho,
    psi::Psi<std::complex<double>> &wfc_g)
{
    if (ib == 0)//once is enough
        ModuleBase::TITLE("IState_Envelope", "set_pw_wfc");

    std::vector<std::complex<double>> Porter(wfc_basis->nrxx);
    // here I refer to v_hartree, but I don't know how to deal with NSPIN=4
    const int nspin0 = (nspin == 2) ? 2 : 1;
    for (int is = 0; is < nspin0; is++)
        for (int ir = 0; ir < wfc_basis->nrxx; ir++)
            Porter[ir] += std::complex<double>(rho[is][ir], 0.0);

    //call FFT
    wfc_basis->real2recip(Porter.data(), &wfc_g(ib,0), ik);
}
