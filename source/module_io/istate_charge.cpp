#include "istate_charge.h"

#include "module_base/blas_connector.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/parallel_common.h"
#include "module_base/scalapack_connector.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/rho_io.h"

IState_Charge::IState_Charge(psi::Psi<double>* psi_gamma_in, Local_Orbital_Charge& loc_in)
    : psi_gamma(psi_gamma_in), loc(&loc_in)
{
}

IState_Charge::~IState_Charge()
{
}

void IState_Charge::begin(Gint_Gamma& gg,
                          elecstate::ElecState* pelec,
                          const ModulePW::PW_Basis* rhopw,
                          const ModulePW::PW_Basis_Big* bigpw)
{
    ModuleBase::TITLE("IState_Charge", "begin");

    std::cout << " Perform |psi(i)|^2 for selected bands." << std::endl;

    if (!GlobalV::GAMMA_ONLY_LOCAL)
    {
        ModuleBase::WARNING_QUIT("IState_Charge::begin", "Only available for GAMMA_ONLY_LOCAL now.");
    }

    int mode = 0;
    if (GlobalV::NBANDS_ISTATE > 0)
        mode = 1;
    else
        mode = 2;

    int fermi_band = 0;
    int bands_below = 0;
    int bands_above = 0;

    // (2) cicle:
    // (2.1) calculate the selected density matrix
    // from wave functions.
    // (2.2) carry out the grid integration to
    // get the charge density.
    this->bands_picked = new int[GlobalV::NBANDS];
    ModuleBase::GlobalFunc::ZEROS(bands_picked, GlobalV::NBANDS);

    // (1)

    // (1.2) read in LOWF_GAMMA.dat
    std::cout << " number of electrons = " << GlobalV::nelec << std::endl;

    // mohan update 2011-03-21
    // if ucell is odd, it's correct,
    // if ucell is even, it's also correct.
    // +1.0e-8 in case like (2.999999999+1)/2
    fermi_band = static_cast<int>((GlobalV::nelec + 1) / 2 + 1.0e-8);
    std::cout << " number of occupied bands = " << fermi_band << std::endl;

    if (mode == 1)
    {
        bands_below = GlobalV::NBANDS_ISTATE;
        bands_above = GlobalV::NBANDS_ISTATE;

        std::cout << " plot band decomposed charge density below fermi surface with " << bands_below << " bands."
                  << std::endl;

        std::cout << " plot band decomposed charge density above fermi surface with " << bands_above << " bands."
                  << std::endl;

        for (int ib = 0; ib < GlobalV::NBANDS; ib++)
        {
            if (ib >= fermi_band - bands_below)
            {
                if (ib < fermi_band + bands_above)
                {
                    bands_picked[ib] = 1;
                }
            }
        }
    }
    else if (mode == 2)
    {
        bool stop = false;
        std::stringstream ss;
        ss << GlobalV::global_out_dir << "istate.info";
        std::cout << " Open the file : " << ss.str() << std::endl;
        if (GlobalV::MY_RANK == 0)
        {
            std::ifstream ifs(ss.str().c_str());
            if (!ifs)
            {
                stop = true;
            }
            else
            {
                // int band_index;
                for (int ib = 0; ib < GlobalV::NBANDS; ++ib)
                {
                    ModuleBase::GlobalFunc::READ_VALUE(ifs, bands_picked[ib]);
                }
            }
        }

#ifdef __MPI
        Parallel_Common::bcast_bool(stop);
        Parallel_Common::bcast_int(bands_picked, GlobalV::NBANDS);
#endif
        if (stop)
        {
            GlobalV::ofs_warning << " Can't find the file : " << ss.str() << std::endl;
            ModuleBase::WARNING_QUIT("IState_Charge::begin", "can't find the istate file.");
        }
    }

    for (int ib = 0; ib < GlobalV::NBANDS; ib++)
    {
        if (bands_picked[ib])
        {
            std::cout << " Perform band decomposed charge density for band " << ib + 1 << std::endl;

            // (1) calculate the density matrix for a partuclar
            // band, whenever it is occupied or not.

#ifdef __MPI
            this->idmatrix(ib, pelec);
#endif
            // (2) zero out of charge density array.
            for (int is = 0; is < GlobalV::NSPIN; is++)
            {
                ModuleBase::GlobalFunc::ZEROS(pelec->charge->rho[is], rhopw->nrxx);
            }

            // (3) calculate charge density for a particular
            // band.
            Gint_inout inout(this->loc->DM, pelec->charge->rho, Gint_Tools::job_type::rho);
            gg.cal_gint(&inout);
            pelec->charge->save_rho_before_sum_band(); // xiaohui add 2014-12-09
            std::stringstream ssc;
            ssc << GlobalV::global_out_dir << "BAND" << ib + 1;
            // 0 means definitely output charge density.
            for (int is = 0; is < GlobalV::NSPIN; is++)
            {
                ssc << "_SPIN" << is << "_CHG.cube";
                const double ef_tmp = pelec->eferm.get_efval(is);
                ModuleIO::write_rho(
#ifdef __MPI
                    bigpw->bz,
                    bigpw->nbz,
                    rhopw->nplane,
                    rhopw->startz_current,
#endif
                    pelec->charge->rho_save[is],
                    is,
                    GlobalV::NSPIN,
                    0,
                    ssc.str(),
                    rhopw->nx,
                    rhopw->ny,
                    rhopw->nz,
                    ef_tmp,
                    &(GlobalC::ucell));
            }
        }
    }

    delete[] bands_picked;
    return;
}

#ifdef __MPI
void IState_Charge::idmatrix(const int& ib, elecstate::ElecState* pelec)
{
    ModuleBase::TITLE("IState_Charge", "idmatrix");

    assert(pelec->wg.nr == GlobalV::NSPIN);
    for (int is = 0; is != GlobalV::NSPIN; ++is)
    {
        std::vector<double> wg_local(this->loc->ParaV->ncol, 0.0);
        const int ib_local = this->loc->ParaV->trace_loc_col[ib];

        int fermi_band = 0;
        fermi_band = static_cast<int>((GlobalV::nelec + 1) / 2 + 1.0e-8);

        if (ib_local >= 0)
        {
            if (ib < fermi_band)
            {
                wg_local[ib_local] = pelec->wg(is, ib);
            }
            else
            {
                wg_local[ib_local] = pelec->wg(is, fermi_band - 1);
            } // unoccupied bands, use occupation of homo
        }

        // wg_wfc(ib,iw) = pelec->wg[ib] * wfc(ib,iw);
        this->psi_gamma->fix_k(is);
        psi::Psi<double> wg_wfc(this->psi_gamma[0], 1);

        for (int ir = 0; ir != wg_wfc.get_nbands(); ++ir)
        {
            BlasConnector::scal(wg_wfc.get_nbasis(), wg_local[ir], wg_wfc.get_pointer() + ir * wg_wfc.get_nbasis(), 1);
        }

        // C++: dm(iw1,iw2) = wfc(ib,iw1).T * wg_wfc(ib,iw2)
        const double one_float = 1.0, zero_float = 0.0;
        const int one_int = 1;
        const char N_char = 'N', T_char = 'T';
        this->loc->dm_gamma.at(is).create(wg_wfc.get_nbands(), wg_wfc.get_nbasis());

        pdgemm_(&N_char,
                &T_char,
                &GlobalV::NLOCAL,
                &GlobalV::NLOCAL,
                &pelec->wg.nc,
                &one_float,
                wg_wfc.get_pointer(),
                &one_int,
                &one_int,
                this->loc->ParaV->desc,
                this->psi_gamma->get_pointer(),
                &one_int,
                &one_int,
                this->loc->ParaV->desc,
                &zero_float,
                this->loc->dm_gamma.at(is).c,
                &one_int,
                &one_int,
                this->loc->ParaV->desc);
    }

    std::cout << " finished calc dm_2d : " << std::endl;

    this->loc->cal_dk_gamma_from_2D_pub();

    std::cout << " finished convert : " << std::endl;
}
#endif
