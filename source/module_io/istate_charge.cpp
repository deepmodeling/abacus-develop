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
                          double** rho,
                          const ModuleBase::matrix& wg,
                          const std::vector<double>& ef_all_spin,
                          const int rhopw_nrxx,
                          const int rhopw_nplane,
                          const int rhopw_startz_current,
                          const int rhopw_nx,
                          const int rhopw_ny,
                          const int rhopw_nz,
                          const int bigpw_bz,
                          const int bigpw_nbz,
                          const bool gamma_only_local,
                          const int nbands_istate,
                          const std::vector<int>& out_band_kb,
                          const int nbands,
                          const double nelec,
                          const int nspin,
                          const int nlocal,
                          const std::string& global_out_dir,
                          const int my_rank,
                          std::ofstream& ofs_warning)
{
    ModuleBase::TITLE("IState_Charge", "begin");

    std::cout << " Perform |psi(i)|^2 for selected bands." << std::endl;

    if (!gamma_only_local)
    {
        ModuleBase::WARNING_QUIT("IState_Charge::begin", "Only available for GAMMA_ONLY_LOCAL now.");
    }

    int mode = 0;
    if (nbands_istate > 0 && static_cast<int>(out_band_kb.size()) == 0)
    {
        mode = 1;
    }
    else if (static_cast<int>(out_band_kb.size()) > 0)
    {
        // If out_band_kb (bands_to_print) is not empty, set mode to 2
        mode = 2;
        std::cout << " Notice: INPUT parameter `nbands_istate` overwritten by `bands_to_print`!" << std::endl;
    }
    else
    {
        mode = 3;
    }

    int fermi_band = 0;
    int bands_below = 0;
    int bands_above = 0;

    // (2) cicle:
    // (2.1) calculate the selected density matrix from wave functions.
    // (2.2) carry out the grid integration to get the charge density.
    this->bands_picked_.resize(nbands);
    ModuleBase::GlobalFunc::ZEROS(bands_picked_.data(), nbands);

    // (1)
    // (1.2) read in WFC_NAO_GAMMA1.dat
    std::cout << " number of electrons = " << nelec << std::endl;

    // mohan update 2011-03-21
    // if ucell is odd, it's correct,
    // if ucell is even, it's also correct.
    // +1.0e-8 in case like (2.999999999+1)/2
    fermi_band = static_cast<int>((nelec + 1) / 2 + 1.0e-8);
    std::cout << " number of occupied bands = " << fermi_band << std::endl;

    if (mode == 1)
    {
        bands_below = nbands_istate;
        bands_above = nbands_istate;

        std::cout << " Plot band decomposed charge density below Fermi surface with " << bands_below << " bands."
                  << std::endl;

        std::cout << " Plot band decomposed charge density above Fermi surface with " << bands_above << " bands."
                  << std::endl;

        for (int ib = 0; ib < nbands; ++ib)
        {
            if (ib >= fermi_band - bands_below)
            {
                if (ib < fermi_band + bands_above)
                {
                    bands_picked_[ib] = 1;
                }
            }
        }
    }
    else if (mode == 2)
    {
        // Check if length of out_band_kb is valid
        if (static_cast<int>(out_band_kb.size()) > nbands)
        {
            ModuleBase::WARNING_QUIT(
                "IState_Charge::begin",
                "The number of bands specified by `bands_to_print` in the INPUT file exceeds `nbands`!");
        }
        // Check if all elements in bands_picked_ are 0 or 1
        for (int value: out_band_kb)
        {
            if (value != 0 && value != 1)
            {
                ModuleBase::WARNING_QUIT(
                    "IState_Charge::begin",
                    "The elements of `bands_to_print` must be either 0 or 1. Invalid values found!");
            }
        }
        // Fill bands_picked_ with values from out_band_kb
        // Remaining bands are already set to 0
        int length = std::min(static_cast<int>(out_band_kb.size()), nbands);
        for (int i = 0; i < length; ++i)
        {
            // out_band_kb rely on function parse_expression from input_conv.cpp
            bands_picked_[i] = out_band_kb[i];
        }

        std::cout << " Plot band decomposed charge density below the Fermi surface: band ";
        for (int i = 0; i + 1 <= fermi_band; ++i)
        {
            if (bands_picked_[i] == 1)
            {
                std::cout << i + 1 << " ";
            }
        }
        std::cout << std::endl;
        std::cout << " Plot band decomposed charge density above the Fermi surface: band ";
        for (int i = fermi_band; i < nbands; ++i)
        {
            if (bands_picked_[i] == 1)
            {
                std::cout << i + 1 << " ";
            }
        }
        std::cout << std::endl;
    }
    else if (mode == 3)
    {
        bool stop = false;
        std::stringstream ss;
        ss << global_out_dir << "istate.info";
        std::cout << " Open the file : " << ss.str() << std::endl;
        if (my_rank == 0)
        {
            std::ifstream ifs(ss.str().c_str());
            if (!ifs)
            {
                stop = true;
            }
            else
            {
                // int band_index;
                for (int ib = 0; ib < nbands; ++ib)
                {
                    ModuleBase::GlobalFunc::READ_VALUE(ifs, bands_picked_[ib]);
                }
            }
        }

#ifdef __MPI
        Parallel_Common::bcast_bool(stop);
        Parallel_Common::bcast_int(bands_picked_.data(), nbands);
#endif
        if (stop)
        {
            ofs_warning << " Can't find the file : " << ss.str() << std::endl;
            ModuleBase::WARNING_QUIT("IState_Charge::begin", "can't find the istate file.");
        }
    }

    for (int ib = 0; ib < nbands; ++ib)
    {
        if (bands_picked_[ib])
        {
            std::cout << " Perform band decomposed charge density for band " << ib + 1 << std::endl;

            // (1) calculate the density matrix for a partuclar
            // band, whenever it is occupied or not.

#ifdef __MPI
            this->idmatrix(ib, nspin, nelec, nlocal, wg);
#endif
            // (2) zero out of charge density array.
            for (int is = 0; is < nspin; ++is)
            {
                ModuleBase::GlobalFunc::ZEROS(rho[is], rhopw_nrxx);
            }

            // (3) calculate charge density for a particular
            // band.
            Gint_inout inout(this->loc->DM, rho, Gint_Tools::job_type::rho);
            gg.cal_gint(&inout);

            // A solution to replace the original implementation of the following code:
            // pelec->charge->save_rho_before_sum_band();
            double** rho_save = new double*[nspin]; // Initialize an array of pointers
            for (int is = 0; is < nspin; is++)
            {
                rho_save[is] = new double[rhopw_nrxx]; // Allocate memory for each internal array
                ModuleBase::GlobalFunc::DCOPY(rho[is], rho_save[is],
                                              rhopw_nrxx); // Copy data after allocation
            }

            std::stringstream ssc;
            ssc << global_out_dir << "BAND" << ib + 1;
            // 0 means definitely output charge density.
            for (int is = 0; is < nspin; ++is)
            {
                ssc << "_SPIN" << is << "_CHG.cube";

                // Use a const vector to store efermi for all spins, replace the original implementation:
                // const double ef_tmp = pelec->eferm.get_efval(is);
                double ef_spin = ef_all_spin[is];
                ModuleIO::write_rho(
#ifdef __MPI
                    bigpw_bz,
                    bigpw_nbz,
                    rhopw_nplane,
                    rhopw_startz_current,
#endif
                    rho_save[is],
                    is,
                    nspin,
                    0,
                    ssc.str(),
                    rhopw_nx,
                    rhopw_ny,
                    rhopw_nz,
                    ef_spin,
                    &(GlobalC::ucell));
            }

            // Release memory of rho_save
            for (int is = 0; is < nspin; is++)
            {
                delete[] rho_save[is]; // Release memory of each internal array
            }
            delete[] rho_save; // Release memory of the array of pointers
        }
    }

    return;
}

#ifdef __MPI
void IState_Charge::idmatrix(const int& ib,
                             const int nspin,
                             const double nelec,
                             const int nlocal,
                             const ModuleBase::matrix& wg)
{
    ModuleBase::TITLE("IState_Charge", "idmatrix");
    assert(wg.nr == nspin);

    int fermi_band = static_cast<int>((nelec + 1) / 2 + 1.0e-8);

    for (int is = 0; is < nspin; ++is)
    {
        std::vector<double> wg_local(this->loc->ParaV->ncol, 0.0);
        const int ib_local = this->loc->ParaV->global2local_col(ib);

        if (ib_local >= 0)
        {
            // For unoccupied bands, use occupation of HOMO
            wg_local[ib_local] = (ib < fermi_band) ? wg(is, ib) : wg(is, fermi_band - 1);
        }

        // wg_wfc(ib,iw) = wg[ib] * wfc(ib,iw);
        this->psi_gamma->fix_k(is);
        psi::Psi<double> wg_wfc(*this->psi_gamma, 1);

        for (int ir = 0; ir < wg_wfc.get_nbands(); ++ir)
        {
            BlasConnector::scal(wg_wfc.get_nbasis(), wg_local[ir], wg_wfc.get_pointer() + ir * wg_wfc.get_nbasis(), 1);
        }

        // dm(iw1,iw2) = wfc(ib,iw1).T * wg_wfc(ib,iw2)
        const double one_float = 1.0, zero_float = 0.0;
        const int one_int = 1;
        const char N_char = 'N', T_char = 'T';

        this->loc->dm_gamma.at(is).create(wg_wfc.get_nbands(), wg_wfc.get_nbasis());

        pdgemm_(&N_char,
                &T_char,
                &nlocal,
                &nlocal,
                &wg.nc,
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

    std::cout << " Finished calculating dm_2d." << std::endl;
    this->loc->cal_dk_gamma_from_2D_pub();
    std::cout << " Finished converting dm_2d to dk_gamma." << std::endl;
}
#endif
