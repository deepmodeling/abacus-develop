#include "elecstate_lcao.h"

#include "math_tools.h"
#include "module_base/timer.h"
#include "src_lcao/grid_technique.h"

namespace elecstate
{
int ElecStateLCAO::out_wfc_lcao = 0;

// for gamma_only(double case) and multi-k(complex<double> case)
template <typename T> void ElecStateLCAO::cal_dm(const ModuleBase::matrix& wg, const psi::Psi<T>& wfc, psi::Psi<T>& dm)
{
    ModuleBase::TITLE("ElecStateLCAO", "cal_dm");

    dm.resize(wfc.get_nk(), this->loc->ParaV->ncol, this->loc->ParaV->nrow);
    const int nbands_local = wfc.get_nbands();
    const int nbasis_local = wfc.get_nbasis();

    // dm = wfc.T * wg * wfc.conj()
    // dm[is](iw1,iw2) = \sum_{ib} wfc[is](ib,iw1).T * wg(is,ib) * wfc[is](ib,iw2).conj()
    for (int ik = 0; ik < wfc.get_nk(); ++ik)
    {
        wfc.fix_k(ik);
        dm.fix_k(ik);
        // wg_wfc(ib,iw) = wg[ib] * wfc(ib,iw);
        psi::Psi<T> wg_wfc(wfc, 1);

        int ib_global = 0;
        for (int ib_local = 0; ib_local < nbands_local; ++ib_local)
        {
            while (ib_local != this->loc->ParaV->trace_loc_col[ib_global])
            {
                ++ib_global;
                if (ib_global >= wg.nc)
                {
                    ModuleBase::WARNING_QUIT("ElecStateLCAO::cal_dm", "please check trace_loc_col!");
                }
            }
            const double wg_local = wg(ik, ib_global);
            T* wg_wfc_pointer = &(wg_wfc(0, ib_local, 0));
            BlasConnector::scal(nbasis_local, wg_local, wg_wfc_pointer, 1);
        }

        // C++: dm(iw1,iw2) = wfc(ib,iw1).T * wg_wfc(ib,iw2)
#ifdef __MPI
        psiMulPsiMpi(wg_wfc, wfc, dm, this->loc->ParaV->desc_wfc, this->loc->ParaV->desc);
#else
        psiMulPsi(wg_wfc, wfc, dm);
#endif
    }

    return;
}

// multi-k case
void ElecStateLCAO::psiToRho(const psi::Psi<std::complex<double>>& psi)
{
    ModuleBase::TITLE("ElecStateLCAO", "psiToRho");
    ModuleBase::timer::tick("ElecStateLCAO", "psiToRho");

    ModuleBase::GlobalFunc::NOTE("Calculate the density matrix.");

    // this part for calculating dm_k in 2d-block format, not used for charge now
    psi::Psi<std::complex<double>> dm_k_2d(psi.get_nk(), psi.get_nbasis(), psi.get_nbasis());

    if (GlobalV::KS_SOLVER == "genelpa" || GlobalV::KS_SOLVER == "scalapack_gvx"
        || GlobalV::KS_SOLVER == "lapack") // Peize Lin test 2019-05-15
    {
        this->cal_dm(this->wg, psi, dm_k_2d);
    }

    // this part for steps:
    // 1. psi_k transform from 2d-block to grid format
    // 2. psi_k_grid -> DM_R
    // 3. DM_R -> rho(r)
    if (GlobalV::KS_SOLVER == "genelpa" || GlobalV::KS_SOLVER == "scalapack_gvx")
    {
        for (int ik = 0; ik < psi.get_nk(); ik++)
        {
            psi.fix_k(ik);
            this->lowf->wfc_2d_to_grid(ElecStateLCAO::out_wfc_lcao, psi.get_pointer(), this->lowf->wfc_k_grid[ik], ik);
        }
    }

    this->loc->cal_dk_k(GlobalC::GridT);
    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        ModuleBase::GlobalFunc::ZEROS(this->charge->rho[is], this->charge->nrxx); // mohan 2009-11-10
    }

    //------------------------------------------------------------
    // calculate the charge density on real space grid.
    //------------------------------------------------------------

    ModuleBase::GlobalFunc::NOTE("Calculate the charge on real space grid!");
    this->uhm->GK.cal_rho_k(this->loc->DM_R);

    this->charge->renormalize_rho();

    ModuleBase::timer::tick("ElecStateLCAO", "psiToRho");
    return;
}

// Gamma_only case
void ElecStateLCAO::psiToRho(const psi::Psi<double>& psi)
{
    ModuleBase::TITLE("ElecStateLCAO", "psiToRho");
    ModuleBase::timer::tick("ElecStateLCAO", "psiToRho");

    this->calculate_weights();

    if (GlobalV::KS_SOLVER == "genelpa" || GlobalV::KS_SOLVER == "scalapack_gvx" || GlobalV::KS_SOLVER == "lapack")
    {
        ModuleBase::timer::tick("ElecStateLCAO", "cal_dm_2d");

        psi::Psi<double> dm_gamma_2d(psi.get_nk(), psi.get_nbasis(), psi.get_nbasis());
        // caution:wfc and dm
        this->cal_dm(this->wg, psi, dm_gamma_2d);

        ModuleBase::timer::tick("ElecStateLCAO", "cal_dm_2d");

        for (int ik = 0; ik < psi.get_nk(); ++ik)
        {
            // for gamma_only case, no convertion occured, just for print.
            if (GlobalV::KS_SOLVER == "genelpa" || GlobalV::KS_SOLVER == "scalapack_gvx")
            {
                psi.fix_k(ik);
                double** wfc_grid = nullptr; // output but not do "2d-to-grid" conversion
                this->lowf->wfc_2d_to_grid(ElecStateLCAO::out_wfc_lcao, psi.get_pointer(), wfc_grid);
            }
            this->loc->dm2dToGrid(dm_gamma_2d, this->loc->DM[ik]); // transform dm_gamma[is].c to this->loc->DM[is]
        }
    }

    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        ModuleBase::GlobalFunc::ZEROS(this->charge->rho[is], this->charge->nrxx); // mohan 2009-11-10
    }

    //------------------------------------------------------------
    // calculate the charge density on real space grid.
    //------------------------------------------------------------
    ModuleBase::GlobalFunc::NOTE("Calculate the charge on real space grid!");
    this->uhm->GG.cal_rho(this->loc->DM);

    this->charge->renormalize_rho();

    ModuleBase::timer::tick("ElecStateLCAO", "psiToRho");
    return;
}

} // namespace elecstate