#include "elecstate.h"

#include "module_base/global_variable.h"
#include "module_base/memory.h"
#include "module_base/tool_title.h"
#include "occupy.h"
#include "module_base/parallel_reduce.h"

namespace elecstate
{

const double* ElecState::getRho(int spin) const
{
    // hamilt::MatrixBlock<double> temp{&(this->charge->rho[spin][0]), 1, this->charge->nrxx}; //
    // this->chr->get_nspin(), this->chr->get_nrxx()};
    return &(this->charge->rho[spin][0]);
}

void ElecState::fixed_weights(const std::vector<double>& ocp_kb)
{

    int num = 0;
    num = this->klist->nks * GlobalV::NBANDS;
    if (num != ocp_kb.size())
    {
        ModuleBase::WARNING_QUIT("ElecState::fixed_weights",
                                 "size of occupation array is wrong , please check ocp_set");
    }

    double num_elec = 0.0;
    for (int i = 0; i < ocp_kb.size(); i++)
    {
        num_elec += ocp_kb[i];
    }
    if (std::abs(num_elec - GlobalV::nelec) > 1.0e-5)
    {
        ModuleBase::WARNING_QUIT("ElecState::fixed_weights",
                                 "total number of occupations is wrong , please check ocp_set");
    }

    for (int ik = 0; ik < this->wg.nr; ik++)
    {
        for (int ib = 0; ib < this->wg.nc; ib++)
        {
            this->wg(ik, ib) = ocp_kb[ik * this->wg.nc + ib];
        }
    }
    this->skip_weights = true;
}

void ElecState::init_nelec_spin()
{
    this->nelec_spin.resize(GlobalV::NSPIN);
    if (GlobalV::NSPIN == 2)
    {
        // in fact, when TWO_EFERMI(nupdown in INPUT is not 0.0), nelec_spin will be fixed.
        this->nelec_spin[0] = (GlobalV::nelec + GlobalV::nupdown) / 2.0;
        this->nelec_spin[1] = (GlobalV::nelec - GlobalV::nupdown) / 2.0;
    }
}

void ElecState::calculate_weights()
{
    ModuleBase::TITLE("ElecState", "calculate_weights");
    if (this->skip_weights)
    {
        return;
    }

    int nbands = this->ekb.nc;
    int nks = this->ekb.nr;

    if (!Occupy::use_gaussian_broadening && !Occupy::fixed_occupations)
    {
        if (GlobalV::TWO_EFERMI)
        {
            Occupy::iweights(nks,
                             this->klist->wk,
                             nbands,
                             this->nelec_spin[0],
                             this->ekb,
                             this->eferm.ef_up,
                             this->wg,
                             0,
                             this->klist->isk);
            Occupy::iweights(nks,
                             this->klist->wk,
                             nbands,
                             this->nelec_spin[1],
                             this->ekb,
                             this->eferm.ef_dw,
                             this->wg,
                             1,
                             this->klist->isk);
            // ef = ( ef_up + ef_dw ) / 2.0_dp need??? mohan add 2012-04-16
        }
        else
        {
            // -1 means don't need to consider spin.
            Occupy::iweights(nks,
                             this->klist->wk,
                             nbands,
                             GlobalV::nelec,
                             this->ekb,
                             this->eferm.ef,
                             this->wg,
                             -1,
                             this->klist->isk);
        }
    }
    else if (Occupy::use_gaussian_broadening)
    {
        if (GlobalV::TWO_EFERMI)
        {
            double demet_up = 0.0;
            double demet_dw = 0.0;
            Occupy::gweights(nks,
                             this->klist->wk,
                             nbands,
                             this->nelec_spin[0],
                             Occupy::gaussian_parameter,
                             Occupy::gaussian_type,
                             this->ekb,
                             this->eferm.ef_up,
                             demet_up,
                             this->wg,
                             0,
                             this->klist->isk);
            Occupy::gweights(nks,
                             this->klist->wk,
                             nbands,
                             this->nelec_spin[1],
                             Occupy::gaussian_parameter,
                             Occupy::gaussian_type,
                             this->ekb,
                             this->eferm.ef_dw,
                             demet_dw,
                             this->wg,
                             1,
                             this->klist->isk);
            this->f_en.demet = demet_up + demet_dw;
        }
        else
        {
            // -1 means is no related to spin.
            Occupy::gweights(nks,
                             this->klist->wk,
                             nbands,
                             GlobalV::nelec,
                             Occupy::gaussian_parameter,
                             Occupy::gaussian_type,
                             this->ekb,
                             this->eferm.ef,
                             this->f_en.demet,
                             this->wg,
                             -1,
                             this->klist->isk);
        }
#ifdef __MPI
        // qianrui fix a bug on 2021-7-21
        Parallel_Reduce::reduce_double_allpool(this->f_en.demet);
#endif
    }
    else if (Occupy::fixed_occupations)
    {
        ModuleBase::WARNING_QUIT("calculate_weights", "other occupations, not implemented");
    }

    return;
}

void ElecState::calEBand()
{
    ModuleBase::TITLE("ElecState", "calEBand");
    // calculate ebands using wg and ekb
    double eband = 0.0;
#ifdef _OPENMP
#pragma omp parallel for collapse(2) reduction(+:eband)
#endif
    for (int ik = 0; ik < this->ekb.nr; ++ik)
    {
        for (int ibnd = 0; ibnd < this->ekb.nc; ibnd++)
        {
            eband += this->ekb(ik, ibnd) * this->wg(ik, ibnd);
        }
    }
    this->f_en.eband = eband;
    if (GlobalV::KPAR != 1 && GlobalV::ESOLVER_TYPE != "sdft")
    {
        //==================================
        // Reduce all the Energy in each cpu
        //==================================
        this->f_en.eband /= GlobalV::NPROC_IN_POOL;
#ifdef __MPI
        Parallel_Reduce::reduce_double_all(this->f_en.eband);
#endif
    }
    return;
}

void ElecState::init_scf(const int istep, const ModuleBase::ComplexMatrix& strucfac)
{
    //---------Charge part-----------------
    // core correction potential.
    this->charge->set_rho_core(strucfac);

    //--------------------------------------------------------------------
    // (2) other effective potentials need charge density,
    // choose charge density from ionic step 0.
    //--------------------------------------------------------------------
    if (istep == 0)
    {
        this->charge->init_rho(this->eferm, strucfac);
    }

    // renormalize the charge density
    this->charge->renormalize_rho();

    //---------Potential part--------------
    this->pot->init_pot(istep, this->charge);
}

void ElecState::init_ks(Charge* chg_in, // pointer for class Charge
                        const K_Vectors* klist_in,
                        int nk_in,
                        ModulePW::PW_Basis* rhopw_in,
                        const ModulePW::PW_Basis_Big* bigpw_in)
{
    this->charge = chg_in;
    this->klist = klist_in;
    this->charge->set_rhopw(rhopw_in);
    this->bigpw = bigpw_in;
    // init nelec_spin with nelec and nupdown
    this->init_nelec_spin();
    // autoset and check GlobalV::NBANDS, nelec_spin is used when NSPIN==2
    this->cal_nbands();
    // initialize ekb and wg
    this->ekb.create(nk_in, GlobalV::NBANDS);
    this->wg.create(nk_in, GlobalV::NBANDS);
}

void ElecState::cal_nbands()
{
    if (GlobalV::ESOLVER_TYPE == "sdft") // qianrui 2021-2-20
    {
        return;
    }
    //=======================================
    // calculate number of bands (setup.f90)
    //=======================================
    double occupied_bands = static_cast<double>(GlobalV::nelec / ModuleBase::DEGSPIN);
    if (GlobalV::LSPINORB == 1)
        occupied_bands = static_cast<double>(GlobalV::nelec);

    if ((occupied_bands - std::floor(occupied_bands)) > 0.0)
    {
        occupied_bands = std::floor(occupied_bands) + 1.0; // mohan fix 2012-04-16
    }

    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "occupied bands", occupied_bands);

    // mohan add 2010-09-04
    // std::cout << "nbands(this-> = " <<GlobalV::NBANDS <<std::endl;
    if (GlobalV::NBANDS == occupied_bands)
    {
        if (Occupy::gauss())
        {
            ModuleBase::WARNING_QUIT("ElecState::cal_nbands", "for smearing, num. of bands > num. of occupied bands");
        }
    }

    if (GlobalV::NBANDS == 0)
    {
        if (GlobalV::NSPIN == 1)
        {
            const int nbands1 = static_cast<int>(occupied_bands) + 10;
            const int nbands2 = static_cast<int>(1.2 * occupied_bands) + 1;
            GlobalV::NBANDS = std::max(nbands1, nbands2);
            if (GlobalV::BASIS_TYPE != "pw")
                GlobalV::NBANDS = std::min(GlobalV::NBANDS, GlobalV::NLOCAL);
        }
        else if (GlobalV::NSPIN == 4)
        {
            const int nbands3 = GlobalV::nelec + 20;
            const int nbands4 = static_cast<int>(1.2 * GlobalV::nelec) + 1;
            GlobalV::NBANDS = std::max(nbands3, nbands4);
            if (GlobalV::BASIS_TYPE != "pw")
                GlobalV::NBANDS = std::min(GlobalV::NBANDS, GlobalV::NLOCAL);
        }
        else if (GlobalV::NSPIN == 2)
        {
            const double max_occ = std::max(this->nelec_spin[0], this->nelec_spin[1]);
            const int nbands3 = static_cast<int>(max_occ) + 11;
            const int nbands4 = static_cast<int>(1.2 * max_occ) + 1;
            GlobalV::NBANDS = std::max(nbands3, nbands4);
            if (GlobalV::BASIS_TYPE != "pw")
                GlobalV::NBANDS = std::min(GlobalV::NBANDS, GlobalV::NLOCAL);
        }
        ModuleBase::GlobalFunc::AUTO_SET("NBANDS", GlobalV::NBANDS);
    }
    // else if ( GlobalV::CALCULATION=="scf" || GlobalV::CALCULATION=="md" || GlobalV::CALCULATION=="relax") //pengfei
    // 2014-10-13
    else
    {
        if (GlobalV::NBANDS < occupied_bands)
            ModuleBase::WARNING_QUIT("unitcell", "Too few bands!");
        if (GlobalV::NSPIN == 2)
        {
            if (GlobalV::NBANDS < this->nelec_spin[0])
            {
                ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "nelec_up", this->nelec_spin[0]);
                ModuleBase::WARNING_QUIT("ElecState::cal_nbands", "Too few spin up bands!");
            }
            if (GlobalV::NBANDS < this->nelec_spin[1])
            {
                ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "nelec_down", this->nelec_spin[1]);
                ModuleBase::WARNING_QUIT("ElecState::cal_nbands", "Too few spin down bands!");
            }
        }
    }

    // mohan update 2021-02-19
    // mohan add 2011-01-5
    if (GlobalV::BASIS_TYPE == "lcao" || GlobalV::BASIS_TYPE == "lcao_in_pw")
    {
        if (GlobalV::NBANDS > GlobalV::NLOCAL)
        {
            ModuleBase::WARNING_QUIT("ElecState::cal_nbandsc", "NLOCAL < NBANDS");
        }
        else
        {
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "NLOCAL", GlobalV::NLOCAL);
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "NBANDS", GlobalV::NBANDS);
        }
    }

    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "NBANDS", GlobalV::NBANDS);
}
} // namespace elecstate
