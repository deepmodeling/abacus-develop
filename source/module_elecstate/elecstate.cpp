#include "elecstate.h"

#include "module_base/global_variable.h"
#include "module_base/tool_title.h"
#include "src_parallel/parallel_reduce.h"
#include "src_pw/global.h"
#include "src_pw/occupy.h"

namespace elecstate
{

const double* ElecState::getRho(int spin) const
{
    // hamilt::MatrixBlock<double> temp{&(this->charge->rho[spin][0]), 1, this->charge->nrxx}; //
    // this->chr->get_nspin(), this->chr->get_nrxx()};
    return &(this->charge->rho[spin][0]);
}

void ElecState::calculate_weights(void)
{
    ModuleBase::TITLE("ElecState", "calculate_weights");

    // for test
    //	std::cout << " gaussian_broadening = " << use_gaussian_broadening << std::endl;
    //	std::cout << " tetrahedron_method = " << use_tetrahedron_method << std::endl;
    //	std::cout << " fixed_occupations = " << fixed_occupations << std::endl;
    double** ekb_tmp = new double*[this->ekb.nr];
    for (int i = 0; i < this->ekb.nr; ++i)
    {
        ekb_tmp[i] = &(this->ekb(i, 0));
    }
    int nbands = this->ekb.nc;

    if (GlobalV::KS_SOLVER == "selinv")
    {
        GlobalV::ofs_running << " Could not calculate occupation." << std::endl;
        return;
    }

    if (!Occupy::use_gaussian_broadening && !Occupy::use_tetrahedron_method && !Occupy::fixed_occupations)
    {
        if (GlobalV::TWO_EFERMI)
        {
            Occupy::iweights(GlobalC::kv.nks,
                             GlobalC::kv.wk,
                             nbands,
                             GlobalC::ucell.magnet.get_nelup(),
                             ekb_tmp,
                             GlobalC::en.ef_up,
                             this->wg,
                             0,
                             GlobalC::kv.isk);
            Occupy::iweights(GlobalC::kv.nks,
                             GlobalC::kv.wk,
                             nbands,
                             GlobalC::ucell.magnet.get_neldw(),
                             ekb_tmp,
                             GlobalC::en.ef_dw,
                             this->wg,
                             1,
                             GlobalC::kv.isk);
            // ef = ( ef_up + ef_dw ) / 2.0_dp need??? mohan add 2012-04-16
        }
        else
        {
            // -1 means don't need to consider spin.
            Occupy::iweights(GlobalC::kv.nks,
                             GlobalC::kv.wk,
                             nbands,
                             GlobalC::CHR.nelec,
                             ekb_tmp,
                             this->ef,
                             this->wg,
                             -1,
                             GlobalC::kv.isk);
        }
    }
    else if (Occupy::use_tetrahedron_method)
    {
        ModuleBase::WARNING_QUIT("calculate_weights", "not implemented yet,coming soon!");
        //		if(my_rank == 0)
        //		{
        //			tweights(GlobalC::kv.nkstot, nspin, nbands, GlobalC::CHR.nelec, ntetra,tetra, GlobalC::wf.et,
        // this->ef, this->wg);
        //		}
    }
    else if (Occupy::use_gaussian_broadening)
    {
        if (GlobalV::TWO_EFERMI)
        {
            double demet_up = 0.0;
            double demet_dw = 0.0;
            Occupy::gweights(GlobalC::kv.nks,
                             GlobalC::kv.wk,
                             nbands,
                             GlobalC::ucell.magnet.get_nelup(),
                             Occupy::gaussian_parameter,
                             Occupy::gaussian_type,
                             ekb_tmp,
                             GlobalC::en.ef_up,
                             demet_up,
                             this->wg,
                             0,
                             GlobalC::kv.isk);
            Occupy::gweights(GlobalC::kv.nks,
                             GlobalC::kv.wk,
                             nbands,
                             GlobalC::ucell.magnet.get_neldw(),
                             Occupy::gaussian_parameter,
                             Occupy::gaussian_type,
                             ekb_tmp,
                             GlobalC::en.ef_dw,
                             demet_dw,
                             this->wg,
                             1,
                             GlobalC::kv.isk);
            GlobalC::en.demet = demet_up + demet_dw;
        }
        else
        {
            // -1 means is no related to spin.
            Occupy::gweights(GlobalC::kv.nks,
                             GlobalC::kv.wk,
                             nbands,
                             GlobalC::CHR.nelec,
                             Occupy::gaussian_parameter,
                             Occupy::gaussian_type,
                             ekb_tmp,
                             this->ef,
                             GlobalC::en.demet,
                             this->wg,
                             -1,
                             GlobalC::kv.isk);
        }

        // qianrui fix a bug on 2021-7-21
        Parallel_Reduce::reduce_double_allpool(GlobalC::en.demet);
    }
    else if (Occupy::fixed_occupations)
    {
        // fix occupations need nelup and neldw.
        // mohan add 2011-04-03
        this->ef = -1.0e+20;
        for (int ik = 0; ik < GlobalC::kv.nks; ik++)
        {
            for (int ibnd = 0; ibnd < nbands; ibnd++)
            {
                if (this->wg(ik, ibnd) > 0.0)
                {
                    this->ef = std::max(this->ef, ekb_tmp[ik][ibnd]);
                }
            }
        }
    }

    if (GlobalV::TWO_EFERMI)
    {
        Parallel_Reduce::gather_max_double_all(GlobalC::en.ef_up);
        Parallel_Reduce::gather_max_double_all(GlobalC::en.ef_dw);
    }
    else
    {
        double ebotom = ekb_tmp[0][0];
        double etop = ekb_tmp[0][0];
        for (int ik = 0; ik < GlobalC::kv.nks; ik++)
        {
            for (int ib = 0; ib < nbands; ib++)
            {
                ebotom = min(ebotom, ekb_tmp[ik][ib]);
                etop = max(etop, ekb_tmp[ik][ib]);
            }
        }

        // parallel
        Parallel_Reduce::gather_max_double_all(this->ef);
        Parallel_Reduce::gather_max_double_all(etop);
        Parallel_Reduce::gather_min_double_all(ebotom);

        // not parallel yet!
        //		OUT(GlobalV::ofs_running,"Top    Energy (eV)", etop * ModuleBase::Ry_to_eV);
        //      OUT(GlobalV::ofs_running,"Fermi  Energy (eV)", this->ef * ModuleBase::Ry_to_eV);
        //		OUT(GlobalV::ofs_running,"Bottom Energy (eV)", ebotom * ModuleBase::Ry_to_eV);
        //		OUT(GlobalV::ofs_running,"Range  Energy (eV)", etop-ebotom * ModuleBase::Ry_to_eV);
    }

    delete[] ekb_tmp;

    return;
}

double ElecState::eBandK(const int& ik)
{
    ModuleBase::TITLE("ElecStatePW", "eBandK");
    double eband_k = 0.0;
    for (int ibnd = 0; ibnd < this->ekb.nc; ibnd++)
    {
        eband_k += this->ekb(ik, ibnd) * this->wg(ik, ibnd);
    }
    return eband_k;
}

void ElecState::print_band(const int& ik, const int& printe, const int& iter)
{
    // check the band energy.
    bool wrong = false;
    int nbands = this->ekb.nc;
    for (int ib = 0; ib < nbands; ++ib)
    {
        if (abs(this->ekb(ik, ib)) > 1.0e10)
        {
            GlobalV::ofs_warning << " ik=" << ik + 1 << " ib=" << ib + 1 << " " << this->ekb(ik, ib) << " Ry"
                                 << std::endl;
            wrong = true;
        }
    }
    if (wrong)
    {
        ModuleBase::WARNING_QUIT("Threshold_Elec::print_eigenvalue", "Eigenvalues are too large!");
    }

    if (GlobalV::MY_RANK == 0)
    {
        // if( GlobalV::DIAGO_TYPE == "selinv" ) xiaohui modify 2013-09-02
        if (GlobalV::KS_SOLVER == "selinv") // xiaohui add 2013-09-02
        {
            GlobalV::ofs_running << " No eigenvalues are available for selected inversion methods." << std::endl;
        }
        else
        {
            if (printe > 0 && ((iter + 1) % printe == 0))
            {
                //	NEW_PART("ENERGY BANDS (Rydberg), (eV)");
                GlobalV::ofs_running << std::setprecision(6);
                GlobalV::ofs_running << " Energy (eV) & Occupations  for spin=" << GlobalV::CURRENT_SPIN + 1
                                     << " K-point=" << ik + 1 << std::endl;
                GlobalV::ofs_running << std::setiosflags(ios::showpoint);
                for (int ib = 0; ib < nbands; ib++)
                {
                    GlobalV::ofs_running << " " << std::setw(6) << ib + 1 << std::setw(15)
                                         << this->ekb(ik, ib) * ModuleBase::Ry_to_eV;
                    // for the first electron iteration, we don't have the energy
                    // spectrum, so we can't get the occupations.
                    GlobalV::ofs_running << std::setw(15) << this->wg(ik, ib);
                    GlobalV::ofs_running << std::endl;
                }
            }
        }
    }
    return;
}

} // namespace elecstate