#include "charge_mixing.h"

#include "module_parameter/parameter.h"
#include "module_base/timer.h"

void Charge_Mixing::allocate_mixing_dmr(const int nnr)
{
    // Note that: we cannot allocate memory for dmr_mdata in set_mixing.
    // since the size of dmr_mdata is given by the size of HContainer.nnr, which is calculated in DensityMatrix::init_DMR().
    // and DensityMatrix::init_DMR() is called in beforescf(). While set_mixing() is called in ESolver_KS::Init().
    ModuleBase::TITLE("Charge_Mixing", "allocate_mixing_dmr");
    ModuleBase::timer::tick("Charge_Mixing", "allocate_mixing_dmr");
    //
    const int dmr_nspin = (PARAM.inp.nspin == 2) ? 2 : 1;
    // allocate memory for dmr_mdata
    if (PARAM.inp.scf_thr_type == 1)
    {
        ModuleBase::WARNING_QUIT("Charge_Mixing", "This Mixing of Density Matrix is not supported for PW basis yet");
    }
    else if (PARAM.inp.scf_thr_type == 2)
    {
        this->mixing->init_mixing_data(this->dmr_mdata, nnr * dmr_nspin, sizeof(double));
    }

    this->dmr_mdata.reset();
    ModuleBase::timer::tick("Charge_Mixing", "allocate_mixing_dmr");

    return;
}

void Charge_Mixing::mix_dmr(elecstate::DensityMatrix<double, double>* DM)
{
    // Notice that DensityMatrix object is a Template class
    ModuleBase::TITLE("Charge_Mixing", "mix_dmr");
    ModuleBase::timer::tick("Charge_Mixing", "mix_dmr");
    //
    std::vector<hamilt::HContainer<double>*> dmr = DM->get_DMR_vector();
    std::vector<std::vector<double>>& dmr_save = DM->get_DMR_save();
    //
    //const int dmr_nspin = (PARAM.inp.nspin == 2) ? 2 : 1;
    double* dmr_in = nullptr;
    double* dmr_out = nullptr;
    if (PARAM.inp.nspin == 1 || PARAM.inp.nspin == 4)
    {
        dmr_in = dmr_save[0].data();
        dmr_out = dmr[0]->get_wrapper();
        this->mixing->push_data(this->dmr_mdata, dmr_in, dmr_out, nullptr, false);    
        this->mixing->mix_data(this->dmr_mdata, dmr_out);
    }
    else if (PARAM.inp.nspin == 2)
    {
        // magnetic density matrix
        double* dmr_mag = nullptr;
        double* dmr_mag_save = nullptr; 
        const int nnr = dmr[0]->get_nnr();
        // allocate dmr_mag[is*nnrx] and dmr_mag_save[is*nnrx]
        dmr_mag = new double[nnr * PARAM.inp.nspin];
        dmr_mag_save = new double[nnr * PARAM.inp.nspin];
        ModuleBase::GlobalFunc::ZEROS(dmr_mag, nnr * PARAM.inp.nspin);
        ModuleBase::GlobalFunc::ZEROS(dmr_mag_save, nnr * PARAM.inp.nspin);
        double* dmr_up;
        double* dmr_down;
        // tranfer dmr into dmr_mag
        dmr_up = dmr[0]->get_wrapper();
        dmr_down = dmr[1]->get_wrapper();
        for (int ir = 0; ir < nnr; ir++)
        {
            dmr_mag[ir] = dmr_up[ir] + dmr_down[ir];
            dmr_mag[ir + nnr] = dmr_up[ir] - dmr_down[ir];
        }
        // tranfer dmr_save into dmr_mag_save
        dmr_up = dmr_save[0].data();
        dmr_down = dmr_save[1].data();
        for (int ir = 0; ir < nnr; ir++)
        {
            dmr_mag_save[ir] = dmr_up[ir] + dmr_down[ir];
            dmr_mag_save[ir + nnr] = dmr_up[ir] - dmr_down[ir];
        }
        //
        dmr_in = dmr_mag_save;
        dmr_out = dmr_mag;
        // no kerker in mixing_dmr
        //auto screen = std::bind(&Charge_Mixing::Kerker_screen_real, this, std::placeholders::_1);
        auto twobeta_mix
            = [this, nnr](double* out, const double* in, const double* sres) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 256)
#endif
            for (int i = 0; i < nnr; ++i)
            {
                out[i] = in[i] + this->mixing_beta * sres[i];
            }
            // magnetism
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 256)
#endif
            for (int i = nnr; i < 2 * nnr; ++i)
            {
                out[i] = in[i] + this->mixing_beta_mag * sres[i];
            }
        };
        this->mixing->push_data(this->dmr_mdata, dmr_in, dmr_out, nullptr, twobeta_mix, false);
        //auto inner_product
        //    = std::bind(&Charge_Mixing::inner_product_real, this, std::placeholders::_1, std::placeholders::_2);
        //this->mixing->cal_coef(this->rho_mdata, inner_product);
        this->mixing->mix_data(this->dmr_mdata, dmr_out);
        // get new dmr from dmr_mag
        dmr_up = dmr[0]->get_wrapper();
        dmr_down = dmr[1]->get_wrapper();
        for (int is = 0; is < PARAM.inp.nspin; is++)
        {
            ModuleBase::GlobalFunc::ZEROS(dmr_up, nnr);
            ModuleBase::GlobalFunc::ZEROS(dmr_down, nnr);
        }
        for (int ir = 0; ir < nnr; ir++)
        {
            dmr_up[ir] = 0.5 * (dmr_mag[ir] + dmr_mag[ir+nnr]);
            dmr_down[ir] = 0.5 * (dmr_mag[ir] - dmr_mag[ir+nnr]);
        }
        // delete
        delete[] dmr_mag;
        delete[] dmr_mag_save;
    }

    ModuleBase::timer::tick("Charge_Mixing", "mix_dmr");

    return;
}

void Charge_Mixing::mix_dmr(elecstate::DensityMatrix<std::complex<double>, double>* DM)
{
    // Notice that DensityMatrix object is a Template class
    ModuleBase::TITLE("Charge_Mixing", "mix_dmr");
    ModuleBase::timer::tick("Charge_Mixing", "mix_dmr");
    //
    std::vector<hamilt::HContainer<double>*> dmr = DM->get_DMR_vector();
    std::vector<std::vector<double>>& dmr_save = DM->get_DMR_save();
    //
    //const int dmr_nspin = (PARAM.inp.nspin == 2) ? 2 : 1;
    double* dmr_in = nullptr;
    double* dmr_out = nullptr;
    if (PARAM.inp.nspin == 1 || PARAM.inp.nspin == 4)
    {
        dmr_in = dmr_save[0].data();
        dmr_out = dmr[0]->get_wrapper();
        this->mixing->push_data(this->dmr_mdata, dmr_in, dmr_out, nullptr, false);    
        this->mixing->mix_data(this->dmr_mdata, dmr_out);
    }
    else if (PARAM.inp.nspin == 2)
    {
        // magnetic density matrix
        double* dmr_mag = nullptr;
        double* dmr_mag_save = nullptr; 
        const int nnr = dmr[0]->get_nnr();
        // allocate dmr_mag[is*nnrx] and dmr_mag_save[is*nnrx]
        dmr_mag = new double[nnr * PARAM.inp.nspin];
        dmr_mag_save = new double[nnr * PARAM.inp.nspin];
        ModuleBase::GlobalFunc::ZEROS(dmr_mag, nnr * PARAM.inp.nspin);
        ModuleBase::GlobalFunc::ZEROS(dmr_mag_save, nnr * PARAM.inp.nspin);
        double* dmr_up;
        double* dmr_down;
        // tranfer dmr into dmr_mag
        dmr_up = dmr[0]->get_wrapper();
        dmr_down = dmr[1]->get_wrapper();
        for (int ir = 0; ir < nnr; ir++)
        {
            dmr_mag[ir] = dmr_up[ir] + dmr_down[ir];
            dmr_mag[ir + nnr] = dmr_up[ir] - dmr_down[ir];
        }
        // tranfer dmr_save into dmr_mag_save
        dmr_up = dmr_save[0].data();
        dmr_down = dmr_save[1].data();
        for (int ir = 0; ir < nnr; ir++)
        {
            dmr_mag_save[ir] = dmr_up[ir] + dmr_down[ir];
            dmr_mag_save[ir + nnr] = dmr_up[ir] - dmr_down[ir];
        }
        //
        dmr_in = dmr_mag_save;
        dmr_out = dmr_mag;
        // no kerker in mixing_dmr
        //auto screen = std::bind(&Charge_Mixing::Kerker_screen_real, this, std::placeholders::_1);
        auto twobeta_mix
            = [this, nnr](double* out, const double* in, const double* sres) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 256)
#endif
            for (int i = 0; i < nnr; ++i)
            {
                out[i] = in[i] + this->mixing_beta * sres[i];
            }
            // magnetism
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 256)
#endif
            for (int i = nnr; i < 2 * nnr; ++i)
            {
                out[i] = in[i] + this->mixing_beta_mag * sres[i];
            }
        };
        this->mixing->push_data(this->dmr_mdata, dmr_in, dmr_out, nullptr, twobeta_mix, false);
        //auto inner_product
        //    = std::bind(&Charge_Mixing::inner_product_real, this, std::placeholders::_1, std::placeholders::_2);
        //this->mixing->cal_coef(this->rho_mdata, inner_product);
        this->mixing->mix_data(this->dmr_mdata, dmr_out);
        // get new dmr from dmr_mag
        dmr_up = dmr[0]->get_wrapper();
        dmr_down = dmr[1]->get_wrapper();
        for (int is = 0; is < PARAM.inp.nspin; is++)
        {
            ModuleBase::GlobalFunc::ZEROS(dmr_up, nnr);
            ModuleBase::GlobalFunc::ZEROS(dmr_down, nnr);
        }
        for (int ir = 0; ir < nnr; ir++)
        {
            dmr_up[ir] = 0.5 * (dmr_mag[ir] + dmr_mag[ir+nnr]);
            dmr_down[ir] = 0.5 * (dmr_mag[ir] - dmr_mag[ir+nnr]);
        }
        // delete
        delete[] dmr_mag;
        delete[] dmr_mag_save;
    }

    ModuleBase::timer::tick("Charge_Mixing", "mix_dmr");

    return;
}