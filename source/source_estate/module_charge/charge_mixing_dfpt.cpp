#include "charge_mixing.h"
#include "source_io/module_parameter/parameter.h"
#include "source_base/timer.h"
#include "source_base/parallel_reduce.h"

void Charge_Mixing::init_mixing_delta_rho(const int length)
{
    if (this->mixing != nullptr)
    {
        this->mixing->init_mixing_data(this->delta_rho_mdata,
                                       length,
                                       sizeof(std::complex<double>));
    }
}

double Charge_Mixing::inner_product_recip_complex(std::complex<double>* rho1,
                                                   std::complex<double>* rho2)
{
    ModuleBase::TITLE("Charge_Mixing", "inner_product_recip_complex");
    ModuleBase::timer::tick("Charge_Mixing", "inner_product_recip_complex");

    const int nspin = PARAM.inp.nspin;
    const int npw = this->rhopw->npw;

    double sum = 0.0;

    // Complex inner product: sum of conj(rho1[i]) * rho2[i]
    // For DFPT, we use a simple inner product without Hartree-like weighting
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : sum)
#endif
    for (int i = 0; i < npw * nspin; ++i)
    {
        sum += (std::conj(rho1[i]) * rho2[i]).real();
    }

#ifdef __MPI
    Parallel_Reduce::reduce_pool(sum);
#endif

    ModuleBase::timer::tick("Charge_Mixing", "inner_product_recip_complex");
    return sum;
}

void Charge_Mixing::Kerker_screen_recip_complex(std::complex<double>* rhog, const int nspin)
{
    // Apply Kerker screening to complex charge density
    // Same logic as Kerker_screen_recip but for complex data
    if (this->mixing_gg0 <= 0.0)
    {
        return;
    }

    const int npw = this->rhopw->npw;
    const double gg0 = this->mixing_gg0;

    for (int is = 0; is < nspin; ++is)
    {
        std::complex<double>* rhog_is = rhog + is * npw;
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 256)
#endif
        for (int ig = 0; ig < npw; ++ig)
        {
            const double gg = this->rhopw->gg[ig];
            const double kerker = gg / (gg + gg0);
            rhog_is[ig] *= kerker;
        }
    }
}

double Charge_Mixing::get_delta_drho(Charge* chr, const int nspin)
{
    ModuleBase::TITLE("Charge_Mixing", "get_delta_drho");
    ModuleBase::timer::tick("Charge_Mixing", "get_delta_drho");

    const int npw = this->rhopw->npw;

    // FFT: real space -> reciprocal space for current and saved delta_rho
    for (int is = 0; is < nspin; ++is)
    {
        this->rhopw->real2recip(chr->delta_rho[is], chr->delta_rhog[is]);
        this->rhopw->real2recip(chr->delta_rho_save[is], chr->delta_rhog_save[is]);
    }

    // Calculate residual: delta_rhog - delta_rhog_save
    std::vector<std::complex<double>> drhog(nspin * npw);

#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 512)
#endif
    for (int is = 0; is < nspin; ++is)
    {
        for (int ig = 0; ig < npw; ++ig)
        {
            drhog[is * npw + ig] = chr->delta_rhog[is][ig] - chr->delta_rhog_save[is][ig];
        }
    }

    // Calculate norm using complex inner product: sum of |drhog|^2 = sum of conj(drhog) * drhog
    double drho = this->inner_product_recip_complex(drhog.data(), drhog.data());

    // Take square root to get L2 norm
    drho = std::sqrt(drho);

    ModuleBase::timer::tick("Charge_Mixing", "get_delta_drho");
    return drho;
}

void Charge_Mixing::mix_delta_rho_recip(Charge* chr, const int nspin)
{
    ModuleBase::TITLE("Charge_Mixing", "mix_delta_rho_recip");
    ModuleBase::timer::tick("Charge_Mixing", "mix_delta_rho_recip");

    const int npw = this->rhopw->npw;

    // Flatten the 2D arrays to 1D for mixing
    // Layout: [spin0: ig0, ig1, ..., ig_{npw-1}], [spin1: ...], ...
    std::vector<std::complex<double>> rhog_in_flat(nspin * npw);
    std::vector<std::complex<double>> rhog_out_flat(nspin * npw);

    for (int is = 0; is < nspin; ++is)
    {
        for (int ig = 0; ig < npw; ++ig)
        {
            rhog_in_flat[is * npw + ig] = chr->delta_rhog_save[is][ig];
            rhog_out_flat[is * npw + ig] = chr->delta_rhog[is][ig];
        }
    }

    // Define inner product using complex conjugate: conj(A) * B
    auto inner_product = std::bind(&Charge_Mixing::inner_product_recip_complex,
                                   this,
                                   std::placeholders::_1,
                                   std::placeholders::_2);

    // Define Kerker screening for complex data
    auto screen = [this, nspin](std::complex<double>* rhog) {
        this->Kerker_screen_recip_complex(rhog, nspin);
    };

    // Push data to mixing object
    this->mixing->push_data(this->delta_rho_mdata,
                            rhog_in_flat.data(),
                            rhog_out_flat.data(),
                            screen,
                            true);

    // Calculate mixing coefficients
    this->mixing->cal_coef(this->delta_rho_mdata, inner_product);

    // Perform mixing
    this->mixing->mix_data(this->delta_rho_mdata, rhog_out_flat.data());

    // Copy back to 2D arrays
    for (int is = 0; is < nspin; ++is)
    {
        for (int ig = 0; ig < npw; ++ig)
        {
            chr->delta_rhog[is][ig] = rhog_out_flat[is * npw + ig];
        }
    }

    ModuleBase::timer::tick("Charge_Mixing", "mix_delta_rho_recip");
}

void Charge_Mixing::mix_delta_rho(Charge* chr, const int nspin)
{
    ModuleBase::TITLE("Charge_Mixing", "mix_delta_rho");
    ModuleBase::timer::tick("Charge_Mixing", "mix_delta_rho");

    const int nrxx = this->rhopw->nrxx;
    const int npw = this->rhopw->npw;

    // Save current delta_rho before mixing
    std::vector<std::complex<double>> delta_rho_tmp(nspin * nrxx);
    for (int is = 0; is < nspin; ++is)
    {
        for (int ir = 0; ir < nrxx; ++ir)
        {
            delta_rho_tmp[is * nrxx + ir] = chr->delta_rho[is][ir];
        }
    }

    // FFT: real space -> reciprocal space
    for (int is = 0; is < nspin; ++is)
    {
        this->rhopw->real2recip(chr->delta_rho_save[is], chr->delta_rhog_save[is]);
        this->rhopw->real2recip(chr->delta_rho[is], chr->delta_rhog[is]);
    }

    // Perform mixing in reciprocal space
    mix_delta_rho_recip(chr, nspin);

    // FFT: reciprocal space -> real space
    for (int is = 0; is < nspin; ++is)
    {
        this->rhopw->recip2real(chr->delta_rhog[is], chr->delta_rho[is]);
    }

    // Update delta_rho_save with the pre-mixing values
    for (int is = 0; is < nspin; ++is)
    {
        for (int ir = 0; ir < nrxx; ++ir)
        {
            chr->delta_rho_save[is][ir] = delta_rho_tmp[is * nrxx + ir];
        }
    }

    ModuleBase::timer::tick("Charge_Mixing", "mix_delta_rho");
}
