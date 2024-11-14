#include "charge_mixing.h"

#include "module_parameter/parameter.h"
#include "module_base/timer.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_base/parallel_reduce.h"

double Charge_Mixing::get_drho(Charge* chr, const double nelec)
{
    ModuleBase::TITLE("Charge_Mixing", "get_drho");
    ModuleBase::timer::tick("Charge_Mixing", "get_drho");
    double drho = 0.0;

    if (PARAM.inp.scf_thr_type == 1)
    {
        for (int is = 0; is < PARAM.inp.nspin; ++is)
        {
            ModuleBase::GlobalFunc::NOTE("Perform FFT on rho(r) to obtain rho(G).");
            chr->rhopw->real2recip(chr->rho[is], chr->rhog[is]);

            ModuleBase::GlobalFunc::NOTE("Perform FFT on rho_save(r) to obtain rho_save(G).");
            chr->rhopw->real2recip(chr->rho_save[is], chr->rhog_save[is]);
        }

        ModuleBase::GlobalFunc::NOTE("Calculate the charge difference between rho(G) and rho_save(G)");
        std::vector<std::complex<double>> drhog(PARAM.inp.nspin * this->rhopw->npw);
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 512)
#endif
        for (int is = 0; is < PARAM.inp.nspin; ++is)
        {
            for (int ig = 0; ig < this->rhopw->npw; ig++)
            {
                drhog[is * rhopw->npw + ig] = chr->rhog[is][ig] - chr->rhog_save[is][ig];
            }
        }

        ModuleBase::GlobalFunc::NOTE("Calculate the norm of the Residual std::vector: < R[rho] | R[rho_save] >");
        drho = this->inner_product_recip_rho(drhog.data(), drhog.data());
    }
    else
    {
        // Note: Maybe it is wrong.
        //       The inner_product_real function (L1-norm) is different from that (L2-norm) in mixing.
        for (int is = 0; is < PARAM.inp.nspin; is++)
        {
            if (is != 0 && is != 3 && PARAM.globalv.domag_z)
            {
                continue;
            }
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : drho)
#endif
            for (int ir = 0; ir < this->rhopw->nrxx; ir++)
            {
                drho += std::abs(chr->rho[is][ir] - chr->rho_save[is][ir]);
            }
        }
#ifdef __MPI
        Parallel_Reduce::reduce_pool(drho);
#endif
        assert(nelec != 0);
        assert(GlobalC::ucell.omega > 0);
        assert(this->rhopw->nxyz > 0);
        drho *= GlobalC::ucell.omega / static_cast<double>(this->rhopw->nxyz);
        drho /= nelec;
    }

    ModuleBase::timer::tick("Charge_Mixing", "get_drho");
    return drho;
}

double Charge_Mixing::get_dkin(Charge* chr, const double nelec)
{
    if (!(XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)) 
    {
        return 0.0;
    };
    ModuleBase::TITLE("Charge_Mixing", "get_dkin");
    ModuleBase::timer::tick("Charge_Mixing", "get_dkin");
    double dkin = 0.0;
    
    // Get dkin from kin_r and kin_r_save for PW and LCAO both, which is different from drho.
    for (int is = 0; is < PARAM.inp.nspin; is++)
    {
        if (is != 0 && is != 3 && PARAM.globalv.domag_z)
        {
            continue;
        }
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : dkin)
#endif
        for (int ir = 0; ir < this->rhopw->nrxx; ir++)
        {
            dkin += std::abs(chr->kin_r[is][ir] - chr->kin_r_save[is][ir]);
        }
    }
#ifdef __MPI
    Parallel_Reduce::reduce_pool(dkin);
#endif
    assert(nelec != 0);
    assert(GlobalC::ucell.omega > 0);
    assert(this->rhopw->nxyz > 0);
    dkin *= GlobalC::ucell.omega / static_cast<double>(this->rhopw->nxyz);
    dkin /= nelec;

    ModuleBase::timer::tick("Charge_Mixing", "get_dkin");
    return dkin;
}

double Charge_Mixing::inner_product_recip_rho(std::complex<double>* rho1, std::complex<double>* rho2)
{
    ModuleBase::TITLE("Charge_Mixing", "inner_product_recip_rho");
    ModuleBase::timer::tick("Charge_Mixing", "inner_product_recip_rho");

    std::complex<double>** rhog1 = new std::complex<double>*[PARAM.inp.nspin];
    std::complex<double>** rhog2 = new std::complex<double>*[PARAM.inp.nspin];
    for (int is = 0; is < PARAM.inp.nspin; is++)
    {
        rhog1[is] = rho1 + is * this->rhopw->npw;
        rhog2[is] = rho2 + is * this->rhopw->npw;
    }

    static const double fac = ModuleBase::e2 * ModuleBase::FOUR_PI / GlobalC::ucell.tpiba2;
    static const double fac2 = ModuleBase::e2 * ModuleBase::FOUR_PI / (ModuleBase::TWO_PI * ModuleBase::TWO_PI);

    double sum = 0.0;

    auto part_of_noncolin = [&]()
    {
        double sum = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : sum)
#endif
        for (int ig = 0; ig < this->rhopw->npw; ++ig)
        {
            if (this->rhopw->gg[ig] < 1e-8) {
                continue;
}
            sum += (conj(rhog1[0][ig]) * rhog2[0][ig]).real() / this->rhopw->gg[ig];
        }
        sum *= fac;
        return sum;
    };

    switch (PARAM.inp.nspin)
    {
    case 1:
        sum += part_of_noncolin();
        break;

    case 2: {
        // (1) First part of density error.
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : sum)
#endif
        for (int ig = 0; ig < this->rhopw->npw; ++ig)
        {
            if (this->rhopw->gg[ig] < 1e-8) {
                continue;
}
            sum += (conj(rhog1[0][ig] + rhog1[1][ig]) * (rhog2[0][ig] + rhog2[1][ig])).real() / this->rhopw->gg[ig];
        }
        sum *= fac;

        if (PARAM.globalv.gamma_only_pw)
        {
            sum *= 2.0;
        }

        // (2) Second part of density error.
        // including |G|=0 term.
        double sum2 = 0.0;

        sum2 += fac2 * (conj(rhog1[0][0] - rhog1[1][0]) * (rhog2[0][0] - rhog2[1][0])).real();

        double mag = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : mag)
#endif
        for (int ig = 0; ig < this->rhopw->npw; ig++)
        {
            mag += (conj(rhog1[0][ig] - rhog1[1][ig]) * (rhog2[0][ig] - rhog2[1][ig])).real();
        }
        mag *= fac2;

        // if(PARAM.globalv.gamma_only_pw);
        if (PARAM.globalv.gamma_only_pw) // Peize Lin delete ; 2020.01.31
        {
            mag *= 2.0;
        }

        // std::cout << " sum=" << sum << " mag=" << mag << std::endl;
        sum2 += mag;
        sum += sum2;
        break;
    }
    case 4:
        // non-collinear spin, added by zhengdy
        if (!PARAM.globalv.domag && !PARAM.globalv.domag_z) {
            sum += part_of_noncolin();
        } else
        {
            // another part with magnetization
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : sum)
#endif
            for (int ig = 0; ig < this->rhopw->npw; ig++)
            {
                if (ig == this->rhopw->ig_gge0) {
                    continue;
}
                sum += (conj(rhog1[0][ig]) * rhog2[0][ig]).real() / this->rhopw->gg[ig];
            }
            sum *= fac;
            const int ig0 = this->rhopw->ig_gge0;
            if (ig0 > 0)
            {
                sum += fac2
                       * ((conj(rhog1[1][ig0]) * rhog2[1][ig0]).real() + (conj(rhog1[2][ig0]) * rhog2[2][ig0]).real()
                          + (conj(rhog1[3][ig0]) * rhog2[3][ig0]).real());
            }
            double fac3 = fac2;
            if (PARAM.globalv.gamma_only_pw)
            {
                fac3 *= 2.0;
            }
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : sum)
#endif
            for (int ig = 0; ig < this->rhopw->npw; ig++)
            {
                if (ig == ig0) {
                    continue;
}
                sum += fac3
                       * ((conj(rhog1[1][ig]) * rhog2[1][ig]).real() + (conj(rhog1[2][ig]) * rhog2[2][ig]).real()
                          + (conj(rhog1[3][ig]) * rhog2[3][ig]).real());
            }
        }
        break;
    }
#ifdef __MPI
    Parallel_Reduce::reduce_pool(sum);
#endif
    ModuleBase::timer::tick("Charge_Mixing", "inner_product_recip_rho");

    sum *= GlobalC::ucell.omega * 0.5;

    delete[] rhog1;
    delete[] rhog2;
    return sum;
}

// a simple inner product, now is not used anywhere. For test only.
double Charge_Mixing::inner_product_recip_simple(std::complex<double>* rho1, std::complex<double>* rho2)
{
    ModuleBase::TITLE("Charge_Mixing", "inner_product_recip_simple");
    ModuleBase::timer::tick("Charge_Mixing", "inner_product_recip_simple");

    double rnorm = 0.0;
    // consider a resize for mixing_angle
    int resize_tmp = 1;
    if (PARAM.inp.nspin == 4 && this->mixing_angle > 0) { resize_tmp = 2;
}
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : rnorm)
#endif
    for (int ig = 0; ig < this->rhopw->npw * PARAM.inp.nspin / resize_tmp; ++ig)
    {
        rnorm += (conj(rho1[ig]) * rho2[ig]).real();
    }
#ifdef __MPI
    Parallel_Reduce::reduce_pool(rnorm);
#endif

    ModuleBase::timer::tick("Charge_Mixing", "inner_product_recip_simple");

    return rnorm;
}

// a Hartree-like inner product
double Charge_Mixing::inner_product_recip_hartree(std::complex<double>* rhog1, std::complex<double>* rhog2)
{
    ModuleBase::TITLE("Charge_Mixing", "inner_product_recip_hartree");
    ModuleBase::timer::tick("Charge_Mixing", "inner_product_recip_hartree");

    static const double fac = ModuleBase::e2 * ModuleBase::FOUR_PI / GlobalC::ucell.tpiba2;
    static const double fac2 = ModuleBase::e2 * ModuleBase::FOUR_PI / (ModuleBase::TWO_PI * ModuleBase::TWO_PI);

    double sum = 0.0;
    const int npw = this->rhopw->npw;

    // a lambda function for summing the charge density
    auto part_of_rho = [&]()
    {
        double sum = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : sum)
#endif
        for (int ig = 0; ig < this->rhopw->npw; ++ig)
        {
            if (this->rhopw->gg[ig] < 1e-8) {
                continue;
}
            sum += (conj(rhog1[ig]) * rhog2[ig]).real() / this->rhopw->gg[ig];
        }
        sum *= fac;
        return sum;
    };
    
    if (PARAM.inp.nspin==1)
    {
        sum += part_of_rho();
    }
    else if (PARAM.inp.nspin==2)
    {
        // charge density part
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : sum)
#endif
        for (int ig = 0; ig < this->rhopw->npw; ++ig)
        {
            if (this->rhopw->gg[ig] < 1e-8) {
                continue;
}
            sum += (conj(rhog1[ig]) * (rhog2[ig])).real() / this->rhopw->gg[ig];
        }
        sum *= fac;

        if (PARAM.globalv.gamma_only_pw)
        {
            sum *= 2.0;
        }

        // (2) Second part of density error.
        // including |G|=0 term.
        double sum2 = 0.0;

        sum2 += fac2 * (conj(rhog1[0 + this->rhopw->npw]) * rhog2[0 + this->rhopw->npw]).real();

        double mag = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : mag)
#endif
        for (int ig = 0; ig < this->rhopw->npw; ig++)
        {
            mag += (conj(rhog1[ig + this->rhopw->npw]) * rhog2[ig + this->rhopw->npw]).real();
        }
        mag *= fac2;

        if (PARAM.globalv.gamma_only_pw)
        {
            mag *= 2.0;
        }

        sum2 += mag;
        sum += sum2;
    }
    else if (PARAM.inp.nspin==4)
    {
        if (!PARAM.globalv.domag && !PARAM.globalv.domag_z)
        {
            sum += part_of_rho();
        }
        else if (this->mixing_angle <= 0)
        {
            // sum for tradtional mixing
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : sum)
#endif
            for (int ig = 0; ig < this->rhopw->npw; ig++)
            {
                if (ig == this->rhopw->ig_gge0) {
                    continue;
}
                sum += (conj(rhog1[ig]) * rhog2[ig]).real() / this->rhopw->gg[ig];
            }
            sum *= fac;
            const int ig0 = this->rhopw->ig_gge0;
            if (ig0 > 0)
            {
                sum += fac2
                       * ((conj(rhog1[ig0 + npw]) * rhog2[ig0 + npw]).real() + (conj(rhog1[ig0 + 2*npw]) * rhog2[ig0 + 2*npw]).real()
                          + (conj(rhog1[ig0 + 3*npw]) * rhog2[ig0 + 3*npw]).real());
            }
            double fac3 = fac2;
            if (PARAM.globalv.gamma_only_pw)
            {
                fac3 *= 2.0;
            }
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : sum)
#endif
            for (int ig = 0; ig < this->rhopw->npw; ig++)
            {
                if (ig == ig0) {
                    continue;
}
                sum += fac3
                       * ((conj(rhog1[ig + npw]) * rhog2[ig + npw]).real() + (conj(rhog1[ig + 2*npw]) * rhog2[ig + 2*npw]).real()
                          + (conj(rhog1[ig + 3*npw]) * rhog2[ig + 3*npw]).real());
            }
        }
        else if (this->mixing_angle > 0)
        {
            // sum for angle mixing
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : sum)
#endif
            for (int ig = 0; ig < this->rhopw->npw; ig++)
            {
                if (ig == this->rhopw->ig_gge0) {
                    continue;
}
                sum += (conj(rhog1[ig]) * rhog2[ig]).real() / this->rhopw->gg[ig];
            }
            sum *= fac;
            const int ig0 = this->rhopw->ig_gge0;
            if (ig0 > 0)
            {
                sum += fac2
                       * ((conj(rhog1[ig0 + this->rhopw->npw]) * rhog2[ig0 + this->rhopw->npw]).real());
            }
            double fac3 = fac2;
            if (PARAM.globalv.gamma_only_pw)
            {
                fac3 *= 2.0;
            }
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : sum)
#endif
            for (int ig = 0; ig < this->rhopw->npw; ig++)
            {
                if (ig == ig0) {
                    continue;
}
                sum += fac3
                       * ((conj(rhog1[ig + this->rhopw->npw]) * rhog2[ig + this->rhopw->npw]).real());
            }
        }
    }
#ifdef __MPI
    Parallel_Reduce::reduce_pool(sum);
#endif

    ModuleBase::timer::tick("Charge_Mixing", "inner_product_recip_hartree");

    sum *= GlobalC::ucell.omega * 0.5;

    return sum;
}

double Charge_Mixing::inner_product_real(double* rho1, double* rho2)
{
    double rnorm = 0.0;
    // consider a resize for mixing_angle
    int resize_tmp = 1;
    if (PARAM.inp.nspin == 4 && this->mixing_angle > 0) { resize_tmp = 2;
}

#ifdef _OPENMP
#pragma omp parallel for reduction(+ : rnorm)
#endif
    for (int ir = 0; ir < this->rhopw->nrxx * PARAM.inp.nspin / resize_tmp; ++ir)
    {
        rnorm += rho1[ir] * rho2[ir];
    }
#ifdef __MPI
    Parallel_Reduce::reduce_pool(rnorm);
#endif
    return rnorm;
}