#include "H_TDDFT_pw.h"

#include "../../input.h"
#include "module_base/constants.h"
#include "module_base/timer.h"
#include "src_lcao/ELEC_evolve.h"

namespace elecstate
{

int H_TDDFT_pw::istep = -1;
//==========================================================
// this function aims to add external time-dependent potential
// (eg: linear potential) used in tddft
// fuxiang add in 2017-05
//==========================================================

void H_TDDFT_pw::cal_fixed_v(double* vl_pseudo)
{
    ModuleBase::TITLE("H_TDDFT_pw", "cal_fixed_v");

    // time evolve
    H_TDDFT_pw::istep++;

    read_parameters();

    // judgement to skip vext
    if (ELEC_evolve::td_vext == 0 || istep > tend || istep < tstart)
    {
        return;
    }
    std::cout << "calculate electric potential" << endl;

    ModuleBase::timer::tick("H_TDDFT_pw", "cal_fixed_v");

    double* vext_space = new double[this->rho_basis_->nrxx];
    double vext_time = 0.0;

    cal_v_space(vext_space);
    vext_time = cal_v_time();

    for (int ir = 0; ir < this->rho_basis_->nrxx; ++ir)
    {
        vl_pseudo[ir] += vext_space[ir] * vext_time;
    }

    delete[] vext_space;

    ModuleBase::timer::tick("H_TDDFT_pw", "cal_fixed_v");
    return;
}

void H_TDDFT_pw::read_parameters()
{
    stype = INPUT.td_stype;

    ttype = INPUT.td_ttype;

    tstart = INPUT.td_tstart;
    tend = INPUT.td_tend;

    // space domain parameters

    // length gauge
    lcut1 = INPUT.td_lcut1;
    lcut2 = INPUT.td_lcut2;

    // time domain parameters

    // Gauss
    gauss_omega = 2 * ModuleBase::PI * ModuleBase::AU_to_FS * INPUT.td_gauss_freq; // time(a.u.)^-1
    gauss_phase = INPUT.td_gauss_phase;
    gauss_sigma2 = INPUT.td_gauss_sigma * INPUT.td_gauss_sigma / (ModuleBase::AU_to_FS * ModuleBase::AU_to_FS);
    gauss_t0 = INPUT.td_gauss_t0;
    gauss_amp = ModuleBase::BOHR_TO_A / ModuleBase::Ry_to_eV * INPUT.td_gauss_amp; // Ry/bohr

    // trapezoid
    trape_omega = 2 * ModuleBase::PI * ModuleBase::AU_to_FS * INPUT.td_trape_freq; // time(a.u.)^-1
    trape_phase = INPUT.td_trape_phase;
    trape_t1 = INPUT.td_trape_t1;
    trape_t2 = INPUT.td_trape_t2;
    trape_t3 = INPUT.td_trape_t3;
    trape_amp = ModuleBase::BOHR_TO_A / ModuleBase::Ry_to_eV * INPUT.td_trape_amp; // Ry/bohr

    // Trigonometric
    trigo_omega1 = 2 * ModuleBase::PI * ModuleBase::AU_to_FS * INPUT.td_trigo_freq1; // time(a.u.)^-1
    trigo_omega2 = 2 * ModuleBase::PI * ModuleBase::AU_to_FS * INPUT.td_trigo_freq1; // time(a.u.)^-1
    trigo_phase1 = INPUT.td_trigo_phase1;
    trigo_phase2 = INPUT.td_trigo_phase2;
    trigo_amp = ModuleBase::BOHR_TO_A / ModuleBase::Ry_to_eV * INPUT.td_trigo_amp; // Ry/bohr

    // Heaviside
    heavi_t0 = INPUT.td_heavi_t0;
    heavi_amp = ModuleBase::BOHR_TO_A / ModuleBase::Ry_to_eV * INPUT.td_heavi_amp; // Ry/bohr

    // HHG
    hhg_amp1 = ModuleBase::BOHR_TO_A / ModuleBase::Ry_to_eV * INPUT.td_hhg_amp1; // Ry/bohr
    hhg_amp2 = ModuleBase::BOHR_TO_A / ModuleBase::Ry_to_eV * INPUT.td_hhg_amp2; // Ry/bohr
    hhg_omega1 = 2 * ModuleBase::PI * ModuleBase::AU_to_FS * INPUT.td_hhg_freq1; // time(a.u.)^-1
    hhg_omega2 = 2 * ModuleBase::PI * ModuleBase::AU_to_FS * INPUT.td_hhg_freq1; // time(a.u.)^-1
    hhg_phase1 = INPUT.td_hhg_phase1;
    hhg_phase2 = INPUT.td_hhg_phase2;
    hhg_t0 = INPUT.td_hhg_t0;
    hhg_sigma2 = INPUT.td_hhg_sigma * INPUT.td_hhg_sigma / (ModuleBase::AU_to_FS * ModuleBase::AU_to_FS);

    return;
}

void H_TDDFT_pw::cal_v_space(double* vext_space)
{
    ModuleBase::TITLE("H_TDDFT_pw", "cal_v_space");
    ModuleBase::timer::tick("H_TDDFT_pw", "cal_v_space");

    switch (stype)
    {
    case 0:
        cal_v_space_length(vext_space);
        break;

    case 1:
        cal_v_space_velocity(vext_space);
        break;

    default:
        std::cout << "space_domain_type of electric field is wrong" << endl;
        break;
    }

    ModuleBase::timer::tick("H_TDDFT_pw", "cal_v_space");
    return;
}

void H_TDDFT_pw::cal_v_space_length(double* vext_space)
{
    ModuleBase::TITLE("H_TDDFT_pw", "cal_v_space_length");
    ModuleBase::timer::tick("H_TDDFT_pw", "cal_v_space_length");

    for (int ir = 0; ir < this->rho_basis_->nrxx; ++ir)
    {
        int i = ir / (this->rho_basis_->ny * this->rho_basis_->nplane);
        int j = ir / this->rho_basis_->nplane - i * this->rho_basis_->ny;
        int k = ir % this->rho_basis_->nplane + this->rho_basis_->startz_current;
        double x = (double)i / this->rho_basis_->nx;
        double y = (double)j / this->rho_basis_->ny;
        double z = (double)k / this->rho_basis_->nz;

        switch (ELEC_evolve::td_vext_dire)
        {
        case 1:
            vext_space[ir] = cal_v_space_length_potential(x);
            break;

        case 2:
            vext_space[ir] = cal_v_space_length_potential(y);
            break;

        case 3:
            vext_space[ir] = cal_v_space_length_potential(z);
            break;

        default:
            std::cout << "direction of electric field is wrong" << endl;
            break;
        }
    }

    ModuleBase::timer::tick("H_TDDFT_pw", "cal_v_space_length");
    return;
}

double H_TDDFT_pw::cal_v_space_length_potential(double i)
{
    double vext_space;
    if (i < this->rho_basis_->nx * lcut1)
    {
        vext_space = ((i / this->rho_basis_->nx - lcut1) / (lcut1 + 1.0 - lcut2) - lcut1) * this->ucell_->lat0;
    }
    else if (i >= this->rho_basis_->nx * lcut1 && i < this->rho_basis_->nx * lcut2)
    {
        vext_space = -i / this->rho_basis_->nx * this->ucell_->lat0;
    }
    else if (i >= this->rho_basis_->nx * lcut2)
    {
        vext_space = ((i / this->rho_basis_->nx - lcut2) / (lcut1 + 1.0 - lcut2) - lcut2) * this->ucell_->lat0;
    }
    return vext_space;
}

void H_TDDFT_pw::cal_v_space_velocity(double* vext_space)
{
    return;
}

double H_TDDFT_pw::cal_v_time()
{
    double vext_time = 0.0;

    switch (ttype)
    {
    case 0:
        vext_time = cal_v_time_Gauss();
        cout << "gauss vtime=" << vext_time << endl;
        break;

    case 1:
        vext_time = cal_v_time_trapezoid();
        cout << "trape vtime=" << vext_time << endl;
        break;

    case 2:
        vext_time = cal_v_time_trigonometric();
        cout << "trigo vtime=" << vext_time << endl;
        break;

    case 3:
        vext_time = cal_v_time_heaviside();
        cout << "heavi vtime=" << vext_time << endl;
        break;

    case 4:
        vext_time = cal_v_time_HHG();
        cout << "hhg vtime=" << vext_time << endl;
        break;

    default:
        std::cout << "time_domain_type of electric field is wrong" << endl;
        break;
    }
    return vext_time;
}

double H_TDDFT_pw::cal_v_time_Gauss()
{
    double vext_time = 0.0;

    double gauss_t = (istep - gauss_t0) * INPUT.mdp.md_dt;
    vext_time = cos(gauss_omega * gauss_t + gauss_phase) * exp(-gauss_t * gauss_t * 0.5 / (gauss_sigma2)) * gauss_amp;

    return vext_time;
}

double H_TDDFT_pw::cal_v_time_trapezoid()
{
    double vext_time = 0.0;

    if (istep < trape_t1)
    {
        vext_time = istep / trape_t1;
    }
    else if (istep < trape_t2)
    {
        vext_time = 1.0;
    }
    else if (istep < trape_t3)
    {
        vext_time = (trape_t3 - istep) / (trape_t3 - trape_t2);
    }

    vext_time = vext_time * trape_amp * cos(trape_omega * istep * INPUT.mdp.md_dt + trape_phase);

    return vext_time;
}

double H_TDDFT_pw::cal_v_time_trigonometric()
{
    double vext_time = 0.0;

    const double timenow = istep * INPUT.mdp.md_dt;

    vext_time = trigo_amp * cos(trigo_omega1 * timenow + trigo_phase1) * sin(trigo_omega2 * timenow + trigo_phase2)
                * sin(trigo_omega2 * timenow + trigo_phase2);

    return vext_time;
}

double H_TDDFT_pw::cal_v_time_heaviside()
{
    if (istep < heavi_t0)
        return heavi_amp;
    else if (istep >= heavi_t0)
        return 0.0;
}

double H_TDDFT_pw::cal_v_time_HHG()
{
    double vext_time = 0.0;

    double hhg_t = (istep - hhg_t0) * INPUT.mdp.md_dt;
    vext_time = (cos(hhg_omega1 * hhg_t + hhg_phase1) * hhg_amp1 + cos(hhg_omega2 * hhg_t + hhg_phase2) * hhg_amp2)
                * exp(-hhg_t * hhg_t * 0.5 / (hhg_sigma2));

    return vext_time;
}

} // namespace elecstate
