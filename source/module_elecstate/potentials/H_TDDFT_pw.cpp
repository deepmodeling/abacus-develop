#include "H_TDDFT_pw.h"

#include "module_io/input.h"
#include "module_base/constants.h"
#include "module_base/timer.h"
#include "module_hamilt_lcao/module_tddft/ELEC_evolve.h"
#include "module_io/input_conv.h"

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

    read_parameters(&INPUT);

    // judgement to skip vext
    if (!ELEC_evolve::td_vext || istep > tend || istep < tstart)
    {
        return;
    }
    std::cout << "calculate electric potential" << endl;

    ModuleBase::timer::tick("H_TDDFT_pw", "cal_fixed_v");

    int count = 0;

    for (auto direc: ELEC_evolve::td_vext_dire_case)
    {
        std::vector<double> vext_space(this->rho_basis_->nrxx, 0.0);
        double vext_time = cal_v_time(ttype[count]); 

        if (ELEC_evolve::out_efield && GlobalV::MY_RANK == 0)
        {
            std::stringstream as;
            as << GlobalV::global_out_dir << "efield_"<<count<<".dat";
            std::ofstream ofs(as.str().c_str(), std::ofstream::app);
            ofs << H_TDDFT_pw::istep*dt*ModuleBase::AU_to_FS << "\t" << vext_time <<endl;
            ofs.close();
        }

        cal_v_space(vext_space, direc); 
        for (size_t ir = 0; ir < this->rho_basis_->nrxx; ++ir)
            vl_pseudo[ir] += vext_space[ir] * vext_time;
        count++;
    }

    ModuleBase::timer::tick("H_TDDFT_pw", "cal_fixed_v");
    return;
}

std::vector<double> H_TDDFT_pw::set_parameters(std::string params, double c)
{
    std::vector<double> params_ori;
    std::vector<double> params_out;
    Input_Conv::parse_expression(params, params_ori);
    for (auto param: params_ori)
        params_out.emplace_back(param * c);

    return params_out;
}

void H_TDDFT_pw::read_parameters(Input *in)
{
    stype = in->td_stype;

    Input_Conv::parse_expression(in->td_ttype, ttype);

    tstart = in->td_tstart;
    tend = in->td_tend;

    dt = in->mdp.md_dt;

    // space domain parameters

    // length gauge
    lcut1 = in->td_lcut1;
    lcut2 = in->td_lcut2;

    // time domain parameters

    // Gauss
    gauss_count = 0;
    gauss_omega = set_parameters(in->td_gauss_freq, 2 * ModuleBase::PI * ModuleBase::AU_to_FS); // time(a.u.)^-1
    gauss_phase = set_parameters(in->td_gauss_phase, 1.0);
    gauss_sigma = set_parameters(in->td_gauss_sigma, 1/ModuleBase::AU_to_FS);
    gauss_t0 = set_parameters(in->td_gauss_t0, 1.0);
    gauss_amp = set_parameters(in->td_gauss_amp, ModuleBase::BOHR_TO_A / ModuleBase::Ry_to_eV); // Ry/bohr

    // trapezoid
    trape_count = 0;
    trape_omega = set_parameters(in->td_trape_freq, 2 * ModuleBase::PI * ModuleBase::AU_to_FS); // time(a.u.)^-1
    trape_phase = set_parameters(in->td_trape_phase, 1.0);
    trape_t1 = set_parameters(in->td_trape_t1, 1.0);
    trape_t2 = set_parameters(in->td_trape_t2, 1.0);
    trape_t3 = set_parameters(in->td_trape_t3, 1.0);
    trape_amp = set_parameters(in->td_trape_amp, ModuleBase::BOHR_TO_A / ModuleBase::Ry_to_eV); // Ry/bohr

    // Trigonometric
    trigo_count = 0;
    trigo_omega1 = set_parameters(in->td_trigo_freq1, 2 * ModuleBase::PI * ModuleBase::AU_to_FS); // time(a.u.)^-1
    trigo_omega2 = set_parameters(in->td_trigo_freq2, 2 * ModuleBase::PI * ModuleBase::AU_to_FS); // time(a.u.)^-1
    trigo_phase1 = set_parameters(in->td_trigo_phase1, 1.0);
    trigo_phase2 = set_parameters(in->td_trigo_phase2, 1.0);
    trigo_amp = set_parameters(in->td_trigo_amp, ModuleBase::BOHR_TO_A / ModuleBase::Ry_to_eV); // Ry/bohr

    // Heaviside
    heavi_count = 0;
    heavi_t0 = set_parameters(in->td_heavi_t0, 1.0);
    heavi_amp = set_parameters(in->td_heavi_amp, ModuleBase::BOHR_TO_A / ModuleBase::Ry_to_eV); // Ry/bohr

    // HHG
    // hhg_count = 0;
    // hhg_amp1 = set_parameters(in->td_hhg_amp1, ModuleBase::BOHR_TO_A / ModuleBase::Ry_to_eV); // Ry/bohr
    // hhg_amp2 = set_parameters(in->td_hhg_amp2, ModuleBase::BOHR_TO_A / ModuleBase::Ry_to_eV); // Ry/bohr
    // hhg_omega1 = set_parameters(in->td_hhg_freq1, 2 * ModuleBase::PI * ModuleBase::AU_to_FS); // time(a.u.)^-1
    // hhg_omega2 = set_parameters(in->td_hhg_freq2, 2 * ModuleBase::PI * ModuleBase::AU_to_FS); // time(a.u.)^-1
    // hhg_phase1 = set_parameters(in->td_hhg_phase1, 1.0);
    // hhg_phase2 = set_parameters(in->td_hhg_phase2, 1.0);
    // hhg_t0 = set_parameters(in->td_hhg_t0, 1.0);
    // hhg_sigma = set_parameters(in->td_hhg_sigma, 1/ModuleBase::AU_to_FS);

    return;
}

void H_TDDFT_pw::cal_v_space(std::vector<double> &vext_space, int direc)
{
    ModuleBase::TITLE("H_TDDFT_pw", "cal_v_space");
    ModuleBase::timer::tick("H_TDDFT_pw", "cal_v_space");

    switch (stype)
    {
    case 0:
        cal_v_space_length(vext_space, direc);
        break;

    case 1:
        cal_v_space_velocity(vext_space, direc);
        break;

    default:
        std::cout << "space_domain_type of electric field is wrong" << endl;
        break;
    }

    ModuleBase::timer::tick("H_TDDFT_pw", "cal_v_space");
    return;
}

void H_TDDFT_pw::cal_v_space_length(std::vector<double> &vext_space, int direc)
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

        switch (direc)
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
    double vext_space=0.0;
    if (i < this->rho_basis_->nx * lcut1)
    {
        vext_space = ((i / this->rho_basis_->nx - lcut1)*(lcut2-lcut1) / (lcut1 + 1.0 - lcut2) - lcut1) * this->ucell_->lat0;
    }
    else if (i >= this->rho_basis_->nx * lcut1 && i < this->rho_basis_->nx * lcut2)
    {
        vext_space = -i / this->rho_basis_->nx * this->ucell_->lat0;
    }
    else if (i >= this->rho_basis_->nx * lcut2)
    {
        vext_space = ((i / this->rho_basis_->nx - lcut2)*(lcut2-lcut1) / (lcut1 + 1.0 - lcut2) - lcut2) * this->ucell_->lat0;
    }
    return vext_space;
}

void H_TDDFT_pw::cal_v_space_velocity(std::vector<double> &vext_space, int direc)
{
    return;
}

double H_TDDFT_pw::cal_v_time(int t_type)
{
    double vext_time = 0.0;

    switch (t_type)
    {
    case 0:
        vext_time = cal_v_time_Gauss();
        break;

    case 2:
        vext_time = cal_v_time_trigonometric();
        break;

    case 3:
        vext_time = cal_v_time_heaviside();
        break;

    // case 4:
    //     vext_time = cal_v_time_HHG();
    //     break;

    default:
        std::cout << "time_domain_type of electric field is wrong" << endl;
        break;
    }
    return vext_time;
}

double H_TDDFT_pw::cal_v_time_Gauss()
{
    double vext_time = 0.0;
    double t0 = *(gauss_t0.begin() + gauss_count);
    double omega = *(gauss_omega.begin() + gauss_count);
    double sigma = *(gauss_sigma.begin() + gauss_count);
    double phase = *(gauss_phase.begin() + gauss_count);
    double amp = *(gauss_amp.begin() + gauss_count);

    double gauss_t = (istep - t0) * dt;
    vext_time = cos(omega * gauss_t + phase) * exp(-gauss_t * gauss_t * 0.5 / (sigma * sigma)) * amp;
    gauss_count++;

    return vext_time;
}

double H_TDDFT_pw::cal_v_time_trapezoid()
{
    double vext_time = 0.0;
    double t1 = *(trape_t1.begin() + trape_count);
    double t2 = *(trape_t2.begin() + trape_count);
    double t3 = *(trape_t3.begin() + trape_count);
    double omega = *(trape_omega.begin() + trape_count);
    double phase = *(trape_phase.begin() + trape_count);
    double amp = *(trape_amp.begin() + trape_count);

    if (istep < t1)
    {
        vext_time = istep / t1;
    }
    else if (istep < t2)
    {
        vext_time = 1.0;
    }
    else if (istep < t3)
    {
        vext_time = (t3 - istep) / (t3 - t2);
    }

    vext_time = vext_time * amp * cos(omega * istep * dt + phase);
    trape_count++;

    return vext_time;
}

double H_TDDFT_pw::cal_v_time_trigonometric()
{
    double vext_time = 0.0;
    double omega1 = *(trigo_omega1.begin() + trigo_count);
    double phase1 = *(trigo_phase1.begin() + trigo_count);
    double omega2 = *(trigo_omega2.begin() + trigo_count);
    double phase2 = *(trigo_phase2.begin() + trigo_count);
    double amp = *(trigo_amp.begin() + trigo_count);

    const double timenow = istep * dt;

    vext_time = amp * cos(omega1 * timenow + phase1) * sin(omega2 * timenow + phase2)
                * sin(omega2 * timenow + phase2);
    trigo_count++;

    return vext_time;
}

double H_TDDFT_pw::cal_v_time_heaviside()
{
    double t0 = *(heavi_t0.begin() + heavi_count);
    double amp = *(heavi_amp.begin() + heavi_count);
    double vext_time = 0.0;
    if (istep < t0)
        vext_time = amp;
    else if (istep >= t0)
        vext_time = 0.0;
    heavi_count++;

    return vext_time;
}

// double H_TDDFT_pw::cal_v_time_HHG()
// {
//     double vext_time = 0.0;
//     double t0 = *(hhg_t0.begin() + hhg_count);
//     double omega1 = *(hhg_omega1.begin() + hhg_count);
//     double phase1 = *(hhg_phase1.begin() + hhg_count);
//     double omega2 = *(hhg_omega2.begin() + hhg_count);
//     double phase2 = *(hhg_phase2.begin() + hhg_count);
//     double amp1 = *(hhg_amp1.begin() + hhg_count);
//     double amp2 = *(hhg_amp2.begin() + hhg_count);
//     double sigma = *(trigo_amp2.begin() + trigo_count);

//     double hhg_t = (istep - t0) * dt;
//     vext_time = (cos(omega1 * hhg_t + phase1) * amp1 + cos(omega2 * hhg_t + phase2) * amp2)
//                 * exp(-hhg_t * hhg_t * 0.5 / (sigma * sigma));
//     hhg_count++;

//     return vext_time;
// }

} // namespace elecstate
