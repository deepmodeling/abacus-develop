#include "H_TDDFT_pw.h"

#include "module_base/constants.h"
#include "module_base/timer.h"
#include "module_hamilt_lcao/module_tddft/evolve_elec.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/input.h"
#include "module_io/input_conv.h"

namespace elecstate
{

int H_TDDFT_pw::istep = -1;

double H_TDDFT_pw::amp;
double H_TDDFT_pw::bmod;
double H_TDDFT_pw::bvec[3];

int H_TDDFT_pw::stype; // 0 : length gauge  1: velocity gauge

std::vector<int> H_TDDFT_pw::ttype;
//  0  Gauss type function.
//  1  trapezoid type function.
//  2  Trigonometric functions, sin^2.
//  3  heaviside function.
//  4  HHG function.

int H_TDDFT_pw::tstart;
int H_TDDFT_pw::tend;
double H_TDDFT_pw::dt;

// space domain parameters

// length gauge
double H_TDDFT_pw::lcut1;
double H_TDDFT_pw::lcut2;

// time domain parameters

// Gauss
int H_TDDFT_pw::gauss_count;
std::vector<double> H_TDDFT_pw::gauss_omega; // time(a.u.)^-1
std::vector<double> H_TDDFT_pw::gauss_phase;
std::vector<double> H_TDDFT_pw::gauss_sigma; // time(a.u.)
std::vector<double> H_TDDFT_pw::gauss_t0;
std::vector<double> H_TDDFT_pw::gauss_amp; // Ry/bohr

// trapezoid
int H_TDDFT_pw::trape_count;
std::vector<double> H_TDDFT_pw::trape_omega; // time(a.u.)^-1
std::vector<double> H_TDDFT_pw::trape_phase;
std::vector<double> H_TDDFT_pw::trape_t1;
std::vector<double> H_TDDFT_pw::trape_t2;
std::vector<double> H_TDDFT_pw::trape_t3;
std::vector<double> H_TDDFT_pw::trape_amp; // Ry/bohr

// Trigonometric
int H_TDDFT_pw::trigo_count;
std::vector<double> H_TDDFT_pw::trigo_omega1; // time(a.u.)^-1
std::vector<double> H_TDDFT_pw::trigo_omega2; // time(a.u.)^-1
std::vector<double> H_TDDFT_pw::trigo_phase1;
std::vector<double> H_TDDFT_pw::trigo_phase2;
std::vector<double> H_TDDFT_pw::trigo_amp; // Ry/bohr

// Heaviside
int H_TDDFT_pw::heavi_count;
std::vector<double> H_TDDFT_pw::heavi_t0;
std::vector<double> H_TDDFT_pw::heavi_amp; // Ry/bohr

void H_TDDFT_pw::cal_fixed_v(double *vl_pseudo)
{
    ModuleBase::TITLE("H_TDDFT_pw", "cal_fixed_v");

    // time evolve
    H_TDDFT_pw::istep++;

    // judgement to skip vext
    if (!module_tddft::Evolve_elec::td_vext || istep > tend || istep < tstart)
    {
        return;
    }
    std::cout << "calculate electric potential" << endl;

    ModuleBase::timer::tick("H_TDDFT_pw", "cal_fixed_v");

    int count = 0;
    gauss_count = 0;
    trape_count = 0;
    trigo_count = 0;
    heavi_count = 0;

    for (auto direc: module_tddft::Evolve_elec::td_vext_dire_case)
    {
        std::vector<double> vext_space(this->rho_basis_->nrxx, 0.0);
        double vext_time = cal_v_time(ttype[count]);

        if (module_tddft::Evolve_elec::out_efield && GlobalV::MY_RANK == 0)
        {
            std::stringstream as;
            as << GlobalV::global_out_dir << "efield_" << count << ".dat";
            std::ofstream ofs(as.str().c_str(), std::ofstream::app);
            ofs << H_TDDFT_pw::istep * dt * ModuleBase::AU_to_FS << "\t"
                << vext_time * ModuleBase::Ry_to_eV / ModuleBase::BOHR_TO_A << endl;
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

    prepare(GlobalC::ucell, direc);

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
            vext_space[ir] = cal_v_space_length_potential(x) / bmod;
            break;

        case 2:
            vext_space[ir] = cal_v_space_length_potential(y) / bmod;
            break;

        case 3:
            vext_space[ir] = cal_v_space_length_potential(z) / bmod;
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
    double vext_space = 0.0;
    if (i < lcut1)
    {
        vext_space = ((i - lcut1) * (lcut2 - lcut1) / (lcut1 + 1.0 - lcut2) - lcut1) * this->ucell_->lat0;
    }
    else if (i >= lcut1 && i < lcut2)
    {
        vext_space = -i * this->ucell_->lat0;
    }
    else if (i >= lcut2)
    {
        vext_space = ((i - lcut2) * (lcut2 - lcut1) / (lcut1 + 1.0 - lcut2) - lcut2) * this->ucell_->lat0;
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

    case 1:
        vext_time = cal_v_time_trapezoid();
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

    vext_time = amp * cos(omega1 * timenow + phase1) * sin(omega2 * timenow + phase2) * sin(omega2 * timenow + phase2);
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

void H_TDDFT_pw::prepare(const UnitCell& cell, int& dir)
{
    if (dir == 1)
    {
        bvec[0] = cell.G.e11;
        bvec[1] = cell.G.e12;
        bvec[2] = cell.G.e13;
    }
    else if (dir == 2)
    {
        bvec[0] = cell.G.e21;
        bvec[1] = cell.G.e22;
        bvec[2] = cell.G.e23;
    }
    else if (dir == 3)
    {
        bvec[0] = cell.G.e31;
        bvec[1] = cell.G.e32;
        bvec[2] = cell.G.e33;
    }
    else
    {
        ModuleBase::WARNING_QUIT("H_TDDFT_pw::prepare", "direction is wrong!");
    }
    bmod = sqrt(pow(bvec[0], 2) + pow(bvec[1], 2) + pow(bvec[2], 2));
}

void H_TDDFT_pw ::compute_force(const UnitCell& cell, ModuleBase::matrix& fe)
{
    int iat = 0;
    for (int it = 0; it < cell.ntype; ++it)
    {
        for (int ia = 0; ia < cell.atoms[it].na; ++ia)
        {
            for (int jj = 0; jj < 3; ++jj)
            {
                fe(iat, jj) = ModuleBase::e2 * amp * cell.atoms[it].ncpp.zv * bvec[jj] / bmod;
            }
            ++iat;
        }
    }
}

} // namespace elecstate
