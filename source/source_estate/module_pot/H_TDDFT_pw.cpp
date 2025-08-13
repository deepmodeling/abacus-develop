#include "H_TDDFT_pw.h"

#include "source_base/constants.h"
#include "source_base/math_integral.h"
#include "source_base/timer.h"
#include "source_io/input_conv.h"
#include "source_io/module_parameter/parameter.h"
#include "source_lcao/module_rt/evolve_elec.h"
#include "source_pw/module_pwdft/global.h"

namespace elecstate
{

int H_TDDFT_pw::istep = -1;
bool H_TDDFT_pw::is_initialized = false;

double H_TDDFT_pw::amp;

// Used for calculating electric field force on ions, summing over directions
vector<double> H_TDDFT_pw::global_vext_time = {0.0, 0.0, 0.0};

int H_TDDFT_pw::stype; // 0 : length gauge  1: velocity gauge

std::vector<int> H_TDDFT_pw::ttype;
//  0: Gaussian type function.
//  1: Trapezoid type function.
//  2: Trigonometric functions, sin^2.
//  3: Heaviside step function.

int H_TDDFT_pw::tstart;
int H_TDDFT_pw::tend;
double H_TDDFT_pw::dt;
// cut dt for integral
double H_TDDFT_pw::dt_int;
int H_TDDFT_pw::istep_int;
// space domain parameters

// length gauge
double H_TDDFT_pw::lcut1;
double H_TDDFT_pw::lcut2;

// velocity gauge
ModuleBase::Vector3<double> H_TDDFT_pw::At;
ModuleBase::Vector3<double> H_TDDFT_pw::At_laststep;
// hybrid gauge
ModuleBase::Vector3<double> H_TDDFT_pw::Et;
// time domain parameters

// Gauss
int H_TDDFT_pw::gauss_count;
std::vector<double> H_TDDFT_pw::gauss_omega; // time(a.u.)^-1
std::vector<double> H_TDDFT_pw::gauss_phase;
std::vector<double> H_TDDFT_pw::gauss_sigma; // time(a.u.)
std::vector<double> H_TDDFT_pw::gauss_t0;
std::vector<double> H_TDDFT_pw::gauss_amp; // Ry/bohr
std::vector<int> H_TDDFT_pw::gauss_ncut;   // cut for integral

// trapezoid
int H_TDDFT_pw::trape_count;
std::vector<double> H_TDDFT_pw::trape_omega; // time(a.u.)^-1
std::vector<double> H_TDDFT_pw::trape_phase;
std::vector<double> H_TDDFT_pw::trape_t1;
std::vector<double> H_TDDFT_pw::trape_t2;
std::vector<double> H_TDDFT_pw::trape_t3;
std::vector<double> H_TDDFT_pw::trape_amp; // Ry/bohr
std::vector<int> H_TDDFT_pw::trape_ncut;   // cut for integral

// Trigonometric
int H_TDDFT_pw::trigo_count;
std::vector<double> H_TDDFT_pw::trigo_omega1; // time(a.u.)^-1
std::vector<double> H_TDDFT_pw::trigo_omega2; // time(a.u.)^-1
std::vector<double> H_TDDFT_pw::trigo_phase1;
std::vector<double> H_TDDFT_pw::trigo_phase2;
std::vector<double> H_TDDFT_pw::trigo_amp; // Ry/bohr
std::vector<int> H_TDDFT_pw::trigo_ncut;   // cut for integral

// Heaviside
int H_TDDFT_pw::heavi_count;
std::vector<double> H_TDDFT_pw::heavi_t0;
std::vector<double> H_TDDFT_pw::heavi_amp; // Ry/bohr

void H_TDDFT_pw::current_step_info(const std::string& file_dir, int& istep)
{
    std::stringstream ssc;
    ssc << file_dir << "Restart_td.txt";
    std::ifstream file(ssc.str().c_str());

    if (!file)
    {
        ModuleBase::WARNING_QUIT("H_TDDFT_pw::current_step_info", "No Restart_td.txt!");
    }

    file >> istep;
    file >> At[0] >> At[1] >> At[2];
    file >> At_laststep[0] >> At_laststep[1] >> At_laststep[2];
    At_laststep = -At_laststep;
    file.close();
}

void H_TDDFT_pw::cal_fixed_v(double* vl_pseudo)
{
    ModuleBase::TITLE("H_TDDFT_pw", "cal_fixed_v");

    // skip if not length gauge
    if (stype != 0)
    {
        return;
    }

    // time evolve
    H_TDDFT_pw::istep++;
    H_TDDFT_pw::istep_int = istep;

    // judgement to skip vext
    if (!PARAM.inp.td_vext || istep > tend || istep < tstart)
    {
        return;
    }

    ModuleBase::timer::tick("H_TDDFT_pw", "cal_fixed_v");

    int count = 0;
    gauss_count = 0;
    trape_count = 0;
    trigo_count = 0;
    heavi_count = 0;

    global_vext_time = {0.0, 0.0, 0.0};

    for (auto direc: PARAM.inp.td_vext_dire)
    {
        std::vector<double> vext_space(this->rho_basis_->nrxx, 0.0);
        double vext_time = cal_v_time(ttype[count], true);

        global_vext_time[direc - 1] += vext_time;

        if (PARAM.inp.out_efield && GlobalV::MY_RANK == 0)
        {
            std::stringstream as;
            as << PARAM.globalv.global_out_dir << "efield_" << count << ".txt";
            std::ofstream ofs(as.str().c_str(), std::ofstream::app);
            ofs << H_TDDFT_pw::istep * dt * ModuleBase::AU_to_FS << "\t"
                << vext_time * ModuleBase::Ry_to_eV / ModuleBase::BOHR_TO_A << std::endl;
            ofs.close();
        }

        cal_v_space(vext_space, direc);
        for (size_t ir = 0; ir < this->rho_basis_->nrxx; ++ir)
        {
            vl_pseudo[ir] += vext_space[ir] * vext_time;
        }
        count++;
    }

    ModuleBase::timer::tick("H_TDDFT_pw", "cal_fixed_v");
    return;
}

void H_TDDFT_pw::cal_v_space(std::vector<double>& vext_space, int direc)
{
    ModuleBase::TITLE("H_TDDFT_pw", "cal_v_space");
    ModuleBase::timer::tick("H_TDDFT_pw", "cal_v_space");

    switch (stype)
    {
    case 0:
        cal_v_space_length(vext_space, direc);
        break;
    default:
        std::cout << "space_domain_type of electric field is wrong" << std::endl;
        break;
    }

    ModuleBase::timer::tick("H_TDDFT_pw", "cal_v_space");
    return;
}

void H_TDDFT_pw::cal_v_space_length(std::vector<double>& vext_space, int direc)
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
            vext_space[ir] = cal_v_space_length_potential(x) * this->ucell_->latvec.e11
                             + cal_v_space_length_potential(y) * this->ucell_->latvec.e21
                             + cal_v_space_length_potential(z) * this->ucell_->latvec.e31;
            break;

        case 2:
            vext_space[ir] = cal_v_space_length_potential(x) * this->ucell_->latvec.e12
                             + cal_v_space_length_potential(y) * this->ucell_->latvec.e22
                             + cal_v_space_length_potential(z) * this->ucell_->latvec.e32;
            break;

        case 3:
            vext_space[ir] = cal_v_space_length_potential(x) * this->ucell_->latvec.e13
                             + cal_v_space_length_potential(y) * this->ucell_->latvec.e23
                             + cal_v_space_length_potential(z) * this->ucell_->latvec.e33;
            break;

        default:
            std::cout << "direction of electric field is wrong" << std::endl;
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
        vext_space = -((i - lcut1) * (lcut2 - lcut1) / (lcut1 + 1.0 - lcut2) - lcut1) * this->ucell_->lat0;
    }
    else if (i >= lcut1 && i < lcut2)
    {
        vext_space = i * this->ucell_->lat0;
    }
    else if (i >= lcut2)
    {
        vext_space = -((i - lcut2) * (lcut2 - lcut1) / (lcut1 + 1.0 - lcut2) - lcut2) * this->ucell_->lat0;
    }
    return vext_space;
}

int H_TDDFT_pw::check_ncut(int t_type)
{
    int ncut = 0;
    switch (t_type)
    {
    case 0:
        ncut = *(gauss_ncut.begin() + gauss_count);
        break;

    case 1:
        ncut = *(trape_ncut.begin() + trape_count);
        break;

    case 2:
        ncut = *(trigo_ncut.begin() + trigo_count);
        break;

    case 3:
        ncut = 2;
        break;

    default:
        std::cout << "time_domain_type of electric field is wrong" << std::endl;
        break;
    }
    return ncut;
}

void H_TDDFT_pw::update_At()
{
    // time evolve
    H_TDDFT_pw::istep++;
    // midpoint rule should be used both in Hamiltonian and here.
    At = At + At_laststep / 2.0;
    At_laststep.set(0.0, 0.0, 0.0);
    Et.set(0.0, 0.0, 0.0);

    // judgement to skip vext
    if (!PARAM.inp.td_vext || istep > tend || istep < tstart)
    {
        return;
    }

    ModuleBase::timer::tick("H_TDDFT_pw", "update_At");

    int count = 0;
    gauss_count = 0;
    trape_count = 0;
    trigo_count = 0;
    heavi_count = 0;
    // parameters for integral
    int ncut = 1;
    bool last = false;
    double out = 0.0;

    for (auto direc: PARAM.inp.td_vext_dire)
    {
        last = false;
        // cut the integral space and initialize relevant parameters
        ncut = check_ncut(ttype[count]);
        istep_int = istep * ncut;
        dt_int = dt / double(ncut);

        // store vext_time for each time point, include the first and last point
        std::vector<double> vext_time(ncut + 1, 0.0); // Use std::vector to manage memory
        for (int i = 0; i <= ncut; i++)
        {
            // if this is the last point, type_count++
            if (i == ncut)
            {
                last = true;
            }
            vext_time[i] = cal_v_time(ttype[count], last);
            istep_int++;
        }
        // Call the Simpson's rule integration using std::vector data
        ModuleBase::Integral::Simpson_Integral(ncut + 1, vext_time.data(), dt_int, out);

        // update At value for its direction
        switch (stype)
        {
        case 1:
            At_laststep[direc - 1] -= out;
            break;
        case 2:
            At_laststep[direc - 1] -= out;
            Et[direc - 1] += vext_time[0];
            break;
        default:
            std::cout << "space_domain_type of electric field is wrong" << std::endl;
            break;
        }

        // output Efield
        if (PARAM.inp.out_efield && GlobalV::MY_RANK == 0)
        {
            std::stringstream as;
            as << PARAM.globalv.global_out_dir << "efield_" << count << ".txt";
            std::ofstream ofs(as.str().c_str(), std::ofstream::app);
            ofs << H_TDDFT_pw::istep * dt * ModuleBase::AU_to_FS << "\t"
                << vext_time[0] * ModuleBase::Ry_to_eV / ModuleBase::BOHR_TO_A << std::endl;
            ofs.close();
        }
        // total count++
        count++;
    }
    At = At + At_laststep / 2.0;

    ModuleBase::timer::tick("H_TDDFT_pw", "update_At");
    return;
}

double H_TDDFT_pw::cal_v_time(int t_type, const bool last)
{
    double vext_time = 0.0;

    switch (t_type)
    {
    case 0:
        vext_time = cal_v_time_Gauss(last);
        break;

    case 1:
        vext_time = cal_v_time_trapezoid(last);
        break;

    case 2:
        vext_time = cal_v_time_trigonometric(last);
        break;

    case 3:
        vext_time = cal_v_time_heaviside(last);
        break;

    default:
        std::cout << "time_domain_type of electric field is wrong" << std::endl;
        break;
    }
    return vext_time;
}

double H_TDDFT_pw::cal_v_time_Gauss(const bool last)
{
    double vext_time = 0.0;
    double t0 = *(gauss_t0.begin() + gauss_count);
    double omega = *(gauss_omega.begin() + gauss_count);
    double sigma = *(gauss_sigma.begin() + gauss_count);
    double phase = *(gauss_phase.begin() + gauss_count);
    double amp = *(gauss_amp.begin() + gauss_count);
    double ncut = *(gauss_ncut.begin() + gauss_count);

    double gauss_t = (istep_int - t0 * ncut) * dt_int;
    vext_time = cos(omega * gauss_t + phase) * exp(-gauss_t * gauss_t * 0.5 / (sigma * sigma)) * amp;
    if (last)
    {
        gauss_count++;
    }

    return vext_time;
}

double H_TDDFT_pw::cal_v_time_trapezoid(const bool last)
{
    double vext_time = 0.0;
    double t1 = *(trape_t1.begin() + trape_count);
    double t2 = *(trape_t2.begin() + trape_count);
    double t3 = *(trape_t3.begin() + trape_count);
    double omega = *(trape_omega.begin() + trape_count);
    double phase = *(trape_phase.begin() + trape_count);
    double amp = *(trape_amp.begin() + trape_count);
    double ncut = *(trape_ncut.begin() + trape_count);

    if (istep < t1)
    {
        vext_time = istep_int / ncut / t1;
    }
    else if (istep < t2)
    {
        vext_time = 1.0;
    }
    else if (istep < t3)
    {
        vext_time = (t3 - istep_int / ncut) / (t3 - t2);
    }

    vext_time = vext_time * amp * cos(omega * istep_int * dt_int + phase);
    if (last)
    {
        trape_count++;
    }

    return vext_time;
}

double H_TDDFT_pw::cal_v_time_trigonometric(const bool last)
{
    double vext_time = 0.0;
    double omega1 = *(trigo_omega1.begin() + trigo_count);
    double phase1 = *(trigo_phase1.begin() + trigo_count);
    double omega2 = *(trigo_omega2.begin() + trigo_count);
    double phase2 = *(trigo_phase2.begin() + trigo_count);
    double amp = *(trigo_amp.begin() + trigo_count);

    const double timenow = istep_int * dt_int;

    vext_time = amp * cos(omega1 * timenow + phase1) * sin(omega2 * timenow + phase2) * sin(omega2 * timenow + phase2);
    if (last)
    {
        trigo_count++;
    }

    return vext_time;
}

double H_TDDFT_pw::cal_v_time_heaviside(const bool last)
{
    double t0 = *(heavi_t0.begin() + heavi_count);
    double amp = *(heavi_amp.begin() + heavi_count);
    double vext_time = 0.0;
    if (istep < t0)
    {
        vext_time = amp;
    }
    else if (istep >= t0)
    {
        vext_time = 0.0;
    }
    if (last)
    {
        heavi_count++;
    }

    return vext_time;
}

void H_TDDFT_pw::compute_force(const UnitCell& cell, ModuleBase::matrix& fe)
{
    int iat = 0;
    for (int it = 0; it < cell.ntype; ++it)
    {
        for (int ia = 0; ia < cell.atoms[it].na; ++ia)
        {
            for (int direc = 0; direc < 3; ++direc)
            {
                // No need to multiply ModuleBase::e2, since the unit of force is Ry/Bohr
                fe(iat, direc) = global_vext_time[direc] * cell.atoms[it].ncpp.zv;
            }
            ++iat;
        }
    }
}

} // namespace elecstate