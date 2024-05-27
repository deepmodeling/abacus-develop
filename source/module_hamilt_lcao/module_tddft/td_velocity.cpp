#include "td_velocity.h"

#include "module_base/timer.h"
#include "module_elecstate/potentials/H_TDDFT_pw.h"

bool TD_Velocity::tddft_velocity = false;
bool TD_Velocity::out_mat_R = false;
bool TD_Velocity::out_vecpot = false;
bool TD_Velocity::init_vecpot_file = false;

TD_Velocity* TD_Velocity::td_vel_op = nullptr;

int TD_Velocity::istep = -1;
int TD_Velocity::max_istep = -1;
std::vector<ModuleBase::Vector3<double>> TD_Velocity::At_from_file;

TD_Velocity::TD_Velocity()
{
    if (init_vecpot_file && istep == -1)
    {
        this->read_cart_At();
    }
    return;
}
TD_Velocity::~TD_Velocity()
{
    this->destroy_HS_R_td_sparse();
}

void TD_Velocity::output_cart_At(const std::string& out_dir)
{
    if (GlobalV::MY_RANK == 0)
    {
        std::string out_file;
        // generate the output file name
        out_file = out_dir + "At.dat";
        std::ofstream ofs;
        // output title
        if (istep == 0)
        {
            ofs.open(out_file.c_str(), std::ofstream::out);
            ofs << std::left << std::setw(8) << "#istep" << std::setw(15) << "A_x" << std::setw(15) << "A_y"
                << std::setw(15) << "A_z" << std::endl;
        }
        else
        {
            ofs.open(out_file.c_str(), std::ofstream::app);
        }
        // output the vector potential
        ofs << std::left << std::setw(8) << istep;
        // divide by 2.0 to get the atomic unit
        for (int i = 0; i < 3; i++)
        {
            ofs << std::scientific << std::setprecision(4) << std::setw(15) << cart_At[i] / 2.0;
        }
        ofs << std::endl;
        ofs.close();
    }
    return;
}

void TD_Velocity::cal_cart_At(const ModuleBase::Vector3<double>& a0,
                              const ModuleBase::Vector3<double>& a1,
                              const ModuleBase::Vector3<double>& a2,
                              const ModuleBase::Vector3<double>& At)
{
    istep++;
    if (init_vecpot_file)
    {
        this->cart_At = At_from_file[istep > max_istep ? max_istep : istep];
    }
    else
    {
        const double l_norm[3] = {a0.norm(), a1.norm(), a2.norm()};
        this->cart_At = a0 * At[0] / l_norm[0] + a1 * At[1] / l_norm[1] + a2 * At[2] / l_norm[2];
    }
    // output the vector potential if needed
    if (out_vecpot == true)
    {
        this->output_cart_At(GlobalV::global_out_dir);
    }
}

void TD_Velocity::read_cart_At(void)
{
    std::string in_file;
    // generate the input file name
    in_file = "At.dat";
    std::ifstream ifs(in_file.c_str());
    // check if the file is exist
    if (!ifs)
    {
        ModuleBase::WARNING_QUIT("TD_Velocity::read_cart_At", "Cannot open Vector potential file!");
    }
    std::string line;
    std::vector<std::string> str_vec;
    // use tmp to skip the istep number
    int tmp = 0;
    while (std::getline(ifs, line))
    {
        // A tmporary vector3 to store the data of this line
        ModuleBase::Vector3<double> At;
        if (line[0] == '#')
        {
            continue;
        }
        std::istringstream iss(line);
        // skip the istep number
        if (!(iss >> tmp))
        {
            ModuleBase::WARNING_QUIT("TD_Velocity::read_cart_At", "Error reading istep!");
        }
        // read the vector potential
        double component = 0;
        // Read three components
        for (int i = 0; i < 3; i++)
        {
            if (!(iss >> component))
            {
                ModuleBase::WARNING_QUIT("TD_Velocity::read_cart_At",
                                         "Error reading component " + std::to_string(i + 1) + " for istep "
                                             + std::to_string(tmp) + "!");
            }
            At[i] = component;
        }
        // unit transform ,change unit into Ry/bohr/e*t_a.u.
        At *= 2.0;
        // add the tmporary vector3 to the vector potential vector
        At_from_file.push_back(At);
    }
    // set the max_istep
    max_istep = At_from_file.size() - 1;
    ifs.close();

    return;
}
void TD_Velocity::destroy_HS_R_td_sparse(void)
{
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>>
        empty_HR_sparse_td_vel_up;
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>>
        empty_HR_sparse_td_vel_down;
    HR_sparse_td_vel[0].swap(empty_HR_sparse_td_vel_up);
    HR_sparse_td_vel[1].swap(empty_HR_sparse_td_vel_down);
}
