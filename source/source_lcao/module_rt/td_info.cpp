#include "td_info.h"

#include "source_estate/module_pot/H_TDDFT_pw.h"
#include "source_io/module_parameter/parameter.h"

bool TD_info::out_mat_R = false;
bool TD_info::out_vecpot = false;
bool TD_info::out_current = false;
bool TD_info::out_current_k = false;
bool TD_info::init_vecpot_file = false;
bool TD_info::evolve_once = false;

TD_info* TD_info::td_vel_op = nullptr;

int TD_info::estep_shift = 0;
int TD_info::istep = -1;
int TD_info::max_istep = -1;
ModuleBase::Vector3<double> TD_info::cart_At;
std::vector<ModuleBase::Vector3<double>> TD_info::At_from_file;

TD_info::TD_info(const UnitCell* ucell_in)
{
    this->ucell = ucell_in;
    if (init_vecpot_file && istep == -1)
    {
        this->read_cart_At();
    }
    //read in restart step
    if(PARAM.inp.mdp.md_restart)
    {
        std::stringstream ssc;
        ssc << PARAM.globalv.global_readin_dir << "Restart_td.txt";
        std::ifstream file(ssc.str().c_str());
        if (!file)
        {
            ModuleBase::WARNING_QUIT("TD_info::TD_info", "No Restart_td.txt!");
        }
        file >> estep_shift;
        //std::cout<<"estep_shift"<<estep_shift<<std::endl;
    }
    this->istep += estep_shift;
    return;
}
TD_info::~TD_info()
{
    if(elecstate::H_TDDFT_pw::stype == 1)
    {
        this->destroy_HS_R_td_sparse();
    }
    for (int dir = 0; dir < 3; dir++)
    {
        if (this->current_term[dir] != nullptr)
        {
            delete this->current_term[dir];
        }
    }
}

void TD_info::output_cart_At(const std::string& out_dir)
{
    if (GlobalV::MY_RANK == 0)
    {
        std::string out_file;
        // generate the output file name
        out_file = out_dir + "At.dat";
        std::ofstream ofs;
        // output title
        if (istep == estep_shift)
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
            ofs << std::scientific << std::setprecision(4) << std::setw(15) << cart_At[i];
        }
        ofs << std::endl;
        ofs.close();
    }
    return;
}

void TD_info::cal_cart_At(const ModuleBase::Vector3<double>& At)
{
    istep++;
    if (init_vecpot_file)
    {
        cart_At = At_from_file[istep > max_istep ? max_istep : istep];
    }
    else
    {
        // transfrom into atomic unit
        cart_At = At / 2.0;
    }
    // output the vector potential if needed
    if (out_vecpot == true)
    {
        this->output_cart_At(PARAM.globalv.global_out_dir);
    }
}

void TD_info::read_cart_At(void)
{
    std::string in_file;
    // generate the input file name
    in_file = "At.dat";
    std::ifstream ifs(in_file.c_str());
    // check if the file is exist
    if (!ifs)
    {
        ModuleBase::WARNING_QUIT("TD_info::read_cart_At", "Cannot open Vector potential file!");
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
            ModuleBase::WARNING_QUIT("TD_info::read_cart_At", "Error reading istep!");
        }
        // read the vector potential
        double component = 0;
        // Read three components
        for (int i = 0; i < 3; i++)
        {
            if (!(iss >> component))
            {
                ModuleBase::WARNING_QUIT("TD_info::read_cart_At",
                                         "Error reading component " + std::to_string(i + 1) + " for istep "
                                             + std::to_string(tmp) + "!");
            }
            At[i] = component;
        }
        // add the tmporary vector3 to the vector potential vector
        At_from_file.push_back(At);
    }
    // set the max_istep
    max_istep = At_from_file.size() - 1;
    ifs.close();

    return;
}
void TD_info::out_restart_info(const int nstep, 
                      const ModuleBase::Vector3<double>& At_current, 
                      const ModuleBase::Vector3<double>& At_laststep)
{
    if (GlobalV::MY_RANK == 0)
    {
        // open file
        std::string outdir = PARAM.globalv.global_out_dir + "Restart_td.txt";
        std::ofstream outFile(outdir);
        if (!outFile) {
            ModuleBase::WARNING_QUIT("out_restart_info", "no Restart_td.txt!");
        }
        // write data
        outFile << nstep << std::endl;
        outFile << At_current[0] << " " << At_current[1] << " " << At_current[2] << std::endl;
        outFile << At_laststep[0] << " " << At_laststep[1] << " " << At_laststep[2] << std::endl;
        outFile.close();
    }
    

    return;
}

void TD_info::initialize_current_term(const hamilt::HContainer<std::complex<double>>* HR,
                                          const Parallel_Orbitals* paraV)
{
    ModuleBase::TITLE("TD_info", "initialize_current_term");
    ModuleBase::timer::tick("TD_info", "initialize_current_term");

    for (int dir = 0; dir < 3; dir++)
    {
        if (this->current_term[dir] == nullptr)
            this->current_term[dir] = new hamilt::HContainer<std::complex<double>>(paraV);
    }

    for (int i = 0; i < HR->size_atom_pairs(); ++i)
    {
        hamilt::AtomPair<std::complex<double>>& tmp = HR->get_atom_pair(i);
        for (int ir = 0; ir < tmp.get_R_size(); ++ir)
        {
            const ModuleBase::Vector3<int> R_index = tmp.get_R_index(ir);
            const int iat1 = tmp.get_atom_i();
            const int iat2 = tmp.get_atom_j();

            hamilt::AtomPair<std::complex<double>> tmp1(iat1, iat2, R_index, paraV);
            for (int dir = 0; dir < 3; dir++)
            {
                this->current_term[dir]->insert_pair(tmp1);
            }
        }
    }
    for (int dir = 0; dir < 3; dir++)
    {
        this->current_term[dir]->allocate(nullptr, true);
    }

    ModuleBase::timer::tick("TD_info", "initialize_current_term");
}

void TD_info::destroy_HS_R_td_sparse(void)
{
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>>
        empty_HR_sparse_td_vel_up;
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>>
        empty_HR_sparse_td_vel_down;
    HR_sparse_td_vel[0].swap(empty_HR_sparse_td_vel_up);
    HR_sparse_td_vel[1].swap(empty_HR_sparse_td_vel_down);
}
