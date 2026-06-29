#include "td_info.h"

#include "source_base/libm/libm.h"
#include "source_estate/module_pot/H_TDDFT_pw.h"
#include "source_io/module_parameter/parameter.h"

bool TD_info::out_mat_R = false;
bool TD_info::out_vecpot = false;
int TD_info::out_current = 0;
bool TD_info::out_current_k = false;
bool TD_info::init_vecpot_file = false;
bool TD_info::evolve_once = false;

TD_info* TD_info::td_vel_op = nullptr;

int TD_info::estep_shift = 0;
int TD_info::istep = -1;
int TD_info::max_istep = -1;
ModuleBase::Vector3<double> TD_info::cart_At;
std::vector<ModuleBase::Vector3<double>> TD_info::At_from_file;

TD_info::TD_info(const UnitCell* ucell_in,const Parallel_Orbitals& pv, const LCAO_Orbitals& orb)
{
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
    if(out_current==2||elecstate::H_TDDFT_pw::stype == 2)
    {
        r_calculator.init(*ucell_in, pv, orb);
    }
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
    for (int dir = 0; dir < 3; dir++)
    {
        if (this->grad_overlap[dir] != nullptr)
        {
            delete this->grad_overlap[dir];
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
    // update hybrid gauge phase
    if(elecstate::H_TDDFT_pw::stype == 2)
    {
        for(const auto& phase_pair : phase_hybrid)
        {
            const ModuleBase::Vector3<int>& r_index = phase_pair.first;
            ModuleBase::Vector3<double> dR = double(r_index.x) * a1 + double(r_index.y) * a2 + double(r_index.z) * a3;
            const double arg_td = cart_At * dR * lat0;
            double sinp, cosp;
            ModuleBase::libm::sincos(arg_td, &sinp, &cosp);
            phase_hybrid[r_index] = std::complex<double>(cosp, sinp);
        }
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
template <typename TR>
void TD_info::initialize_phase_hybrid(const UnitCell& ucell, const hamilt::HContainer<TR>* hR)
{
    this->a1 = ucell.a1;
    this->a2 = ucell.a2;
    this->a3 = ucell.a3;
    this->lat0 = ucell.lat0;
    for (int i = 0; i < hR->size_atom_pairs(); ++i)
    {
        hamilt::AtomPair<TR>& tmp = hR->get_atom_pair(i);
        for(int ir = 0;ir < tmp.get_R_size(); ++ir )
        {
            const ModuleBase::Vector3<int> r_index = tmp.get_R_index(ir);
            if(phase_hybrid.count(r_index))continue;


            ModuleBase::Vector3<double> dR = double(r_index.x) * a1 + double(r_index.y) * a2 + double(r_index.z) * a3;
            const double arg_td = cart_At * dR * lat0;
            double sinp, cosp;
            ModuleBase::libm::sincos(arg_td, &sinp, &cosp);
            phase_hybrid[r_index] = std::complex<double>(cosp, sinp);
        }
    }
}
void TD_info::initialize_current_term(const hamilt::HContainer<std::complex<double>>* HR,
                                          const Parallel_Orbitals* paraV)
{
    ModuleBase::TITLE("TD_info", "initialize_current_term");
    ModuleBase::timer::start("TD_info", "initialize_current_term");

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

    ModuleBase::timer::end("TD_info", "initialize_current_term");
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

void TD_info::calculate_grad_overlap(const Parallel_Orbitals& paraV,
                                     const UnitCell& ucell,
                                     const Grid_Driver& GridD,
                                     const std::vector<double>& orb_cutoff,
                                     const TwoCenterIntegrator* intor)
{
    ModuleBase::TITLE("TD_info", "calculate_grad_overlap");
    ModuleBase::timer::start("TD_info", "calculate_grad_overlap");
    for (int dir=0;dir<3;dir++)
    {
        if (this->grad_overlap[dir] != nullptr)
        {
            delete this->grad_overlap[dir];
        }
        this->grad_overlap[dir] = new hamilt::HContainer<double>(&paraV);
    }
    for (int iat1 = 0; iat1 < ucell.nat; iat1++)
    {
        auto tau1 = ucell.get_tau(iat1);
        int T1=0;
        int I1=0;
        ucell.iat2iait(iat1, &I1, &T1);
        AdjacentAtomInfo adjs;
        GridD.Find_atom(ucell, tau1, T1, I1, &adjs);
        for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
        {
            const int T2 = adjs.ntype[ad];
            const int I2 = adjs.natom[ad];
            int iat2 = ucell.itia2iat(T2, I2);
            if (paraV.get_nrow_atom(iat1) <= 0 || paraV.get_ncol_atom(iat2) <= 0)
            {
                continue;
            }
            const ModuleBase::Vector3<int>& R_index = adjs.box[ad];
            // choose the real adjacent atoms
            // Note: the distance of atoms should less than the cutoff radius,
            // When equal, the theoretical value of matrix element is zero,
            // but the calculated value is not zero due to the numerical error, which would lead to result changes.
            if (ucell.cal_dtau(iat1, iat2, R_index).norm() * ucell.lat0
                >= orb_cutoff[T1] + orb_cutoff[T2])
            {
                continue;
            }
            hamilt::AtomPair<double> tmp(iat1, iat2, R_index, &paraV);
            for (int dir=0;dir<3;dir++)
            {
                this->grad_overlap[dir]->insert_pair(tmp);
            }
        }
    }
    // allocate the memory of BaseMatrix in grad_overlap, and set the new values to zero
    for (int dir=0;dir<3;dir++)
    {
        this->grad_overlap[dir]->allocate(nullptr, true);
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int iap = 0; iap < this->grad_overlap[0]->size_atom_pairs(); ++iap)
    {
        hamilt::AtomPair<double>& tmp = this->grad_overlap[0]->get_atom_pair(iap);
        const int iat1 = tmp.get_atom_i();
        const int iat2 = tmp.get_atom_j();

        for (int iR = 0; iR < tmp.get_R_size(); ++iR)
        {
            const ModuleBase::Vector3<int> R_index = tmp.get_R_index(iR);
            auto dtau = ucell.cal_dtau(iat1, iat2, R_index);

            double* p_data[3] = {nullptr, nullptr, nullptr};
            for (int i = 0; i < 3; i++)
            {
                p_data[i] = this->grad_overlap[i]->get_atom_pair(iap).get_pointer(iR);
            }
            // ---------------------------------------------
            // get info of orbitals of atom1 and atom2 from ucell
            // ---------------------------------------------
            int T1=0;
            int I1=0;
            ucell.iat2iait(iat1, &I1, &T1);
            int T2=0;
            int I2=0;
            ucell.iat2iait(iat2, &I2, &T2);
            Atom& atom1 = ucell.atoms[T1];
            Atom& atom2 = ucell.atoms[T2];

            // npol is the number of polarizations,
            // 1 for non-magnetic (one Hamiltonian matrix only has spin-up or spin-down),
            // 2 for magnetic (one Hamiltonian matrix has both spin-up and spin-down)
            const int npol = ucell.get_npol();

            const int* iw2l1 = atom1.iw2l.data();
            const int* iw2n1 = atom1.iw2n.data();
            const int* iw2m1 = atom1.iw2m.data();
            const int* iw2l2 = atom2.iw2l.data();
            const int* iw2n2 = atom2.iw2n.data();
            const int* iw2m2 = atom2.iw2m.data();

            // ---------------------------------------------
            // calculate the overlap matrix for each pair of orbitals
            // ---------------------------------------------
            double olm[3] = {0, 0, 0};
            auto row_indexes = paraV.get_indexes_row(iat1);
            auto col_indexes = paraV.get_indexes_col(iat2);
            const int step_trace = col_indexes.size() + 1;
            for (int iw1l = 0; iw1l < row_indexes.size(); iw1l += npol)
            {
                const int iw1 = row_indexes[iw1l] / npol;
                const int L1 = iw2l1[iw1];
                const int N1 = iw2n1[iw1];
                const int m1 = iw2m1[iw1];

                // convert m (0,1,...2l) to M (-l, -l+1, ..., l-1, l)
                int M1 = (m1 % 2 == 0) ? -m1 / 2 : (m1 + 1) / 2;

                for (int iw2l = 0; iw2l < col_indexes.size(); iw2l += npol)
                {
                    const int iw2 = col_indexes[iw2l] / npol;
                    const int L2 = iw2l2[iw2];
                    const int N2 = iw2n2[iw2];
                    const int m2 = iw2m2[iw2];

                    // convert m (0,1,...2l) to M (-l, -l+1, ..., l-1, l)
                    int M2 = (m2 % 2 == 0) ? -m2 / 2 : (m2 + 1) / 2;
                    intor->calculate(T1, L1, N1, M1, T2, L2, N2, M2, dtau * ucell.lat0, nullptr, olm);
                    for (int dir = 0; dir < 3; dir++)
                    {
                        for (int ipol = 0; ipol < npol; ipol++)
                        {
                            p_data[dir][ipol * step_trace] += olm[dir];
                        }
                        p_data[dir] += npol;
                    }
                }
                for (int dir = 0; dir < 3; dir++)
                {
                    p_data[dir] += (npol - 1) * col_indexes.size();
                }
                
            }
        }
    }
    ModuleBase::timer::end("TD_info", "calculate_grad_overlap");
}
template
void TD_info::initialize_phase_hybrid<std::complex<double>>(const UnitCell& ucell, const hamilt::HContainer<std::complex<double>>* hR);
template
void TD_info::initialize_phase_hybrid<double>(const UnitCell& ucell, const hamilt::HContainer<double>* hR);