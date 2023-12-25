#include "module_io/to_qo.h"
#include "module_base/libm/libm.h"

void toQO::unwrap_unitcell(UnitCell* p_ucell)
{
    p_ucell_ = p_ucell;
    ntype_ = p_ucell->ntype;
    std::for_each(p_ucell->atoms, p_ucell->atoms + p_ucell->ntype, [this](Atom& atom){
        symbols_.push_back(atom.label);
        charges_.push_back(atom.ncpp.zv);
    });
    nmax_.resize(ntype_);
    na_.resize(ntype_);
    for(int itype = 0; itype < ntype_; itype++)
    {
        if(strategy_ != "energy")
        {
            nmax_[itype] = atom_database_.principle_quantum_number[symbols_[itype]];
        }
        else
        {
            nmax_[itype] = atom_database_.atom_Z[symbols_[itype]];
        }
        na_[itype] = p_ucell_->atoms[itype].na;
    }
}

template <typename T>
void toQO::eliminate_duplicate_vector3(std::vector<ModuleBase::Vector3<T>> &v)
{
    std::vector<std::vector<T>> v_;
    // convert vector3 to vector
    for(int i = 0; i < v.size(); i++)
    {
        v_.push_back(std::vector<T>{v[i].x, v[i].y, v[i].z});
    }
    std::sort(v_.begin(), v_.end());
    v_.erase(std::unique(v_.begin(), v_.end()), v_.end());
    v.clear();
    v.resize(v_.size());
    for(int i = 0; i < v_.size(); i++)
    {
        v[i] = ModuleBase::Vector3<T>(v_[i][0], v_[i][1], v_[i][2]);
    }
}
template void toQO::eliminate_duplicate_vector3<int>(std::vector<ModuleBase::Vector3<int>> &v);

void toQO::deallocate_ovlp(const bool is_R)
{
    if(nchi_ == 0 || nphi_ == 0)
    {
        ModuleBase::WARNING_QUIT("toQO::deallocate_ovlp", "nchi_ or nphi_ is zero, which means not properly initialized.");
    }
    if(is_R)
    {
        int nR = save_mem_? 1 : nR_;
        for(int iR = 0; iR < nR; iR++)
        {
            for(int i = 0; i < ovlp_R_[iR].size(); i++)
            {
                ovlp_R_[iR][i].clear();
            }
            ovlp_R_[iR].clear();
        }
        ovlp_R_.clear();
    }
    else
    {
        for(int i = 0; i < ovlp_k_.size(); i++)
        {
            ovlp_k_[i].clear();
        }
        ovlp_k_.clear();
    }
}

void toQO::allocate_ovlp(const bool is_R)
{

    if(nchi_ == 0 || nphi_ == 0)
    {
        ModuleBase::WARNING_QUIT("toQO::allocate_ovlp", "nchi_ or nphi_ is zero, which means not properly initialized.");
    }
    if(is_R)
    {
        int nR = save_mem_? 1 : nR_;
        // allocate memory for ovlp_R_
        for(int iR = 0; iR < nR; iR++)
        {
            RealMatrix matrix;
            for(int i = 0; i < nchi_; i++)
            {
                std::vector<double> row;
                for(int j = 0; j < nphi_; j++)
                {
                    row.push_back(0.0);
                }
                matrix.push_back(row);
            }
            ovlp_R_.push_back(matrix);
        }
    }
    else
    {
        // allocate memory for ovlp_k_
        for(int i = 0; i < nchi_; i++)
        {
            std::vector<std::complex<double>> row;
            for(int j = 0; j < nphi_; j++)
            {
                row.push_back(std::complex<double>(0.0, 0.0));
            }
            ovlp_k_.push_back(row);
        }
    }
}

void toQO::clean_up()
{
    if(nchi_ == 0 || nphi_ == 0)
    {
        ModuleBase::WARNING_QUIT("toQO::allocate", "nchi_ or nphi_ is zero, which means not properly initialized.");
    }
    deallocate_ovlp(true);
    deallocate_ovlp(false);
    allocate_ovlp(true);
    allocate_ovlp(false);
}

void toQO::zero_out_ovlps(const bool is_R)
{
    if(nchi_ == 0 || nphi_ == 0)
    {
        ModuleBase::WARNING_QUIT("toQO::zero_out_ovlps", "nchi_ or nphi_ is zero, which means not properly initialized.");
    }
    if(is_R)
    {
        int nR = save_mem_? 1 : nR_;
        for(int iR = 0; iR < nR; iR++)
        {
            for(int i = 0; i < nchi_; i++)
            {
                for(int j = 0; j < nphi_; j++)
                {
                    ovlp_R_[iR][i][j] = 0.0;
                }
            }
        }
    }
    else
    {
        for(int i = 0; i < nchi_; i++)
        {
            for(int j = 0; j < nphi_; j++)
            {
                ovlp_k_[i][j] = std::complex<double>(0.0, 0.0);
            }
        }
    }
}

std::vector<ModuleBase::Vector3<int>> toQO::scan_supercell_for_atom(int it, int ia, int start_it, int start_ia)
{
    std::vector<ModuleBase::Vector3<int>> n1n2n3;
    // cutoff radius of numerical atomic orbital of atom itia
    double rcut_i = ao_->rcut_max();
    if(rcut_i > 10)
    {
        #ifdef __MPI
        if(GlobalV::MY_RANK == 0)
        {
        #endif
        std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl
                  << "! Warning: rcut_i of atom in type " << it << " and index " << ia << " is larger than 10 bohr: " << std::fixed << rcut_i << " ." << std::endl
                  << "! This value has been larger than maximal cutoff radius of numerical orbitals, " << std::endl
                  << "! will brings about high computational cost and make k <-> R transform" << std::endl
                  << "! possibly not reversible. Suggest to try other qo_basis." << std::endl
                  << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
        #ifdef __MPI
        }
        #endif
    }
    // lattice vectors
    for(int itype = start_it; itype < p_ucell_->ntype; itype++)
    {
        for(int iatom = start_ia; iatom < p_ucell_->atoms[itype].na; iatom++)
        {
            double rcut_j = nao_->rcut_max();
            ModuleBase::Vector3<double> rij = p_ucell_->atoms[itype].tau[iatom] - p_ucell_->atoms[it].tau[ia]; // in unit lat0?
            int n1 = 0; int n2 = 0; int n3 = 0;
            // calculate the sup of n1, n2, n3
            // rcut_i, j in bohr! a1, a2 and a3 are in lat0, so multiply with lat0
            // int n1max = int(std::ceil((rcut_i + rcut_j)/p_ucell_->a1.norm()/p_ucell_->lat0));
            // int n2max = int(std::ceil((rcut_i + rcut_j)/p_ucell_->a2.norm()/p_ucell_->lat0));
            // int n3max = int(std::ceil((rcut_i + rcut_j)/p_ucell_->a3.norm()/p_ucell_->lat0));
            ModuleBase::Vector3<double> a1_in_Bohr = p_ucell_->a1 * p_ucell_->lat0;
            ModuleBase::Vector3<double> a2_in_Bohr = p_ucell_->a2 * p_ucell_->lat0;
            ModuleBase::Vector3<double> a3_in_Bohr = p_ucell_->a3 * p_ucell_->lat0;
            double rcut_ij = rcut_i + rcut_j;
            std::vector<int> n1n2n3_max = rcut_to_supercell_index(rcut_ij, a1_in_Bohr, a2_in_Bohr, a3_in_Bohr);
            // scan n1, n2, n3
            for(int n1 = -n1n2n3_max[0]; n1 <= n1n2n3_max[0]; n1++)
            {
                for(int n2 = -n1n2n3_max[1]; n2 <= n1n2n3_max[1]; n2++)
                {
                    for(int n3 = -n1n2n3_max[2]; n3 <= n1n2n3_max[2]; n3++)
                    {
                        double f = norm2_rij_supercell(rij, n1, n2, n3);
                        if(f < std::pow(rcut_i + rcut_j, 2))
                        {
                            n1n2n3.push_back(ModuleBase::Vector3<int>(n1, n2, n3));
                        }
                    }
                }
            }
        }
    }
    eliminate_duplicate_vector3<int>(n1n2n3);
    return n1n2n3;
}

double cosine_between_vector3(ModuleBase::Vector3<double> v1, ModuleBase::Vector3<double> v2)
{
    double f = v1 * v2;
    f /= v1.norm();
    f /= v2.norm();
    return f;
}

std::vector<int> toQO::rcut_to_supercell_index(double rcut, ModuleBase::Vector3<double> a, ModuleBase::Vector3<double> b, ModuleBase::Vector3<double> c)
{
    double fab = std::sqrt(1-std::pow(cosine_between_vector3(a, b), 2));
    double fac = std::sqrt(1-std::pow(cosine_between_vector3(a, c), 2));
    double fbc = std::sqrt(1-std::pow(cosine_between_vector3(b ,c), 2));
    double fa = std::min(fab, fac);
    double fb = std::min(fab, fbc);
    double fc = std::min(fac, fbc);
    int n1max = int(std::ceil(rcut/a.norm()/fa));
    int n2max = int(std::ceil(rcut/b.norm()/fb));
    int n3max = int(std::ceil(rcut/c.norm()/fc));
    std::vector<int> n1n2n3 = {n1max, n2max, n3max};
    return n1n2n3;
}

double toQO::norm2_rij_supercell(ModuleBase::Vector3<double> rij, int n1, int n2, int n3)
{
    double f = rij * rij;
    f += n1*n1*(p_ucell_->a1*p_ucell_->a1);
    f += n2*n2*(p_ucell_->a2*p_ucell_->a2);
    f += n3*n3*(p_ucell_->a3*p_ucell_->a3);
    f += 2*n1*n2*(p_ucell_->a1*p_ucell_->a2);
    f += 2*n1*n3*(p_ucell_->a1*p_ucell_->a3);
    f += 2*n2*n3*(p_ucell_->a2*p_ucell_->a3);
    f += 2*n1*(p_ucell_->a1*rij);
    f += 2*n2*(p_ucell_->a2*rij);
    f += 2*n3*(p_ucell_->a3*rij);
    f += rij*rij;
    return f;
}

void toQO::scan_supercell()
{
    std::vector<ModuleBase::Vector3<int>> n1n2n3_overall;
    for(int it = 0; it < p_ucell_->ntype; it++)
    {
        for(int ia = 0; ia < p_ucell_->atoms[it].na; ia++)
        {
            std::vector<ModuleBase::Vector3<int>> n1n2n3 = scan_supercell_for_atom(it, ia);
            n1n2n3_overall.insert(n1n2n3_overall.end(), n1n2n3.begin(), n1n2n3.end());
        }
    }
    // delete duplicates
    eliminate_duplicate_vector3<int>(n1n2n3_overall);

    supercells_ = n1n2n3_overall;
    nR_ = supercells_.size();
}

ModuleBase::Vector3<double> toQO::cal_two_center_vector(ModuleBase::Vector3<double> rij,
                                                        ModuleBase::Vector3<int> R)
{
    ModuleBase::Vector3<double> Rij;
    Rij.x = rij.x + R.x * p_ucell_->a1.x 
                  + R.y * p_ucell_->a2.x 
                  + R.z * p_ucell_->a3.x;
    Rij.y = rij.y + R.x * p_ucell_->a1.y 
                  + R.y * p_ucell_->a2.y 
                  + R.z * p_ucell_->a3.y;
    Rij.z = rij.z + R.x * p_ucell_->a1.z 
                  + R.y * p_ucell_->a2.z 
                  + R.z * p_ucell_->a3.z;
    return Rij;
}

void toQO::fold_ovlp_R(int ik)
{
    if(save_mem_) // exception handling
    {
        ModuleBase::WARNING_QUIT("toQO::folding_ovlp_R", "save_mem_ is true, which means ovlp_R_ has only one S(R).");
    }
    int nrow = ovlp_R_[0].size();
    int ncol = ovlp_R_[0][0].size();
    // clean up
    zero_out_ovlps(false);
    // calculate
    for(int iR = 0; iR < supercells_.size(); iR++)
    {
        ModuleBase::Vector3<double> R(double(supercells_[iR].x), double(supercells_[iR].y), double(supercells_[iR].z));
        double arg = (kvecs_d_[ik] * R) * ModuleBase::TWO_PI;
        double sinp, cosp;
        ModuleBase::libm::sincos(arg, &sinp, &cosp);
        std::complex<double> phase = std::complex<double>(cosp, sinp);
        for(int i = 0; i < nrow; i++)
        {
            for(int j = 0; j < ncol; j++)
            {
                ovlp_k_[i][j] += phase * ovlp_R_[iR][i][j];
            }
        }
    }
}

void toQO::append_ovlp_R_eiRk(int ik, int iR)
{
    // save memory mode, so ovlp_R_ has only one S(R)
    if(ovlp_R_.size() > 1)
    {
        ModuleBase::WARNING_QUIT("toQO::append_ovlp_R_eiRk", "ovlp_R_ has more than one S(R).");
    }
    int nrow = ovlp_R_[0].size();
    int ncol = ovlp_R_[0][0].size();

    // calculate
    ModuleBase::Vector3<double> R(double(supercells_[iR].x), double(supercells_[iR].y), double(supercells_[iR].z));
    double arg = (kvecs_d_[ik] * R) * ModuleBase::TWO_PI;
    double sinp, cosp;
    ModuleBase::libm::sincos(arg, &sinp, &cosp);
    std::complex<double> phase = std::complex<double>(cosp, sinp);
    for(int i = 0; i < nrow; i++)
    {
        for(int j = 0; j < ncol; j++)
        {
            ovlp_k_[i][j] += phase * ovlp_R_[0][i][j];
        }
    }
}

template <typename T>
void toQO::write_ovlp(const std::vector<std::vector<T>> &ovlp, const int& ik)
{
    std::string filename = "QO_ovlp_" + std::to_string(ik) + ".dat";
    std::ofstream ofs(GlobalV::global_out_dir + filename);
    if(!ofs.is_open())
    {
        ModuleBase::WARNING_QUIT("toQO::write_ovlp", "can not open file: " + filename);
    }
    // first write kpoint coordinate
    ofs << "KPOINT_COORDINATE: " << std::setw(22) << std::setprecision(14) << std::right << std::scientific << kvecs_d_[ik].x << " "
                                 << std::setw(22) << std::setprecision(14) << std::right << std::scientific << kvecs_d_[ik].y << " "
                                 << std::setw(22) << std::setprecision(14) << std::right << std::scientific << kvecs_d_[ik].z << std::endl;
    for(int i = 0; i < ovlp.size(); i++)
    {
        for(int j = 0; j < ovlp[i].size(); j++)
        {
            ofs << std::setw(22) << std::setprecision(14) << std::right << std::scientific << ovlp[i][j] << " ";
        }
        ofs << std::endl;
    }
    ofs.close();
}
template void toQO::write_ovlp<double>(const std::vector<std::vector<double>>& ovlp, const int& ik);
template void toQO::write_ovlp<std::complex<double>>(const std::vector<std::vector<std::complex<double>>>& ovlp, const int& ik);

void toQO::write_supercells()
{
    std::ofstream ofs(GlobalV::global_out_dir + "QO_supercells.dat");
    if(!ofs.is_open())
    {
        ModuleBase::WARNING_QUIT("toQO::write_supercells", "can not open file: QO_supercells.dat");
    }
    for(int i = 0; i < supercells_.size(); i++)
    {
        ofs << std::setw(5) << std::right << supercells_[i].x << " "
            << std::setw(5) << std::right << supercells_[i].y << " "
            << std::setw(5) << std::right << supercells_[i].z << std::endl;
    }
    ofs.close();
}