#include "module_io/to_qo.h"
#include "module_base/libm/libm.h"

void toQO::unwrap_unitcell(UnitCell* p_ucell)
{
    p_ucell_ = p_ucell;
    ntype_ = p_ucell->ntype;
    std::for_each(p_ucell->atoms, p_ucell->atoms + p_ucell->ntype, [this](Atom& atom){
        symbols_.push_back(atom.ncpp.psd);
        na_.push_back(atom.na);
    });
    nmax_.resize(ntype_);
    charges_.resize(ntype_);
    for(int itype = 0; itype < ntype_; itype++)
    {
        std::cout << "type " << itype << " " << symbols_[itype] << " strategy: " << strategies_[itype] << std::endl;
        nmax_[itype] = (strategies_[itype].substr(0, 6) != "energy")? atom_database_.principle_quantum_number[symbols_[itype]]: atom_database_.atom_Z[symbols_[itype]];
        charges_[itype] = atom_database_.atom_Z[symbols_[itype]];
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

void toQO::allocate_ovlp(const bool& is_R)
{
    if(is_R) ovlpR_.resize(nchi_ * nphi_, 0.0);
    else ovlpk_.resize(nchi_ * nphi_, std::complex<double>(0.0, 0.0));
}

void toQO::deallocate_ovlp(const bool& is_R)
{
    if(is_R) ovlpR_.clear();
    else ovlpk_.clear();
}

void toQO::zero_out_ovlps(const bool& is_R)
{
    if(is_R) std::fill(ovlpR_.begin(), ovlpR_.end(), 0.0);
    else std::fill(ovlpk_.begin(), ovlpk_.end(), std::complex<double>(0.0, 0.0));
}

std::vector<ModuleBase::Vector3<int>> toQO::scan_supercell_for_atom(int it, int ia, int start_it, int start_ia)
{
    std::vector<ModuleBase::Vector3<int>> n1n2n3;
    // cutoff radius of numerical atomic orbital of atom itia
    double rcut_i = ao_->rcut_max(it);
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
            double rcut_j = nao_->rcut_max(itype);
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

void toQO::append_ovlpR_eiRk(int ik, int iR)
{
    // calculate
    ModuleBase::Vector3<double> R(double(supercells_[iR].x), double(supercells_[iR].y), double(supercells_[iR].z));
    double arg = (kvecs_d_[ik] * R) * ModuleBase::TWO_PI;
    double sinp, cosp;
    ModuleBase::libm::sincos(arg, &sinp, &cosp);
    std::complex<double> phase = std::complex<double>(cosp, sinp);
    std::transform(ovlpR_.begin(), ovlpR_.end(), ovlpk_.begin(), ovlpk_.begin(),
        [&](const double& ovlpR, const std::complex<double>& ovlpk) {
            return ovlpR * phase + ovlpk;
        });
}

// template function definition
// 20240310 make compatible with R space matrices
template <typename T>
void toQO::write_ovlp(const std::string& dir,
                      const std::vector<T> &ovlp, 
                      const int& nrows,
                      const int& ncols,
                      const bool& is_R,
                      const int& i)
{
    std::string filename = is_R? "QO_ovlpR_" + std::to_string(i) + ".dat": "QO_ovlp_" + std::to_string(i) + ".dat";
    std::ofstream ofs(dir + filename);
    if(!ofs.is_open())
    {
        ModuleBase::WARNING_QUIT("toQO::write_ovlp", "can not open file: " + filename);
    }
    if(is_R)
    {
        ofs << "SUPERCELL_COORDINATE: " << std::setw(5) << std::right << supercells_[i].x << " "
                                        << std::setw(5) << std::right << supercells_[i].y << " "
                                        << std::setw(5) << std::right << supercells_[i].z << std::endl;
    }
    else
    {
        ofs << "KPOINT_COORDINATE: " << std::setw(22) << std::setprecision(14) << std::right << std::scientific << kvecs_d_[i].x << " "
                                     << std::setw(22) << std::setprecision(14) << std::right << std::scientific << kvecs_d_[i].y << " "
                                     << std::setw(22) << std::setprecision(14) << std::right << std::scientific << kvecs_d_[i].z << std::endl;
    }
    for(int irow = 0; irow < nrows; irow++)
    {
        for(int icol = 0; icol < ncols; icol++)
        {
            ofs << std::setw(22) << std::setprecision(14) << std::right << std::scientific << ovlp[irow*nrows+icol] << " ";
        }
        ofs << std::endl;
    }
    ofs.close();
}
// explicit instantiation
template void toQO::write_ovlp<double>(const std::string& dir,
                                       const std::vector<double>& ovlp, 
                                       const int& nrows,
                                       const int& ncols,
                                       const bool& is_R,
                                       const int& ik);
template void toQO::write_ovlp<std::complex<double>>(const std::string& dir,
                                                     const std::vector<std::complex<double>>& ovlp, 
                                                     const int& nrows,
                                                     const int& ncols,
                                                     const bool& is_R,
                                                     const int& ik);
// a free function to convert string storing C++ std::complex to std::complex
// format: (real,imag), both part in scientific format
std::complex<double> str2complex(const std::string& str)
{
    std::string real_str, imag_str;
    int i = 1; // skip '('
    while(str[i] != ',') real_str += str[i]; i++;
    i++; // skip ','
    while(str[i] != ')') imag_str += str[i]; i++;
    return std::complex<double>(std::stod(real_str), std::stod(imag_str));
}
// complete I/O of QO module
void toQO::read_ovlp(const std::string& dir,
                     const int& nrows,
                     const int& ncols,
                     const bool& is_R,
                     const int& ik)
{
    std::string filename = is_R? "QO_ovlpR_" + std::to_string(ik) + ".dat": "QO_ovlp_" + std::to_string(ik) + ".dat";
    std::ifstream ifs(dir + filename);
    if(!ifs.is_open())
    {
        ModuleBase::WARNING_QUIT("toQO::read_ovlp", "can not open file: " + filename);
    }
    // read header
    std::string line;
    std::getline(ifs, line);
    // read ovlp
    int inum = 0;
    while(ifs.good())
    {
        std::string anum;
        ifs >> anum;
        if(is_R) ovlpR_[inum] = std::stod(anum);
        else ovlpk_[inum] = str2complex(anum);
        inum++;
    }
}

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

#include "../module_base/parallel_common.h"
void toQO::mpi_plan()
{
    
}

void toQO::radialcollection_indexing(const RadialCollection& radcol,
                                     const std::vector<int>& natoms,
                                     std::map<std::tuple<int,int,int,int,int>,int>& index_map,
                                     std::map<int,std::tuple<int,int,int,int,int>>& index_map_reverse)
{
    // in RadialCollection, radials are stored type by type and actually not relevant with exact atom index,
    // so the number of atom of each type is external information.
    // the map should be like: (itype, iatom, l, izeta, m) -> index and the reverse one is index -> (itype, iatom, l, izeta, m)
    int index = 0;
    for(int itype = 0; itype < radcol.ntype(); itype++)
    {
        for(int iatom = 0; iatom < natoms[itype]; iatom++)
        {
            for(int l = 0; l <= radcol.lmax(itype); l++)
            {
                // here should be an orbital_filter operation
                // temporary choice
                if((!orbital_filter(l, strategies_[itype]))&&(qo_basis_ == "pswfc")) continue;
                std::vector<int> ms;
                for(int m_abs = 0; m_abs <= l; m_abs++)
                {
                    ms.push_back(m_abs);
                    if(m_abs != 0) ms.push_back(-m_abs);
                }
                for(int izeta = 0; izeta < radcol.nzeta(itype, l); izeta++)
                {
                    for(int m: ms)
                    {
                        index_map[std::make_tuple(itype, iatom, l, izeta, m)] = index;
                        index_map_reverse[index] = std::make_tuple(itype, iatom, l, izeta, m);
                        index++;
                    }
                }
            }
        }
    }
}
