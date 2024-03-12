#include "module_io/to_qo.h"
#include "module_base/libm/libm.h"
#ifdef __MPI
#include "../module_base/parallel_common.h"
#endif
void toQO::read_structures(const UnitCell* p_ucell, 
                           const std::vector<ModuleBase::Vector3<double>>& kvecs_d,
                           const int& rank,
                           const int& nranks)
{
    // assume p_ucell will be totally available for all processors if MPI is enabled
    p_ucell_ = p_ucell;

    // atom related properties
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

    // k point related properties
    kvecs_d_ = kvecs_d;
    nks_ = kvecs_d.size();
    nks_tot_ = nks_;

    // scatter k points to all ranks if MPI is enabled
#ifdef __MPI
    // scatter kvecs_d_ to all ranks
    std::vector<std::vector<int>> nks_divided(nranks); // indiced by iproc, then list of indices of kvecs_d_
    if(rank == 0)
    {
        int nks = nks_tot_;
        int nks_perrank = nks / nranks;
        int nks_remain = nks % nranks;
        
        int start_ik = 0;
        for(int i = 0; i < nranks; i++)
        {
            int nks_this_rank = nks_perrank + int(i < nks_remain);
            std::vector<int> nks_this_rank_indices;
            for(int j = 0; j < nks_this_rank; j++)
            {
                nks_this_rank_indices.push_back(start_ik + j);
            }
            // for those ranks with less k points, pad with -1
            if((i >= nks_remain)&&(nks_remain > 0)) nks_this_rank_indices.push_back(-1);
            start_ik += nks_this_rank;
            nks_divided[i] = nks_this_rank_indices;
        }
    }
    for(int i = 0; i < nranks; i++)
    {
        int nks_dim;
        if(iproc_ == 0) nks_dim = nks_divided[i].size();
        Parallel_Common::bcast_int(nks_dim);
        if(iproc_ != 0) nks_divided[i].resize(nks_dim);
        Parallel_Common::bcast_int(nks_divided[i].data(), nks_dim);
    }
    //bcast_stdvector_ofvector3double(kvecs_d_); // because kvecs_d is already broadcasted in the main program
    iks_.clear();
    for(int i = 0; i < nks_divided[rank].size(); i++)
    {
        if(nks_divided[rank][i] != -1) iks_.push_back(nks_divided[rank][i]); // we want it is explicitly -1
        else iks_.push_back(-1);
    }
    nks_ = iks_.size();
#endif
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

std::vector<ModuleBase::Vector3<int>> toQO::scan_supercell_for_atom(int it, int ia, int start_it, int start_ia)
{
    std::vector<ModuleBase::Vector3<int>> n1n2n3;
    // cutoff radius of numerical atomic orbital of atom itia
    double rcut_i = ao_->rcut_max(it);
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

void toQO::scan_supercell(const int& rank, const int& nranks)
{
    if(rank == 0)
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
        nR_tot_ = nR_;
        
        write_supercells();
    }
    /*-------------------------------------------------------------------------------------------*/
#ifdef __MPI // scatter supercells_ to all ranks
    Parallel_Common::bcast_int(nR_);
    Parallel_Common::bcast_int(nR_tot_);
    bcast_stdvector_ofvector3int(supercells_);
    // scatter
    std::vector<std::vector<int>> nR_divided(nranks);  // indiced by iproc, then list of indices of supercells_
    if(rank == 0)
    {
        // divide nR into std::vector of std::vector<int>, each std::vector<int> is a chunk of indices of supercells_
        int nRs = nR_;
        int nRs_perrank = nRs / nranks;
        int nRs_remain = nRs % nranks;

        int start_iR = 0;
        for(int i = 0; i < nranks; i++)
        {
            int nR_this_rank = nRs_perrank + int(i < nRs_remain);
            std::vector<int> nR_this_rank_indices;
            for(int j = 0; j < nR_this_rank; j++)
            {
                nR_this_rank_indices.push_back(start_iR + j);
            }
            start_iR += nR_this_rank;
            nR_divided[i] = nR_this_rank_indices;
        }
    }
    for(int i = 0; i < nranks; i++)
    {
        int nR_dim;
        if(rank == 0) nR_dim = nR_divided[i].size();
        Parallel_Common::bcast_int(nR_dim);
        if(rank != 0) nR_divided[i].resize(nR_dim);
        Parallel_Common::bcast_int(nR_divided[i].data(), nR_dim);
    }
    iRs_.clear();
    for(int i = 0; i < nR_divided[rank].size(); i++) iRs_.push_back(nR_divided[rank][i]);
    nR_ = iRs_.size();
#endif
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
    // add all values of ovlpR_ to ovlpk_ with multiplication of phase
    for(int i = 0; i < nchi_ * nphi_; i++) ovlpk_[i] += ovlpR_[i] * phase;
}

void toQO::write_supercells()
{
    std::ofstream ofs(out_dir_ + "QO_supercells.dat");
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

