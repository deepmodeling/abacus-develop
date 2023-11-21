#include "module_io/to_qo.h"

void toQO::eliminate_duplicate_vector3(std::vector<ModuleBase::Vector3<int>> &v)
{
    std::vector<std::vector<int>> v_;
    for(int i = 0; i < v.size(); i++)
    {
        v_.push_back(std::vector<int>{v[i].x, v[i].y, v[i].z});
    }
    std::sort(v_.begin(), v_.end());
    v_.erase(std::unique(v_.begin(), v_.end()), v_.end());
    v.clear();
    v.resize(v_.size());
    for(int i = 0; i < v_.size(); i++)
    {
        v[i] = ModuleBase::Vector3<int>(v_[i][0], v_[i][1], v_[i][2]);
    }
}

void toQO::allocate_ovlps()
{
    if(nchi_ == 0 || nphi_ == 0)
    {
        ModuleBase::WARNING_QUIT("toQO::allocate_ovlps", "nchi_ or nphi_ is zero, which means not properly initialized.");
    }
    // deallocate memory for ovlp_ao_nao_R_ and ovlp_ao_nao_k_
    for(int iR = 0; iR < ovlp_ao_nao_R_.size(); iR++)
    {
        for(int i = 0; i < ovlp_ao_nao_R_[iR].size(); i++)
        {
            ovlp_ao_nao_R_[iR][i].clear();
        }
        ovlp_ao_nao_R_[iR].clear();
    }
    ovlp_ao_nao_R_.clear();
    for(int ik = 0; ik < ovlp_ao_nao_k_.size(); ik++)
    {
        for(int i = 0; i < ovlp_ao_nao_k_[ik].size(); i++)
        {
            ovlp_ao_nao_k_[ik][i].clear();
        }
        ovlp_ao_nao_k_[ik].clear();
    }
    ovlp_ao_nao_k_.clear();
    // allocate memory for ovlp_ao_nao_R_
    for(int iR = 0; iR < nR_; iR++)
    {
        std::vector<std::vector<std::complex<double>>> matrix;
        for(int i = 0; i < nchi_; i++)
        {
            std::vector<std::complex<double>> row;
            for(int j = 0; j < nphi_; j++)
            {
                row.push_back(std::complex<double>(0.0, 0.0));
            }
            matrix.push_back(row);
        }
        ovlp_ao_nao_R_.push_back(matrix);
    }
    // allocate memory for ovlp_ao_nao_k_
    for(int ik = 0; ik < nkpts_; ik++)
    {
        std::vector<std::vector<std::complex<double>>> matrix;
        for(int i = 0; i < nchi_; i++)
        {
            std::vector<std::complex<double>> row;
            for(int j = 0; j < nphi_; j++)
            {
                row.push_back(std::complex<double>(0.0, 0.0));
            }
            matrix.push_back(row);
        }
        ovlp_ao_nao_k_.push_back(matrix);
    }
}

void toQO::zero_out_ovlps(const bool is_R)
{
    if(nchi_ == 0 || nphi_ == 0)
    {
        ModuleBase::WARNING_QUIT("toQO::zero_out_ovlps", "nchi_ or nphi_ is zero, which means not properly initialized.");
    }
    if(is_R)
    {
        for(int iR = 0; iR < supercells_.size(); iR++)
        {
            for(int i = 0; i < nchi_; i++)
            {
                for(int j = 0; j < nphi_; j++)
                {
                    ovlp_ao_nao_R_[iR][i][j] = std::complex<double>(0.0, 0.0);
                }
            }
        }
    }
    else
    {
        for(int ik = 0; ik < nkpts_; ik++)
        {
            for(int i = 0; i < nchi_; i++)
            {
                for(int j = 0; j < nphi_; j++)
                {
                    ovlp_ao_nao_k_[ik][i][j] = std::complex<double>(0.0, 0.0);
                }
            }
        }
    }
}

std::vector<ModuleBase::Vector3<int>> toQO::scan_supercell_for_atom(int it, int ia, int start_it, int start_ia)
{
    std::vector<ModuleBase::Vector3<int>> n1n2n3;
    // cutoff radius of numerical atomic orbital of atom itia
    double rcut_i = 10.0; //WARNING! one should get rcut_i of AO here
    // lattice vectors
    for(int itype = start_it; itype < p_ucell_->ntype; itype++)
    {
        for(int iatom = start_ia; iatom < p_ucell_->atoms[itype].na; iatom++)
        {
            double rcut_j = 10.0; //WARNING! one should get rcut_j of NAO here
            ModuleBase::Vector3<double> rij = p_ucell_->atoms[itype].tau[iatom] - p_ucell_->atoms[it].tau[ia];
            int n1 = 0; int n2 = 0; int n3 = 0;
            // calculate the sup of n1, n2, n3
            int n1max = (rcut_i + rcut_j)/p_ucell_->a1.norm();
            int n2max = (rcut_i + rcut_j)/p_ucell_->a2.norm();
            int n3max = (rcut_i + rcut_j)/p_ucell_->a3.norm();
            // scan n1, n2, n3
            for(int n1 = -n1max; n1 <= n1max; n1++)
            {
                for(int n2 = -n2max; n2 <= n2max; n2++)
                {
                    for(int n3 = -n3max; n3 <= n3max; n3++)
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
    eliminate_duplicate_vector3(n1n2n3);
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
    eliminate_duplicate_vector3(n1n2n3_overall);
    supercells_ = n1n2n3_overall;
}

std::vector<std::vector<std::complex<double>>> toQO::folding_ovlp_R(ModuleBase::Vector3<double> kvec_c)
{
    std::vector<std::vector<std::complex<double>>> matrices_k;
    int nrow = ovlp_ao_nao_R_[0].size();
    int ncol = ovlp_ao_nao_R_[0][0].size();
    // initialize matrices_k with zeros
    for(int i = 0; i < nrow; i++)
    {
        std::vector<std::complex<double>> row;
        for(int j = 0; j < ncol; j++)
        {
            row.push_back(std::complex<double>(0.0, 0.0));
        }
        matrices_k.push_back(row);
    }
    for(int iR = 0; iR < supercells_.size(); iR++)
    {
        ModuleBase::Vector3<double> R = p_ucell_->a1 * double(supercells_[iR][0]) + p_ucell_->a2 * double(supercells_[iR][1]) + p_ucell_->a3 * double(supercells_[iR][2]);
        double arg = kvec_c * R;
        std::complex<double> phase = std::exp(std::complex<double>(0.0, 1.0) * arg);
        for(int i = 0; i < nrow; i++)
        {
            for(int j = 0; j < ncol; j++)
            {
                matrices_k[i][j] += phase * ovlp_ao_nao_R_[iR][i][j];
            }
        }
    }
    return matrices_k;
}
