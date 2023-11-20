#include "module_io/to_qo.h"
#include "module_base/vector3.h"

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
    std::sort(n1n2n3_overall.begin(), n1n2n3_overall.end());
    n1n2n3_overall.erase(std::unique(n1n2n3_overall.begin(), n1n2n3_overall.end()), n1n2n3_overall.end());
    supercells_ = n1n2n3_overall;
}

Matrix toQO::folding_matrixR(ModuleBase::Vector3<double> kvec_c)
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

using Matrix = std::vector<std::vector<std::complex<double>>>;