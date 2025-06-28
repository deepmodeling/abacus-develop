#include <map>
#include "gint_dvlocal.h"
#include "phi_operator.h"
#include "source_base/parallel_reduce.h"

namespace ModuleGint
{

void Gint_dvlocal::cal_dvlocal()
{
    ModuleBase::TITLE("Gint", "cal_gint_dvlocal");
    ModuleBase::timer::tick("Gint", "cal_gint_dvlocal");
    init_hr_gint_();
    cal_hr_gint_();
    ModuleBase::timer::tick("Gint", "cal_gint_dvlocal");
}

void Gint_dvlocal::init_hr_gint_()
{
    pvdpRx = gint_info_->get_hr<double>();
    pvdpRy = gint_info_->get_hr<double>();
    pvdpRz = gint_info_->get_hr<double>();
}

void Gint_dvlocal::cal_hr_gint_()
{
#pragma omp parallel
    {
        PhiOperator phi_op;
        std::vector<double> phi;
        std::vector<double> phi_vldr3;
        std::vector<double> dphi_x;
        std::vector<double> dphi_y;
        std::vector<double> dphi_z;
#pragma omp for schedule(dynamic)
        for(const auto& biggrid: gint_info_->get_biggrids())
        {
            if(biggrid->get_atoms().empty())
            {
                continue;
            }
            phi_op.set_bgrid(biggrid);
            const int phi_len = phi_op.get_rows() * phi_op.get_cols();
            phi.resize(phi_len);
            phi_vldr3.resize(phi_len);
            dphi_x.resize(phi_len);
            dphi_y.resize(phi_len);
            dphi_z.resize(phi_len);
            phi_op.set_phi_dphi(phi.data(), dphi_x.data(), dphi_y.data(), dphi_z.data());
            phi_op.phi_mul_vldr3(vr_eff_, dr3_, phi.data(), phi_vldr3.data());
            phi_op.phi_mul_phi(phi_vldr3.data(), dphi_x.data(), pvdpRx, PhiOperator::Triangular_Matrix::Upper);
            phi_op.phi_mul_phi(phi_vldr3.data(), dphi_y.data(), pvdpRy, PhiOperator::Triangular_Matrix::Upper);
            phi_op.phi_mul_phi(phi_vldr3.data(), dphi_z.data(), pvdpRz, PhiOperator::Triangular_Matrix::Upper);
        }
    }
}

void Gint_dvlocal::cal_dvlocal_R_sparseMatrix(
    const int nspin,
    const int cspin,
    const int nlocal,
    const double sparse_thr, 
    const Parallel_Orbitals& pv,
    const UnitCell& ucell,
    const Grid_Driver& gdriver,
    LCAO_HS_Arrays& hs_arrays)
{
    ModuleBase::TITLE("Gint", "cal_dvlocal_R_sparseMatrix");
    ModuleBase::timer::tick("Gint", "cal_dvlocal_R_sparseMatrix");
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> pvdpRx_sparseMatrix;
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> pvdpRy_sparseMatrix;
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> pvdpRz_sparseMatrix;
    
    double temp_value_double;

    Vec3d tau1, dtau;
    for (int iap = 0; iap < pvdpRx.size_atom_pairs(); iap++)
    {
        const auto& ap = pvdpRx.get_atom_pair(iap);
        const int iat1 = ap.get_atom_i();
        const int iat2 = ap.get_atom_j();
        const int it1 = ucell.iat2it[iat1];
        const int it2 = ucell.iat2it[iat2];
        const Atom* atom1 = &ucell.atoms[it1];
        const Atom* atom2 = &ucell.atoms[it2];
        const int start1 = ucell.itiaiw2iwt(it1, ucell.iat2ia[iat1], 0);
        const int start2 = ucell.itiaiw2iwt(it2, ucell.iat2ia[iat2], 0);

        for (int ir = 0; ir < ap.get_R_size(); ir++)
        {
            const ModuleBase::Vector3<int> R = ap.get_R_index(ir);
            Abfs::Vector3_Order<int> dR(R.x, R.y, R.z);
            double* p_pvdpRx = pvdpRx.get_atom_pair(iap).get_pointer(ir);
            double* p_pvdpRy = pvdpRy.get_atom_pair(iap).get_pointer(ir);
            double* p_pvdpRz = pvdpRz.get_atom_pair(iap).get_pointer(ir);

            for (int iw = 0; iw < atom1->nw * npol_; iw++)
            {
                for (int iw2 = 0; iw2 < atom2->nw * npol_; iw2++)
                {
                    const int nw = atom2->nw;
                    const int mug0 = iw / npol_;
                    const int nug0 = iw2 / npol_;
                    const int iw_nowg = mug0 * nw + nug0;
                    
                    double temp_value = p_pvdpRx[iw_nowg];
                    if (std::abs(temp_value) > sparse_thr)
                    {
                        pvdpRx_sparseMatrix[dR][start1 + iw][start2 + iw2] = temp_value;
                    }
                    temp_value = p_pvdpRy[iw_nowg];
                    if (std::abs(temp_value) > sparse_thr)
                    {
                        pvdpRy_sparseMatrix[dR][start1 + iw][start2 + iw2] = temp_value;
                    }
                    temp_value = p_pvdpRz[iw_nowg];
                    if (std::abs(temp_value) > sparse_thr)
                    {
                        pvdpRz_sparseMatrix[dR][start1 + iw][start2 + iw2] = temp_value;
                    }
                }
            }
        }
    }
    distribute_pvdpR_sparseMatrix(cspin, 0, nlocal, sparse_thr, pvdpRx_sparseMatrix, pv, hs_arrays);
    distribute_pvdpR_sparseMatrix(cspin, 1, nlocal, sparse_thr, pvdpRy_sparseMatrix, pv, hs_arrays);
    distribute_pvdpR_sparseMatrix(cspin, 2, nlocal, sparse_thr, pvdpRz_sparseMatrix, pv, hs_arrays);
    ModuleBase::timer::tick("Gint", "cal_dvlocal_R_sparseMatrix");
}


void Gint_dvlocal::distribute_pvdpR_sparseMatrix(
    const int cspin,
    const int dim,
    const int nlocal,
    const double sparse_threshold,
    const std::map<Abfs::Vector3_Order<int>,
                    std::map<size_t, std::map<size_t, double>>>&
        pvdpR_sparseMatrix,
    const Parallel_Orbitals& pv,
    LCAO_HS_Arrays& hs_arrays)
{
    int total_R_num = hs_arrays.all_R_coor.size();
    std::vector<int> nonzero_num(total_R_num);
    std::vector<int> minus_nonzero_num(total_R_num);
    int count = 0;
    for (const auto& R_coor: hs_arrays.all_R_coor)
    {
        auto iter = pvdpR_sparseMatrix.find(R_coor);
        if (iter != pvdpR_sparseMatrix.end())
        {
            for (auto& row_loop: iter->second)
            {
                nonzero_num[count] += row_loop.second.size();
            }
        }

        auto minus_R_coor = -1 * R_coor;

        iter = pvdpR_sparseMatrix.find(minus_R_coor);
        if (iter != pvdpR_sparseMatrix.end())
        {
            for (auto& row_loop: iter->second)
            {
                minus_nonzero_num[count] += row_loop.second.size();
            }
        }
        count++;
    }

    Parallel_Reduce::reduce_all(nonzero_num.data(), total_R_num);
    Parallel_Reduce::reduce_all(minus_nonzero_num.data(), total_R_num);

    std::vector<double> tmp(nlocal);
    count = 0;

    const std::vector<int>& trace_lo = gint_info_->get_trace_lo();
    for (const auto& R_coor: hs_arrays.all_R_coor)
    {
        if (nonzero_num[count] != 0 || minus_nonzero_num[count] != 0)
        {
            auto minus_R_coor = -1 * R_coor;

            for (int row = 0; row < nlocal; ++row)
            {
                tmp.assign(tmp.size(), 0);

                auto iter = pvdpR_sparseMatrix.find(R_coor);
                if (iter != pvdpR_sparseMatrix.end())
                {

                    if (trace_lo[row] >= 0)
                    {
                        auto row_iter = iter->second.find(row);
                        if (row_iter != iter->second.end())
                        {
                            for (auto& value: row_iter->second)
                            {
                                tmp[value.first] = value.second;
                            }
                        }
                    }
                }

                auto minus_R_iter = pvdpR_sparseMatrix.find(minus_R_coor);
                if (minus_R_iter != pvdpR_sparseMatrix.end())
                {
                    for (int col = 0; col < row; ++col)
                    {
                        if (trace_lo[col] >= 0)
                        {
                            auto row_iter = minus_R_iter->second.find(col);
                            if (row_iter != minus_R_iter->second.end())
                            {
                                auto col_iter = row_iter->second.find(row);
                                if (col_iter != row_iter->second.end())
                                {
                                    tmp[col] = col_iter->second;
                                }
                            }
                        }
                    }
                }

                Parallel_Reduce::reduce_pool(tmp.data(), nlocal);

                if (pv.global2local_row(row) >= 0)
                {
                    for (int col = 0; col < nlocal; ++col)
                    {
                        if (pv.global2local_col(col) >= 0)
                        {
                            if (std::abs(tmp[col]) > sparse_threshold)
                            {
                                if (dim == 0)
                                {
                                    double& value = hs_arrays.dHRx_sparse[cspin][R_coor][row][col];
                                    value += tmp[col];
                                    if (std::abs(value) <= sparse_threshold)
                                    {
                                        hs_arrays.dHRx_sparse[cspin][R_coor][row].erase(col);
                                    }
                                }
                                if (dim == 1)
                                {
                                    double& value = hs_arrays.dHRy_sparse[cspin][R_coor][row][col];
                                    value += tmp[col];
                                    if (std::abs(value) <= sparse_threshold)
                                    {
                                        hs_arrays.dHRy_sparse[cspin][R_coor][row].erase(col);
                                    }
                                }
                                if (dim == 2)
                                {
                                    double& value = hs_arrays.dHRz_sparse[cspin][R_coor][row][col];
                                    value += tmp[col];
                                    if (std::abs(value) <= sparse_threshold)
                                    {
                                        hs_arrays.dHRz_sparse[cspin][R_coor][row].erase(col);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        count++;
    }
}
}