#include "pos_op_calc.h"

#include "source_base/timer.h"
#include "source_base/ylm.h"

PosOpCalc::PosOpCalc(const PosOpBasis& basis) : basis_(basis)
{
}

ModuleBase::Vector3<double> PosOpCalc::pos_matrix(const ModuleBase::Vector3<double>& R1,
                                                   int T1, int L1, int m1, int N1,
                                                   const ModuleBase::Vector3<double>& R2,
                                                   int T2, int L2, int m2, int N2) const
{
    ModuleBase::Vector3<double> origin_point(0.0, 0.0, 0.0);
    double factor = sqrt(ModuleBase::FOUR_PI / 3.0);
    const ModuleBase::Vector3<double>& distance = R2 - R1;

    double overlap_o = basis_.get_center2_orb11().at(T1).at(T2).at(L1).at(N1).at(L2).at(N2).cal_overlap(origin_point,
                                                                                                        distance,
                                                                                                        m1,
                                                                                                        m2);

    double overlap_x = -1 * factor * basis_.get_center2_orb21_r().at(T1).at(T2).at(L1).at(N1).at(L2).at(N2).cal_overlap(origin_point,
                                                                                                                        distance,
                                                                                                                        m1,
                                                                                                                        1,
                                                                                                                        m2); // m =  1

    double overlap_y = -1 * factor * basis_.get_center2_orb21_r().at(T1).at(T2).at(L1).at(N1).at(L2).at(N2).cal_overlap(origin_point,
                                                                                                                        distance,
                                                                                                                        m1,
                                                                                                                        2,
                                                                                                                        m2); // m = -1

    double overlap_z = factor * basis_.get_center2_orb21_r().at(T1).at(T2).at(L1).at(N1).at(L2).at(N2).cal_overlap(origin_point,
                                                                                                                  distance,
                                                                                                                  m1,
                                                                                                                  0,
                                                                                                                  m2); // m =  0

    ModuleBase::Vector3<double> temp_prp = ModuleBase::Vector3<double>(overlap_x, overlap_y, overlap_z) + R1 * overlap_o;

    return temp_prp;
}

ModuleBase::Vector3<double> PosOpCalc::pos_grad_matrix(const ModuleBase::Vector3<double>& R1,
                                                        int T1, int L1, int m1, int N1,
                                                        const ModuleBase::Vector3<double>& R2,
                                                        int T2, int L2, int m2, int N2,
                                                        const ModuleBase::Vector3<double>& Efield,
                                                        const ModuleBase::Vector3<double>& dR) const
{
    ModuleBase::Vector3<double> origin_point(0.0, 0.0, 0.0);
    double factor = sqrt(ModuleBase::FOUR_PI / 3.0);
    const ModuleBase::Vector3<double>& distance = R2 - R1;

    ModuleBase::Vector3<double> grad_o = basis_.get_center2_orb11().at(T1).at(T2).at(L1).at(N1).at(L2).at(N2).cal_grad_overlap(origin_point, distance, m1, m2);

    ModuleBase::Vector3<double> grad_rx = -1 * factor
                                          * basis_.get_center2_orb21_r().at(T1).at(T2).at(L1).at(N1).at(L2).at(N2).cal_grad_overlap(origin_point,
                                                                                                                                      distance,
                                                                                                                                      m1,
                                                                                                                                      1,
                                                                                                                                      m2); // m =  1

    ModuleBase::Vector3<double> grad_ry = -1 * factor
                                          * basis_.get_center2_orb21_r().at(T1).at(T2).at(L1).at(N1).at(L2).at(N2).cal_grad_overlap(origin_point,
                                                                                                                                      distance,
                                                                                                                                      m1,
                                                                                                                                      2,
                                                                                                                                      m2); // m = -1

    ModuleBase::Vector3<double> grad_rz = factor
                                          * basis_.get_center2_orb21_r().at(T1).at(T2).at(L1).at(N1).at(L2).at(N2).cal_grad_overlap(origin_point,
                                                                                                                                      distance,
                                                                                                                                      m1,
                                                                                                                                      0,
                                                                                                                                      m2); // m =  0

    ModuleBase::Vector3<double> temp_prp = Efield[0] * grad_rx + Efield[1] * grad_ry + Efield[2] * grad_rz + (Efield * (R1 - dR)) * grad_o;

    return temp_prp;
}

void PosOpCalc::pos_beta_matrix(const UnitCell& ucell,
                                 std::vector<std::vector<double>>& nlm,
                                 const ModuleBase::Vector3<double>& R1,
                                 int T1, int L1, int m1, int N1,
                                 const ModuleBase::Vector3<double>& R2,
                                 int T2) const
{
    ModuleBase::Vector3<double> origin_point(0.0, 0.0, 0.0);
    double factor = sqrt(ModuleBase::FOUR_PI / 3.0);
    const ModuleBase::Vector3<double>& distance = R2 - R1;
    const NonlocalInfoBase& infoNL_ = *ucell.infoNL;
    const int nproj = infoNL_.get_nproj(T2);
    nlm.resize(4);
    if (nproj == 0)
    {
        for (int i = 0; i < 4; i++)
        {
            nlm[i].resize(1);
        }
        return;
    }

    int natomwfc = 0;
    for (int ip = 0; ip < nproj; ip++)
    {
        const int L2 = infoNL_.get_proj_L(T2, ip); // mohan add 2021-05-07
        natomwfc += 2 * L2 + 1;
    }
    for (int i = 0; i < 4; i++)
    {
        nlm[i].resize(natomwfc);
    }
    int index = 0;
    for (int ip = 0; ip < nproj; ip++)
    {
        const int L2 = infoNL_.get_proj_L(T2, ip);
        for (int m2 = 0; m2 < 2 * L2 + 1; m2++)
        {
            double overlap_o = basis_.get_center2_orb11_nonlocal().at(T1).at(T2).at(L1).at(N1).at(ip).cal_overlap(origin_point, distance, m1, m2);

            double overlap_x = -1 * factor
                               * basis_.get_center2_orb21_r_nonlocal().at(T1).at(T2).at(L1).at(N1).at(ip).cal_overlap(origin_point,
                                                                                                                      distance,
                                                                                                                      m1,
                                                                                                                      1,
                                                                                                                      m2); // m =  1

            double overlap_y = -1 * factor
                               * basis_.get_center2_orb21_r_nonlocal().at(T1).at(T2).at(L1).at(N1).at(ip).cal_overlap(origin_point,
                                                                                                                      distance,
                                                                                                                      m1,
                                                                                                                      2,
                                                                                                                      m2); // m = -1

            double overlap_z = factor
                               * basis_.get_center2_orb21_r_nonlocal().at(T1).at(T2).at(L1).at(N1).at(ip).cal_overlap(origin_point,
                                                                                                                      distance,
                                                                                                                      m1,
                                                                                                                      0,
                                                                                                                      m2); // m =  0

            nlm[0][index] = overlap_o;
            nlm[1][index] = overlap_x + (R1 * overlap_o).x;
            nlm[2][index] = overlap_y + (R1 * overlap_o).y;
            nlm[3][index] = overlap_z + (R1 * overlap_o).z;
            index++;
        }
    }
}
