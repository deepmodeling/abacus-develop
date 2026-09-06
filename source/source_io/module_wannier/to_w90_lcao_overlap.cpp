#include "to_w90_lcao.h"

#include "source_lcao/module_ri/abfs_vector3_order.h"
#include "source_io/module_parameter/parameter.h"
#include "fr_overlap.h"
#include "source_base/math_integral.h"
#include "source_base/math_polyint.h"
#include "source_base/math_sphbes.h"
#include "source_base/math_ylmreal.h"
#include "source_base/parallel_reduce.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_hamilt/module_hcontainer/atom_pair.h"

#include <fstream>
#include <functional>

#ifdef __LCAO

void toW90_LCAO::cal_orbA_overlap_R(const UnitCell& ucell)
{
    int row = this->ParaV->get_row_size();
    int R_num = R_coor_car.size();
    int global_ir = 0;

    psi_psiA_R.resize(row);
    for (int ir = 0; ir < row; ir++)
    {
        psi_psiA_R[ir].resize(num_wannier);

        for (int ic = 0; ic < num_wannier; ++ic)
        {
            psi_psiA_R[ir][ic].resize(R_num);
        }
    }

    double bs2, bs3, bs6, bs12;
    bs2 = 1.0 / sqrt(2.0);
    bs3 = 1.0 / sqrt(3.0);
    bs6 = 1.0 / sqrt(6.0);
    bs12 = 1.0 / sqrt(12.0);

    for (int ir = 0; ir < row; ir++)
    {
        global_ir = ParaV->local2global_row(ir);
        int orb_index_row = global_ir / npol_;

        int it1 = iw2it[orb_index_row];
        int ia1 = iw2ia[orb_index_row];
        int im1 = iw2im[orb_index_row];

        for (int wannier_index = 0; wannier_index < num_wannier; wannier_index++)
        {
            if (L[wannier_index] >= 0)
            {
                for (int iR = 0; iR < R_num; iR++)
                {
                    ModuleBase::Vector3<double> R_car = R_coor_car[iR];
                    ModuleBase::Vector3<double> orb_center
                        = (ucell.atoms[it1].tau[ia1] + R_car) * ucell.lat0;
                    ModuleBase::Vector3<double> project_orb_center = R_centre[wannier_index] * ucell.lat0;

                    double overlap_o
                        = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(0).cal_overlap(orb_center,
                                                                                                   project_orb_center,
                                                                                                   im1,
                                                                                                   m[wannier_index]);
                    psi_psiA_R[ir][wannier_index][iR] = overlap_o;
                }
            }
            else
            {
                if (L[wannier_index] == -1)
                {
                    double tmp_bs2 = 0;
                    if (m[wannier_index] == 0) {
                        tmp_bs2 = bs2;
}
                    if (m[wannier_index] == -1) {
                        tmp_bs2 = -bs2;
}

                    for (int iR = 0; iR < R_num; iR++)
                    {
                        ModuleBase::Vector3<double> R_car = R_coor_car[iR];
                        ModuleBase::Vector3<double> orb_center
                            = (ucell.atoms[it1].tau[ia1] + R_car) * ucell.lat0;
                        ModuleBase::Vector3<double> project_orb_center = R_centre[wannier_index] * ucell.lat0;

                        double overlap_s = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(0).cal_overlap(
                            orb_center,
                            project_orb_center,
                            im1,
                            0);
                        double overlap_px = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(1).cal_overlap(
                            orb_center,
                            project_orb_center,
                            im1,
                            1);
                        psi_psiA_R[ir][wannier_index][iR] = bs2 * overlap_s + tmp_bs2 * overlap_px;
                    }
                }

                if (L[wannier_index] == -2)
                {
                    if (m[wannier_index] == 0 || m[wannier_index] == 1)
                    {
                        double tmp_bs2 = bs2;
                        if (m[wannier_index] == -1) {
                            tmp_bs2 = -bs2;
}

                        for (int iR = 0; iR < R_num; iR++)
                        {
                            ModuleBase::Vector3<double> R_car = R_coor_car[iR];
                            ModuleBase::Vector3<double> orb_center
                                = (ucell.atoms[it1].tau[ia1] + R_car) * ucell.lat0;
                            ModuleBase::Vector3<double> project_orb_center
                                = R_centre[wannier_index] * ucell.lat0;

                            double overlap_s = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(0).cal_overlap(
                                orb_center,
                                project_orb_center,
                                im1,
                                0);
                            double overlap_px = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index]
                                                    .at(1)
                                                    .cal_overlap(orb_center, project_orb_center, im1, 1);
                            double overlap_py = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index]
                                                    .at(1)
                                                    .cal_overlap(orb_center, project_orb_center, im1, 2);
                            psi_psiA_R[ir][wannier_index][iR]
                                = bs3 * overlap_s - bs6 * overlap_px + tmp_bs2 * overlap_py;
                        }
                    }
                    else if (m[wannier_index] == 2)
                    {
                        for (int iR = 0; iR < R_num; iR++)
                        {
                            ModuleBase::Vector3<double> R_car = R_coor_car[iR];
                            ModuleBase::Vector3<double> orb_center
                                = (ucell.atoms[it1].tau[ia1] + R_car) * ucell.lat0;
                            ModuleBase::Vector3<double> project_orb_center
                                = R_centre[wannier_index] * ucell.lat0;

                            double overlap_s = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(0).cal_overlap(
                                orb_center,
                                project_orb_center,
                                im1,
                                0);
                            double overlap_px = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index]
                                                    .at(1)
                                                    .cal_overlap(orb_center, project_orb_center, im1, 1);
                            psi_psiA_R[ir][wannier_index][iR] = bs3 * overlap_s + 2.0 * bs6 * overlap_px;
                        }
                    }
                }

                if (L[wannier_index] == -3)
                {
                    double m_px = 1.0;
                    double m_py = 1.0;
                    double m_pz = 1.0;

                    if (m[wannier_index] == 1)
                    {
                        m_py = -1.0;
                        m_pz = -1.0;
                    }
                    else if (m[wannier_index] == 2)
                    {
                        m_px = -1.0;
                        m_pz = -1.0;
                    }
                    else if (m[wannier_index] == 3)
                    {
                        m_px = -1.0;
                        m_py = -1.0;
                    }

                    for (int iR = 0; iR < R_num; iR++)
                    {
                        ModuleBase::Vector3<double> R_car = R_coor_car[iR];
                        ModuleBase::Vector3<double> orb_center
                            = (ucell.atoms[it1].tau[ia1] + R_car) * ucell.lat0;
                        ModuleBase::Vector3<double> project_orb_center = R_centre[wannier_index] * ucell.lat0;

                        double overlap_s = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(0).cal_overlap(
                            orb_center,
                            project_orb_center,
                            im1,
                            0);
                        double overlap_px = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(1).cal_overlap(
                            orb_center,
                            project_orb_center,
                            im1,
                            1);
                        double overlap_py = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(1).cal_overlap(
                            orb_center,
                            project_orb_center,
                            im1,
                            2);
                        double overlap_pz = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(1).cal_overlap(
                            orb_center,
                            project_orb_center,
                            im1,
                            0);
                        psi_psiA_R[ir][wannier_index][iR]
                            = 0.5 * (overlap_s + m_px * overlap_px + m_py * overlap_py + m_pz * overlap_pz);
                    }
                }

                if (L[wannier_index] == -4)
                {
                    if (m[wannier_index] == 0 || m[wannier_index] == 1)
                    {
                        double tmp_bs2 = bs2;
                        if (m[wannier_index] == -1) {
                            tmp_bs2 = -bs2;
}

                        for (int iR = 0; iR < R_num; iR++)
                        {
                            ModuleBase::Vector3<double> R_car = R_coor_car[iR];
                            ModuleBase::Vector3<double> orb_center
                                = (ucell.atoms[it1].tau[ia1] + R_car) * ucell.lat0;
                            ModuleBase::Vector3<double> project_orb_center
                                = R_centre[wannier_index] * ucell.lat0;

                            double overlap_s = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(0).cal_overlap(
                                orb_center,
                                project_orb_center,
                                im1,
                                0);
                            double overlap_px = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index]
                                                    .at(1)
                                                    .cal_overlap(orb_center, project_orb_center, im1, 1);
                            double overlap_py = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index]
                                                    .at(1)
                                                    .cal_overlap(orb_center, project_orb_center, im1, 2);
                            psi_psiA_R[ir][wannier_index][iR]
                                = bs3 * overlap_s - bs6 * overlap_px + tmp_bs2 * overlap_py;
                        }
                    }
                    else if (m[wannier_index] == 2)
                    {
                        for (int iR = 0; iR < R_num; iR++)
                        {
                            ModuleBase::Vector3<double> R_car = R_coor_car[iR];
                            ModuleBase::Vector3<double> orb_center
                                = (ucell.atoms[it1].tau[ia1] + R_car) * ucell.lat0;
                            ModuleBase::Vector3<double> project_orb_center
                                = R_centre[wannier_index] * ucell.lat0;

                            double overlap_s = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(0).cal_overlap(
                                orb_center,
                                project_orb_center,
                                im1,
                                0);
                            double overlap_px = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index]
                                                    .at(1)
                                                    .cal_overlap(orb_center, project_orb_center, im1, 1);
                            psi_psiA_R[ir][wannier_index][iR] = bs3 * overlap_s + 2.0 * bs6 * overlap_px;
                        }
                    }
                    else if (m[wannier_index] == 3 || m[wannier_index] == 4)
                    {
                        double m_pz = 1.0;
                        if (m[wannier_index] == 4) {
                            m_pz = -1.0;
}

                        for (int iR = 0; iR < R_num; iR++)
                        {
                            ModuleBase::Vector3<double> R_car = R_coor_car[iR];
                            ModuleBase::Vector3<double> orb_center
                                = (ucell.atoms[it1].tau[ia1] + R_car) * ucell.lat0;
                            ModuleBase::Vector3<double> project_orb_center
                                = R_centre[wannier_index] * ucell.lat0;

                            double overlap_pz = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index]
                                                    .at(1)
                                                    .cal_overlap(orb_center, project_orb_center, im1, 0);
                            double overlap_dz2 = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index]
                                                     .at(2)
                                                     .cal_overlap(orb_center, project_orb_center, im1, 0);
                            psi_psiA_R[ir][wannier_index][iR] = bs2 * (m_pz * overlap_pz + overlap_dz2);
                        }
                    }
                }

                if (L[wannier_index] == -5)
                {
                    if (m[wannier_index] == 0 || m[wannier_index] == 1)
                    {
                        double tmp_bs2 = -bs2;
                        double tmp_bs12 = -bs12;
                        double tmp_d = 0.5;

                        if (m[wannier_index] == 1)
                        {
                            tmp_bs2 = bs2;
                        }

                        for (int iR = 0; iR < R_num; iR++)
                        {
                            ModuleBase::Vector3<double> R_car = R_coor_car[iR];
                            ModuleBase::Vector3<double> orb_center
                                = (ucell.atoms[it1].tau[ia1] + R_car) * ucell.lat0;
                            ModuleBase::Vector3<double> project_orb_center
                                = R_centre[wannier_index] * ucell.lat0;

                            double overlap_s = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(0).cal_overlap(
                                orb_center,
                                project_orb_center,
                                im1,
                                0);
                            double overlap_px = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index]
                                                    .at(1)
                                                    .cal_overlap(orb_center, project_orb_center, im1, 1);
                            double overlap_dz2 = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index]
                                                     .at(2)
                                                     .cal_overlap(orb_center, project_orb_center, im1, 0);
                            double overlap_dx2_y2 = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index]
                                                        .at(2)
                                                        .cal_overlap(orb_center, project_orb_center, im1, 3);
                            psi_psiA_R[ir][wannier_index][iR] = bs6 * overlap_s + tmp_bs2 * overlap_px
                                                                + tmp_bs12 * overlap_dz2 + tmp_d * overlap_dx2_y2;
                        }
                    }
                    else if (m[wannier_index] == 2 || m[wannier_index] == 3)
                    {
                        double tmp_bs2 = -bs2;
                        double tmp_bs12 = -bs12;
                        double tmp_d = -0.5;

                        if (m[wannier_index] == 3)
                        {
                            tmp_bs2 = bs2;
                        }

                        for (int iR = 0; iR < R_num; iR++)
                        {
                            ModuleBase::Vector3<double> R_car = R_coor_car[iR];
                            ModuleBase::Vector3<double> orb_center
                                = (ucell.atoms[it1].tau[ia1] + R_car) * ucell.lat0;
                            ModuleBase::Vector3<double> project_orb_center
                                = R_centre[wannier_index] * ucell.lat0;

                            double overlap_s = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(0).cal_overlap(
                                orb_center,
                                project_orb_center,
                                im1,
                                0);
                            double overlap_py = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index]
                                                    .at(1)
                                                    .cal_overlap(orb_center, project_orb_center, im1, 2);
                            double overlap_dz2 = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index]
                                                     .at(2)
                                                     .cal_overlap(orb_center, project_orb_center, im1, 0);
                            double overlap_dx2_y2 = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index]
                                                        .at(2)
                                                        .cal_overlap(orb_center, project_orb_center, im1, 3);
                            psi_psiA_R[ir][wannier_index][iR] = bs6 * overlap_s + tmp_bs2 * overlap_py
                                                                + tmp_bs12 * overlap_dz2 + tmp_d * overlap_dx2_y2;
                        }
                    }
                    else if (m[wannier_index] == 4 || m[wannier_index] == 5)
                    {
                        double tmp_pz = -1.0;

                        if (m[wannier_index] == 5) {
                            tmp_pz = 1.0;
}

                        for (int iR = 0; iR < R_num; iR++)
                        {
                            ModuleBase::Vector3<double> R_car = R_coor_car[iR];
                            ModuleBase::Vector3<double> orb_center
                                = (ucell.atoms[it1].tau[ia1] + R_car) * ucell.lat0;
                            ModuleBase::Vector3<double> project_orb_center
                                = R_centre[wannier_index] * ucell.lat0;

                            double overlap_s = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(0).cal_overlap(
                                orb_center,
                                project_orb_center,
                                im1,
                                0);
                            double overlap_pz = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index]
                                                    .at(1)
                                                    .cal_overlap(orb_center, project_orb_center, im1, 0);
                            double overlap_dz2 = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index]
                                                     .at(2)
                                                     .cal_overlap(orb_center, project_orb_center, im1, 0);
                            psi_psiA_R[ir][wannier_index][iR]
                                = bs6 * overlap_s + tmp_pz * bs2 * overlap_pz + bs3 * overlap_dz2;
                        }
                    }
                }
            }
        }
    }
}

void toW90_LCAO::unkdotA(const K_Vectors& kv,
                               const psi::Psi<std::complex<double>>& psi_in,
                               const int& ik,
                               ModuleBase::ComplexMatrix& Amn)
{
    Amn.create(num_bands, num_wannier);

    int row = this->ParaV->get_row_size();
    int index_band = -1;
    int R_num = R_coor_car.size();
    if (nspin_ != 4)
    {
        for (int ib = 0; ib < nbands_; ib++)
        {
            if (exclude_bands.count(ib)) {
                continue;
}
            index_band++;

            int ic = this->ParaV->global2local_col(ib);
            if (ic >= 0)
            {
                for (int wannier_index = 0; wannier_index < num_wannier; wannier_index++)
                {
                    for (int iR = 0; iR < R_num; iR++)
                    {
                        ModuleBase::Vector3<double> R = R_coor_car[iR];
                        double kRn = -1.0 * kv.kvec_c[ik] * R * ModuleBase::TWO_PI;
                        std::complex<double> kRn_phase(cos(kRn), sin(kRn));
                        std::complex<double> tmp(0.0, 0.0);

                        for (int ir = 0; ir < row; ir++)
                        {
                            tmp += std::conj(psi_in(ik, ic, ir)) * psi_psiA_R[ir][wannier_index][iR];
                        }

                        Amn(index_band, wannier_index) += kRn_phase * tmp;
                    }
                }
            }
        }
    }
    else
    {
        for (int ib = 0; ib < nbands_; ib++)
        {
            if (exclude_bands.count(ib)) {
                continue;
}
            index_band++;

            int ic = this->ParaV->global2local_col(ib);
            if (ic >= 0)
            {
                for (int wannier_index = 0; wannier_index < num_wannier; wannier_index++)
                {
                    for (int iR = 0; iR < R_num; iR++)
                    {
                        ModuleBase::Vector3<double> R = R_coor_car[iR];
                        double kRn = kv.kvec_c[ik] * R * ModuleBase::TWO_PI;
                        std::complex<double> kRn_phase(cos(kRn), sin(kRn));
                        std::complex<double> tmp(0.0, 0.0);
                        int half_row = row / 2;
                        for (int ir = 0; ir < half_row; ir++)
                        {
                            tmp += up_con[wannier_index] * std::conj(psi_in(ik, ic, 2 * ir))
                                       * psi_psiA_R[2 * ir][wannier_index][iR]
                                   + dn_con[wannier_index] * std::conj(psi_in(ik, ic, 2 * ir + 1))
                                         * psi_psiA_R[2 * ir + 1][wannier_index][iR];
                        }

                        Amn(index_band, wannier_index) += kRn_phase * tmp;
                    }
                }
            }
        }
    }

    Parallel_Reduce::reduce_all(Amn.c, Amn.size);
}

#endif
