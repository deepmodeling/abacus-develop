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

void toW90_LCAO::initialize_orb_table(const UnitCell& ucell)
{
#ifdef __LCAO
    const int ntype = orb_.get_ntype();
    int lmax_orb = -1;
    for (int it = 0; it < ntype; it++)
    {
        lmax_orb = std::max(lmax_orb, orb_.Phi[it].getLmax());
    }
    const double dr = orb_.get_dR();
    const double dk = orb_.get_dk();
    const int kmesh = orb_.get_kmesh() * 4 + 1;
    int Rmesh = static_cast<int>(orb_.get_Rmax() / dr) + 4;
    Rmesh += 1 - Rmesh % 2;

    const int Lmax = lmax_orb + 1;
    const int Lmax_used = 2 * lmax_orb + 1;
    Center2_Orb::init_Table_Spherical_Bessel(Lmax_used,
                                             dr,
                                             dk,
                                             kmesh,
                                             Rmesh,
                                             psb_);
    ModuleBase::Ylm::set_coefficients();
    mgt_.init_Gaunt_CH(Lmax);
    mgt_.init_Gaunt(Lmax);
#endif
}

void toW90_LCAO::set_R_coor(const UnitCell& ucell, const Grid_Driver& gd)
{
    int R_minX = int(-gd.getGlayerX_minus());
    int R_minY = int(-gd.getGlayerY_minus());
    int R_minZ = int(-gd.getGlayerZ_minus());

    int R_x = gd.getGlayerX() + gd.getGlayerX_minus();
    int R_y = gd.getGlayerY() + gd.getGlayerY_minus();
    int R_z = gd.getGlayerZ() + gd.getGlayerZ_minus();

    int R_num = R_x * R_y * R_z;
    R_coor_car.resize(R_num);

    int count = 0;
    for (int ix = 0; ix < R_x; ix++)
    {
        for (int iy = 0; iy < R_y; iy++)
        {
            for (int iz = 0; iz < R_z; iz++)
            {
                ModuleBase::Vector3<double> tmpR(ix + R_minX, iy + R_minY, iz + R_minZ);
                R_coor_car[count] = tmpR * ucell.latvec;
                count++;
            }
        }
    }
}

void toW90_LCAO::count_delta_k(const UnitCell& ucell, const K_Vectors& kv)
{
    std::set<Coordinate_3D> delta_k_all_tmp;
    for (int ik = 0; ik < cal_num_kpts; ik++)
    {
        for (int ib = 0; ib < nntot; ib++)
        {
            int ikb = nnlist[ik][ib];
            ModuleBase::Vector3<double> G = nncell[ik][ib];

            int cal_ik = ik + start_k_index;
            int cal_ikb = ikb + start_k_index;

            ModuleBase::Vector3<double> ik_car = kv.kvec_c[ik];
            ModuleBase::Vector3<double> ikb_car = kv.kvec_c[ikb] + G * ucell.G;
            Abfs::Vector3_Order<double> dk = (ikb_car - ik_car) * ucell.tpiba;
            Coordinate_3D temp_dk(dk.x, dk.y, dk.z);
            delta_k_all_tmp.insert(temp_dk);
        }
    }

    delta_k_all.resize(delta_k_all_tmp.size());

    int index = 0;
    for (auto& delta_k: delta_k_all_tmp)
    {
        delta_k_all_index[delta_k] = index;
        delta_k_all[index] = delta_k;
        index++;
    }
}

void toW90_LCAO::unkdotkb(const UnitCell& ucell,
                                const K_Vectors& kv,
                                const psi::Psi<std::complex<double>>& psi_in,
                                const int& ik,
                                const int& ikb,
                                const ModuleBase::Vector3<double> G,
                                ModuleBase::ComplexMatrix& Mmn)
{
    Mmn.create(num_bands, num_bands);

    int row = this->ParaV->get_row_size();
    int col = this->ParaV->get_col_size();
    int nloc = row * col;
    std::vector<std::complex<double>> midmatrix(nloc);

    int R_num = R_coor_car.size();
    ModuleBase::Vector3<double> ik_car = kv.kvec_c[ik];
    ModuleBase::Vector3<double> ikb_car = kv.kvec_c[ikb] + G * ucell.G;
    Abfs::Vector3_Order<double> dk = (ikb_car - ik_car) * ucell.tpiba;
    Coordinate_3D temp_dk(dk.x, dk.y, dk.z);
    int delta_k_index = delta_k_all_index[temp_dk];

    hamilt::HContainer<std::complex<double>>* tmp_FR_container = FR[delta_k_index].get_FR_pointer();
    auto row_indexes = ParaV->get_indexes_row();
    auto col_indexes = ParaV->get_indexes_col();

    for (int iap = 0; iap < tmp_FR_container->size_atom_pairs(); ++iap)
    {
        int atom_i = tmp_FR_container->get_atom_pair(iap).get_atom_i();
        int atom_j = tmp_FR_container->get_atom_pair(iap).get_atom_j();
        int start_i = ParaV->atom_begin_row[atom_i];
        int start_j = ParaV->atom_begin_col[atom_j];
        int row_size = ParaV->get_nrow_atom(atom_i);
        int col_size = ParaV->get_ncol_atom(atom_j);
        for (int iR = 0; iR < tmp_FR_container->get_atom_pair(iap).get_R_size(); ++iR)
        {
            auto& matrix = tmp_FR_container->get_atom_pair(iap).get_HR_values(iR);
            const ModuleBase::Vector3<int> r_index = tmp_FR_container->get_atom_pair(iap).get_R_index(iR);
            ModuleBase::Vector3<double> dR
                = ModuleBase::Vector3<double>(r_index.x, r_index.y, r_index.z) * ucell.latvec;
            double phase = ikb_car * dR * ModuleBase::TWO_PI;
            std::complex<double> kRn_phase = std::exp(ModuleBase::IMAG_UNIT * phase);
            for (int i = 0; i < row_size; ++i)
            {
                int mu = row_indexes[start_i + i];
                int ir = ParaV->global2local_row(mu);
                for (int j = 0; j < col_size; ++j)
                {
                    int nu = col_indexes[start_j + j];
                    int ic = ParaV->global2local_col(nu);
                    int index = ic * row + ir;
                    midmatrix[index] += kRn_phase * matrix.get_value(i, j);
                }
            }
        }
    }

    char transa = 'C';
    char transb = 'N';
    int Bands = nbands_;
    int nlocal = PARAM.globalv.nlocal;
    std::complex<double> alpha = {1.0, 0.0}, beta = {0.0, 0.0};
    int one = 1;

    std::vector<std::complex<double>> C_matrix(nloc);
    std::vector<std::complex<double>> out_matrix(nloc);

#ifdef __MPI
    ScalapackConnector::gemm(transa,
            transb,
            Bands,
            nlocal,
            nlocal,
            alpha,
            &psi_in(ik, 0, 0),
            one,
            one,
            this->ParaV->desc,
            midmatrix.data(),
            one,
            one,
            this->ParaV->desc,
            beta,
            C_matrix.data(),
            one,
            one,
            this->ParaV->desc);

    ScalapackConnector::gemm(transb,
            transb,
            Bands,
            Bands,
            nlocal,
            alpha,
            C_matrix.data(),
            one,
            one,
            this->ParaV->desc,
            &psi_in(ikb, 0, 0),
            one,
            one,
            this->ParaV->desc,
            beta,
            out_matrix.data(),
            one,
            one,
            this->ParaV->desc);
#endif

    int count_m = -1;
    for (int m = 0; m < nbands_; m++)
    {
        if (exclude_bands.count(m)) {
            continue;
}
        count_m++;

        int ir = this->ParaV->global2local_row(m);
        if (ir >= 0)
        {
            int count_n = -1;
            for (int n = 0; n < nbands_; n++)
            {
                if (exclude_bands.count(n)) {
                    continue;
}
                count_n++;

                int ic = this->ParaV->global2local_col(n);
                if (ic >= 0)
                {
                    int index = ic * row + ir;
                    Mmn(count_m, count_n) = out_matrix[index];
                }
            }
        }
    }

    Parallel_Reduce::reduce_all(Mmn.c, num_bands * num_bands);
}

void toW90_LCAO::produce_basis_orb()
{
    int mat_Nr = orb_.Phi[0].PhiLN(0, 0).getNr();
    int count_Nr = 0;

    orbs.resize(orb_.get_ntype());
    for (int T = 0; T < orb_.get_ntype(); ++T)
    {
        count_Nr = orb_.Phi[T].PhiLN(0, 0).getNr();
        if (count_Nr > mat_Nr)
        {
            mat_Nr = count_Nr;
            orb_r_ntype = T;
        }

        orbs[T].resize(orb_.Phi[T].getLmax() + 1);
        for (int L = 0; L <= orb_.Phi[T].getLmax(); ++L)
        {
            orbs[T][L].resize(orb_.Phi[T].getNchi(L));
            for (int N = 0; N < orb_.Phi[T].getNchi(L); ++N)
            {
                const auto& orb_origin = orb_.Phi[T].PhiLN(L, N);
                orbs[T][L][N].set_orbital_info(orb_origin.getLabel(),
                                               orb_origin.getType(),
                                               orb_origin.getL(),
                                               orb_origin.getChi(),
                                               orb_origin.getNr(),
                                               orb_origin.getRab(),
                                               orb_origin.getRadial(),
                                               Numerical_Orbital_Lm::Psi_Type::Psi,
                                               orb_origin.getPsi(),
                                               static_cast<int>(orb_origin.getNk() * kmesh_times) | 1,
                                               orb_origin.getDk(),
                                               orb_origin.getDruniform(),
                                               false,
                                               true,
                                               PARAM.inp.cal_force);
            }
        }
    }
}

void toW90_LCAO::produce_trial_in_lcao()
{
    A_orbs.resize(num_wannier);

    std::vector<double> r(mesh_r);
    std::vector<double> rab(mesh_r);

    for (int ir = 0; ir < mesh_r; ir++)
    {
        rab[ir] = dr;
        r[ir] = ir * dr;
    }

    const auto& orb_origin = orb_.Phi[orb_r_ntype].PhiLN(0, 0);

    std::vector<double> psi(mesh_r);
    std::vector<double> psir(mesh_r);
    std::vector<double> inner(mesh_r);
    for (int i = 0; i < num_wannier; i++)
    {
        double alfa32 = pow(alfa[i], 3.0 / 2.0);
        double alfa_new = alfa[i];
        int wannier_index = i;

        if (rvalue[i] == 1)
        {
            for (int ir = 0; ir < mesh_r; ir++)
            {
                psi[ir] = 2.0 * alfa32 * exp(-alfa_new * r[ir]);
            }
        }

        if (rvalue[i] == 2)
        {
            for (int ir = 0; ir < mesh_r; ir++)
            {
                psi[ir] = 1.0 / sqrt(8.0) * alfa32 * (2.0 - alfa_new * r[ir]) * exp(-alfa_new * r[ir] * 0.5);
            }
        }

        if (rvalue[i] == 3)
        {
            for (int ir = 0; ir < mesh_r; ir++)
            {
                psi[ir] = sqrt(4.0 / 27.0) * alfa32
                          * (1.0 - 2.0 / 3.0 * alfa_new * r[ir] + 2.0 / 27.0 * pow(alfa_new, 2.0) * r[ir] * r[ir])
                          * exp(-alfa_new * r[ir] * 1.0 / 3.0);
            }
        }

        for (int ir = 0; ir < mesh_r; ir++)
        {
            psir[ir] = psi[ir] * r[ir];
        }

        // renormalize radial wave functions
        for (int ir = 0; ir < mesh_r; ir++)
        {
            inner[ir] = psir[ir] * psir[ir];
        }
        double unit = 0.0;
        ModuleBase::Integral::Simpson_Integral(mesh_r, inner.data(), rab.data(), unit);

        for (int ir = 0; ir < mesh_r; ir++)
        {
            psi[ir] /= sqrt(unit);
        }

        if (L[i] >= 0)
        {
            A_orbs[i].resize(1);

            A_orbs[i][0].set_orbital_info(orb_origin.getLabel(),
                                          orb_origin.getType(),
                                          L[i],
                                          1,
                                          mesh_r,
                                          rab.data(),
                                          r.data(),
                                          Numerical_Orbital_Lm::Psi_Type::Psi,
                                          psi.data(),
                                          static_cast<int>(orb_origin.getNk() * kmesh_times) | 1,
                                          orb_origin.getDk(),
                                          orb_origin.getDruniform(),
                                          false,
                                          true,
                                          PARAM.inp.cal_force);
        }
        else
        {
            int tmp_size = 0;

            if (L[i] == -1 || L[i] == -2 || L[i] == -3) {
                tmp_size = 2;
}

            if (L[i] == -4 || L[i] == -5) {
                tmp_size = 3;
}

            A_orbs[i].resize(tmp_size);

            for (int tmp_L = 0; tmp_L < tmp_size; tmp_L++)
            {
                A_orbs[i][tmp_L].set_orbital_info(orb_origin.getLabel(),
                                                  orb_origin.getType(),
                                                  tmp_L,
                                                  1,
                                                  mesh_r,
                                                  rab.data(),
                                                  r.data(),
                                                  Numerical_Orbital_Lm::Psi_Type::Psi,
                                                  psi.data(),
                                                  static_cast<int>(orb_origin.getNk() * kmesh_times) | 1,
                                                  orb_origin.getDk(),
                                                  orb_origin.getDruniform(),
                                                  false,
                                                  true,
                                                  PARAM.inp.cal_force);
            }
        }
    }
}

void toW90_LCAO::construct_overlap_table_project()
{
    int row = this->ParaV->get_row_size();
    int global_ir = 0;

    for (int ir = 0; ir < row; ir += npol_)
    {
        global_ir = ParaV->local2global_row(ir);
        int orb_index_row = global_ir / npol_;

        for (int wannier_index = 0; wannier_index < num_wannier; wannier_index++)
        {
            if (L[wannier_index] >= 0)
            {
                center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].insert(std::make_pair(
                    0,
                    Center2_Orb::Orb11(orbs[iw2it[orb_index_row]][iw2iL[orb_index_row]][iw2iN[orb_index_row]],
                                       A_orbs[wannier_index][0],
                                       psb_,
                                       mgt_)));
            }
            else
            {
                int tmp_size = 0;

                if (L[wannier_index] == -1 || L[wannier_index] == -2 || L[wannier_index] == -3) {
                    tmp_size = 2;
}

                if (L[wannier_index] == -4 || L[wannier_index] == -5) {
                    tmp_size = 3;
}

                for (int tmp_L = 0; tmp_L < tmp_size; tmp_L++)
                {
                    center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].insert(std::make_pair(
                        tmp_L,
                        Center2_Orb::Orb11(orbs[iw2it[orb_index_row]][iw2iL[orb_index_row]][iw2iN[orb_index_row]],
                                           A_orbs[wannier_index][tmp_L],
                                           psb_,
                                           mgt_)));
                }
            }
        }
    }

    for (auto& co1: center2_orb11_A)
    {
        for (auto& co2: co1.second)
        {
            for (auto& co3: co2.second)
            {
                co3.second.init_radial_table();
            }
        }
    }
}

#endif
