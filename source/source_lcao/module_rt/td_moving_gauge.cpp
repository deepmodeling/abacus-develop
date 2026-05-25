#include "td_moving_gauge.h"

#include "source_base/global_function.h"
#include "source_base/libm/libm.h" // sincos

namespace module_rt
{

TD_MovingGauge::~TD_MovingGauge()
{
    for (int i = 0; i < nat_; ++i)
    {
        delete DR_x_[i];
        delete DR_y_[i];
        delete DR_z_[i];
    }
}

template <typename T_sR>
void TD_MovingGauge::init_DR(const hamilt::HContainer<T_sR>* sR_template,
                             const UnitCell* ucell,
                             const Parallel_Orbitals* paraV,
                             TwoCenterIntegrator* intor)
{
    nat_ = ucell->nat;
    DR_x_.resize(nat_, nullptr);
    DR_y_.resize(nat_, nullptr);
    DR_z_.resize(nat_, nullptr);

    // 1. Allocate an HContainer for each atom
    for (int K = 0; K < nat_; ++K)
    {
        DR_x_[K] = new hamilt::HContainer<double>(paraV);
        DR_y_[K] = new hamilt::HContainer<double>(paraV);
        DR_z_[K] = new hamilt::HContainer<double>(paraV);
    }

    // 2. Construct the sparsity pattern based on sR_template, only allocate terms where delta_{JK} is non-zero
    for (int iap = 0; iap < sR_template->size_atom_pairs(); ++iap)
    {
        const auto& ap = sR_template->get_atom_pair(iap);
        int iat1 = ap.get_atom_i();
        int iat2 = ap.get_atom_j(); // target ket atom J

        hamilt::AtomPair<double> ap_x(iat1, iat2, paraV);
        hamilt::AtomPair<double> ap_y(iat1, iat2, paraV);
        hamilt::AtomPair<double> ap_z(iat1, iat2, paraV);

        for (int ir = 0; ir < ap.get_R_size(); ++ir)
        {
            auto R_idx = ap.get_R_index(ir);
            ap_x.get_HR_values(R_idx.x, R_idx.y, R_idx.z);
            ap_y.get_HR_values(R_idx.x, R_idx.y, R_idx.z);
            ap_z.get_HR_values(R_idx.x, R_idx.y, R_idx.z);
        }

        // Only insert this pair into the container corresponding to atom iat2
        DR_x_[iat2]->insert_pair(ap_x);
        DR_y_[iat2]->insert_pair(ap_y);
        DR_z_[iat2]->insert_pair(ap_z);
    }

    // 3. Allocate memory for all DR_[K] containers
    for (int K = 0; K < nat_; ++K)
    {
        DR_x_[K]->allocate(nullptr, true);
        DR_y_[K]->allocate(nullptr, true);
        DR_z_[K]->allocate(nullptr, true);
    }

    // 4. Calculate and fill the R-space derivatives of the overlap matrix
    int npol = ucell->get_npol();
    for (int iap = 0; iap < sR_template->size_atom_pairs(); ++iap)
    {
        const auto& ap = sR_template->get_atom_pair(iap);
        int iat1 = ap.get_atom_i();
        int iat2 = ap.get_atom_j();

        int T1 = ucell->iat2it[iat1];
        int T2 = ucell->iat2it[iat2];
        const Atom& atom1 = ucell->atoms[T1];
        const Atom& atom2 = ucell->atoms[T2];

        auto row_indexes = paraV->get_indexes_row(iat1);
        auto col_indexes = paraV->get_indexes_col(iat2);

        for (int ir = 0; ir < ap.get_R_size(); ++ir)
        {
            auto R_idx = ap.get_R_index(ir);
            ModuleBase::Vector3<double> dtau = ucell->cal_dtau(iat1, iat2, R_idx);

            int R_arr[3] = {R_idx.x, R_idx.y, R_idx.z};
            double* dx = DR_x_[iat2]->data(iat1, iat2, R_arr);
            double* dy = DR_y_[iat2]->data(iat1, iat2, R_arr);
            double* dz = DR_z_[iat2]->data(iat1, iat2, R_arr);

            for (int iw1l = 0; iw1l < row_indexes.size(); iw1l += npol)
            {
                const int iw1 = row_indexes[iw1l] / npol;
                int L1 = atom1.iw2l[iw1];
                int N1 = atom1.iw2n[iw1];
                int m1 = atom1.iw2m[iw1];
                int M1 = (m1 % 2 == 0) ? -m1 / 2 : (m1 + 1) / 2;

                for (int iw2l = 0; iw2l < col_indexes.size(); iw2l += npol)
                {
                    const int iw2 = col_indexes[iw2l] / npol;
                    int L2 = atom2.iw2l[iw2];
                    int N2 = atom2.iw2n[iw2];
                    int m2 = atom2.iw2m[iw2];
                    int M2 = (m2 % 2 == 0) ? -m2 / 2 : (m2 + 1) / 2;

                    double olm[4] = {0.0, 0.0, 0.0, 0.0};
                    // out stores the integral value in olm[0], grad_out stores the gradient in olm[1] to olm[3]
                    intor->calculate(T1, L1, N1, M1, T2, L2, N2, M2, dtau * ucell->lat0, &olm[0], &olm[1]);

                    // Handle the spin dimension (the overlap is diagonal in spin space)
                    for (int is1 = 0; is1 < npol; ++is1)
                    {
                        for (int is2 = 0; is2 < npol; ++is2)
                        {
                            int r_offset = iw1l + is1;
                            int c_offset = iw2l + is2;
                            int linear_idx = r_offset * col_indexes.size() + c_offset;

                            if (is1 == is2)
                            {
                                dx[linear_idx] = olm[1];
                                dy[linear_idx] = olm[2];
                                dz[linear_idx] = olm[3];
                            }
                            else
                            {
                                dx[linear_idx] = 0.0;
                                dy[linear_idx] = 0.0;
                                dz[linear_idx] = 0.0;
                            }
                        }
                    }
                }
            }
        }
    }
}

template <typename T_sR>
void TD_MovingGauge::update_DR(const hamilt::HContainer<T_sR>* sR_template,
                               const UnitCell* ucell,
                               const Parallel_Orbitals* paraV,
                               TwoCenterIntegrator* intor)
{
    // Update the R-space derivatives of the overlap matrix with the current atomic positions
    int npol = ucell->get_npol();
    for (int iap = 0; iap < sR_template->size_atom_pairs(); ++iap)
    {
        const auto& ap = sR_template->get_atom_pair(iap);
        int iat1 = ap.get_atom_i();
        int iat2 = ap.get_atom_j();

        int T1 = ucell->iat2it[iat1];
        int T2 = ucell->iat2it[iat2];
        const Atom& atom1 = ucell->atoms[T1];
        const Atom& atom2 = ucell->atoms[T2];

        auto row_indexes = paraV->get_indexes_row(iat1);
        auto col_indexes = paraV->get_indexes_col(iat2);

        for (int ir = 0; ir < ap.get_R_size(); ++ir)
        {
            auto R_idx = ap.get_R_index(ir);
            ModuleBase::Vector3<double> dtau = ucell->cal_dtau(iat1, iat2, R_idx);

            int R_arr[3] = {R_idx.x, R_idx.y, R_idx.z};
            double* dx = DR_x_[iat2]->data(iat1, iat2, R_arr);
            double* dy = DR_y_[iat2]->data(iat1, iat2, R_arr);
            double* dz = DR_z_[iat2]->data(iat1, iat2, R_arr);

            for (int iw1l = 0; iw1l < row_indexes.size(); iw1l += npol)
            {
                const int iw1 = row_indexes[iw1l] / npol;
                int L1 = atom1.iw2l[iw1];
                int N1 = atom1.iw2n[iw1];
                int m1 = atom1.iw2m[iw1];
                int M1 = (m1 % 2 == 0) ? -m1 / 2 : (m1 + 1) / 2;

                for (int iw2l = 0; iw2l < col_indexes.size(); iw2l += npol)
                {
                    const int iw2 = col_indexes[iw2l] / npol;
                    int L2 = atom2.iw2l[iw2];
                    int N2 = atom2.iw2n[iw2];
                    int m2 = atom2.iw2m[iw2];
                    int M2 = (m2 % 2 == 0) ? -m2 / 2 : (m2 + 1) / 2;

                    double olm[4] = {0.0, 0.0, 0.0, 0.0};
                    intor->calculate(T1, L1, N1, M1, T2, L2, N2, M2, dtau * ucell->lat0, &olm[0], &olm[1]);

                    for (int is1 = 0; is1 < npol; ++is1)
                    {
                        for (int is2 = 0; is2 < npol; ++is2)
                        {
                            int r_offset = iw1l + is1;
                            int c_offset = iw2l + is2;
                            int linear_idx = r_offset * col_indexes.size() + c_offset;

                            if (is1 == is2)
                            {
                                dx[linear_idx] = olm[1];
                                dy[linear_idx] = olm[2];
                                dz[linear_idx] = olm[3];
                            }
                            else
                            {
                                dx[linear_idx] = 0.0;
                                dy[linear_idx] = 0.0;
                                dz[linear_idx] = 0.0;
                            }
                        }
                    }
                }
            }
        }
    }
}

template <typename TK>
void TD_MovingGauge::get_D_k(int K, const ModuleBase::Vector3<double>& kvec_d, TK* Dk_x, TK* Dk_y, TK* Dk_z, int hk_ld)
    const
{
    hamilt::folding_HR(*(DR_x_[K]), Dk_x, kvec_d, hk_ld, 1);
    hamilt::folding_HR(*(DR_y_[K]), Dk_y, kvec_d, hk_ld, 1);
    hamilt::folding_HR(*(DR_z_[K]), Dk_z, kvec_d, hk_ld, 1);
}

template <typename TK>
void TD_MovingGauge::get_P_k(const UnitCell* ucell,
                             const ModuleBase::Vector3<double>& kvec_d,
                             TK* P_k,
                             int matrix_size,
                             int hk_ld) const
{
    std::vector<TK> Dk_x(matrix_size, TK(0.0, 0.0));
    std::vector<TK> Dk_y(matrix_size, TK(0.0, 0.0));
    std::vector<TK> Dk_z(matrix_size, TK(0.0, 0.0));

    for (int K = 0; K < nat_; ++K)
    {
        std::fill(Dk_x.begin(), Dk_x.end(), TK(0.0, 0.0));
        std::fill(Dk_y.begin(), Dk_y.end(), TK(0.0, 0.0));
        std::fill(Dk_z.begin(), Dk_z.end(), TK(0.0, 0.0));

        this->get_D_k(K, kvec_d, Dk_x.data(), Dk_y.data(), Dk_z.data(), hk_ld);

        // Obtain the real-time velocity of atom K from the UnitCell (in Hartree atomic units)
        int it = ucell->iat2it[K];
        int ia = ucell->iat2ia[K];
        double vx = ucell->atoms[it].vel[ia].x;
        double vy = ucell->atoms[it].vel[ia].y;
        double vz = ucell->atoms[it].vel[ia].z;

        // Construct the coefficients: P = -i * v * D
        // Unit conversion: Hartree a.u. to Rydberg a.u. requires multiplying
        TK coef_x(0.0, -2.0 * vx);
        TK coef_y(0.0, -2.0 * vy);
        TK coef_z(0.0, -2.0 * vz);

        // Accumulate the contribution from atom K to the P_k matrix
        for (int i = 0; i < matrix_size; ++i)
        {
            P_k[i] += coef_x * Dk_x[i] + coef_y * Dk_y[i] + coef_z * Dk_z[i];
        }
    }
}

template void TD_MovingGauge::init_DR<double>(const hamilt::HContainer<double>* sR_template,
                                              const UnitCell* ucell,
                                              const Parallel_Orbitals* paraV,
                                              TwoCenterIntegrator* intor);

template void TD_MovingGauge::init_DR<std::complex<double>>(const hamilt::HContainer<std::complex<double>>* sR_template,
                                                            const UnitCell* ucell,
                                                            const Parallel_Orbitals* paraV,
                                                            TwoCenterIntegrator* intor);

template void TD_MovingGauge::update_DR<double>(const hamilt::HContainer<double>* sR_template,
                                                const UnitCell* ucell,
                                                const Parallel_Orbitals* paraV,
                                                TwoCenterIntegrator* intor);

template void TD_MovingGauge::update_DR<std::complex<double>>(
    const hamilt::HContainer<std::complex<double>>* sR_template,
    const UnitCell* ucell,
    const Parallel_Orbitals* paraV,
    TwoCenterIntegrator* intor);

template void TD_MovingGauge::get_D_k<std::complex<double>>(int K,
                                                            const ModuleBase::Vector3<double>& kvec_d,
                                                            std::complex<double>* Dk_x,
                                                            std::complex<double>* Dk_y,
                                                            std::complex<double>* Dk_z,
                                                            int hk_ld) const;

template void TD_MovingGauge::get_P_k<std::complex<double>>(const UnitCell* ucell,
                                                            const ModuleBase::Vector3<double>& kvec_d,
                                                            std::complex<double>* P_k,
                                                            int matrix_size,
                                                            int hk_ld) const;

} // namespace module_rt
