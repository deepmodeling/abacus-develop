#include "td_moving_gauge.h"

#include "source_base/basis_flattened_indexing.h"
#include "source_base/parallel_blkcyc.h"

#include <cstddef>

namespace module_rt
{

namespace
{

struct LocalBasis
{
    const BasisLabel* label;
    int offset;
};

struct LocalBasisBlock
{
    std::vector<LocalBasis> spatial;
    int size = 0;
};

enum class BasisAxis
{
    Row,
    Column
};

LocalBasisBlock local_basis_block(const BasisFlattenedIndexing& basis,
                                  const BlockCyclicDistribution& layout,
                                  const int iat,
                                  const BasisAxis axis)
{
    LocalBasisBlock block;
    const BasisRange basis_range = basis.range(iat);
    for (std::size_t ibas = basis_range.ibas0; ibas < basis_range.ibas1; ++ibas)
    {
        const int iloc = axis == BasisAxis::Row
                             ? layout.irloc(static_cast<int>(ibas))
                             : layout.icloc(static_cast<int>(ibas));
        if (iloc < 0)
        {
            continue;
        }
        const BasisLabel& label = basis.label(ibas);
        if (label.spin == 0)
        {
            const LocalBasis local_basis = {&label, block.size};
            block.spatial.push_back(local_basis);
        }
        ++block.size;
    }
    return block;
}

} // namespace

template <typename T_sR>
void TD_MovingGauge::init_DR(const hamilt::HContainer<T_sR>* sR_template,
                             const UnitCell* ucell,
                             const BasisFlattenedIndexing& basis,
                             const BlockCyclicDistribution& layout,
                             TwoCenterIntegrator* intor)
{
    nat_ = ucell->nat;
    dr_x_.resize(nat_);
    dr_y_.resize(nat_);
    dr_z_.resize(nat_);

    // 1. Create an HContainer wrapper for each atom.
    for (int iat = 0; iat < nat_; ++iat)
    {
        dr_x_[iat].reset(new hamilt::HContainer<double>(basis, layout, nullptr));
        dr_y_[iat].reset(new hamilt::HContainer<double>(basis, layout, nullptr));
        dr_z_[iat].reset(new hamilt::HContainer<double>(basis, layout, nullptr));
    }

    // 2. Construct the sparsity pattern based on sR_template, only allocate terms where delta_{JK} is non-zero
    for (int iap = 0; iap < sR_template->size_atom_pairs(); ++iap)
    {
        const auto& ap = sR_template->get_atom_pair(iap);
        const int iat1 = ap.get_atom_i();
        const int iat2 = ap.get_atom_j();

        for (int ir = 0; ir < ap.get_R_size(); ++ir)
        {
            const ModuleBase::Vector3<int> R_index = ap.get_R_index(ir);
            dr_x_[iat2]->insert_pair(iat1, iat2, R_index);
            dr_y_[iat2]->insert_pair(iat1, iat2, R_index);
            dr_z_[iat2]->insert_pair(iat1, iat2, R_index);
        }
    }

    // 3. Allocate memory for all dr_[iat] containers.
    for (int iat = 0; iat < nat_; ++iat)
    {
        dr_x_[iat]->allocate(nullptr, true);
        dr_y_[iat]->allocate(nullptr, true);
        dr_z_[iat]->allocate(nullptr, true);
    }

    // 4. Calculate and fill the R-space derivatives of the overlap matrix
    this->fill_DR(sR_template, ucell, basis, layout, intor);
}

template <typename T_sR>
void TD_MovingGauge::update_DR(const hamilt::HContainer<T_sR>* sR_template,
                               const UnitCell* ucell,
                               const BasisFlattenedIndexing& basis,
                               const BlockCyclicDistribution& layout,
                               TwoCenterIntegrator* intor)
{
    this->fill_DR(sR_template, ucell, basis, layout, intor);
}

template <typename T_sR>
void TD_MovingGauge::fill_DR(const hamilt::HContainer<T_sR>* sR_template,
                             const UnitCell* ucell,
                             const BasisFlattenedIndexing& basis,
                             const BlockCyclicDistribution& layout,
                             TwoCenterIntegrator* intor)
{
    const int npol = basis.npol();
    std::vector<LocalBasisBlock> rows(basis.nat());
    std::vector<LocalBasisBlock> columns(basis.nat());
    for (int iat = 0; iat < basis.nat(); ++iat)
    {
        rows[iat] = local_basis_block(basis, layout, iat, BasisAxis::Row);
        columns[iat] = local_basis_block(basis, layout, iat, BasisAxis::Column);
    }

    for (int iap = 0; iap < sR_template->size_atom_pairs(); ++iap)
    {
        const auto& ap = sR_template->get_atom_pair(iap);
        const int iat1 = ap.get_atom_i();
        const int iat2 = ap.get_atom_j();
        const LocalBasisBlock& atom_rows = rows[iat1];
        const LocalBasisBlock& atom_columns = columns[iat2];

        for (int ir = 0; ir < ap.get_R_size(); ++ir)
        {
            const ModuleBase::Vector3<int> R_index = ap.get_R_index(ir);
            const ModuleBase::Vector3<double> dtau = ucell->cal_dtau(iat1, iat2, R_index);
            int R_arr[3] = {R_index.x, R_index.y, R_index.z};
            double* dx = dr_x_[iat2]->data(iat1, iat2, R_arr);
            double* dy = dr_y_[iat2]->data(iat1, iat2, R_arr);
            double* dz = dr_z_[iat2]->data(iat1, iat2, R_arr);

            for (const LocalBasis& row : atom_rows.spatial)
            {
                for (const LocalBasis& column : atom_columns.spatial)
                {
                    double olm[4] = {0.0, 0.0, 0.0, 0.0};
                    intor->calculate(row.label->itype,
                                     row.label->l,
                                     row.label->izeta,
                                     row.label->m,
                                     column.label->itype,
                                     column.label->l,
                                     column.label->izeta,
                                     column.label->m,
                                     dtau * ucell->lat0,
                                     &olm[0],
                                     &olm[1]);

                    for (int spin_row = 0; spin_row < npol; ++spin_row)
                    {
                        for (int spin_column = 0; spin_column < npol; ++spin_column)
                        {
                            const int local_row = row.offset + spin_row;
                            const int local_column = column.offset + spin_column;
                            const int linear_index = local_row * atom_columns.size + local_column;

                            if (spin_row == spin_column)
                            {
                                dx[linear_index] = olm[1];
                                dy[linear_index] = olm[2];
                                dz[linear_index] = olm[3];
                            }
                            else
                            {
                                dx[linear_index] = 0.0;
                                dy[linear_index] = 0.0;
                                dz[linear_index] = 0.0;
                            }
                        }
                    }
                }
            }
        }
    }
}

template <typename TK>
void TD_MovingGauge::get_D_k(const int iat,
                             const ModuleBase::Vector3<double>& kvec_d,
                             TK* Dk_x,
                             TK* Dk_y,
                             TK* Dk_z,
                             const int hk_ld) const
{
    hamilt::folding_HR(*(dr_x_[iat]), Dk_x, kvec_d, hk_ld, 1);
    hamilt::folding_HR(*(dr_y_[iat]), Dk_y, kvec_d, hk_ld, 1);
    hamilt::folding_HR(*(dr_z_[iat]), Dk_z, kvec_d, hk_ld, 1);
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

    for (int iat = 0; iat < nat_; ++iat)
    {
        std::fill(Dk_x.begin(), Dk_x.end(), TK(0.0, 0.0));
        std::fill(Dk_y.begin(), Dk_y.end(), TK(0.0, 0.0));
        std::fill(Dk_z.begin(), Dk_z.end(), TK(0.0, 0.0));

        this->get_D_k(iat, kvec_d, Dk_x.data(), Dk_y.data(), Dk_z.data(), hk_ld);

        // Obtain the real-time velocity of atom iat from the UnitCell (in Hartree atomic units).
        const int it = ucell->iat2it[iat];
        const int ia = ucell->iat2ia[iat];
        const double vx = ucell->atoms[it].vel[ia].x;
        const double vy = ucell->atoms[it].vel[ia].y;
        const double vz = ucell->atoms[it].vel[ia].z;

        // Construct the coefficients: P = -i * v * D
        // Unit conversion: Hartree a.u. to Rydberg a.u. requires multiplying
        const TK coef_x(0.0, -2.0 * vx);
        const TK coef_y(0.0, -2.0 * vy);
        const TK coef_z(0.0, -2.0 * vz);

        // Accumulate the contribution from atom iat to the P_k matrix.
        for (int i = 0; i < matrix_size; ++i)
        {
            P_k[i] += coef_x * Dk_x[i] + coef_y * Dk_y[i] + coef_z * Dk_z[i];
        }
    }
}

template void TD_MovingGauge::init_DR<double>(const hamilt::HContainer<double>* sR_template,
                                              const UnitCell* ucell,
                                              const BasisFlattenedIndexing& basis,
                                              const BlockCyclicDistribution& layout,
                                              TwoCenterIntegrator* intor);

template void TD_MovingGauge::init_DR<std::complex<double>>(const hamilt::HContainer<std::complex<double>>* sR_template,
                                                            const UnitCell* ucell,
                                                            const BasisFlattenedIndexing& basis,
                                                            const BlockCyclicDistribution& layout,
                                                            TwoCenterIntegrator* intor);

template void TD_MovingGauge::update_DR<double>(const hamilt::HContainer<double>* sR_template,
                                                const UnitCell* ucell,
                                                const BasisFlattenedIndexing& basis,
                                                const BlockCyclicDistribution& layout,
                                                TwoCenterIntegrator* intor);

template void TD_MovingGauge::update_DR<std::complex<double>>(
    const hamilt::HContainer<std::complex<double>>* sR_template,
    const UnitCell* ucell,
    const BasisFlattenedIndexing& basis,
    const BlockCyclicDistribution& layout,
    TwoCenterIntegrator* intor);

template void TD_MovingGauge::get_D_k<std::complex<double>>(int iat,
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
