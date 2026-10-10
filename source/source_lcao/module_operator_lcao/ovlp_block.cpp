#include "ovlp_block.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_basis/module_nao/two_center_integrator.h"
#include "source_cell/unitcell.h"
#include <complex>

namespace hamilt
{
template <typename TR>
void cal_overlap_block(const UnitCell& cell, const TwoCenterIntegrator& integrator,
    const int atom_i, const int atom_j, const Parallel_Orbitals& distribution,
    const ModuleBase::Vector3<double>& displacement, TR* values)
{
    const int type_i = cell.iat2it[atom_i];
    const int type_j = cell.iat2it[atom_j];
    const Atom& bra = cell.atoms[type_i];
    const Atom& ket = cell.atoms[type_j];
    const int npol = cell.get_npol();
    const auto rows = distribution.get_indexes_row(atom_i);
    const auto cols = distribution.get_indexes_col(atom_j);
    const int column_count = cols.size();
    const int spin_stride = column_count + 1;
    for (std::size_t row = 0; row < rows.size(); row += npol)
    {
        const int mu = rows[row] / npol;
        const int l1 = bra.iw2l[mu];
        const int n1 = bra.iw2n[mu];
        const int bra_m = bra.iw2m[mu];
        const int m1 = bra_m % 2 == 0 ? -bra_m / 2 : (bra_m + 1) / 2;
        for (std::size_t col = 0; col < cols.size(); col += npol)
        {
            const int nu = cols[col] / npol;
            const int l2 = ket.iw2l[nu];
            const int n2 = ket.iw2n[nu];
            const int ket_m = ket.iw2m[nu];
            const int m2 = ket_m % 2 == 0 ? -ket_m / 2 : (ket_m + 1) / 2;
            double integral = 0.0;
            integrator.calculate(type_i, l1, n1, m1, type_j, l2, n2, m2, displacement, &integral);
            const std::size_t offset = row * column_count + col;
            for (int spin = 0; spin < npol; ++spin)
            {
                values[offset + spin * spin_stride] += integral;
            }
        }
    }
}

template void cal_overlap_block(const UnitCell&, const TwoCenterIntegrator&, int, int,
    const Parallel_Orbitals&, const ModuleBase::Vector3<double>&, double*);
template void cal_overlap_block(const UnitCell&, const TwoCenterIntegrator&, int, int,
    const Parallel_Orbitals&, const ModuleBase::Vector3<double>&, std::complex<double>*);
}
