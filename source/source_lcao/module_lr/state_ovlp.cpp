#include "state_ovlp.h"
#include "source_base/matrix3.h"
#include "source_base/module_external/blas_connector.h"
#include <cmath>
#include <complex>
#include <stdexcept>
#include <type_traits>
#include <algorithm>
#include "source_base/parallel_2d.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_cell/unitcell.h"
#include "source_lcao/module_operator_lcao/ovlp_block.h"
#include "source_hamilt/module_hcontainer/hcontainer_funcs.h"
#ifdef __MPI
#include "source_base/module_external/scalapack_connector.h"
#endif

namespace LR
{
std::vector<ModuleBase::Vector3<int>> cross_image_indices(
    const ModuleBase::Vector3<double>& displacement,
    const ModuleBase::Matrix3& lattice, const double cutoff)
{
    const double volume = lattice.Det();
    if (cutoff <= 0.0 || std::abs(volume) < 1e-14)
    {
        throw std::invalid_argument("LR root overlap: invalid lattice or orbital cutoff");
    }
    const ModuleBase::Matrix3 inverse = lattice.Inverse();
    const ModuleBase::Vector3<double> fractional = displacement * inverse;
    const ModuleBase::Vector3<double> reciprocal_x(inverse.e11, inverse.e21, inverse.e31);
    const ModuleBase::Vector3<double> reciprocal_y(inverse.e12, inverse.e22, inverse.e32);
    const ModuleBase::Vector3<double> reciprocal_z(inverse.e13, inverse.e23, inverse.e33);
    const double extent_x = cutoff * reciprocal_x.norm();
    const double extent_y = cutoff * reciprocal_y.norm();
    const double extent_z = cutoff * reciprocal_z.norm();
    const double lower_x = -fractional.x - extent_x;
    const double upper_x = -fractional.x + extent_x;
    const double lower_y = -fractional.y - extent_y;
    const double upper_y = -fractional.y + extent_y;
    const double lower_z = -fractional.z - extent_z;
    const double upper_z = -fractional.z + extent_z;
    const int xmin = static_cast<int>(std::ceil(lower_x));
    const int xmax = static_cast<int>(std::floor(upper_x));
    const int ymin = static_cast<int>(std::ceil(lower_y));
    const int ymax = static_cast<int>(std::floor(upper_y));
    const int zmin = static_cast<int>(std::ceil(lower_z));
    const int zmax = static_cast<int>(std::floor(upper_z));
    std::vector<ModuleBase::Vector3<int>> images;
    for (int x = xmin; x <= xmax; ++x)
    {
        for (int y = ymin; y <= ymax; ++y)
        {
            for (int z = zmin; z <= zmax; ++z)
            {
                const ModuleBase::Vector3<int> image(x, y, z);
                const ModuleBase::Vector3<double> shifted = displacement + image * lattice;
                if (shifted.norm() < cutoff)
                {
                    images.push_back(image);
                }
            }
        }
    }
    return images;
}

// Build ordered old-bra/new-ket S(R) blocks on the existing AO 2D distribution.
// Fourier folding preserves their nonsymmetry; Gamma is simply k = 0.
std::vector<double> cross_ao_overlap(const UnitCell& cell,
    const std::vector<ModuleBase::Vector3<double>>& old_positions,
    const std::vector<double>& cutoff, const TwoCenterIntegrator& integrator,
    const Parallel_Orbitals& pmat)
{
    if (cell.get_npol() != 1)
    {
        throw std::invalid_argument("LR root overlap requires collinear spin");
    }
    ModuleBase::Matrix3 lattice = cell.latvec;
    lattice *= cell.lat0;
    hamilt::HContainer<double> overlap_r(&pmat);
    for (int i = 0; i < cell.nat; ++i)
    {
        const int type_i = cell.iat2it[i];
        for (int j = 0; j < cell.nat; ++j)
        {
            if (pmat.is_invalid_atom_pair(i, j))
            {
                continue;
            }
            const int type_j = cell.iat2it[j];
            const auto displacement = (cell.get_tau(j) - old_positions[i]) * cell.lat0;
            const double radius = cutoff[type_i] + cutoff[type_j];
            const auto indices = cross_image_indices(displacement, lattice, radius);
            for (const auto& index : indices)
            {
                const hamilt::AtomPair<double> pair(i, j, index, &pmat);
                overlap_r.insert_pair(pair);
            }
        }
    }
    overlap_r.allocate(nullptr, true);
    for (int pair_index = 0; pair_index < overlap_r.size_atom_pairs(); ++pair_index)
    {
        auto& pair = overlap_r.get_atom_pair(pair_index);
        const int i = pair.get_atom_i();
        const int j = pair.get_atom_j();
        const auto displacement = (cell.get_tau(j) - old_positions[i]) * cell.lat0;
        for (int ir = 0; ir < pair.get_R_size(); ++ir)
        {
            const auto index = pair.get_R_index(ir);
            const auto shifted = displacement + index * lattice;
            double* const values = pair.get_pointer(ir);
            hamilt::cal_overlap_block(cell, integrator, i, j, pmat, shifted, values);
        }
    }
    // Use the complex Fourier overload even at Gamma: it folds every retained R
    // without fix_gamma(), so the same S(R) can later be folded at nonzero k.
    const int local_count = pmat.get_local_size();
    const std::complex<double> zero(0.0, 0.0);
    std::vector<std::complex<double>> local_overlap(local_count, zero);
    const ModuleBase::Vector3<double> gamma(0.0, 0.0, 0.0);
    const int leading_dimension = pmat.get_row_size();
    if (local_count > 0)
    {
        hamilt::folding_HR(overlap_r, local_overlap.data(), gamma, leading_dimension, 1);
    }
    std::vector<double> overlap(local_count, 0.0);
    for (int index = 0; index < local_count; ++index)
    {
        overlap[index] = local_overlap[index].real();
    }
    return overlap;
}

namespace
{
void check_grid(const Parallel_2D& a, const Parallel_2D& b)
{
#ifdef __MPI
    if (a.blacs_ctxt != b.blacs_ctxt || a.get_block_size() != b.get_block_size())
    {
        throw std::invalid_argument("Root overlap matrices require the same BLACS grid and block size");
    }
#endif
}
}

template <typename T>
std::vector<T> mo_overlap_dist(const std::vector<T>& old_coefficients,
    const std::vector<T>& ao, const std::vector<T>& current_coefficients,
    const Parallel_2D& pc, const Parallel_2D& ps, const Parallel_2D& po, const int first_band)
{
    check_grid(pc, ps);
    check_grid(pc, po);
    const int naos = pc.get_global_row_size();
    const int bands = po.get_global_row_size();
    if (ps.get_global_row_size() != naos || ps.get_global_col_size() != naos
        || bands <= 0 || first_band < 0 || first_band + bands > pc.get_global_col_size()
        || po.get_global_col_size() != bands
        || old_coefficients.size() < pc.get_local_size()
        || current_coefficients.size() < pc.get_local_size() || ao.size() < ps.get_local_size())
    {
        throw std::invalid_argument("Root overlap matrix dimensions do not match their distributions");
    }
    const T zero = T(0);
    const T one = T(1);
    // Align the selected old/new coefficient columns before either GEMM.
    // Only this occupied or virtual block is formed, never the full MO overlap.
    Parallel_2D coefficients;
#ifdef __MPI
    const int block_size = pc.get_block_size();
    coefficients.set(naos, bands, block_size, pc.blacs_ctxt);
#else
    coefficients.set_serial(naos, bands);
#endif
    const int coefficient_count = std::max<int>(1, coefficients.get_local_size());
    const int overlap_count = std::max<int>(1, po.get_local_size());
    std::vector<T> coefficient_block(coefficient_count, zero);
    std::vector<T> intermediate(coefficient_count, zero);
    std::vector<T> overlap(overlap_count, zero);
    const T* const old_data = old_coefficients.empty() ? &zero : old_coefficients.data();
    const T* const new_data = current_coefficients.empty() ? &zero : current_coefficients.data();
    const T* const ao_data = ao.empty() ? &zero : ao.data();
    const char adjoint = std::is_same<T, double>::value ? 'T' : 'C';
#ifdef __MPI
    const int first = 1;
    const int source_column = first_band + 1;
    // GEMR2D's legacy wrapper is non-const but only reads these sources.
    T* const old_source = const_cast<T*>(old_data);
    T* const new_source = const_cast<T*>(new_data);
    int* const source_desc = const_cast<int*>(pc.desc);
    Cpxgemr2d(naos, bands, new_source, first, source_column, source_desc,
        coefficient_block.data(), first, first, coefficients.desc, pc.blacs_ctxt);
    // Empty local blocks join all redistribution and PBLAS calls.
    ScalapackConnector::gemm('N', 'N', naos, bands, naos, one,
        ao_data, first, first, ps.desc, coefficient_block.data(), first, first, coefficients.desc,
        zero, intermediate.data(), first, first, coefficients.desc);
    // Reuse the coefficient buffer once S*C_current has been evaluated.
    Cpxgemr2d(naos, bands, old_source, first, source_column, source_desc,
        coefficient_block.data(), first, first, coefficients.desc, pc.blacs_ctxt);
    ScalapackConnector::gemm(adjoint, 'N', bands, bands, naos, one,
        coefficient_block.data(), first, first, coefficients.desc,
        intermediate.data(), first, first, coefficients.desc,
        zero, overlap.data(), first, first, po.desc);
#else
    const int source_offset = first_band * naos;
    const int selected_count = naos * bands;
    std::copy_n(new_data + source_offset, selected_count, coefficient_block.data());
    BlasConnector::gemm_cm('N', 'N', naos, bands, naos, one,
        ao_data, naos, coefficient_block.data(), naos, zero, intermediate.data(), naos);
    std::copy_n(old_data + source_offset, selected_count, coefficient_block.data());
    BlasConnector::gemm_cm(adjoint, 'N', bands, bands, naos, one,
        coefficient_block.data(), naos, intermediate.data(), naos, zero, overlap.data(), bands);
#endif
    overlap.resize(po.get_local_size());
    return overlap;
}

template <typename T>
std::vector<T> project_reference_dist(const std::vector<T>& previous,
    const std::vector<T>& occupied_overlap, const std::vector<T>& virtual_overlap,
    const Parallel_2D& pi, const Parallel_2D& pv, const Parallel_2D& px)
{
    check_grid(pi, px);
    check_grid(pv, px);
    const int nocc = px.get_global_col_size();
    const int nvirt = px.get_global_row_size();
    if (pi.get_global_row_size() != nocc || pi.get_global_col_size() != nocc
        || pv.get_global_row_size() != nvirt || pv.get_global_col_size() != nvirt
        || previous.size() != px.get_local_size()
        || occupied_overlap.size() != pi.get_local_size()
        || virtual_overlap.size() != pv.get_local_size())
    {
        throw std::invalid_argument("Root projection requires matching windows and local matrix sizes");
    }
    const T zero = T(0);
    const T one = T(1);
    const int local_count = std::max<int>(1, px.get_local_size());
    std::vector<T> intermediate(local_count, zero);
    std::vector<T> projected(local_count, zero);
    const T* const old_data = previous.empty() ? &zero : previous.data();
    const T* const occupied_data = occupied_overlap.empty() ? &zero : occupied_overlap.data();
    const T* const virtual_data = virtual_overlap.empty() ? &zero : virtual_overlap.data();
    const char adjoint = std::is_same<T, double>::value ? 'T' : 'C';
#ifdef __MPI
    const int first = 1;
    ScalapackConnector::gemm(adjoint, 'N', nvirt, nocc, nvirt, one,
        virtual_data, first, first, pv.desc,
        old_data, first, first, px.desc, zero,
        intermediate.data(), first, first, px.desc);
    ScalapackConnector::gemm('N', 'N', nvirt, nocc, nocc, one,
        intermediate.data(), first, first, px.desc,
        occupied_data, first, first, pi.desc, zero,
        projected.data(), first, first, px.desc);
#else
    BlasConnector::gemm_cm(adjoint, 'N', nvirt, nocc, nvirt, one,
        virtual_data, nvirt, old_data, nvirt, zero, intermediate.data(), nvirt);
    BlasConnector::gemm_cm('N', 'N', nvirt, nocc, nocc, one,
        intermediate.data(), nvirt, occupied_data, nocc, zero, projected.data(), nvirt);
#endif
    projected.resize(px.get_local_size());
    return projected;
}

template std::vector<double> mo_overlap_dist(const std::vector<double>&,
    const std::vector<double>&, const std::vector<double>&,
    const Parallel_2D&, const Parallel_2D&, const Parallel_2D&, int);
template std::vector<std::complex<double>> mo_overlap_dist(
    const std::vector<std::complex<double>>&, const std::vector<std::complex<double>>&,
    const std::vector<std::complex<double>>&, const Parallel_2D&, const Parallel_2D&, const Parallel_2D&, int);
template std::vector<double> project_reference_dist(const std::vector<double>&,
    const std::vector<double>&, const std::vector<double>&, const Parallel_2D&, const Parallel_2D&, const Parallel_2D&);
template std::vector<std::complex<double>> project_reference_dist(
    const std::vector<std::complex<double>>&, const std::vector<std::complex<double>>&,
    const std::vector<std::complex<double>>&, const Parallel_2D&, const Parallel_2D&, const Parallel_2D&);

}
