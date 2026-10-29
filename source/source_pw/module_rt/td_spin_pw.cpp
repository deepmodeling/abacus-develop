#include "source_pw/module_rt/td_spin_pw.h"

#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_pw/module_proj/onsite_proj.h"
#include "source_pw/module_rt/td_pw.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>

namespace pw
{
namespace
{
double image_dist(const ModuleBase::Vector3<double>& delta, const ModuleBase::Matrix3& lattice, const std::array<double, 3>& inverse_norm)
{
    ModuleBase::Vector3<double> nearest;
    for (int a = 0; a < 3; ++a)
    {
        nearest[a] = delta[a] - std::round(delta[a]);
    }
    double best = (nearest * lattice).norm2();
    const double radius = std::sqrt(best);
    int lower[3];
    int upper[3];
    for (int a = 0; a < 3; ++a)
    {
        // Every better image lies inside these bounds, even for highly skew cells.
        const double extent = radius * inverse_norm[a] + 1e-12;
        lower[a] = static_cast<int>(std::ceil(-delta[a] - extent));
        upper[a] = static_cast<int>(std::floor(-delta[a] + extent));
    }
    for (int i = lower[0]; i <= upper[0]; ++i)
    {
        for (int j = lower[1]; j <= upper[1]; ++j)
        {
            for (int k = lower[2]; k <= upper[2]; ++k)
            {
                const ModuleBase::Vector3<double> image(delta.x + i, delta.y + j, delta.z + k);
                const double distance_sq = (image * lattice).norm2();
                best = std::min(best, distance_sq);
            }
        }
    }
    return best;
}
} // namespace

std::vector<int> td_voronoi_pw(const UnitCell& cell, const ModulePW::PW_Basis& grid)
{
    const ModuleBase::Matrix3 inverse = cell.latvec.Inverse();
    const std::array<double, 3> inverse_norm{
        {std::sqrt(inverse.e11 * inverse.e11 + inverse.e21 * inverse.e21 + inverse.e31 * inverse.e31),
         std::sqrt(inverse.e12 * inverse.e12 + inverse.e22 * inverse.e22 + inverse.e32 * inverse.e32),
         std::sqrt(inverse.e13 * inverse.e13 + inverse.e23 * inverse.e23 + inverse.e33 * inverse.e33)}};
    std::vector<ModuleBase::Vector3<double>> atoms;
    for (int it = 0; it < cell.ntype; ++it)
    {
        for (int ia = 0; ia < cell.atoms[it].na; ++ia)
        {
            atoms.push_back(cell.atoms[it].taud[ia]);
        }
    }
    std::vector<int> labels(grid.nrxx, 0);
    for (int ir = 0; ir < grid.nrxx; ++ir)
    {
        const int ix = ir / (grid.ny * grid.nplane);
        const int iy = ir / grid.nplane % grid.ny;
        const int iz = ir % grid.nplane + grid.startz_current;
        const ModuleBase::Vector3<double> point(double(ix) / grid.nx, double(iy) / grid.ny, double(iz) / grid.nz);
        double best = std::numeric_limits<double>::infinity();
        for (int atom = 0; atom < cell.nat; ++atom)
        {
            const ModuleBase::Vector3<double> delta = point - atoms[atom];
            const double distance = image_dist(delta, cell.latvec, inverse_norm);
            const double tolerance = 64 * std::numeric_limits<double>::epsilon() * std::max(1.0, distance);
            if (atom == 0 || distance < best - tolerance)
            {
                best = distance;
                labels[ir] = atom;
            }
        }
    }
    return labels;
}

template <typename Real, typename Device>
void TDSpinPW<Real, Device>::initialize(const UnitCell& cell,
                                        const ModulePW::PW_Basis& grid,
                                        const ModulePW::PW_Basis_K& basis,
                                        const K_Vectors& kv,
                                        Structure_Factor& sf,
                                        pseudopot_cell_vnl& pp,
                                        const std::string& mag_partition,
                                        const std::string& orbital_dir,
                                        double onsite_radius,
                                        bool calc_torque)
{
    if (cell_ != nullptr || cell.get_npol() != 2)
    {
        ModuleBase::WARNING_QUIT("TDSpinPW", "Initialize once with a spinor basis.");
    }
    if (mag_partition != "none" && mag_partition != "voronoi" && mag_partition != "onsite")
    {
        ModuleBase::WARNING_QUIT("TDSpinPW", "Unknown magnetic partition.");
    }
    cell_ = &cell;
    grid_ = &grid;
    basis_ = &basis;
    kv_ = &kv;
    sf_ = &sf;
    pp_ = &pp;
    mag_partition_ = mag_partition;
    orbital_dir_ = orbital_dir;
    onsite_radius_ = onsite_radius;
    calc_torque_ = calc_torque;
}

template <typename Real, typename Device>
void TDSpinPW<Real, Device>::sync_positions()
{
    bool changed = pos_cache_.size() != static_cast<size_t>(cell_->nat);
    int atom = 0;
    for (int it = 0; !changed && it < cell_->ntype; ++it)
    {
        for (int ia = 0; !changed && ia < cell_->atoms[it].na; ++ia, ++atom)
        {
            changed = !(pos_cache_[atom] == cell_->atoms[it].taud[ia]);
        }
    }
    if (!changed)
    {
        return;
    }
    // Atomic phases and nearest-atom labels change even at fixed grid, k point and vector potential.
    labels_.clear();
    if (onsite_projector_)
    {
        onsite_projector_->invalidate_proj();
    }
    pos_cache_.clear();
    for (int it = 0; it < cell_->ntype; ++it)
    {
        for (int ia = 0; ia < cell_->atoms[it].na; ++ia)
        {
            pos_cache_.push_back(cell_->atoms[it].taud[ia]);
        }
    }
}

template <typename Real, typename Device>
void TDSpinPW<Real, Device>::density(const elecstate::ElecState& state, TDSpinResultPW* result)
{
    const UnitCell& cell = *cell_;
    const ModulePW::PW_Basis& grid = *grid_;
    const bool voronoi = mag_partition_ == "voronoi";
    if (voronoi && labels_.size() != static_cast<size_t>(grid.nrxx))
    {
        labels_ = td_voronoi_pw(cell, grid);
    }
    const double volume = cell.omega / grid.nxyz;
    const double* const* rho = state.charge->rho;
    const ModuleBase::matrix& potential = state.pot->get_eff_v();
    double totals[8] = {};
    if (voronoi)
    {
        result->atom_mag.assign(3 * cell.nat, 0.0);
    }
    for (int ir = 0; ir < grid.nrxx; ++ir)
    {
        double norm = 0.0;
        for (int a = 0; a < 3; ++a)
        {
            const double m = rho[a + 1][ir];
            totals[a] += m * volume;
            norm += m * m;
            if (voronoi)
            {
                result->atom_mag[3 * labels_[ir] + a] += m * volume;
            }
            if (calc_torque_)
            {
                // The effective spin field is in Ry; the angular-momentum torque here uses Hartree time.
                const int b = (a + 1) % 3;
                const int c = (a + 2) % 3;
                totals[5 + a] += 0.5 * volume * (potential(b + 1, ir) * rho[c + 1][ir] - potential(c + 1, ir) * rho[b + 1][ir]);
            }
        }
        totals[3] += std::sqrt(norm) * volume;
        totals[4] += rho[0][ir] * volume;
    }
    // Density already contains all k pools. Reduce only the grid partition.
    Parallel_Reduce::reduce_pool(totals, 8);
    if (voronoi)
    {
        const int count = result->atom_mag.size();
        Parallel_Reduce::reduce_pool(result->atom_mag.data(), count);
    }
    for (int a = 0; a < 3; ++a)
    {
        result->mag[a] = totals[a];
        result->xc[a] = totals[5 + a];
    }
    result->mag_abs = totals[3];
    result->electrons = totals[4];
}

template <typename Real, typename Device>
void TDSpinPW<Real, Device>::soc(const psi::Psi<Complex, Device>& psi,
                                 const ModuleBase::matrix& weights,
                                 const ModuleBase::Vector3<double>& a,
                                 TDSpinResultPW* result)
{
    const UnitCell& cell = *cell_;
    const ModulePW::PW_Basis_K& basis = *basis_;
    pseudopot_cell_vnl& pp = *pp_;
    result->atom_soc.assign(3 * cell.nat, 0.0);
    if (pp.nkb == 0)
    {
        return;
    }
    using CtDevice = typename ct::PsiToContainer<Device>::type;
    const ct::DeviceType device = ct::DeviceTypeToEnum<CtDevice>::value;
    const int nbands = psi.get_nbands();
    // Projection kernels count spin-component columns, not physical bands.
    const int ncols = 2 * nbands;
    const int ld_psi = psi.get_nbasis();
    const std::int64_t size = std::max<std::int64_t>(1, static_cast<std::int64_t>(pp.nkb) * basis.npwk_max);
    // Empty local PW partitions still need a correctly typed buffer before projector reductions.
    if (beta_.NumElements() < size || beta_.data_type() != ct::DataTypeToEnum<Complex>::value || beta_.device_type() != device)
    {
        beta_ = ct::Tensor(ct::DataTypeToEnum<Complex>::value, device, {size});
    }
    Complex* beta = beta_.template data<Complex>();
    const size_t count = static_cast<size_t>(pp.nkb) * ncols;
    projections_.resize(count);
    contracted_.resize(count);
    Device* context = nullptr;
    for (int ik = 0; ik < basis.nks; ++ik)
    {
        psi.fix_k(ik);
        const int npw = basis.npwk[ik];
        if (npw > 0)
        {
            pp.getvnl(context, cell, ik, a, beta);
        }
        nonlocal_.project_contract(cell, pp, basis, 0, ncols, ld_psi, 2, npw, beta, psi.get_pointer());
        if (basis.poolrank != 0)
        {
            // B and DB already contain the full pool sum; count each k-point contribution only once.
            continue;
        }
        base_device::memory::synchronize_memory_op<Complex, base_device::DEVICE_CPU, Device>()(projections_.data(),
                                                                                               nonlocal_.projections(),
                                                                                               count);
        base_device::memory::synchronize_memory_op<Complex, base_device::DEVICE_CPU, Device>()(contracted_.data(),
                                                                                               nonlocal_.contracted(),
                                                                                               count);
        int offset = 0;
        int atom = 0;
        const Complex imaginary(0, 1);
        for (int it = 0; it < cell.ntype; ++it)
        {
            const int nh = cell.atoms[it].ncpp.nh;
            for (int ia = 0; ia < cell.atoms[it].na; ++ia, ++atom)
            {
                for (int p = offset; p < offset + nh; ++p)
                {
                    for (int band = 0; band < nbands; ++band)
                    {
                        const Complex up = std::conj(projections_[(2 * band) * pp.nkb + p]);
                        const Complex down = std::conj(projections_[(2 * band + 1) * pp.nkb + p]);
                        const Complex du = contracted_[p * ncols + 2 * band];
                        const Complex dd = contracted_[p * ncols + 2 * band + 1];
                        const double weight = 0.5 * weights(ik, band);
                        result->atom_soc[3 * atom] += weight * std::imag(up * dd + down * du);
                        result->atom_soc[3 * atom + 1] += weight * std::imag(-imaginary * up * dd + imaginary * down * du);
                        result->atom_soc[3 * atom + 2] += weight * std::imag(up * du - down * dd);
                    }
                }
                offset += nh;
            }
        }
    }
    const int values = result->atom_soc.size();
    Parallel_Reduce::reduce_all(result->atom_soc.data(), values);
    for (int atom = 0; atom < cell.nat; ++atom)
    {
        for (int a = 0; a < 3; ++a)
        {
            result->soc[a] += result->atom_soc[3 * atom + a];
        }
    }
}

template <typename Real, typename Device>
void TDSpinPW<Real, Device>::onsite(const psi::Psi<Complex, Device>& psi,
                                    const elecstate::ElecState& state,
                                    const ModuleBase::Vector3<double>& a,
                                    TDSpinResultPW* result)
{
    const UnitCell& cell = *cell_;
    const ModulePW::PW_Basis_K& basis = *basis_;
    const K_Vectors& kv = *kv_;
    Structure_Factor& sf = *sf_;
    // Onsite owns its atomic result; the independently integrated total moment remains unchanged.
    result->atom_mag.assign(3 * cell.nat, 0.0);
    const double dq = pp_->table_dq();
    using OnsiteValue = std::complex<double>;
    const int nbands = psi.get_nbands();
    const int ld_psi = psi.get_nbasis();
    if (!onsite_psi_)
    {
        const int npw = basis.npwk[0];
        onsite_psi_.reset(new psi::Psi<OnsiteValue, Device>(1, nbands, ld_psi, npw, true));
        onsite_projector_.reset(new projectors::OnsiteProjector<double, Device>());
        q_bound_ = td_momentum_bound(basis, cell.tpiba);
        const int nq = static_cast<int>(std::ceil(q_bound_ / dq)) + 4;
        onsite_projector_->init(orbital_dir_, &cell, *onsite_psi_, kv, basis, sf, onsite_radius_, nq, dq, state.wg, state.ekb);
    }
    projectors::OnsiteProjector<double, Device>* projector = onsite_projector_.get();
    if (onsite_psi_->get_nbands() != nbands || onsite_psi_->get_nbasis() != ld_psi)
    {
        ModuleBase::WARNING_QUIT("TDSpinPW", "The onsite wavefunction layout changed after initialization.");
    }
    const double bound = q_bound_ + a.norm();
    projector->ensure_range(bound, dq);
    const int nkb = projector->get_tot_nproj();
    const int npwx = ld_psi / 2;
    const size_t size = static_cast<size_t>(nbands) * ld_psi;
    for (int ik = 0; ik < basis.nks; ++ik)
    {
        psi.fix_k(ik);
        OnsiteValue* converted = onsite_psi_->get_pointer();
        if (size > 0)
        {
            // Onsite projections use double precision, so float spinors require a value conversion, not a pointer cast.
            base_device::memory::cast_memory_op<OnsiteValue, Complex, Device, Device>()(converted, psi.get_pointer(), size);
        }
        const int columns = 2 * nbands;
        projector->overlap_td(ik, columns, converted, npwx, a);
        if (basis.poolrank != 0)
        {
            // Overlaps are replicated after the pool reduction; count each occupied orbital only once.
            continue;
        }
        const OnsiteValue* b = projector->get_h_becp();
        int offset = 0;
        for (int atom = 0; atom < cell.nat; ++atom)
        {
            const int nh = projector->get_nh(atom);
            for (int p = offset; p < offset + nh; ++p)
            {
                for (int band = 0; band < nbands; ++band)
                {
                    const OnsiteValue up = b[2 * band * nkb + p];
                    const OnsiteValue down = b[(2 * band + 1) * nkb + p];
                    const double weight = state.wg(ik, band);
                    const OnsiteValue cross = std::conj(up) * down;
                    result->atom_mag[3 * atom] += 2 * weight * cross.real();
                    result->atom_mag[3 * atom + 1] += 2 * weight * cross.imag();
                    result->atom_mag[3 * atom + 2] += weight * (std::norm(up) - std::norm(down));
                }
            }
            offset += nh;
        }
    }
    const int count = result->atom_mag.size();
    Parallel_Reduce::reduce_all(result->atom_mag.data(), count);
}

template <typename Real, typename Device>
TDSpinResultPW TDSpinPW<Real, Device>::calculate(const psi::Psi<Complex, Device>& psi,
                                                 const elecstate::ElecState& state,
                                                 const ModuleBase::Vector3<double>& a)
{
    ModuleBase::timer::start("TDSpinPW", "calculate");
    if (cell_ == nullptr)
    {
        ModuleBase::WARNING_QUIT("TDSpinPW", "Initialize before sampling observables.");
    }
    if (psi.get_npol() != 2 || psi.get_nbasis() != 2 * basis_->npwk_max)
    {
        ModuleBase::WARNING_QUIT("TDSpinPW", "Magnetic observables require the initialized spinor layout.");
    }
    sync_positions();
    const int old_k = psi.get_current_k();
    TDSpinResultPW result;
    density(state, &result);
    if (mag_partition_ == "onsite")
    {
        onsite(psi, state, a, &result);
    }
    if (calc_torque_)
    {
        soc(psi, state.wg, a, &result);
    }
    // Observable evaluation must preserve the caller's k-point cursor.
    psi.fix_k(old_k);
    ModuleBase::timer::end("TDSpinPW", "calculate");
    return result;
}

template class TDSpinPW<float, base_device::DEVICE_CPU>;
template class TDSpinPW<double, base_device::DEVICE_CPU>;
#if defined(__CUDA) || defined(__ROCM)
template class TDSpinPW<float, base_device::DEVICE_GPU>;
template class TDSpinPW<double, base_device::DEVICE_GPU>;
#endif
} // namespace pw
