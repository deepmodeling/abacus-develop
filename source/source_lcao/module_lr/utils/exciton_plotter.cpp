#include "exciton_plotter.h"

#include "source_base/constants.h"
#include "source_base/matrix3.h"
#include "source_base/module_external/blas_connector.h"
#include "source_base/parallel_global.h"
#include "source_base/ylm.h"
#include "source_cell/atom_spec.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <utility>

namespace
{
using Vec3 = ModuleBase::Vector3<double>;

struct SliceGeometry
{
    std::string plane;
    double slice_pos = 0.0;
    Vec3 a_bohr, b_bohr, c_bohr;
    Vec3 u_vec, v_vec, perp_offset;
    int u_axis = 0, v_axis = 1;
    int nk_u = 1, nk_v = 1;
    int u_start_cells = 0, u_end_cells = 1;
    int v_start_cells = 0, v_end_cells = 1;
    double u_start = 0.0, u_end = 1.0;
    double v_start = 0.0, v_end = 1.0;
    int res = 1, nu = 2, nv = 2;
};

SliceGeometry make_slice_geometry(const UnitCell& ucell,
                                  const K_Vectors& kv,
                                  const std::string& plane,
                                  const double slice_pos,
                                  const int npoints,
                                  const std::vector<int>& range)
{
    if (range.size() != 4 || range[0] >= range[1] || range[2] >= range[3])
    {
        ModuleBase::WARNING_QUIT("ExcitonPlotter", "Invalid slice range: expected ustart < uend and vstart < vend.");
    }

    SliceGeometry geom;
    geom.plane = plane;
    geom.slice_pos = slice_pos;
    geom.a_bohr = ucell.a1 * ucell.lat0;
    geom.b_bohr = ucell.a2 * ucell.lat0;
    geom.c_bohr = ucell.a3 * ucell.lat0;

    if (plane == "ab")
    {
        geom.u_vec = geom.a_bohr;
        geom.v_vec = geom.b_bohr;
        geom.perp_offset = geom.c_bohr * (slice_pos / geom.c_bohr.norm());
        geom.u_axis = 0;
        geom.v_axis = 1;
        geom.nk_u = std::max(1, kv.nmp[0]);
        geom.nk_v = std::max(1, kv.nmp[1]);
    }
    else if (plane == "bc")
    {
        geom.u_vec = geom.b_bohr;
        geom.v_vec = geom.c_bohr;
        geom.perp_offset = geom.a_bohr * (slice_pos / geom.a_bohr.norm());
        geom.u_axis = 1;
        geom.v_axis = 2;
        geom.nk_u = std::max(1, kv.nmp[1]);
        geom.nk_v = std::max(1, kv.nmp[2]);
    }
    else if (plane == "ca")
    {
        geom.u_vec = geom.c_bohr;
        geom.v_vec = geom.a_bohr;
        geom.perp_offset = geom.b_bohr * (slice_pos / geom.b_bohr.norm());
        geom.u_axis = 2;
        geom.v_axis = 0;
        geom.nk_u = std::max(1, kv.nmp[2]);
        geom.nk_v = std::max(1, kv.nmp[0]);
    }
    else
    {
        ModuleBase::WARNING_QUIT("ExcitonPlotter", "Unknown slice plane: " + plane + ". Use ab, bc, or ca.");
    }

    geom.u_start_cells = range[0];
    geom.u_end_cells = range[1];
    geom.v_start_cells = range[2];
    geom.v_end_cells = range[3];
    geom.u_start = static_cast<double>(geom.u_start_cells);
    geom.u_end = static_cast<double>(geom.u_end_cells);
    geom.v_start = static_cast<double>(geom.v_start_cells);
    geom.v_end = static_cast<double>(geom.v_end_cells);

    const int total_u = geom.u_end_cells - geom.u_start_cells;
    const int total_v = geom.v_end_cells - geom.v_start_cells;
    const int max_cells = std::max(total_u, total_v);
    geom.res = std::max(1, npoints / std::max(1, max_cells));
    geom.nu = total_u * geom.res + 1;
    geom.nv = total_v * geom.res + 1;
    return geom;
}

double bloch_phase_arg(const Vec3& kvec_d, const int n1, const int n2, const int n3)
{
    return ModuleBase::TWO_PI * (kvec_d.x * n1 + kvec_d.y * n2 + kvec_d.z * n3);
}

double bloch_phase_arg_axis(const Vec3& kvec_d, const int axis, const int ncell)
{
    return ModuleBase::TWO_PI * kvec_d[axis] * ncell;
}

std::complex<double> to_complex(const double value)
{
    return std::complex<double>(value, 0.0);
}

std::complex<double> to_complex(const std::complex<double>& value)
{
    return value;
}

char adjoint_flag(const double*)
{
    return 'T';
}

char adjoint_flag(const std::complex<double>*)
{
    return 'C';
}

template <typename TK>
double density_from_dmk(const TK* dmk, const std::vector<std::complex<double>>& phi, const int naos)
{
    std::complex<double> rho(0.0, 0.0);
    for (int mu = 0; mu < naos; ++mu)
    {
        std::complex<double> row_sum(0.0, 0.0);
        for (int nu = 0; nu < naos; ++nu)
        {
            // cal_effective_dmk_* stores D^T for DensityMatrix; recover D(mu, nu)
            // when contracting the matrix directly in real space.
            row_sum += to_complex(dmk[nu + mu * naos]) * std::conj(phi[nu]);
        }
        rho += phi[mu] * row_sum;
    }
    return std::max(0.0, rho.real());
}

void write_slice_data(const UnitCell& ucell,
                      const SliceGeometry& geom,
                      const std::vector<double>& density,
                      const std::string& filename,
                      const int istate,
                      const double energy_ry,
                      const std::string& density_kind,
                      const bool has_hole,
                      const std::array<double, 3>& r_h_fix)
{
    if (density.size() != static_cast<std::size_t>(geom.nu * geom.nv))
    {
        ModuleBase::WARNING_QUIT("ExcitonPlotter", "Invalid slice density dimensions.");
    }

    std::ofstream ofs(filename);
    if (!ofs)
    {
        ModuleBase::WARNING_QUIT("ExcitonPlotter", "Cannot open slice file: " + filename);
    }

    ofs << std::setprecision(12);
    ofs << "# ABACUS_EXCITON_SLICE 1\n";
    ofs << "# state " << istate << " energy_Ry " << energy_ry << " density_kind " << density_kind << "\n";
    const std::string fixed_particle
        = density_kind == "conditional_hole" ? "electron" : (density_kind == "conditional_elec" ? "hole" : "none");
    ofs << "# fixed_particle " << fixed_particle;
    if (has_hole)
    {
        ofs << " position_bohr " << r_h_fix[0] << " " << r_h_fix[1] << " " << r_h_fix[2];
    }
    ofs << "\n";
    ofs << "# plane " << geom.plane << " slice_pos_bohr " << geom.slice_pos << "\n";
    ofs << "# grid " << geom.nu << " " << geom.nv << " u_range " << geom.u_start << " " << geom.u_end << " v_range "
        << geom.v_start << " " << geom.v_end << "\n";
    ofs << "# bvk " << geom.nk_u << " " << geom.nk_v << "\n";
    ofs << "# lattice_bohr\n";
    ofs << "# a " << geom.a_bohr.x << " " << geom.a_bohr.y << " " << geom.a_bohr.z << "\n";
    ofs << "# b " << geom.b_bohr.x << " " << geom.b_bohr.y << " " << geom.b_bohr.z << "\n";
    ofs << "# c " << geom.c_bohr.x << " " << geom.c_bohr.y << " " << geom.c_bohr.z << "\n";
    ofs << "# axes_bohr\n";
    ofs << "# u " << geom.u_vec.x << " " << geom.u_vec.y << " " << geom.u_vec.z << "\n";
    ofs << "# v " << geom.v_vec.x << " " << geom.v_vec.y << " " << geom.v_vec.z << "\n";
    ofs << "# atoms " << ucell.nat << "\n";
    for (int iat = 0; iat < ucell.nat; ++iat)
    {
        const int it = ucell.iat2it[iat];
        const Vec3 tau = ucell.get_tau(iat) * ucell.lat0; // tau in Bohr
        ofs << "# " << ucell.atoms[it].label << " " << tau.x << " " << tau.y << " " << tau.z << "\n";
    }
    ofs << "# data\n";
    ofs << std::scientific << std::setprecision(8);
    for (int iu = 0; iu < geom.nu; ++iu)
    {
        for (int iv = 0; iv < geom.nv; ++iv)
        {
            ofs << density[iu * geom.nv + iv] << (iv + 1 < geom.nv ? " " : "");
        }
        ofs << "\n";
    }
}
} // namespace

// ============================================================================
// OrbitalEvaluator non-template method implementations
// ============================================================================

void OrbitalEvaluator::init(const LCAO_Orbitals& orb, const UnitCell& ucell)
{
    atypes_.resize(ucell.ntype);
    for (int it = 0; it < ucell.ntype; ++it)
    {
        const Numerical_Orbital& no = orb.Phi[it];
        const Atom& atom = ucell.atoms[it];
        auto& td = atypes_[it];

        td.nw = atom.nw;
        td.nwl = atom.nwl;
        td.rcut = no.getRcut();
        td.dr_uniform = no.PhiLN(0, 0).dr_uniform;
        td.radial_blocks.reserve(atom.nw);

        for (int iw = 0; iw < atom.nw; ++iw)
        {
            if (atom.iw2_new[iw])
            {
                const int l = atom.iw2l[iw];
                const int n = atom.iw2n[iw];
                const auto& phi_ln = no.PhiLN(l, n);

                AtomTypeData::RadialBlock block;
                block.begin_iw = iw;
                block.size = 2 * l + 1;
                block.ylm_begin = atom.iw2_ylm[iw];
                block.psi_uniform = phi_ln.psi_uniform.data();
                block.dpsi_uniform = phi_ln.dpsi_uniform.data();
                td.radial_blocks.push_back(block);
            }
        }
    }
}

/// @brief Replicates the cubic Hermite interpolation from GintAtom::set_phi()
void OrbitalEvaluator::eval_phi(const int it, const ModuleBase::Vector3<double>& dr, double* phi_out) const
{
    const auto& td = atypes_[it];
    const double dist = dr.norm() < 1e-9 ? 1e-9 : dr.norm();
    if (dist > td.rcut)
    {
        ModuleBase::GlobalFunc::ZEROS(phi_out, td.nw);
        return;
    }

    std::vector<double> ylma;
    ModuleBase::Ylm::sph_harm(td.nwl, dr.x / dist, dr.y / dist, dr.z / dist, ylma);

    const double position = dist / td.dr_uniform;
    const int ip = static_cast<int>(position);
    const double dx = position - ip;
    const double dx2 = dx * dx;
    const double dx3 = dx2 * dx;
    const double c3 = 3.0 * dx2 - 2.0 * dx3;
    const double c1 = 1.0 - c3;
    const double c2 = (dx - 2.0 * dx2 + dx3) * td.dr_uniform;
    const double c4 = (dx3 - dx2) * td.dr_uniform;

    const auto* blocks = td.radial_blocks.data();
    const int num_blocks = td.radial_blocks.size();
    for (int ib = 0; ib < num_blocks; ++ib)
    {
        const auto& block = blocks[ib];
        const double* psi_uniform = block.psi_uniform;
        const double* dpsi_uniform = block.dpsi_uniform;
        const double psi
            = c1 * psi_uniform[ip] + c2 * dpsi_uniform[ip] + c3 * psi_uniform[ip + 1] + c4 * dpsi_uniform[ip + 1];

        const int begin_iw = block.begin_iw;
        const int end_iw = begin_iw + block.size;
        int idx_lm = block.ylm_begin;
        for (int iw = begin_iw; iw < end_iw; ++iw, ++idx_lm)
        {
            phi_out[iw] = psi * ylma[idx_lm];
        }
    }
}

// =============================================================================
// Evaluate the Bloch-summed numerical atomic orbital at position r:
//   phi^B_mu(r, k) = Sum_R exp(i * 2pi * kvec_d · (n1,n2,n3)) * phi_mu(r - tau - R)
// where kvec_d is the direct-coordinate k-vector. R runs over all lattice vectors
// within the orbital cutoff radius R = n1*a + n2*b + n3*c.
// =============================================================================
void OrbitalEvaluator::eval_phi_all_bloch(const ModuleBase::Vector3<double>& r,
                                          const UnitCell& ucell,
                                          std::complex<double>* phi_bloch,
                                          const ModuleBase::Vector3<double>& kvec_d) const
{
    // Cell vectors in Bohr for real-space distance checks
    const double lat0 = ucell.lat0;
    ModuleBase::Vector3<double> a = ucell.a1 * lat0;
    ModuleBase::Vector3<double> b = ucell.a2 * lat0;
    ModuleBase::Vector3<double> c = ucell.a3 * lat0;

    // Determine how many lattice translations to check in each direction.
    // nR = ceil(rcut / min(|a|,|b|,|c|)) + 1 ensures all overlapping orbitals are captured.
    double max_rcut = 0.0;
    for (const auto& td: atypes_)
        max_rcut = std::max(max_rcut, td.rcut);
    const int nR = static_cast<int>(std::ceil(max_rcut / std::min({a.norm(), b.norm(), c.norm()}))) + 1;

    std::vector<double> phi_tmp;
    int iw_global = 0;

    // For each atom, iterate over neighboring cell images
    for (int iat = 0; iat < ucell.nat; ++iat)
    {
        const int it = ucell.iat2it[iat];
        const int ia = ucell.iat2ia[iat];
        const ModuleBase::Vector3<double> tau = ucell.atoms[it].tau[ia] * lat0;
        const int nw = atypes_[it].nw;
        const double rcut = atypes_[it].rcut;
        phi_tmp.resize(nw);

        // Triple loop over lattice vectors R = n1*a + n2*b + n3*c
        for (int n1 = -nR; n1 <= nR; ++n1)
        {
            for (int n2 = -nR; n2 <= nR; ++n2)
            {
                for (int n3 = -nR; n3 <= nR; ++n3)
                {
                    ModuleBase::Vector3<double> R = a * (double)n1 + b * (double)n2 + c * (double)n3;
                    ModuleBase::Vector3<double> dr = r - tau - R;
                    if (dr.norm() > rcut)
                        continue; // skip if orbital does not reach

                    double kdotR = bloch_phase_arg(kvec_d, n1, n2, n3);
                    std::complex<double> exp_ikR(std::cos(kdotR), std::sin(kdotR));

                    // Evaluate real-space orbital and accumulate with Bloch phase
                    eval_phi(it, dr, phi_tmp.data());
                    for (int iw = 0; iw < nw; ++iw)
                        phi_bloch[iw_global + iw] += exp_ikR * phi_tmp[iw];
                }
            }
        }
        iw_global += nw;
    }
}

namespace LR_Util
{

template <typename T>
std::vector<container::Tensor> ExcitonPlotter<T>::cal_effective_dmk_hole(const int istate)
{
    ModuleBase::TITLE("ExcitonPlotter", "cal_effective_dmk_hole");
    assert(this->nspin_x == 1);
    const int offset_b = istate * this->ldim;
    const int nocc = this->nocc[0];
    const int nvirt = this->nvirt[0];
    const auto data_type = container::DataTypeToEnum<T>::value;
    std::vector<container::Tensor> result(this->nk,
                                          container::Tensor(data_type, DEV::CpuDevice, {this->naos, this->naos}));
    const T alpha(1.0);
    const T beta(0.0);
    char adjoint = adjoint_flag(static_cast<const T*>(nullptr));
    char normal = 'N';
    char transpose = 'T';

    for (int ik = 0; ik < this->nk; ++ik)
    {
        this->psi_ks_vec[0].fix_k(ik);
        T* const dmk = result[ik].data<T>();
        ModuleBase::GlobalFunc::ZEROS(dmk, this->naos * this->naos);
        const T* const xk = this->X + offset_b + ik * (this->ldim / this->nk);

        // Step 1: M_k = X_k^H * X_k → nocc x nocc Hermitian
        container::Tensor M_k(data_type, DEV::CpuDevice, {nocc, nocc});
        T* const M_k_data = M_k.data<T>();
        ModuleBase::GlobalFunc::ZEROS(M_k_data, nocc * nocc);
        BlasConnector::gemm_cm(adjoint, normal, nocc, nocc, nvirt,
                               alpha, xk, nvirt, xk, nvirt,
                               beta, M_k_data, nocc);

        // Step 2: temp = M_k * C_occ^H  →  nocc x naos
        container::Tensor temp(data_type, DEV::CpuDevice, {nocc, this->naos});
        T* const temp_data = temp.data<T>();
        BlasConnector::gemm_cm(normal, adjoint, nocc, this->naos, nocc,
                               alpha, M_k_data, nocc,
                               this->psi_ks_vec[0].get_pointer(0), this->naos,
                               beta, temp_data, nocc);

        // Step 3: dmk = temp^T * C_occ^T = (C_occ * M_k * C_occ^H)^T.
        // DensityMatrix expects the AO indices in this transposed order.
        BlasConnector::gemm_cm(transpose, transpose, this->naos, this->naos, nocc,
                               alpha, temp_data, nocc,
                               this->psi_ks_vec[0].get_pointer(0), this->naos,
                               beta, dmk, this->naos);
    }
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "cal effective dmk hole");
    return result;
}

template <typename T>
std::vector<container::Tensor> ExcitonPlotter<T>::cal_effective_dmk_elec(const int istate)
{
    ModuleBase::TITLE("ExcitonPlotter", "cal_effective_dmk_elec");
    assert(this->nspin_x == 1);
    const int offset_b = istate * this->ldim;
    const int nocc = this->nocc[0];
    const int nvirt = this->nvirt[0];
    const auto data_type = container::DataTypeToEnum<T>::value;
    std::vector<container::Tensor> result(this->nk,
                                          container::Tensor(data_type, DEV::CpuDevice, {this->naos, this->naos}));
    const T alpha(1.0);
    const T beta(0.0);
    char adjoint = adjoint_flag(static_cast<const T*>(nullptr));
    char normal = 'N';
    char transpose = 'T';

    for (int ik = 0; ik < this->nk; ++ik)
    {
        this->psi_ks_vec[0].fix_k(ik);
        T* const dmk = result[ik].data<T>();
        ModuleBase::GlobalFunc::ZEROS(dmk, this->naos * this->naos);
        const T* const xk = this->X + offset_b + ik * (this->ldim / this->nk);

        // Step 1: N_k = X_k * X_k^H  →  nvirt x nvirt Hermitian
        container::Tensor N_k(data_type, DEV::CpuDevice, {nvirt, nvirt});
        T* const N_k_data = N_k.data<T>();
        ModuleBase::GlobalFunc::ZEROS(N_k_data, nvirt * nvirt);
        BlasConnector::gemm_cm(normal, adjoint, nvirt, nvirt, nocc,
                               alpha, xk, nvirt, xk, nvirt,
                               beta, N_k_data, nvirt);

        // Step 2: temp = N_k * C_virt^H  →  nvirt x naos
        container::Tensor temp(data_type, DEV::CpuDevice, {nvirt, this->naos});
        T* const temp_data = temp.data<T>();
        BlasConnector::gemm_cm(normal, adjoint, nvirt, this->naos, nvirt,
                               alpha, N_k_data, nvirt,
                               this->psi_ks_vec[0].get_pointer(nocc), this->naos,
                               beta, temp_data, nvirt);

        // Step 3: dmk = temp^T * C_virt^T = (C_virt * N_k * C_virt^H)^T.
        // DensityMatrix expects the AO indices in this transposed order.
        BlasConnector::gemm_cm(transpose, transpose, this->naos, this->naos, nvirt,
                               alpha, temp_data, nvirt,
                               this->psi_ks_vec[0].get_pointer(nocc), this->naos,
                               beta, dmk, this->naos);
    }
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "cal effective dmk elec");
    return result;
}

template <typename T>
void ExcitonPlotter<T>::plot_average_density(const int istate, const std::string& type)
{
    ModuleBase::TITLE("ExcitonPlotter", "plot_average_density");
    if (type != "hole" && type != "elec")
    {
        ModuleBase::WARNING_QUIT("ExcitonPlotter", "Unknown average density type: " + type + ". Use hole or elec.");
    }
    const auto dmk = type == "hole" ? cal_effective_dmk_hole(istate) : cal_effective_dmk_elec(istate);
    elecstate::DensityMatrix<T, double> dm(&this->pmat, this->nspin_x, this->kv.kvec_d, this->nk);
    for (int ik = 0; ik < this->nk; ++ik)
    {
        dm.set_DMK_pointer(ik, dmk[ik].template data<T>());
    }
    LR_Util::initialize_DMR(dm, this->pmat, this->ucell, this->gd_, this->orb_cutoff_);
    dm.cal_DMR();

    double** rho_result = nullptr;
    LR_Util::_allocate_2order_nested_ptr(rho_result, this->nspin_x, this->rho_basis.nrxx);
    for (int is = 0; is < this->nspin_x; ++is)
    {
        ModuleBase::GlobalFunc::ZEROS(rho_result[is], this->rho_basis.nrxx);
    }
    ModuleGint::cal_gint_rho(dm.get_DMR_vector(), this->nspin_x, rho_result, false);

    for (int is = 0; is < this->nspin_x; ++is)
    {
        const std::string filename = this->output_dir_ + "Exciton_avg_" + type + "_state" + std::to_string(istate)
                                     + "_spin" + std::to_string(is) + ".cube";
        ModuleIO::write_vdata_palgrid(this->Pgrid,
                                      rho_result[is],
                                      is,
                                      this->nspin_x,
                                      0,
                                      filename,
                                      this->eig[istate],
                                      &this->ucell,
                                      8,
                                      0,
                                      false, /*two_fermi*/
                                      false);
    }
    LR_Util::_deallocate_2order_nested_ptr(rho_result, this->nspin_x);
}

template <typename T>
void ExcitonPlotter<T>::plot_average_slice(const int istate,
                                           const std::string& type,
                                           const std::string& plane,
                                           const double slice_pos,
                                           const int npoints,
                                           const std::vector<int>& range)
{
    ModuleBase::TITLE("ExcitonPlotter", "plot_average_slice");
    assert(this->nspin_x == 1);
    if (!this->orb_)
    {
        ModuleBase::WARNING_QUIT("ExcitonPlotter",
                                 "plot_average_slice requires LCAO_Orbitals; pass orb to constructor.");
    }
    if (type != "hole" && type != "elec")
    {
        ModuleBase::WARNING_QUIT("ExcitonPlotter", "Unknown average density type: " + type + ". Use hole or elec.");
    }

    const auto dmk = type == "hole" ? cal_effective_dmk_hole(istate) : cal_effective_dmk_elec(istate);
    const SliceGeometry geom = make_slice_geometry(this->ucell, this->kv, plane, slice_pos, npoints, range);
    const double du = (geom.u_end - geom.u_start) / (geom.nu - 1);
    const double dv = (geom.v_end - geom.v_start) / (geom.nv - 1);
    const int cell_res = geom.res;

    std::vector<double> density(geom.nu * geom.nv, 0.0);
    std::vector<double> rho_cache(cell_res * cell_res, 0.0);
    std::vector<bool> cache_valid(cell_res * cell_res, false);
    std::vector<std::complex<double>> phi_bloch(this->naos);

    for (int iu = 0; iu < geom.nu; ++iu)
    {
        const double u_frac = geom.u_start + du * iu;
        const double u_cell = u_frac - std::floor(u_frac);
        int uc_idx = static_cast<int>(std::round(u_cell * cell_res));
        if (uc_idx >= cell_res)
        {
            uc_idx -= cell_res;
        }

        for (int iv = 0; iv < geom.nv; ++iv)
        {
            const double v_frac = geom.v_start + dv * iv;
            const double v_cell = v_frac - std::floor(v_frac);
            int vc_idx = static_cast<int>(std::round(v_cell * cell_res));
            if (vc_idx >= cell_res)
            {
                vc_idx -= cell_res;
            }

            const int cache_idx = uc_idx * cell_res + vc_idx;
            if (!cache_valid[cache_idx])
            {
                const Vec3 position = geom.perp_offset + geom.u_vec * u_cell + geom.v_vec * v_cell;
                double rho = 0.0;
                for (int ik = 0; ik < this->nk; ++ik)
                {
                    std::fill(phi_bloch.begin(), phi_bloch.end(), std::complex<double>(0.0, 0.0));
                    const Vec3 kvec_d
                        = ik < static_cast<int>(this->kv.kvec_d.size()) ? this->kv.kvec_d[ik] : Vec3(0.0, 0.0, 0.0);
                    this->orb_eval_.eval_phi_all_bloch(position, this->ucell, phi_bloch.data(), kvec_d);
                    rho += density_from_dmk(dmk[ik].template data<T>(), phi_bloch, this->naos);
                }
                rho_cache[cache_idx] = rho;
                cache_valid[cache_idx] = true;
            }
            density[iu * geom.nv + iv] = rho_cache[cache_idx];
        }
    }

    const std::string filename
        = this->output_dir_ + "Exciton_avg_" + type + "_slice_state" + std::to_string(istate) + ".dat";
    write_slice_data(this->ucell,
                     geom,
                     density,
                     filename,
                     istate,
                     this->eig[istate],
                     "average_" + type,
                     false,
                     {0.0, 0.0, 0.0});
    std::cout << "Average " << type << " density slice written to " << filename << " (" << geom.nu << "x" << geom.nv
              << ", range [" << geom.u_start_cells << ", " << geom.u_end_cells << "] x [" << geom.v_start_cells
              << ", " << geom.v_end_cells << "], " << geom.res << " pts/cell)" << std::endl;
}

template <typename T>
std::vector<std::vector<std::complex<double>>> ExcitonPlotter<T>::build_conditional_coefficients(
    const int istate,
    const std::array<double, 3>& r_fix_in,
    const bool plot_electron)
{
    using Complex = std::complex<double>;
    const Vec3 r_fix(r_fix_in[0], r_fix_in[1], r_fix_in[2]);
    const int nocc = this->nocc[0];
    const int nvirt = this->nvirt[0];
    const int nmix = plot_electron ? nvirt : nocc;
    const int offset_b = istate * this->ldim;

    std::vector<std::vector<Complex>> mixing(this->nk, std::vector<Complex>(nmix, Complex(0.0, 0.0)));
    for (int ik = 0; ik < this->nk; ++ik)
    {
        const int x_start = ik * nocc * nvirt;
        if (plot_electron)
        {
            for (int io = 0; io < nocc; ++io)
            {
                const Complex fixed_wfc = std::conj(
                    this->orb_eval_
                        .eval_wfc_bloch<T>(r_fix, ik, io, this->psi_ks_vec[0], this->ucell, this->kv.kvec_d[ik]));
                for (int iv = 0; iv < nvirt; ++iv)
                {
                    mixing[ik][iv] += this->X[offset_b + x_start + iv + io * nvirt] * fixed_wfc;
                }
            }
        }
        else
        {
            for (int iv = 0; iv < nvirt; ++iv)
            {
                const Complex fixed_wfc = this->orb_eval_.eval_wfc_bloch<T>(r_fix,
                                                                            ik,
                                                                            nocc + iv,
                                                                            this->psi_ks_vec[0],
                                                                            this->ucell,
                                                                            this->kv.kvec_d[ik]);
                for (int io = 0; io < nocc; ++io)
                {
                    mixing[ik][io] += this->X[offset_b + x_start + iv + io * nvirt] * fixed_wfc;
                }
            }
        }
    }

    std::vector<std::vector<Complex>> coefficients(this->nk, std::vector<Complex>(this->naos, Complex(0.0, 0.0)));
    for (int ik = 0; ik < this->nk; ++ik)
    {
        this->psi_ks_vec[0].fix_k(ik);
        const int band_start = plot_electron ? nocc : 0;
        if (plot_electron)
        {
            for (int ib = 0; ib < nmix; ++ib)
            {
                const T* const band = this->psi_ks_vec[0].get_pointer(band_start + ib);
                for (int mu = 0; mu < this->naos; ++mu)
                {
                    coefficients[ik][mu] += band[mu] * mixing[ik][ib];
                }
            }
        }
        else
        {
            for (int ib = 0; ib < nmix; ++ib)
            {
                const T* const band = this->psi_ks_vec[0].get_pointer(band_start + ib);
                for (int mu = 0; mu < this->naos; ++mu)
                {
                    coefficients[ik][mu] += std::conj(band[mu]) * mixing[ik][ib];
                }
            }
        }
    }
    return coefficients;
}

template <typename T>
void ExcitonPlotter<T>::plot_cond_slice(const int istate,
                                        const std::array<double, 3>& r_fix_in,
                                        const std::string& plane,
                                        const double slice_pos,
                                        const int npoints,
                                        const std::vector<int>& range,
                                        const std::string& type)
{
    using Complex = std::complex<double>;
    ModuleBase::TITLE("ExcitonPlotter", "plot_cond_slice");
    assert(this->nspin_x == 1);
    if (!this->orb_)
    {
        ModuleBase::WARNING_QUIT("ExcitonPlotter", "plot_cond_slice requires LCAO_Orbitals; pass orb to constructor.");
    }
    const bool plot_electron = type == "elec";
    if (!plot_electron && type != "hole")
    {
        ModuleBase::WARNING_QUIT("ExcitonPlotter", "Unknown conditional slice type: " + type + ". Use elec or hole.");
    }
    const auto coefficients = build_conditional_coefficients(istate, r_fix_in, plot_electron);
    const SliceGeometry geom = make_slice_geometry(this->ucell, this->kv, plane, slice_pos, npoints, range);
    const double du = (geom.u_end - geom.u_start) / (geom.nu - 1);
    const double dv = (geom.v_end - geom.v_start) / (geom.nv - 1);

    std::vector<double> phase_u(this->nk);
    std::vector<double> phase_v(this->nk);
    for (int ik = 0; ik < this->nk; ++ik)
    {
        phase_u[ik] = bloch_phase_arg_axis(this->kv.kvec_d[ik], geom.u_axis, 1);
        phase_v[ik] = bloch_phase_arg_axis(this->kv.kvec_d[ik], geom.v_axis, 1);
    }

    const int cell_res = geom.res;
    std::vector<double> density(geom.nu * geom.nv, 0.0);
    std::vector<Complex> phi_bloch(this->naos);
    std::vector<std::vector<Complex>> psi_cache(cell_res * cell_res, std::vector<Complex>(this->nk, Complex(0.0, 0.0)));
    std::vector<bool> cache_valid(cell_res * cell_res, false);

    for (int iu = 0; iu < geom.nu; ++iu)
    {
        const double u_frac = geom.u_start + du * iu;
        const double u_cell = u_frac - std::floor(u_frac);
        const int u_translation = static_cast<int>(std::floor(u_frac));
        int uc_idx = static_cast<int>(std::round(u_cell * cell_res));
        if (uc_idx >= cell_res)
        {
            uc_idx -= cell_res;
        }

        for (int iv = 0; iv < geom.nv; ++iv)
        {
            const double v_frac = geom.v_start + dv * iv;
            const double v_cell = v_frac - std::floor(v_frac);
            const int v_translation = static_cast<int>(std::floor(v_frac));
            int vc_idx = static_cast<int>(std::round(v_cell * cell_res));
            if (vc_idx >= cell_res)
            {
                vc_idx -= cell_res;
            }

            const int cache_index = uc_idx * cell_res + vc_idx;
            if (!cache_valid[cache_index])
            {
                const Vec3 home_position = geom.perp_offset + geom.u_vec * u_cell + geom.v_vec * v_cell;
                for (int ik = 0; ik < this->nk; ++ik)
                {
                    std::fill(phi_bloch.begin(), phi_bloch.end(), Complex(0.0, 0.0));
                    this->orb_eval_.eval_phi_all_bloch(home_position,
                                                       this->ucell,
                                                       phi_bloch.data(),
                                                       this->kv.kvec_d[ik]);
                    Complex wfc_k(0.0, 0.0);
                    if (plot_electron)
                    {
                        for (int mu = 0; mu < this->naos; ++mu)
                        {
                            wfc_k += coefficients[ik][mu] * phi_bloch[mu];
                        }
                    }
                    else
                    {
                        for (int mu = 0; mu < this->naos; ++mu)
                        {
                            wfc_k += coefficients[ik][mu] * std::conj(phi_bloch[mu]);
                        }
                    }
                    psi_cache[cache_index][ik] = wfc_k;
                }
                cache_valid[cache_index] = true;
            }

            Complex conditional_wfc(0.0, 0.0);
            for (int ik = 0; ik < this->nk; ++ik)
            {
                double phase_arg = phase_u[ik] * u_translation + phase_v[ik] * v_translation;
                if (!plot_electron)
                {
                    phase_arg = -phase_arg;
                }
                const Complex phase(std::cos(phase_arg), std::sin(phase_arg));
                conditional_wfc += phase * psi_cache[cache_index][ik];
            }
            density[iu * geom.nv + iv] = std::norm(conditional_wfc);
        }
    }

    std::cout << "Bloch-sum evaluation complete (cached " << cell_res << "x" << cell_res << " home-cell positions)."
              << std::endl;
    const std::string filename
        = this->output_dir_ + "Exciton_cond_" + type + "_slice_state" + std::to_string(istate) + ".dat";
    write_slice_data(this->ucell,
                     geom,
                     density,
                     filename,
                     istate,
                     this->eig[istate],
                     "conditional_" + type,
                     true,
                     r_fix_in);
    std::cout << "Conditional density slice written to " << filename << " (" << geom.nu << "x" << geom.nv << ", range "
              << "[" << geom.u_start_cells << ", " << geom.u_end_cells << "] x [" << geom.v_start_cells << ", "
              << geom.v_end_cells << "], " << geom.res << " pts/cell)" << std::endl;
}

template class ExcitonPlotter<double>;
template class ExcitonPlotter<std::complex<double>>;

} // namespace LR_Util
