// The KB-projector construction of DFPT_Pert (radial_vq, real_ylm,
// grad_real_ylm, build_vkb, build_vkb_dk) lives in dfpt_pert_vkb.cpp
// and the nonlocal first/second-order potentials (dVnl_dtau,
// apply_d2vnl) in dfpt_pert_nl.cpp.

#include "dfpt_pert.h"

#include "source_base/constants.h"
#include "source_base/global_function.h"
#include "source_base/math_integral.h"
#include "source_base/truncated_func.h"
#include "source_cell/atom_pseudo.h"
#include "source_cell/atom_spec.h"
#include "source_pw/module_pwdft/stru_fac.h"

#include <cmath>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

#include "source_base/timer.h"
#include "source_base/tool_title.h"

namespace ModuleDFPT
{

DFPT_Pert::DFPT_Pert()
{
}

DFPT_Pert::~DFPT_Pert()
{
}

void DFPT_Pert::init(UnitCell& ucell, ModulePW::PW_Basis* pw_rho, ModulePW::PW_Basis_K* pw_wfc, Structure_Factor& sf)
{
    ModuleBase::TITLE("DFPT_Pert", "init");
    ModuleBase::timer::start("DFPT_Pert", "init");
    ucell_ = &ucell;
    pw_rho_ = pw_rho;
    pw_wfc_ = pw_wfc;
    sf_ = &sf;
    ModuleBase::timer::end("DFPT_Pert", "init");
}

void DFPT_Pert::atom_index(int atom_idx, int& it, int& ia) const
{
    ModuleBase::TITLE("DFPT_Pert", "atom_index");
    ModuleBase::timer::start("DFPT_Pert", "atom_index");
    it = 0;
    ia = atom_idx;
    for (int it_type = 0; it_type < ucell_->ntype; ++it_type)
    {
        if (ia < ucell_->atoms[it_type].na)
        {
            it = it_type;
            ModuleBase::timer::end("DFPT_Pert", "atom_index");
            return;
        }
        ia -= ucell_->atoms[it_type].na;
    }
    // out of range: leave it/ia at the last type / last picture and let the
    // caller guard; dV requests with invalid indices simply produce nothing.
    ia = -1;
    ModuleBase::timer::end("DFPT_Pert", "atom_index");
}

void DFPT_Pert::rho_gvec(int ig, ModuleBase::Vector3<double>& gcar) const
{
    ModuleBase::TITLE("DFPT_Pert", "rho_gvec");
    ModuleBase::timer::start("DFPT_Pert", "rho_gvec");
    const int isz = pw_rho_->ig2isz[ig];
    int iz = isz % pw_rho_->nz;
    const int is = isz / pw_rho_->nz;
    const int ixy = pw_rho_->is2fftixy[is];
    int ix = ixy / pw_rho_->fftny;
    int iy = ixy % pw_rho_->fftny;
    if (ix >= int(pw_rho_->nx / 2) + 1)
    {
        ix -= pw_rho_->nx;
    }
    if (iy >= int(pw_rho_->ny / 2) + 1)
    {
        iy -= pw_rho_->ny;
    }
    if (iz >= int(pw_rho_->nz / 2) + 1)
    {
        iz -= pw_rho_->nz;
    }
    gcar = ModuleBase::Vector3<double>(ix, iy, iz) * ucell_->G;
    ModuleBase::timer::end("DFPT_Pert", "rho_gvec");
}

double DFPT_Pert::vloc_at_g(int it, double g2) const
{
    ModuleBase::TITLE("DFPT_Pert", "vloc_at_g");
    ModuleBase::timer::start("DFPT_Pert", "vloc_at_g");
    // g2 is the squared magnitude in bohr^-2 units.
    const double g_zero_tol = 1.0e-8; ///< empirical parameter: |G| floor (bohr^-1) for the G=0 radial integral
    const Atom* atom = &ucell_->atoms[it];
    const double zv = atom->ncpp.zv;
    if (atom->coulomb_potential)
    {
        // analytic Coulomb local potential (vl_pw.cpp::vloc_coulomb)
        ModuleBase::timer::end("DFPT_Pert", "vloc_at_g");
        return -zv * ModuleBase::e2 * ModuleBase::FOUR_PI / ucell_->omega / g2;
    }
    // numeric pseudopotential: mirror vl_pw.cpp::vloc_of_g at the requested
    // magnitude instead of interpolating the precomputed shell table. This
    // keeps the rho-grid kernel consistent with the ground-state local
    // potential for every magnitude |Delta+q|.
    const int msh = atom->ncpp.msh;
    const double fac = zv * ModuleBase::e2;
    std::vector<double> aux(msh);
    const double g = std::sqrt(g2);
    if (g < g_zero_tol)
    {
        double v0 = 0.0;
        for (int ir = 0; ir < msh; ++ir)
        {
            aux[ir] = atom->ncpp.r[ir] * (atom->ncpp.r[ir] * atom->ncpp.vloc_at[ir] + fac);
        }
        ModuleBase::Integral::Simpson_Integral(msh, aux.data(), atom->ncpp.rab.data(), v0);
        ModuleBase::timer::end("DFPT_Pert", "vloc_at_g");
        return v0 * ModuleBase::FOUR_PI / ucell_->omega;
    }
    for (int ir = 0; ir < msh; ++ir)
    {
        aux[ir] = (atom->ncpp.r[ir] * atom->ncpp.vloc_at[ir] + fac * std::erf(atom->ncpp.r[ir]))
                  * std::sin(g * atom->ncpp.r[ir]) / g;
    }
    double v = 0.0;
    ModuleBase::Integral::Simpson_Integral(msh, aux.data(), atom->ncpp.rab.data(), v);
    // erf(r)-compensating gaussian subtraction (same as vloc_of_g)
    v -= fac * ModuleBase::truncated_exp(-g2 * 0.25) / g2;
    ModuleBase::timer::end("DFPT_Pert", "vloc_at_g");
    return v * ModuleBase::FOUR_PI / ucell_->omega;
}

void DFPT_Pert::dVloc_dtau(int atom_idx,
                           int dir,
                           const ModuleBase::Vector3<double>& q,
                           std::vector<std::complex<double>>& dv)
{
    ModuleBase::TITLE("DFPT_Pert", "dVloc_dtau");
    ModuleBase::timer::start("DFPT_Pert", "dVloc_dtau");
    if (pw_rho_ == nullptr || pw_rho_->gamma_only)
    {
        ModuleBase::WARNING_QUIT("DFPT_Pert::dVloc_dtau",
                                 "DFPT requires a complex (gamma_only=false) real-space basis.");
    }
    int it = 0;
    int ia = 0;
    atom_index(atom_idx, it, ia);
    if (ia < 0)
    {
        ModuleBase::timer::end("DFPT_Pert", "dVloc_dtau");
        return;
    }
    const ModuleBase::Vector3<double>& tau = ucell_->atoms[it].tau[ia];
    const int npw = pw_rho_->npw;
    dv.assign(npw, std::complex<double>(0.0, 0.0));
    ModuleBase::Vector3<double> gcar;
    const double w2_floor = 1.0e-12; ///< empirical parameter: |Delta+q|^2 zero-shell guard (2*pi/lat0 units)
    for (int ig = 0; ig < npw; ++ig)
    {
        rho_gvec(ig, gcar);
        const ModuleBase::Vector3<double> w = gcar + q; // Delta + q, 2*pi/lat0 units
        const double w2 = w * w;
        // the Delta == -q component carries no displacement gradient (constant
        // potential shift) and is dropped, consistently with the G=0 handling
        // of the ground-state local potential.
        if (w2 < w2_floor)
        {
            continue;
        }
        const double g_bohr2 = w2 * ucell_->tpiba2;
        const double vloc = vloc_at_g(it, g_bohr2);
        // GS structure-factor convention (stru_fac.cpp: ci_tpi =
        // NEG_IMAG_UNIT * 2pi): exp(-i 2pi (g.tau)), with g in 1/lat0 units
        // and tau in lat0 units; 2pi/lat0 = tpiba only multiplies the
        // magnitude (vl_pw.cpp: qnorm = |g| * tpiba).
        const double arg = -ModuleBase::TWO_PI * (w * tau);
        const std::complex<double> phase(std::cos(arg), std::sin(arg));
        // dV_loc / d tau_direction = -i g_dir * Vloc * exp(-i (Delta+q).tau)
        const std::complex<double> iw_dir = std::complex<double>(0.0, -1.0) * (ucell_->tpiba * w[dir]);
        dv[ig] = iw_dir * vloc * phase;
    }
    ModuleBase::timer::end("DFPT_Pert", "dVloc_dtau");
}

void DFPT_Pert::build_dv(int q_idx, int atom_idx, int dir, DFPT_PW_Data& data)
{
    ModuleBase::TITLE("DFPT_Pert", "build_dv");
    ModuleBase::timer::start("DFPT_Pert", "build_dv");
    // the local first-order potential is assembled on the rho grid in reciprocal
    // space (reciprocal coefficients indexed by the rho-basis ig), then brought
    // to the shared real-space grid where apply_dv performs the convolution.
    if (pw_rho_ == nullptr)
    {
        ModuleBase::timer::end("DFPT_Pert", "build_dv");
        return;
    }
    const ModuleBase::Vector3<double> q_cart = data.get_qvec(q_idx) * ucell_->G;
    std::vector<std::complex<double>> dv_recip;
    dVloc_dtau(atom_idx, dir, q_cart, dv_recip);
    data.set_dv_recip_c(q_idx, 0, dv_recip);

    std::vector<std::complex<double>> dv_real(pw_rho_->nrxx);
    pw_rho_->recip2real(dv_recip.data(), dv_real.data());
    data.set_dv_rc(q_idx, 0, dv_real);
    data.set_pert_atom(atom_idx);
    data.set_pert_dir(dir);

    // DFT+U perturbation reservation (U0): append the first-order Hubbard
    // potential dV_U when a DFT+U provider is wired. Physical implementation
    // lands in C1 (frozen projector term) and C3 (occupation response).
    if (data.with_u())
    {
        build_dv_u(q_idx, atom_idx, dir, data);
    }
    ModuleBase::timer::end("DFPT_Pert", "build_dv");
}

void DFPT_Pert::real_space_dv(int q_idx,
                              int k_idx,
                              const psi::Psi<std::complex<double>>& psi,
                              DFPT_PW_Data& data,
                              const DFPT_KQ_Basis& kq,
                              std::vector<std::vector<std::complex<double>>>& dv_psi) const
{
    ModuleBase::TITLE("DFPT_Pert", "real_space_dv");
    ModuleBase::timer::start("DFPT_Pert", "real_space_dv");
    const std::vector<std::complex<double>> dv_rc = data.get_dv_rc(q_idx, 0);
    if (dv_rc.empty() || dv_rc.size() != static_cast<size_t>(pw_rho_->nrxx))
    {
        ModuleBase::timer::end("DFPT_Pert", "real_space_dv");
        return;
    }
    apply_vr_core(k_idx, dv_rc, psi, kq, dv_psi);
    ModuleBase::timer::end("DFPT_Pert", "real_space_dv");
}

void DFPT_Pert::apply_vr(int q_idx,
                         int k_idx,
                         const std::vector<std::complex<double>>& v_rc,
                         const psi::Psi<std::complex<double>>& psi,
                         const ModuleBase::Vector3<double>& q_cart,
                         std::vector<std::vector<std::complex<double>>>& dv_psi) const
{
    ModuleBase::TITLE("DFPT_Pert", "apply_vr");
    ModuleBase::timer::start("DFPT_Pert", "apply_vr");
    (void)q_idx;
    if (pw_rho_ == nullptr || pw_wfc_ == nullptr || v_rc.size() != static_cast<size_t>(pw_rho_->nrxx))
    {
        dv_psi.clear();
        ModuleBase::timer::end("DFPT_Pert", "apply_vr");
        return;
    }
    DFPT_KQ_Basis kq;
    kq.init(pw_wfc_, pw_rho_, q_cart, k_idx);
    apply_vr_core(k_idx, v_rc, psi, kq, dv_psi);
    ModuleBase::timer::end("DFPT_Pert", "apply_vr");
}

void DFPT_Pert::apply_vr_core(int k_idx,
                              const std::vector<std::complex<double>>& v_rc,
                              const psi::Psi<std::complex<double>>& psi,
                              const DFPT_KQ_Basis& kq,
                              std::vector<std::vector<std::complex<double>>>& dv_psi) const
{
    ModuleBase::TITLE("DFPT_Pert", "apply_vr_core");
    ModuleBase::timer::start("DFPT_Pert", "apply_vr_core");
    const int nbands = psi.get_nbands();
    const int npwk_kq = kq.get_npwk();
    std::vector<std::complex<double>> u_r(pw_rho_->nrxx);
    std::vector<std::complex<double>> d_r(pw_rho_->nrxx);
    std::vector<std::complex<double>> d_recip(pw_rho_->npw);
    dv_psi.assign(nbands, std::vector<std::complex<double>>(npwk_kq, std::complex<double>(0.0, 0.0)));
    for (int iband = 0; iband < nbands; ++iband)
    {
        pw_wfc_->recip2real(&psi(k_idx, iband, 0), u_r.data(), k_idx);
        for (int ir = 0; ir < pw_rho_->nrxx; ++ir)
        {
            d_r[ir] = u_r[ir] * v_rc[ir];
        }
        pw_rho_->real2recip(d_r.data(), d_recip.data());
        std::vector<std::complex<double>> dpsi(npwk_kq, std::complex<double>(0.0, 0.0));
        for (int igl = 0; igl < npwk_kq; ++igl)
        {
            const int ig_rho = kq.get_ig_rho(igl);
            if (ig_rho >= 0)
            {
                dpsi[igl] = d_recip[ig_rho];
            }
        }
        dv_psi[iband] = dpsi;
    }
    ModuleBase::timer::end("DFPT_Pert", "apply_vr_core");
}

void DFPT_Pert::apply_dv(int q_idx, int k_idx, const psi::Psi<std::complex<double>>& psi, DFPT_PW_Data& data)
{
    ModuleBase::TITLE("DFPT_Pert", "apply_dv");
    ModuleBase::timer::start("DFPT_Pert", "apply_dv");
    const int atom_idx = data.get_pert_atom();
    const int dir = data.get_pert_dir();
    const ModuleBase::Vector3<double> q_cart = data.get_qvec(q_idx) * ucell_->G;

    DFPT_KQ_Basis kq;
    kq.init(pw_wfc_, pw_rho_, q_cart, k_idx);

    const int nbands = psi.get_nbands();
    std::vector<std::vector<std::complex<double>>> dv_psi(nbands);

    // local contribution: dVloc(r) * psi on the shared real-space grid
    real_space_dv(q_idx, k_idx, psi, data, kq, dv_psi);

    // nonlocal contribution: dVnl/dtau |psi> (per displaced atom)
    std::vector<std::vector<std::complex<double>>> dv_psi_nl;
    dVnl_dtau(atom_idx, dir, q_cart, psi, k_idx, dv_psi_nl);
    if (dv_psi_nl.size() == static_cast<size_t>(nbands))
    {
        for (int iband = 0; iband < nbands; ++iband)
        {
            if (dv_psi[iband].size() != dv_psi_nl[iband].size())
            {
                continue;
            }
            for (size_t i = 0; i < dv_psi[iband].size(); ++i)
            {
                dv_psi[iband][i] += dv_psi_nl[iband][i];
            }
        }
    }

    for (int iband = 0; iband < nbands; ++iband)
    {
        data.set_dpsi(q_idx, k_idx, iband, dv_psi[iband]);
    }
    ModuleBase::timer::end("DFPT_Pert", "apply_dv");
}

void DFPT_Pert::build_dv_u(int q_idx, int atom_idx, int dir, DFPT_PW_Data& data)
{
    ModuleBase::TITLE("DFPT_Pert", "build_dv_u");
    ModuleBase::timer::start("DFPT_Pert", "build_dv_u");
    // C1 frozen term of the first-order Hubbard potential:
    //   |dphi(k+q)/dtau> V_eff <phi(k)|psi> + adjoint
    // The provider is only usable when its occupation matrices are
    // initialized (u_active()); DFPT_PW::init additionally rejects a wired
    // provider outright (every U hook is a no-op U0 reservation), so this
    // guard is defense in depth; the diamond DFT+U test (C7) will exercise
    // this path once OnsiteProjector integration on the DFPT k+q basis is
    // finalized.
    if (!data.u_active())
    {
        return;
    }
    (void)q_idx;
    (void)atom_idx;
    (void)dir;
    // TODO(C7): Onsite_Proj_tools on the DFPT k+q basis; occupation response
    // (|phi(k+q)> U(diag*delta - docc) <phi(k)|psi>) lands in C3 after docc.
}

void DFPT_Pert::d2vloc_r(int atom_idx, int da, int db, std::vector<std::complex<double>>& dv2_r) const
{
    ModuleBase::TITLE("DFPT_Pert", "d2vloc_r");
    ModuleBase::timer::start("DFPT_Pert", "d2vloc_r");
    if (pw_rho_ == nullptr)
    {
        ModuleBase::timer::end("DFPT_Pert", "d2vloc_r");
        return;
    }
    int it = 0;
    int ia = 0;
    atom_index(atom_idx, it, ia);
    if (ia < 0)
    {
        dv2_r.clear();
        ModuleBase::timer::end("DFPT_Pert", "d2vloc_r");
        return;
    }
    const ModuleBase::Vector3<double>& tau = ucell_->atoms[it].tau[ia];
    const int npw = pw_rho_->npw;
    std::vector<std::complex<double>> dv2_recip(npw, std::complex<double>(0.0, 0.0));
    ModuleBase::Vector3<double> gcar;
    const double w2_floor = 1.0e-12; ///< empirical parameter: |G|^2 zero-shell guard (2*pi/lat0 units)
    for (int ig = 0; ig < npw; ++ig)
    {
        rho_gvec(ig, gcar);
        // QE ground truth (dynmat_us.f90): the mixed (+q,-q) second-order
        // local potential is the integer-G, q-independent kernel
        // -tpiba^2 G_da G_db vloc(|G|) exp(-i G.tau); the (+q,-q) dressings
        // collapse the carrier to 0 for every q, not only when 2q is
        // reciprocal.
        const ModuleBase::Vector3<double> w = gcar;
        const double w2 = w * w;
        if (w2 < w2_floor)
        {
            continue;
        }
        const double vloc = vloc_at_g(it, w2 * ucell_->tpiba2);
        const double arg = -ModuleBase::TWO_PI * (w * tau);
        const std::complex<double> phase(std::cos(arg), std::sin(arg));
        // (-i g_da)(-i g_db) = -g_da g_db
        dv2_recip[ig] = -(ucell_->tpiba * w[da]) * (ucell_->tpiba * w[db]) * vloc * phase;
    }
    dv2_r.assign(pw_rho_->nrxx, std::complex<double>(0.0, 0.0));
    pw_rho_->recip2real(dv2_recip.data(), dv2_r.data());
    ModuleBase::timer::end("DFPT_Pert", "d2vloc_r");
}

void DFPT_Pert::build_efield(const ModuleBase::Vector3<double>& field, DFPT_PW_Data& data)
{
    ModuleBase::TITLE("DFPT_Pert", "build_efield");
    ModuleBase::timer::start("DFPT_Pert", "build_efield");
    // first-order electric-field potential: delta V(r) = - r . E (q=0 limit,
    // position operator in the periodic cell). Computed directly on the shared
    // real-space grid. Only relevant for the Q0 dielectric response (C6).
    if (pw_rho_ == nullptr)
    {
        ModuleBase::timer::end("DFPT_Pert", "build_efield");
        return;
    }
    if (pw_rho_->gamma_only)
    {
        ModuleBase::WARNING_QUIT("DFPT_Pert::build_efield",
                                 "DFPT requires a complex (gamma_only=false) real-space basis.");
    }
    std::vector<std::complex<double>> dv_real(pw_rho_->nrxx, std::complex<double>(0.0, 0.0));
    const ModuleBase::Matrix3& latvec = ucell_->latvec;
    const double lat0 = ucell_->lat0;
    // shared real-space grid layout (serial pool): ir = (ix*ny + iy)*nz + iz,
    // i.e. z runs fastest (verified against the impulse response of the FFT).
    for (int ir = 0; ir < pw_rho_->nrxx; ++ir)
    {
        const int iz = ir % pw_rho_->nz;
        const int rem = ir / pw_rho_->nz;
        const int iy = rem % pw_rho_->ny;
        const int ix = rem / pw_rho_->ny;
        const double fx = ix / static_cast<double>(pw_rho_->nx);
        const double fy = iy / static_cast<double>(pw_rho_->ny);
        const double fz = iz / static_cast<double>(pw_rho_->nz);
        ModuleBase::Vector3<double> r;
        r.x = (fx * latvec.e11 + fy * latvec.e12 + fz * latvec.e13) * lat0;
        r.y = (fx * latvec.e21 + fy * latvec.e22 + fz * latvec.e23) * lat0;
        r.z = (fx * latvec.e31 + fy * latvec.e32 + fz * latvec.e33) * lat0;
        dv_real[ir] = -(field * r); // -e r.E (e absorbed in field convention)
    }
    data.set_dv_rc(0, 0, dv_real);
    ModuleBase::timer::end("DFPT_Pert", "build_efield");
}

} // namespace ModuleDFPT
