#include "dfpt_q0.h"

#include "source_base/timer.h"
#include "source_base/tool_title.h"

#include <cmath>
#include <complex>
#include <vector>

namespace ModuleDFPT
{

namespace
{

/// element accessor for ModuleBase::Matrix3 (row i, column j); the public
/// interface only exposes the named e11..e33 members
inline double me(const ModuleBase::Matrix3& m, int i, int j)
{
    switch (3 * i + j)
    {
    case 0:
        return m.e11;
    case 1:
        return m.e12;
    case 2:
        return m.e13;
    case 3:
        return m.e21;
    case 4:
        return m.e22;
    case 5:
        return m.e23;
    case 6:
        return m.e31;
    case 7:
        return m.e32;
    default:
        return m.e33;
    }
}

/// folded fractional equality with the lattice periodicity absorbed
inline bool folded_equal(double a, double b, double tol)
{
    const double d = std::abs(a - b);
    return d < tol || std::abs(d - 1.0) < tol;
}

/// true if kp differs from every already-folded star member
bool is_new_member(const std::vector<ModuleBase::Vector3<double>>& kfolds, const ModuleBase::Vector3<double>& kp)
{
    const double kfold_tol = 1.0e-5; ///< empirical parameter: folded fractional-coordinate match tolerance
    for (size_t im = 0; im < kfolds.size(); ++im)
    {
        if (folded_equal(kp.x, kfolds[im].x, kfold_tol) && folded_equal(kp.y, kfolds[im].y, kfold_tol)
            && folded_equal(kp.z, kfolds[im].z, kfold_tol))
        {
            return false;
        }
    }
    return true;
}

/// atom image map iat -> image atom of the j-th symmetry operation
/// (direct-space gmatrix/gtrans pair; species map onto themselves).
/// Returns false when some atom has no same-species image, i.e. the
/// operation set is inconsistent with the structure.
bool make_atom_map(const UnitCell& ucell, const ModuleSymmetry::Symmetry& symm, int j, std::vector<int>& atom_map)
{
    const int nat = ucell.nat;
    const double tau_match_tol = 1.0e-4; ///< empirical parameter: direct-space atom-image match tolerance
    atom_map.assign(nat, -1);
    for (int iat = 0; iat < nat; ++iat)
    {
        const int it = ucell.iat2it[iat];
        const int ia = ucell.iat2ia[iat];
        ModuleBase::Vector3<double> tp = ucell.atoms[it].taud[ia] * symm.gmatrix[j] + symm.gtrans[j];
        tp.x -= std::floor(tp.x);
        tp.y -= std::floor(tp.y);
        tp.z -= std::floor(tp.z);
        for (int jat = 0; jat < nat; ++jat)
        {
            if (ucell.iat2it[jat] != it)
            {
                continue; // a species maps onto itself
            }
            const int ja = ucell.iat2ia[jat];
            const ModuleBase::Vector3<double>& tq = ucell.atoms[it].taud[ja];
            if (folded_equal(tp.x, tq.x, tau_match_tol) && folded_equal(tp.y, tq.y, tau_match_tol)
                && folded_equal(tp.z, tq.z, tau_match_tol))
            {
                atom_map[iat] = jat;
                break;
            }
        }
        if (atom_map[iat] < 0)
        {
            return false;
        }
    }
    return true;
}

/// [ik][v][ig] wave-function / response coefficients over the k mesh
typedef std::vector<std::vector<std::vector<std::complex<double>>>> KMeshCoeffs;

/// wg-weighted partial chi_k[ik](a, b) = sum_occ Re<Y^a|dpsi^E,b> at
/// every stored k (QE dielec.f90: eps -= 4*(4pi/Omega)*wk*Re<dvpsi|dpsi>)
void accumulate_chi_eps(const ModuleBase::matrix& wg,
                        const std::vector<KMeshCoeffs>& yr,
                        const std::vector<KMeshCoeffs>& de,
                        int nk,
                        std::vector<ModuleBase::matrix>& chi_k)
{
    const int nbands = wg.nc;
    for (int ik = 0; ik < nk; ++ik)
    {
        for (int v = 0; v < nbands; ++v)
        {
            if (!dfpt_band_occupied(wg, ik, v))
            {
                continue; // empty
            }
            for (int a = 0; a < 3; ++a)
            {
                const int npw = static_cast<int>(yr[a][ik][v].size());
                if (npw <= 0)
                {
                    continue;
                }
                for (int b = 0; b < 3; ++b)
                {
                    if (static_cast<int>(de[b][ik][v].size()) != npw)
                    {
                        continue;
                    }
                    std::complex<double> dot(0.0, 0.0);
                    for (int ig = 0; ig < npw; ++ig)
                    {
                        dot += std::conj(yr[a][ik][v][ig]) * de[b][ik][v][ig];
                    }
                    chi_k[ik](a, b) += wg(ik, v) * dot.real();
                }
            }
        }
    }
}

/// star average of the partial tensors: eps += (1/n_star) R chi(k) R^T
void star_average_eps(const std::vector<std::vector<DFPT_Q0::StarMember>>& stars,
                      const std::vector<ModuleBase::matrix>& chi_k,
                      ModuleBase::matrix& eps)
{
    const int nk = static_cast<int>(stars.size());
    for (int ik = 0; ik < nk; ++ik)
    {
        const double inv_nstar = 1.0 / static_cast<double>(stars[ik].size());
        for (size_t im = 0; im < stars[ik].size(); ++im)
        {
            double rot[9];
            DFPT_Q0::rotate_tensor(stars[ik][im].cart, chi_k[ik], rot);
            for (int a = 0; a < 3; ++a)
            {
                for (int b = 0; b < 3; ++b)
                {
                    eps(a, b) += inv_nstar * rot[3 * a + b];
                }
            }
        }
    }
}

/// wg-weighted partial chi_k[ik](a, idir) of one atom displacement
/// direction, paired with the position legs Y^a of the q = 0 mesh
void accumulate_disp_chi(const ModuleBase::matrix& wg,
                         const std::vector<KMeshCoeffs>& yr,
                         const KMeshCoeffs& disp,
                         int idir,
                         int nbasis,
                         std::vector<ModuleBase::matrix>& chi_k)
{
    const int nk = static_cast<int>(disp.size());
    const int nbands = wg.nc;
    for (int ik = 0; ik < nk; ++ik)
    {
        for (int v = 0; v < nbands; ++v)
        {
            if (!dfpt_band_occupied(wg, ik, v))
            {
                continue; // empty
            }
            const int npw = static_cast<int>(disp[ik][v].size());
            if (npw <= 0 || npw > nbasis)
            {
                continue; // unsolved slot or inconsistent basis
            }
            // <dpsi^kappa_idir(scf)|Y^a_v> per field direction
            for (int a = 0; a < 3; ++a)
            {
                if (static_cast<int>(yr[a][ik][v].size()) != npw)
                {
                    continue;
                }
                std::complex<double> dot(0.0, 0.0);
                for (int ig = 0; ig < npw; ++ig)
                {
                    dot += std::conj(disp[ik][v][ig]) * yr[a][ik][v][ig];
                }
                chi_k[ik](a, idir) += wg(ik, v) * dot.real();
            }
        }
    }
}

/// star-credit the atom-resolved partials: the partial at member Rk is
/// R chi(k) R^T and is credited to the image atom R(iat)
void credit_star_born(const std::vector<std::vector<DFPT_Q0::StarMember>>& stars,
                      const std::vector<ModuleBase::matrix>& chi_k,
                      int iat,
                      std::vector<ModuleBase::matrix>& zacc)
{
    const int nk = static_cast<int>(stars.size());
    for (int ik = 0; ik < nk; ++ik)
    {
        const double inv_nstar = 1.0 / static_cast<double>(stars[ik].size());
        for (size_t im = 0; im < stars[ik].size(); ++im)
        {
            const DFPT_Q0::StarMember& mem = stars[ik][im];
            const int jat = (mem.atom_map.empty()) ? iat : mem.atom_map[iat];
            double rot[9];
            DFPT_Q0::rotate_tensor(mem.cart, chi_k[ik], rot);
            for (int a = 0; a < 3; ++a)
            {
                for (int d = 0; d < 3; ++d)
                {
                    zacc[jat](a, d) += inv_nstar * rot[3 * a + d];
                }
            }
        }
    }
}

/// assemble Z*_iat = -2 Zacc + Z_ion on the diagonal and stash per atom
void set_born_charges(const UnitCell& ucell, const std::vector<ModuleBase::matrix>& zacc, DFPT_PW_Data& data)
{
    const int nat = ucell.nat;
    for (int iat = 0; iat < nat; ++iat)
    {
        ModuleBase::matrix zstar(3, 3, true);
        for (int a = 0; a < 3; ++a)
        {
            for (int d = 0; d < 3; ++d)
            {
                zstar(a, d) = -2.0 * zacc[iat](a, d);
            }
        }
        // ionic rigid-ion charge on the diagonal (a == b directions)
        const int it = ucell.iat2it[iat];
        const double zion = ucell.atoms[it].ncpp.zv;
        for (int d = 0; d < 3; ++d)
        {
            zstar(d, d) += zion;
        }
        data.set_born(iat, zstar);
    }
}

} // namespace

void DFPT_Q0::init(UnitCell& ucell, ModulePW::PW_Basis* pw_rho, ModulePW::PW_Basis_K* pw_wfc, DFPT_Pert* pert)
{
    ModuleBase::TITLE("DFPT_Q0", "init");
    ModuleBase::timer::start("DFPT_Q0", "init");
    ucell_ = &ucell;
    pw_rho_ = pw_rho;
    pw_wfc_ = pw_wfc;
    pert_ = pert;
    stars_.clear();
    ModuleBase::timer::end("DFPT_Q0", "init");
}

void DFPT_Q0::build_stars(int nk)
{
    ModuleBase::TITLE("DFPT_Q0", "build_stars");
    ModuleBase::timer::start("DFPT_Q0", "build_stars");
    // every k starts with the identity member (also the permanent fallback)
    stars_.assign(nk, std::vector<StarMember>(1, StarMember()));
    if (ucell_ == nullptr || pw_wfc_ == nullptr || pw_wfc_->kvec_d == nullptr)
    {
        ModuleBase::timer::end("DFPT_Q0", "build_stars");
        return;
    }
    const ModuleSymmetry::Symmetry& symm = ucell_->symm;
    if (symm.nrotk <= 0)
    {
        // no point-group analysis (symmetry off / unreduced mesh): the
        // stored list is already the full mesh, identity members only
        ModuleBase::timer::end("DFPT_Q0", "build_stars");
        return;
    }
    const int nat = ucell_->nat;
    std::vector<ModuleBase::Vector3<double>> kfolds;
    for (int ik = 0; ik < nk; ++ik)
    {
        kfolds.clear();
        // the pre-filled identity member owns the folded k itself
        ModuleBase::Vector3<double> k0 = pw_wfc_->kvec_d[ik];
        k0.x -= std::round(k0.x);
        k0.y -= std::round(k0.y);
        k0.z -= std::round(k0.z);
        kfolds.push_back(k0);
        for (int j = 0; j < symm.nrotk; ++j)
        {
            ModuleBase::Vector3<double> kp = pw_wfc_->kvec_d[ik] * symm.kgmatrix[j];
            // fold to [-0.5, 0.5): star members are grid points, the
            // folded coordinates identify the distinct mesh points
            kp.x -= std::round(kp.x);
            kp.y -= std::round(kp.y);
            kp.z -= std::round(kp.z);
            if (!is_new_member(kfolds, kp))
            {
                continue;
            }
            kfolds.push_back(kp);
            StarMember mem;
            // cartesian form of the same operation: k_frac' = k_frac * K,
            // k_cart = k_frac * G, hence k_cart' = k_cart * (G^-1 K G). That
            // product is the row-convention operator; rotate_tensor applies
            // the column form chi' = R chi R^T, so store the transpose
            const ModuleBase::Matrix3 krow = ucell_->G.Inverse() * symm.kgmatrix[j] * ucell_->G;
            mem.cart = ModuleBase::Matrix3(krow.e11,
                                           krow.e21,
                                           krow.e31,
                                           krow.e12,
                                           krow.e22,
                                           krow.e32,
                                           krow.e13,
                                           krow.e23,
                                           krow.e33);
            if (!make_atom_map(*ucell_, symm, j, mem.atom_map))
            {
                // inconsistent operation set: fall back to identity-only
                // stars for every k (the unreduced-sum behavior)
                stars_.assign(nk, std::vector<StarMember>(1, StarMember()));
                ModuleBase::timer::end("DFPT_Q0", "build_stars");
                return;
            }
            stars_[ik].push_back(mem);
        }
    }
    ModuleBase::timer::end("DFPT_Q0", "build_stars");
}

void DFPT_Q0::rotate_tensor(const ModuleBase::Matrix3& r, const ModuleBase::matrix& chi, double (&chi_rot)[9])
{
    ModuleBase::TITLE("DFPT_Q0", "rotate_tensor");
    ModuleBase::timer::start("DFPT_Q0", "rotate_tensor");
    for (int a = 0; a < 3; ++a)
    {
        for (int b = 0; b < 3; ++b)
        {
            double s = 0.0;
            for (int ap = 0; ap < 3; ++ap)
            {
                for (int bp = 0; bp < 3; ++bp)
                {
                    s += me(r, a, ap) * me(r, b, bp) * chi(ap, bp);
                }
            }
            chi_rot[3 * a + b] = s;
        }
    }
    ModuleBase::timer::end("DFPT_Q0", "rotate_tensor");
}

void DFPT_Q0::compute_eps(const ModuleBase::matrix& wg, DFPT_PW_Data& data)
{
    ModuleBase::TITLE("DFPT_Q0", "compute_eps");
    ModuleBase::timer::start("DFPT_Q0", "compute_eps");
    if (ucell_ == nullptr)
    {
        ModuleBase::timer::end("DFPT_Q0", "compute_eps");
        return;
    }
    const int nk = wg.nr;

    // bare position legs Y^a and converged E-field responses dpsi^E,b of
    // the q = 0 mesh (DFPT_PW::solve_pos_resp / solve_efield_resp)
    std::vector<KMeshCoeffs> yr(3);
    std::vector<KMeshCoeffs> de(3);
    for (int a = 0; a < 3; ++a)
    {
        yr[a] = data.get_pos_resp(a);
        de[a] = data.get_dpsi_efield(a);
        if (static_cast<int>(yr[a].size()) != nk || static_cast<int>(de[a].size()) != nk)
        {
            ModuleBase::timer::end("DFPT_Q0", "compute_eps");
            return; // responses not solved: nothing to accumulate
        }
    }

    build_stars(nk);
    std::vector<ModuleBase::matrix> chi_k(nk, ModuleBase::matrix(3, 3, true));
    accumulate_chi_eps(wg, yr, de, nk, chi_k);
    ModuleBase::matrix eps(3, 3, true);
    // wg carries the full k weight (star size included) times the spin
    // factor 2, so the star-averaged partials sum to the complete
    // Brillouin-zone average: no extra 1/nk normalization
    star_average_eps(stars_, chi_k, eps);
    for (int a = 0; a < 3; ++a)
    {
        for (int b = 0; b < 3; ++b)
        {
            // 16 pi / Omega: QE dielec.f90 form eps = 1 - 4*(4pi/Omega)*wk*
            // Re<Y^a|dpsi^E,b> (validated against QE 7.2 Si to 0.06%:
            // 23.6825 here vs 23.6685 QE)
            eps(a, b) *= -16.0 * ModuleBase::PI / ucell_->omega;
            if (a == b)
            {
                eps(a, b) += 1.0;
            }
        }
    }
    data.set_dielectric(eps);
    ModuleBase::timer::end("DFPT_Q0", "compute_eps");
}

void DFPT_Q0::compute_born(const psi::Psi<std::complex<double>>& psi,
                           const ModuleBase::matrix& wg,
                           const ModuleBase::matrix& eig,
                           DFPT_PW_Data& data)
{
    ModuleBase::TITLE("DFPT_Q0", "compute_born");
    ModuleBase::timer::start("DFPT_Q0", "compute_born");
    if (ucell_ == nullptr)
    {
        ModuleBase::timer::end("DFPT_Q0", "compute_born");
        return;
    }
    const int nk = psi.get_nk();
    const int nbasis = psi.get_nbasis();
    (void)eig;

    // solved position legs Y^a_{k,v} = P_c x_a|psi_v> of the q = 0 mesh
    // (DFPT_PW::solve_pos_resp stashes them per direction)
    std::vector<KMeshCoeffs> yr(3);
    for (int a = 0; a < 3; ++a)
    {
        yr[a] = data.get_pos_resp(a);
        if (static_cast<int>(yr[a].size()) != nk)
        {
            ModuleBase::timer::end("DFPT_Q0", "compute_born");
            return; // position responses not solved: nothing to accumulate
        }
    }

    build_stars(nk);
    // star-rotated electronic partials, credited to the image atom under
    // each star member: zacc[kappa](a, idir)
    const int nat = ucell_->nat;
    std::vector<ModuleBase::matrix> zacc(nat, ModuleBase::matrix(3, 3, true));

    for (int iat = 0; iat < nat; ++iat)
    {
        // wg-weighted partial chi_k[ik](a, idir) of THIS atom at every k
        std::vector<ModuleBase::matrix> chi_k(nk, ModuleBase::matrix(3, 3, true));
        for (int idir = 0; idir < 3; ++idir)
        {
            // converged screened displacement response dpsi(scf)/du of this
            // mode, stashed by solve_displacement before compute_born runs
            const KMeshCoeffs disp = data.get_dpsi_disp(iat, idir);
            if (static_cast<int>(disp.size()) != nk)
            {
                continue;
            }
            accumulate_disp_chi(wg, yr, disp, idir, nbasis, chi_k);
        }
        credit_star_born(stars_, chi_k, iat, zacc);
    }

    set_born_charges(*ucell_, zacc, data);
    ModuleBase::timer::end("DFPT_Q0", "compute_born");
}

void DFPT_Q0::compute_q0_response(DFPT_PW_Data& data)
{
    ModuleBase::TITLE("DFPT_Q0", "compute_q0_response");
    ModuleBase::timer::start("DFPT_Q0", "compute_q0_response");
    // DFT+U reservation (U0): V_U is nonlocal (onsite projector), so the
    // position operator does NOT commute with the DFT+U potential. The
    // [r, V_U] commutator term must be handled separately in addition to
    // the occupation-matrix response (docc) when u_active() runs; this is
    // the hardest DFT+U piece and is deferred with the Plus_U wiring.
    (void)data;
    ModuleBase::timer::end("DFPT_Q0", "compute_q0_response");
}

} // namespace ModuleDFPT
