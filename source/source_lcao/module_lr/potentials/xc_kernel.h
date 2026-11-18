#ifndef ABACUS_SOURCE_LCAO_MODULE_LR_POTENTIALS_XC_KERNEL_H
#define ABACUS_SOURCE_LCAO_MODULE_LR_POTENTIALS_XC_KERNEL_H

#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/unitcell.h"
#include "source_base/parallel_grid.h"
#include "source_estate/module_charge/charge.h"
#include <stdexcept>
#include <memory>
#define CREF(x) const std::vector<double>& x = x##_
namespace LR
{
    /// libxc component layout for the spin-POLARIZED case, shared by the closed- and open-shell
    /// $g^{xc}$ code. Spin: 0=u, 1=d. Sigma: 0=uu, 1=ud, 2=dd.
    /// libxc stores only the unique (unordered) index combinations, which is where all the
    /// multiplicities in the contractions come from.
    namespace libxc_idx
    {
        /// unordered sigma PAIR -> v2sigma2 (6) and the per-rho block of v3rhosigma2
        constexpr int p2[3][3] = { {0,1,2},{1,3,4},{2,4,5} };
        /// unordered sigma TRIPLE -> v3sigma3 (10):
        /// (000)(001)(002)(011)(012)(022)(111)(112)(122)(222)
        constexpr int p3[3][3][3] = {
            { {0,1,2},{1,3,4},{2,4,5} },
            { {1,3,4},{3,6,7},{4,7,8} },
            { {2,4,5},{4,7,8},{5,8,9} } };
        /// v3rho3 (4): (uuu,uud,udd,ddd) -- indexed by the number of d's
        inline constexpr int r3(const int s0, const int s1, const int s2)
        {
            return s0 + s1 + s2;
        }
        /// v2rho2 (3): (uu,ud,dd)
        inline constexpr int r2(const int s0, const int s1)
        {
            return s0 + s1;
        }
        /// v2rhosigma (6): [rho u,d] x [sigma uu,ud,dd]
        inline constexpr int rs(const int s, const int a)
        {
            return 3 * s + a;
        }
        /// v3rho2sigma (9): [rho-pair uu,ud,dd] x [sigma uu,ud,dd]
        inline constexpr int r2s(const int s0, const int s1, const int a)
        {
            return 3 * (s0 + s1) + a;
        }
        /// v3rhosigma2 (12): [rho u,d] x [unordered sigma pair]
        inline constexpr int rs2(const int s, const int a, const int b)
        {
            return 6 * s + p2[a][b];
        }
        /// $\partial\sigma_a/\partial\nabla\rho_\tau = \theta^\tau_a\,\nabla\rho_{c(\tau,a)}$
        /// tau=u: (uu -> 2 grad rho_u, ud -> 1 grad rho_d, dd -> 0)
        /// tau=d: (uu -> 0,            ud -> 1 grad rho_u, dd -> 2 grad rho_d)
        constexpr double theta[2][3] = { {2., 1., 0.}, {0., 1., 2.} };
        /// which density-gradient channel goes with (tau, a); -1 where theta vanishes
        constexpr int chan[2][3] = { {0, 1, -1}, {-1, 0, 1} };
    }

    /// Derivative of the local energy's sigma-channel partial derivative along a density
    /// perturbation. Input coefficients are the packed LibXC arrays for ONE grid point;
    /// density_response has two components and sigma_response has (uu, ud, dd).
    /// The half-open spin range selects one input spin or both. Only structurally active
    /// sigma channels are read, avoiding artificial zero-times-NaN contributions.
    /// Used by the unrestricted second-order potential contraction.
    inline double sigma_potential_derivative(const double* rho_sigma, const double* sigma_sigma,
        const int channel, const double* density_response, const double* sigma_response,
        const int first_spin, const int end_spin)
    {
        double response = 0.;
        for (int spin = first_spin; spin < end_spin; ++spin)
        {
            response += rho_sigma[libxc_idx::rs(spin, channel)] * density_response[spin];
        }
        for (int a = first_spin; a < end_spin + 1; ++a)
        {
            response += sigma_sigma[libxc_idx::p2[channel][a]] * sigma_response[a];
        }
        return response;
    }

    /// @brief Calculate the exchange-correlation (XC) kernel ($f_{xc}=\delta^2E_xc/\delta\rho^2$) and store its components.
    class KernelXC
    {
        using Tvec = std::vector<double>;
    public:
        /// Which spin combinations of the $g^{xc}$ coefficient set (`GxcCoef`) to build.
        enum GxcSpin { NoGxc = 0, Singlet = 1, Triplet = 2, BothSpins = 3 };
        /// GS response starts from raw sigma; non-finite HSE coefficients trigger a reported fallback.
        enum class Purpose { LinearResponse, GroundStateResponse };
        KernelXC(const ModulePW::PW_Basis& rho_basis,
            const UnitCell& ucell,
            const Charge& chg_gs,
            const Parallel_Grid& pgrid,
            const int& nspin,
            const std::string& kernel_name,
            const std::vector<std::string>& lr_init_xc_kernel,
            const bool openshell = false,
            const int gxc_spin = GxcSpin::NoGxc);
        /// Explicit response role; GS response requests only second derivatives.
        /// Density and grid are borrowed and never modified; their spin layouts must agree.
        KernelXC(Purpose purpose, const ModulePW::PW_Basis& rho_basis,
            const UnitCell& ucell, const Charge& chg_gs, const Parallel_Grid& pgrid,
            const int& nspin, const std::string& kernel_name,
            const std::vector<std::string>& lr_init_xc_kernel, bool openshell, int gxc_spin);
        ~KernelXC() {}

        /// Metadata belongs to this coefficient set, independent of later global XC changes.
        int nspin() const
        {
            return nspin_;
        }
        int functional_type() const
        {
            return functional_type_;
        }
        const std::string& name() const
        {
            return kernel_name_;
        }

        /// Construction metadata: true only after a collective HSE GS-response fallback.
        /// This flags an approximate orbital Hessian; the SCF potential is unchanged.
        bool gs_sigma_fallback_used() const { return gs_sigma_fallback_used_; }

        // const references
        CREF(vrho);CREF(vsigma); CREF(v2rho2); CREF(v2rhosigma); CREF(v2sigma2);
        CREF(v3rho3);
        CREF(v3rho2sigma);
        CREF(v3rhosigma2);
        CREF(v3sigma3);
        /// Scalar mixed-density/gradient and gradient/gradient GGA response coefficients.
        /// The common ground-state gradient is supplied by the potential contraction.
        /// These same arrays supply the transition-gradient terms of gxc; no duplicate cache.
        struct ClosedShellGga
        {
            Tvec rho_sigma;
            Tvec sigma_sigma;
        };
        const ClosedShellGga& closed_shell_gga(const bool triplet) const
        {
            if (openshell_)
            {
                throw std::logic_error("Closed-shell GGA coefficients requested for an open-shell kernel");
            }
            return (nspin_ == 1 || !triplet) ? closed_shell_gga_s_ : closed_shell_gga_t_;
        }

        /// @brief The coefficients of
        ///     $v^{(2)}(r)=\iint dr'dr''\,g^{xc}(r,r',r'')\rho^1(r')\rho^1(r'')$
        /// for ONE spin combination, stored exactly as they appear in the final formula:
        ///     $v^{(2)} = a_{s^2}s^2 + a_{st}\,s\,t + a_{t^2}t^2 + a_q\,q
        ///        - \nabla\cdot[\,(e_{s^2}s^2 + e_{st}\,s\,t + e_{t^2}t^2 + e_q\,q)\,\nabla\rho
        ///          + (c_s\,s + c_t\,t)\,\nabla\rho^1\,]$
        /// with $s=\rho^1$, $t=\nabla\rho\cdot\nabla\rho^1$, $q=\nabla\rho^1\cdot\nabla\rho^1$.
        ///
        /// The transition-gradient scalars c_s and c_t reuse ClosedShellGga (factor 2 for
        /// nspin=2, factor 8 for singlet-normalized nspin=1). The other factors are folded
        /// in here, so `PotGradXCLR::cal_v_eff`
        /// is a literal transcription of the formula with no arithmetic of its own, and nspin=1,
        /// singlet and triplet all run through the same code. For LDA only `a_s2` is filled.
        /// $v^{(2)} = A - \nabla\cdot E$
        ///
        /// NOTE the four $e$ are *scalars*, with the common $\nabla\rho^{gs}$ factored out of the sum.
        /// Should an open-shell version ever need $\nabla\rho_u\ne\nabla\rho_d$ under the same
        /// divergence, this factorization no longer holds and they must go back to `Vector3`.
        struct GxcCoef
        {
            /// Availability is independent of the number of points assigned to this rank.
            bool built = false;
            // Local part A.
            std::vector<double> a_s2;
            std::vector<double> a_st;
            std::vector<double> a_t2;
            std::vector<double> a_q;
            // Divergence coefficients, all multiplying the ground-state density gradient.
            std::vector<double> e_s2;
            std::vector<double> e_st;
            std::vector<double> e_t2;
            std::vector<double> e_q;
        };
        /// nspin=1 has no singlet/triplet distinction, so it always returns the one set that is built.
        /// Throws instead of handing back an empty set when the requested combination was not
        /// requested at construction -- silently returning zeros would look like a physics bug.
        const GxcCoef& gxc(const bool triplet) const
        {
            const GxcCoef& ret = (nspin_ == 1 || !triplet) ? gxc_s_ : gxc_t_;
            if (!ret.built)
            {
                throw std::runtime_error(
                    "KernelXC: the " + std::string(triplet ? "triplet" : "singlet")
                    + " g^xc coefficients were not built; pass the matching `GxcSpin` flag to the constructor.");
            }
            return ret;
        }

        /// Componentwise derivative of prepared sigma; empty means identity.
        /// Borrowed for the kernel lifetime, point-major: one component at nspin=1,
        /// or (uu, ud, dd) at nspin=2. Only implicit coefficient derivatives use it.
        const Tvec& sigma_jacobian() const { return sigma_jacobian_; }

        const bool& openshell = openshell_;
        const std::vector<std::vector<ModuleBase::Vector3<double>>>& drho_gs = drho_gs_;
        /// Whether THIS kernel (built for its own functional name, which may differ from the
        /// ground state's `dft_functional` in a cross-functional run such as TDLDA@PBE) needs
        /// the GGA gradient terms. `drho_gs_` is only ever filled when this kernel's own `is_gga`
        /// was true at construction (see `f_xc_libxc`), so its emptiness is a reliable per-kernel
        /// proxy -- unlike the global `XC_Functional::get_func_type()`, which reflects the
        /// ground state's functional and disagrees with this kernel whenever the two differ.
        bool is_gga() const
        {
            return !this->drho_gs_.empty();
        }

    private:
#ifdef __LIBXC
        /// @brief Calculate the XC kernel using libxc.
        void f_xc_libxc(const int& nspin, const double& omega, const double& tpiba, const double* const* const rho_gs, const double* const rho_core = nullptr);
        /// calculate the input rho, grad rho, and sigma for libxc
        void get_rho_drho_sigma(const int& nspin,
            const double& tpiba,
            const double* const* const rho_gs,
            const double* const rho_core,
            const bool& is_gga,
            std::vector<double>& rho,
            std::vector<std::vector<ModuleBase::Vector3<double>>>& gradrho,
            std::vector<double>& sigma);
#endif
        // See https://libxc.gitlab.io/manual/libxc-5.1.x/ for the naming convention of the following members.
        // std::map<std::string, std::vector<double>> kernel_set_; // [kernel_type][nrxx][nspin]

        // ================================== XC kernels ============================================
        std::vector<double> vrho_;
        std::vector<double> vsigma_;
        std::vector<double> v2rho2_;
        std::vector<double> v2rhosigma_;
        std::vector<double> v2sigma2_;
        // std::map<std::string, std::vector<ModuleBase::Vector3<double>>> grad_kernel_set_;// [kernel_type][nrxx][nspin],  intermediate terms for GGA
        // for nspin=1, gga
        std::vector<std::vector<ModuleBase::Vector3<double>>> drho_gs_; ///< $\nabla\rho$ of the ground state, size: nspin*nrxx
        ClosedShellGga closed_shell_gga_s_;
        ClosedShellGga closed_shell_gga_t_;

        // ================================== XC kernels ============================================
        // ================================== XC kernel Gradiants ====================================
        Tvec v3rho3_;
        Tvec v3rho2sigma_;
        Tvec v3rhosigma2_;
        Tvec v3sigma3_;

        // The two spin combinations of $v^{(2)}$'s coefficients (see `GxcCoef` above).
        // `gxc_t_` stays empty for nspin=1, where there is no triplet.
        GxcCoef gxc_s_;
        GxcCoef gxc_t_;
        /// @brief Fill `dst` for one spin combination. All the spin algebralives here, driven by the weight
        /// vectors that distinguish singlet from triplet -- the two differ only in those weights.
        void build_gxc_coef(GxcCoef& dst, const bool triplet, const int& nspin, const bool& is_gga);
        const int nspin_;
        const std::string kernel_name_;
        int functional_type_ = 0;
        const int gxc_spin_ = GxcSpin::NoGxc;   ///< which `GxcCoef` sets to build, see `GxcSpin`
        const Purpose purpose_;
        Tvec sigma_jacobian_;
        bool gs_sigma_fallback_used_ = false;
        // ================================== XC kernel Gradiants ====================================
        const ModulePW::PW_Basis& rho_basis_;
        const bool openshell_ = false;
    };

    /// Explicitly evaluate a coefficient set for one functional and ground-state density.
    /// Sharing potentials must use compatible spin layouts and derivative requirements.
    std::shared_ptr<const KernelXC> build_xc_kernel(const std::string& kernel_name,
        const ModulePW::PW_Basis& rho_basis, const UnitCell& ucell, const Charge& chg_gs,
        const Parallel_Grid& pgrid, const int nspin, const bool openshell, const int gxc_spin,
        const std::vector<std::string>& lr_init_xc_kernel);
    /// Build immutable GS response from raw coefficients, retrying non-finite HSE with the LR floor.
    /// Uses the supplied GS functional, density and grid, requesting no third derivatives.
    std::shared_ptr<const KernelXC> build_gs_xc_kernel(const std::string& kernel_name,
        const ModulePW::PW_Basis& rho_basis, const UnitCell& ucell, const Charge& chg_gs,
        const Parallel_Grid& pgrid, int nspin, bool openshell,
        const std::vector<std::string>& lr_init_xc_kernel);

}

#endif // ABACUS_SOURCE_LCAO_MODULE_LR_POTENTIALS_XC_KERNEL_H
