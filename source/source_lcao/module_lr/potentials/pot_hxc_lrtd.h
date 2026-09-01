#pragma once
#include "source_estate/module_pot/h_hartree_pw.h"
#include "xc_kernel.h"
#include "pot_lr_base.h"
#include <unordered_map>
#include <memory>

namespace LR
{
    class PotHxcLR : public PotLRBase
    {
    public:
        /// S1: the nspin=1 singlet LR kernel, 2*K^Hartree + 2*K^xc.
        ///     The factor 2 on *both* terms is what makes it equal to S2_singlet. 
        /// S1_gs: the nspin=1 counterpart of S2_gs, i.e. the *ground-state* Hxc kernel
        ///     K^Hartree + K^xc = S1 / 2. This is what `pot_hxc_gs` needs at nspin=1.
        /// S2_singlet: 2*K^Hartree + K^xc_{upup} + K^xc_{updown}
        /// S2_triplet: K^xc_{upup} - K^xc_{updown}
        /// S2_updown: K^Hartree + (K^xc_{upup}, K^xc_{updown},  K^xc_{downup} or K^xc_{downdown}), according to `ispin_op` (for spin-polarized systems)
        /// S2_gs: the nspin=2 *ground-state* Hxc kernel K^Hartree + (K^xc_{upup} + K^xc_{updown})/2 = S2_singlet / 2.
        ///     Used for `pot_hxc_gs` in LR gradients, where the convention is `K_Hxc(singlet) = 2 * pot_hxc_gs`.
        ///     The 1/2 on the xc part is not a convention but the chain rule: the derivative is
        ///     taken w.r.t. the *total* density matrix, and $\partial v_u/\partial\rho =
        ///     (f_{uu}+f_{ud})/2$ because $\rho_u=\rho_d=\rho/2$. The Hartree part needs no halving.
        ///     Do NOT use S1 here when nspin=2: `KernelXC` is built with `PARAM.inp.nspin`, so the
        ///     kernel arrays carry 3 spin components per grid point while the S1 integrand indexes
        ///     them as if there were 1.
        enum SpinType { S1 = 0, S2_singlet = 1, S2_triplet = 2, S2_updown = 3, S2_gs = 4, S1_gs = 5 };
        /// XCType here is to determin the method of integration from kernel to potential, not the way calculating the kernel
        enum XCType { None = 0, LDA = 1, GGA = 2, HYB_GGA = 4 };
        /// constructor building exchange-correlation kernel
        PotHxcLR(const std::string& xc_kernel, const ModulePW::PW_Basis& rho_basis,
            const UnitCell& ucell, const Charge& chg_gs/*ground state*/, const Parallel_Grid& pgrid,
            const SpinType& st = SpinType::S1, const std::vector<std::string>& lr_init_xc_kernel = { "default" },
            const int gxc_spin = KernelXC::GxcSpin::NoGxc);
        /// Constructor taking an already-built kernel. Several `PotHxcLR` can share the same* $f^{xc}$ arrays.
        /// The caller is responsible for the ordering: `KernelXC` calls `XC_Functional::set_xc_type`,
        /// and this constructor reads the resulting global `get_func_type()`, so build the kernel
        /// immediately before the potentials that use it.
        PotHxcLR(std::shared_ptr<const KernelXC> kernel, const std::string& xc_kernel,
            const ModulePW::PW_Basis& rho_basis, const UnitCell& ucell, const int nrxx,
            const SpinType& st = SpinType::S1);
        ~PotHxcLR() {}
        virtual void cal_v_eff(double** rho, const UnitCell& ucell, ModuleBase::matrix& v_eff, const std::vector<int>& ispin_op = { 0,0 })  const override;

        /// Build a kernel that can be shared by several `PotHxcLR` (see the constructor above).
        /// `openshell` must match what every sharing potential would have passed, i.e.
        /// `st == SpinType::S2_updown`; `gxc_spin` must cover every combination they will ask for.
        static std::shared_ptr<const KernelXC> make_kernel(const std::string& xc_kernel,
            const ModulePW::PW_Basis& rho_basis, const UnitCell& ucell, const Charge& chg_gs,
            const Parallel_Grid& pgrid, const bool openshell, const int gxc_spin,
            const std::vector<std::string>& lr_init_xc_kernel = { "default" });

        const KernelXC& xc_kernel_components() const { return *xc_kernel_components_; }
    private:
        std::unique_ptr<elecstate::PotHartree> pot_hartree_;
        /// different components of local and semi-local xc kernels:
        /// LDA: v2rho2
        /// GGA: v2rho2, v2rhosigma, v2sigma2
        /// meta-GGA: v2rho2, v2rhosigma, v2sigma2, v2rholap, v2rhotau, v2sigmalap, v2sigmatau, v2laptau, v2lap2, v2tau2
        /// To allow different potential objects sharing the same kernel.
        const std::shared_ptr<const KernelXC> xc_kernel_components_;
        const std::string xc_kernel_;
        const SpinType spin_type_ = SpinType::S1;
        XCType xc_type_ = XCType::None;

        // enum class as key for unordered_map is not supported in C++11 sometimes
        // https://github.com/llvm/llvm-project/issues/49601
        // struct SpinHash { std::size_t operator()(const SpinType& s) const { return std::hash<int>()(static_cast<int>(s)); } };
        using Tfunc = std::function<void(const double* const /**<[in] rho*/,
            ModuleBase::matrix& /**<[out] v_eff */,
            const std::vector<int>& ispin_op) >;
        // std::unordered_map<SpinType, Tfunc, SpinHash> kernel_to_potential_;
        std::map<SpinType, Tfunc> kernel_to_potential_;

        void set_integral_func(const SpinType& s, const XCType& xc);

        // ---- scratch buffers ------------------------------------------------------------
        // `cal_v_eff` used to allocate (and value-initialize) every temporary on every call:
        // at a 200^3 grid that is ~450 MB of memset per call, all of it overwritten before it
        // is read, plus the page faults of freshly mapped pages. They now come from a pool
        // that grows on demand and is SHARED by every `PotHxcLR` (a calculation holds three of
        // them -- singlet, triplet and the ground-state kernel -- and giving each its own copy
        // cost ~0.9 GB on the bigger grids). Sharing is safe because `cal_v_eff` is only ever
        // entered from a single thread: all the OpenMP lives inside the loops it calls.
        struct Scratch
        {
            std::vector<ModuleBase::Vector3<double>> drho;   ///< $\nabla\rho^X$
            std::vector<ModuleBase::Vector3<double>> gdot;   ///< integrand of the divergence
            std::vector<double> vxc;                         ///< accumulated $v^{xc}$
            std::vector<std::complex<double>> rhog;          ///< $\rho^X(G)$, npw
            std::vector<std::complex<double>> vg;            ///< G-space scratch, nmaxgr
            void alloc(const int nrxx, const int npw, const int nmaxgr, const bool gga);
        };
        static Scratch& scratch();
        /// Hartree potential of the transition density, built from `scratch().rhog` and
        /// accumulated into `v_eff` scaled by `factor`. Replaces `H_Hartree_pw::v_hartree`,
        /// which for the LR use case (a) re-did the forward FFT of $\rho^X$ that the GGA branch
        /// needs anyway, (b) accumulated a Hartree "energy" of the transition density and pushed
        /// it through `Parallel_Reduce::reduce_pool` on every call -- a per-call collective whose
        /// result is never read, and which clobbers the global `H_Hartree_pw::hartree_energy` --
        /// and (c) returned a full `matrix` by value, to which `v_eff += 2 * (...)` then added a
        /// second temporary of the same size.
        void add_v_hartree(const UnitCell& ucell, ModuleBase::matrix& v_eff, const double factor) const;
        // ---- pre-contracted spin combinations ------------------------------------------
        // At nspin=2 the singlet/triplet integrands read `v2rho2[3ir] +- v2rho2[3ir+1]` and
        // `2*vsigma[3ir] +- vsigma[3ir+1]` at every grid point of every call. Both combinations
        // depend only on the ground state, so they are formed once here. This also turns two
        // strided reads into one contiguous one, which is where most of the gain is.
        // Only the combination this potential's own `spin_type_` needs is built (one scalar
        // array each, so 8 B/point, and nothing at all for nspin=1 or the open-shell branch),
        // which is why they live here rather than in the shared `KernelXC`.
        mutable std::vector<double> v2rho2_comb_;
        mutable std::vector<double> vsigma_comb_;   ///< GGA only
        void build_spin_combos(const bool gga) const;
    };

} // namespace LR
