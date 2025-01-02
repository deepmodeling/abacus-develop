#pragma once
#include "module_cell/klist.h"
#include "module_lr/utils/gint_template.h"
#include "module_psi/psi.h"
#include "module_elecstate/module_dm/density_matrix.h"
#include "module_lr/utils/lr_util.h"
#include "module_basis/module_nao/two_center_bundle.h"
#include "module_hamilt_lcao/module_tddft/td_current.h"
namespace LR
{
    template<typename T>
    class LR_Spectrum
    {
    public:
        LR_Spectrum(const int& nspin_global, const int& naos, const std::vector<int>& nocc, const std::vector<int>& nvirt,
            typename TGint<T>::type* gint, const ModulePW::PW_Basis& rho_basis, psi::Psi<T>& psi_ks_in,
            const UnitCell& ucell, const K_Vectors& kv_in, const Grid_Driver& gd, const std::vector<double>& orb_cutoff,
            const TwoCenterBundle& two_center_bundle_,
            const std::vector<Parallel_2D>& pX_in, const Parallel_2D& pc_in, const Parallel_Orbitals& pmat_in,
            const double* eig, const T* X, const int& nstate, const bool& openshell,
            const std::string& gauge = "length") :
            nspin_x(openshell ? 2 : 1), naos(naos), nocc(nocc), nvirt(nvirt), nk(kv_in.get_nks() / nspin_global),
            gint(gint), rho_basis(rho_basis), ucell(ucell), kv(kv_in), gd_(gd),
            orb_cutoff_(orb_cutoff), two_center_bundle_(two_center_bundle_),
            pX(pX_in), pc(pc_in), pmat(pmat_in),
            eig(eig), X(X), nstate(nstate),
            ldim(nk* (nspin_x == 2 ? pX_in[0].get_local_size() + pX_in[1].get_local_size() : pX_in[0].get_local_size())),
            gdim(nk* std::inner_product(nocc.begin(), nocc.end(), nvirt.begin(), 0))
        {
            for (int is = 0;is < nspin_global;++is) { psi_ks.emplace_back(LR_Util::get_psi_spin(psi_ks_in, is, nk)); }
            gauge == "velocity" ? this->cal_transition_dipoles_velocity() : this->cal_transition_dipoles_length();
            this->oscillator_strength();
        };
        /// @brief calculate the optical absorption spectrum with $Im[1/[(w+i\eta)^2-\Omega_S^2]]$
        void optical_absorption_method1(const std::vector<double>& freq, const double eta);
        /// @brief calculate the optical absorption spectrum with lorentzian delta function
        void optical_absorption_method2(const std::vector<double>& freq, const double eta);
        /// @brief print out the transition dipole moment and the main contributions to the transition amplitude
        void transition_analysis(const std::string& spintype);

        //========================================== test functions ==============================================
        /// @brief write transition dipole
        void write_transition_dipole(const std::string& filename);
        /// @brief calculate transition dipole in velocity gauge using ks eigenvalues instead of excitation energies
        void test_transition_dipoles_velocity_ks(const double* const ks_eig);
        //======================================================================================================
    private:
        /// $$2/3\Omega\sum_{ia\sigma} |\braket{\psi_{i}|\mathbf{r}|\psi_{a}} |^2\int \rho_{\alpha\beta}(\mathbf{r}) \mathbf{r} d\mathbf{r}$$
        void oscillator_strength();
        /// calculate the transition dipole of state S in length gauge: $\sum_{iak}X^S_{iak}<ik|r|ak>$
        ModuleBase::Vector3<T> cal_transition_dipole_istate_length(const int istate);
        /// calculate the transition dipole of all states in length gauge
        void cal_transition_dipoles_length();
        /// calculate the transition dipole of state S in velocity gauge: $i(\sum_{iak}X^S_{iak}<ik|v|ak>)/\Omega_S$
        ModuleBase::Vector3<T> cal_transition_dipole_istate_velocity_R(const int istate, const TD_current& vR);
        ModuleBase::Vector3<T> cal_transition_dipole_istate_velocity_k(const int istate, const TD_current& vR);
        /// calculate the transition dipole of all states in velocity gauge
        void cal_transition_dipoles_velocity();
        double cal_mean_squared_dipole(ModuleBase::Vector3<T> dipole);
        /// calculate the transition density matrix
        elecstate::DensityMatrix<T, T> cal_transition_density_matrix(const int istate, const T* X_in = nullptr, const bool need_R = true);
        const int nspin_x = 1;   ///< 1 for singlet/triplet, 2 for updown(openshell)
        const int naos = 1;
        const std::vector<int>& nocc;
        const std::vector<int>& nvirt;
        const int nk = 1;
        const int nstate = 1;
        const int ldim = 1;///< local leading dimension of X, or the data size of each state
        const int gdim = 1;///< global leading dimension of X
        const double ana_thr = 0.3;     ///< {abs(X) > thr} will appear in the transition analysis log
        const double* eig;
        const T* X;
        const K_Vectors& kv;
        std::vector<psi::Psi<T>> psi_ks;
        const std::vector<Parallel_2D>& pX;
        const Parallel_2D& pc;
        const Parallel_Orbitals& pmat;
        typename TGint<T>::type* gint = nullptr;
        const ModulePW::PW_Basis& rho_basis;
        const Grid_Driver& gd_;
        const UnitCell& ucell;
        const std::vector<double>& orb_cutoff_;
        const TwoCenterBundle& two_center_bundle_;

        void cal_gint_rho(double** rho, const int& nrxx);
        std::map<std::string, int> get_pair_info(const int i); ///< given the index in X, return its ispin, ik, iocc, ivirt

        std::vector<ModuleBase::Vector3<T>> transition_dipole_;   ///< $\braket{ \psi_{i} | \mathbf{r} | \psi_{a} }$
        std::vector<double> mean_squared_transition_dipole_;    /// $|dipole|^2/3$, atomic unit (Hartree)
        std::vector<double> oscillator_strength_;///< $2/3\Omega |\sum_{ia\sigma} \braket{\psi_{i}|\mathbf{r}|\psi_{a}} |^2$, atomic unit (Hartree)
    };
}
