#ifndef MODULEHAMILTPW_ONSITEPROJECTOR_H
#define MODULEHAMILTPW_ONSITEPROJECTOR_H
#include "module_base/module_device/device.h"
#include "module_hsolver/kernels/math_kernel_op.h"
#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"
#include "module_basis/module_pw/pw_basis_k.h"
#include "module_hamilt_pw/hamilt_pwdft/radial_proj.h"
#include "module_psi/psi.h"
#include "module_hamilt_pw/hamilt_pwdft/onsite_proj_tools.h"

#include <string>
#include <vector>
#include <complex>
namespace projectors
{
    template <typename T, typename Device>
    class OnsiteProjector
    {
        public:

        /**
         * @brief initialize the radial projector for real-space projection involving operators
         * 
         * @param orbital_dir You know what it is
         * @param orb_files You know what it is
         * @param nproj # of projectors for each type defined in UnitCell, can be zero
         * @param lproj angular momentum for each projector
         * @param iproj index of zeta function that each projector generated from
         * @param onsite_r onsite-radius for all valid projectors
         * @param rgrid [out] the radial grid shared by all projectors
         * @param projs [out] projectors indexed by `iproj`
         * @param it2iproj [out] for each type, the projector index (across all types)
         */
        void init_proj(const std::string& orbital_dir,
                       const std::vector<std::string>& orb_files,
                       const std::vector<int>& nproj,           // for each type, the number of projectors
                       const std::vector<int>& lproj,           // angular momentum of projectors within the type (l of zeta function)
                       const std::vector<int>& iproj,           // index of projectors within the type (izeta)
                       const std::vector<double>& onsite_r); // for each type, the projector index (across all types)

        /**
         * @brief calculate the onsite projectors in reciprocal space(|G+K>) for all atoms
         */
        void tabulate_atomic(const int ik, const char grad = 'n');
        
        void overlap_proj_psi(
                    const int npm,
                    const std::complex<double>* ppsi
                    );
        void read_abacus_orb(std::ifstream& ifs,
                            std::string& elem,
                            double& ecut,
                            int& nr,
                            double& dr,
                            std::vector<int>& nzeta,
                            std::vector<std::vector<double>>& radials,
                            const int rank = 0);
        /// @brief static access to this class instance
        static OnsiteProjector<T, Device>* get_instance();
        void init(const std::string& orbital_dir,
                    const UnitCell* ucell_in,
                    const psi::Psi<std::complex<T>, Device>& psi,
                    const K_Vectors& kv,
                    const ModulePW::PW_Basis_K& pw_basis,             // level1: the plane wave basis, need ik
                    Structure_Factor& sf,                              // level2: the structure factor calculator
                    const double onsite_radius,
                    const int nq,
                    const double dq,
                    const ModuleBase::matrix& wg,
                    const ModuleBase::matrix& ekb);
        
        /// @brief calculate and print the occupations of all lm orbitals
        void cal_occupations(const psi::Psi<std::complex<T>, Device>* psi, const ModuleBase::matrix& wg_in);

        int get_size_becp() const { return size_becp; }
        std::complex<double>* get_becp() const { return becp; }
        std::complex<double>* get_h_becp() const { return h_becp; }
        std::complex<double>* get_tab_atomic() const { return tab_atomic_; }
        int get_tot_nproj() const { return tot_nproj; }
        int get_npw() const { return npw_; }
        int get_npwx() const { return npwx_; }
        const int& get_nh(int iat) const { return iat_nh[iat]; }

        hamilt::Onsite_Proj_tools<T, Device>* get_fs_tools() const { return fs_tools; }

        private:
        OnsiteProjector(){};
        ~OnsiteProjector();

        Device* ctx = {};
        base_device::DEVICE_CPU* cpu_ctx = {};
        base_device::AbacusDevice_t device = {};
        static OnsiteProjector<T, Device> *instance;

        hamilt::Onsite_Proj_tools<T, Device>* fs_tools = nullptr;

        std::complex<double>* tab_atomic_ = nullptr;
        std::complex<double>* becp = nullptr;  // nbands * nkb
        // save becp in CPU memory, only used when Device is GPU
        std::complex<double>* h_becp;

        int size_becp = 0;
        int size_vproj = 0;
        int tot_nproj = 0;
        int npw_ = 0;
        int npwx_ = 0;
        int ik_ = 0;
        std::vector<std::vector<int>> it2ia;
        std::vector<double> rgrid;
        std::vector<std::vector<double>> projs;
        std::vector<std::vector<int>> it2iproj;
        std::vector<int> lproj;
        std::vector<int> iat_nh;

        const UnitCell* ucell = nullptr;

        const ModulePW::PW_Basis_K* pw_basis_ = nullptr;             // level1: the plane wave basis, need ik
        Structure_Factor* sf_ = nullptr;                             // level2: the structure factor calculator
        int ntype = 0;

        RadialProjection::RadialProjector rp_;
        std::vector<int> irow2it_;
        std::vector<int> irow2iproj_;
        std::vector<int> irow2m_;
        std::map<std::tuple<int, int, int, int>, int> itiaiprojm2irow_;

        ModuleBase::realArray tab;
        ModuleBase::matrix nhtol;

        bool initialed = false;

        /// @brief rename the operators for CPU/GPU device
        using gemm_op = hsolver::gemm_op<std::complex<T>, Device>;

        using resmem_complex_op = base_device::memory::resize_memory_op<std::complex<T>, Device>;
        using resmem_complex_h_op = base_device::memory::resize_memory_op<std::complex<T>, base_device::DEVICE_CPU>;
        using setmem_complex_op = base_device::memory::set_memory_op<std::complex<T>, Device>;
        using delmem_complex_op = base_device::memory::delete_memory_op<std::complex<T>, Device>;
        using delmem_complex_h_op = base_device::memory::delete_memory_op<std::complex<T>, base_device::DEVICE_CPU>;
        using syncmem_complex_h2d_op
            = base_device::memory::synchronize_memory_op<std::complex<T>, Device, base_device::DEVICE_CPU>;
        using syncmem_complex_d2h_op
            = base_device::memory::synchronize_memory_op<std::complex<T>, base_device::DEVICE_CPU, Device>;

        using resmem_var_op = base_device::memory::resize_memory_op<T, Device>;
        using resmem_var_h_op = base_device::memory::resize_memory_op<T, base_device::DEVICE_CPU>;
        using setmem_var_op = base_device::memory::set_memory_op<T, Device>;
        using delmem_var_op = base_device::memory::delete_memory_op<T, Device>;
        using delmem_var_h_op = base_device::memory::delete_memory_op<T, base_device::DEVICE_CPU>;
        using syncmem_var_h2d_op = base_device::memory::synchronize_memory_op<T, Device, base_device::DEVICE_CPU>;
        using syncmem_var_d2h_op = base_device::memory::synchronize_memory_op<T, base_device::DEVICE_CPU, Device>;

        using resmem_int_op = base_device::memory::resize_memory_op<int, Device>;
        using delmem_int_op = base_device::memory::delete_memory_op<int, Device>;
        using syncmem_int_h2d_op = base_device::memory::synchronize_memory_op<int, Device, base_device::DEVICE_CPU>;
    };
}// namespace projectors

#endif // MODULEHAMILTPW_ONSITEPROJECTOR_H