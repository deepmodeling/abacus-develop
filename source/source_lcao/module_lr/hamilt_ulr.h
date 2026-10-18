#ifndef ABACUS_SOURCE_LCAO_MODULE_LR_HAMILT_ULR_H
#define ABACUS_SOURCE_LCAO_MODULE_LR_HAMILT_ULR_H

#include "source_psi/psi.h"
#include "source_hamilt/operator.h"
#include <functional>
#include <memory>
#include <string>
#include <vector>
class UnitCell;
class Grid_Driver;
class K_Vectors;
class Parallel_2D;
class Parallel_Orbitals;
namespace ModuleBase
{
class matrix;
}
namespace module_dm
{
template <typename TK, typename TR>
class DensityMatrix;
}
#ifdef __EXX
template<typename T> class Exx_LRI;
#endif
namespace LR
{
    class PotHxcLR;
    /// Unristricted TDDFT (TDA) for open-shell systems
    /// The A matrix is diveded by 4 blocks: uu, ud, du, dd
    template<typename T>
    class HamiltULR
    {
    public:
        HamiltULR(const std::string& xc_kernel,
            const int& nspin,
            const int& naos,
            const std::vector<int>& nocc,   ///< {up, down}
            const std::vector<int>& nvirt,   ///< {up, down}
            const UnitCell& ucell_in,
            const std::vector<double>& orb_cutoff,
            Grid_Driver& gd_in,
            const psi::Psi<T>& psi_ks_in,
            const ModuleBase::matrix& eig_ks,
#ifdef __EXX
            std::weak_ptr<Exx_LRI<T>> exx_lri_in,
            const double& exx_alpha,
#endif
            std::vector<std::shared_ptr<PotHxcLR>>& pot_in,
            const K_Vectors& kv_in,
            const std::vector<Parallel_2D>& pX_in,   ///< {up, down}
            const Parallel_2D& pc_in,
            const Parallel_Orbitals& pmat_in);
        ~HamiltULR();
        void hPsi(const T* const psi_in, T* const hpsi, const int ld_psi, const int& nband) const;
        std::vector<T> matrix() const;

    private:
        const std::vector<int>& nocc;
        const std::vector<int>& nvirt;

        const std::vector<Parallel_2D>& pX;


        int nk = 1;
        int ldim = 1;
        int gdim = 1;

        /// 4 operator lists: uu, ud, du, dd
        std::vector<hamilt::Operator<T>*> ops;

        /// transition density matrix in AO representation
        /// Hxc only: size=1, calculate on the same address for each bands
        /// Hxc+Exx: size=nbands, store the result of each bands for common use
        std::unique_ptr<module_dm::DensityMatrix<T, T>> DM_trans;

        std::function<void(const int&, const T* const)> cal_dm_trans;
        const bool tdm_sym = false;     ///< whether to symmetrize the transition density matrix
    };
}

#endif // ABACUS_SOURCE_LCAO_MODULE_LR_HAMILT_ULR_H
