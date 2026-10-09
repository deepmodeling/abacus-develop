#ifndef HAMILT_LCAO_H
#define HAMILT_LCAO_H

#include "source_basis/module_nao/two_center_bundle.h"
#include "source_cell/klist.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_hamilt/hamilt.h"
#include "source_hamilt/hs_matrix_k.h"
#include "source_hamilt/module_hcontainer/hcontainer.h"

#include <memory>
#include <vector>

// elecstate::Potential forward declaration, full definition in potential_new.h (moved to .cpp)
namespace elecstate { class Potential; }

// module_dm::DensityMatrix forward declaration, full definition in density_matrix.h (moved to .cpp)
namespace module_dm { template <typename TK, typename TR> class DensityMatrix; }

// Setup_DeePKS forward declaration, full definition in setup_deepks.h (moved to .cpp)
template <typename TK> class Setup_DeePKS;
// Plus_U_Base forward declaration, full definition in source_pw/module_pwdft/dftu_base.h
class Plus_U_Base;

// Exx_NAO forward declaration, full definition in setup_exx.h (moved to .cpp)
template <typename TK> class Exx_NAO;

/// Exx_Info forward declaration, full definition in exx_info.h
struct Exx_Info;

// Input_para forward declaration, full definition in input_parameter.h
struct Input_para;

namespace hamilt
{

// OperatorLCAO forward declaration, full definition in
// module_operator_lcao/operator_lcao.h (moved to .cpp)
template <typename TK, typename TR> class OperatorLCAO;

// template first for type of k space H matrix elements
// template second for type of temporary matrix, 
// gamma_only fix-gamma-matrix + S-gamma, 
// multi-k fix-Real + S-Real
template <typename TK, typename TR>
class HamiltLCAO : public Hamilt<TK>
{
  public:

    /**
     * @brief Constructor of Hamiltonian for LCAO base
     * HR and SR will be allocated with Operators
     */
    HamiltLCAO(const UnitCell& ucell,
               const Grid_Driver& grid_d,
               const Parallel_Orbitals* paraV,
               elecstate::Potential* pot_in,
               const K_Vectors& kv_in,
               const TwoCenterBundle& two_center_bundle,
               const LCAO_Orbitals& orb,
               module_dm::DensityMatrix<TK, double>* DM_in,
               Plus_U_Base* p_dftu, // mohan add 2025-11-05
               Setup_DeePKS<TK> &deepks,
               const int istep,
               Exx_NAO<TK> &exx_nao,
               const Exx_Info& exx_info,
               const Input_para& inp,
               const bool load_exx_flag);

    /**
     * @brief Constructor of vacuum Operators, only HR and SR will be initialed as empty HContainer
     */
    HamiltLCAO(const UnitCell& ucell,
               const Grid_Driver& grid_d,
               const Parallel_Orbitals* paraV,
               const K_Vectors& kv_in,
               const TwoCenterIntegrator& intor_overlap_orb,
               const std::vector<double>& orb_cutoff);

    ~HamiltLCAO()
    {
        delete this->ops;
    }

    /// get pointer of Operator<TK> ops
    Operator<TK>*& getOperator();

    /// get H(k) pointer
    TK* getHk() const
    {
        return this->hsk->get_hk();
    }

    /// get S(k) pointer
    TK* getSk() const
    {
        return this->hsk->get_sk();
    }

    /// get HR pointer of *this->hR, which is a HContainer<TR> and contains H(R)
    HContainer<TR>* getHR()
    {
        return this->hR.get();
    }
    const HContainer<TR>* getHR() const
    {
        return this->hR.get();
    }

    /// get SR pointer of *this->sR, which is a HContainer<TR> and contains S(R)
    HContainer<TR>* getSR()
    {
        return this->sR.get();
    }
    const HContainer<TR>* getSR() const
    {
        return this->sR.get();
    }

#ifdef __MLALGO
    /// get V_delta_R pointer of *this->V_delta_R, which is a HContainer<TR> and contains V_delta(R)
    HContainer<TR>*& get_V_delta_R()
    {
        return this->V_delta_R;
    }
#endif

    /// get hRS2 buffer for NSPIN=2 case (spin-up in first half, spin-down in second half)
    std::vector<TR>& getHRS2() { return this->hRS2; }

    /// Get HR as a vector of HContainer pointers (one per spin).
    /// For nspin=2, returns pointers to internally managed per-spin wrappers over hRS2.
    /// Returned pointers are owned by this class; caller must NOT delete them.
    std::vector<HContainer<TR>*> getHR_vector();

    /// refresh the status of HR
    void refresh(bool yes) override;

    // for target K point, update consequence of hPsi() and matrix()
    void updateHk(const int ik) override;

    /**
     * @brief special for LCAO, update SK only
     *
     * @param ik target K point
     * @param kvec_d: direct coordinates of k-points
     * @param hk_type 0: SK is row-major, 1: SK is collumn-major
     * @return void
     */
    void updateSk(const int ik, const int hk_type);

    // core function: return H(k) and S(k) matrixs for direct solving eigenvalues.
    // not used in PW base
    void matrix(MatrixBlock<TK>& hk_in, MatrixBlock<TK>& sk_in) override;

  private:

    const K_Vectors* kv = nullptr;

    //! Real space Hamiltonian H(R), where R is the Bravis lattice vector
    std::unique_ptr<HContainer<TR>> hR;

    //! Real space overlap matrix S(R), where R is the Bravis lattice vector
    std::unique_ptr<HContainer<TR>> sR;

#ifdef __MLALGO
    HContainer<TR>* V_delta_R = nullptr;
#endif

    //! Hamiltonian and overlap matrices for a specific k point
    std::unique_ptr<HS_Matrix_K<TK>> hsk;

    // special case for NSPIN=2 , data of HR should be separated into two parts
    // save them in this->hRS2;
    std::vector<TR> hRS2;

    /// Per-spin HContainer wrappers for nspin=2 (owned by this class).
    /// Rebuilt by getHR_vector() whenever hRS2 is resized.
    std::unique_ptr<HContainer<TR>> hr_spin_up_;
    std::unique_ptr<HContainer<TR>> hr_spin_dn_;

    int refresh_times = 1;

    //! current_spin for NSPIN=2 case 
    //! 0: Hamiltonian for spin up, 
    //! 1: Hamiltonian for spin down
    int current_spin = 0;

    //! snapshot of inp.nspin taken at construction; avoids PARAM dependency
    int nspin = 1;

    //! snapshot of inp.vl_in_h taken at construction; avoids PARAM dependency
    bool vl_in_h = true;

    //! cached downcast of this->ops to OperatorLCAO, filled on first use
    //! to avoid repeating dynamic_cast in updateHk/refresh
    OperatorLCAO<TK, TR>* ops_lcao_ = nullptr;

    /// get this->ops downcast to OperatorLCAO<TK, TR>*, cached in ops_lcao_
    OperatorLCAO<TK, TR>* getOperatorLCAO();
};

} // namespace hamilt

#endif
