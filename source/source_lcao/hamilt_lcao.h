#ifndef HAMILT_LCAO_H 
#define HAMILT_LCAO_H 

#include "source_basis/module_nao/two_center_bundle.h"
#include "source_cell/klist.h"
#include "source_cell/module_neighbor/sltk_atom_arrange.h"
#include "source_estate/module_dm/density_matrix.h"
#include "source_estate/module_pot/potential_new.h"
#include "source_hamilt/hamilt.h"
#include "source_lcao/hs_matrix_k.hpp"
#include "source_lcao/module_gint/gint_gamma.h"
#include "source_lcao/module_gint/gint_k.h"
#include "source_lcao/module_hcontainer/hcontainer.h"

#include <vector>

#ifdef __MLALGO
#include "source_lcao/module_deepks/LCAO_deepks.h"
#endif

#ifdef __EXX
#include "source_lcao/module_ri/Exx_LRI.h"
#endif
namespace hamilt
{

// template first for type of k space H matrix elements
// template second for type of temporary matrix, 
// gamma_only fix-gamma-matrix + S-gamma, 
// multi-k fix-Real + S-Real
template <typename TK, typename TR>
class HamiltLCAO : public Hamilt<TK>
{
  public:


    using TAC = std::pair<int, std::array<int, 3>>;


    /**
     * @brief Constructor of Hamiltonian for LCAO base
     * HR and SR will be allocated with Operators
     */
    HamiltLCAO(Gint_Gamma* GG_in,
               Gint_k* GK_in,
               const UnitCell& ucell,
               const Grid_Driver& grid_d,
			   const Parallel_Orbitals* paraV,
			   elecstate::Potential* pot_in,
			   const K_Vectors& kv_in,
			   const TwoCenterBundle& two_center_bundle,
               const LCAO_Orbitals& orb,
               elecstate::DensityMatrix<TK, double>* DM_in
#ifdef __MLALGO
               ,
               LCAO_Deepks<TK>* ld_in
#endif
#ifdef __EXX
               ,
               const int istep,
               int* exx_two_level_step = nullptr,
               std::vector<std::map<int, std::map<TAC, RI::Tensor<double>>>>* Hexxd = nullptr,
               std::vector<std::map<int, std::map<TAC, RI::Tensor<std::complex<double>>>>>* Hexxc = nullptr
#endif
    );

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
        if (this->ops != nullptr)
        {
            delete this->ops;
        }
        delete this->hR;
        delete this->sR;
        delete this->hsk;
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

    int get_size_hsk() const
    {
        return this->hsk->get_size();
    }

    /// get HR pointer of *this->hR, which is a HContainer<TR> and contains H(R)
    HContainer<TR>*& getHR()
    {
        return this->hR;
    }

    /// get SR pointer of *this->sR, which is a HContainer<TR> and contains S(R)
    HContainer<TR>*& getSR()
    {
        return this->sR;
    }

#ifdef __MLALGO
    /// get V_delta_R pointer of *this->V_delta_R, which is a HContainer<TR> and contains V_delta(R)
    HContainer<TR>*& get_V_delta_R()
    {
        return this->V_delta_R;
    }
#endif

    /// refresh the status of HR
    void refresh() override;

    // for target K point, update consequence of hPsi() and matrix()
    virtual void updateHk(const int ik) override;

    /**
     * @brief special for LCAO, update SK only
     *
     * @param ik target K point
     * @param kvec_d: direct coordinates of k-points
     * @param hk_type 0: SK is row-major, 1: SK is collumn-major
     * @return void
     */
	void updateSk(const int ik, const int hk_type = 0);

    // core function: return H(k) and S(k) matrixs for direct solving eigenvalues.
    // not used in PW base
    void matrix(MatrixBlock<TK>& hk_in, MatrixBlock<TK>& sk_in) override;

  private:

    const K_Vectors* kv = nullptr;

    //! Real space Hamiltonian H(R), where R is the Bravis lattice vector
    HContainer<TR>* hR = nullptr;

    //! Real space overlap matrix S(R), where R is the Bravis lattice vector
    HContainer<TR>* sR = nullptr;

#ifdef __MLALGO
    HContainer<TR>* V_delta_R = nullptr;
#endif

    //! Hamiltonian and overlap matrices for a specific k point
    HS_Matrix_K<TK>* hsk = nullptr;

    // special case for NSPIN=2 , data of HR should be separated into two parts
    // save them in this->hRS2;
    std::vector<TR> hRS2;

    int refresh_times = 1;

    //! current_spin for NSPIN=2 case 
    //! 0: Hamiltonian for spin up, 
    //! 1: Hamiltonian for spin down
    int current_spin = 0;

    const int istep = 0;
};

} // namespace hamilt

#endif
