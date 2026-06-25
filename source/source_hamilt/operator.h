#ifndef OPERATOR_H
#define OPERATOR_H

#include "source_psi/psi.h"

#include <complex>

namespace hamilt
{

enum class calculation_type
{
    no,
    pw_ekinetic,
    pw_nonlocal,
    pw_veff,
    pw_meta,
    pw_onsite,
    pw_exx,
    lcao_overlap,
    lcao_fixed,
    lcao_gint,
    lcao_deepks,
    lcao_exx,
    lcao_dftu,
    lcao_sc_lambda,
    lcao_tddft_periodic,
};

// Basic class for operator module,
// it is designed for "O|psi>" and "<psi|O|psi>"
// Operator "O" might have several different types, which should be calculated one by one.
// In basic class , function add() is designed for combine all operators together with a chain.
template <typename T, typename Device = base_device::DEVICE_CPU>
class Operator
{
  public:
    Operator();
    virtual ~Operator();

    /// @brief input type of hPsi: (psi_input, range, hpsi_pointer)
    /// @details
    /// hpsi_info bundles the input and output of hPsi():
    /// - std::get<0>: pointer to psi::Psi<T, Device> (input wavefunction)
    /// - std::get<1>: psi::Range specifying which bands to operate on
    /// - std::get<2>: T* pointer to output hpsi buffer (can equal psi_in for in-place)
    typedef std::tuple<const psi::Psi<T, Device>*, const psi::Range, T*> hpsi_info;

    /// @brief Core hot-path: compute H|psi> by traversing the operator chain.
    /// @details
    /// This is the central computational kernel of ABACUS, called O(n_bands * n_iter * n_kpoints)
    /// times during a single SCF calculation. It accounts for 18-25% of total runtime.
    ///
    /// Algorithm:
    /// 1. Unwrap the input hpsi_info to get psi, band range, and output buffer
    /// 2. Allocate or reuse the temporary hpsi buffer via get_hpsi()
    /// 3. Call act() on the first operator node (is_first_node=true) -- this node zeros hpsi
    /// 4. Iterate through next_op linked list, calling act() on each subsequent node (is_first_node=false)
    ///    Each node accumulates its contribution: hpsi += O|psi>
    /// 5. If in_place mode, copy temporary hpsi back to the caller-provided buffer
    /// 6. Return wrapped hpsi_info for downstream use
    ///
    /// The operator chain typically includes (in order):
    /// Ekinetic → Veff → Nonlocal → Meta → OnsiteProj
    ///
    /// @param input hpsi_info tuple: (psi_input, band_range, hpsi_output_pointer)
    ///   - psi_input: the wavefunction Psi object
    ///   - band_range: which bands to compute (range_1 to range_2 inclusive)
    ///   - hpsi_output_pointer: pre-allocated buffer for H|psi> result.
    ///     If equal to psi_input pointer, in_place mode is used (temporary buffer allocated internally).
    /// @return hpsi_info containing (internal_hpsi, range, caller_hpsi_pointer)
    /// @note This function is performance-critical. The operator chain traversal is the innermost
    ///       loop of iterative diagonalization methods (CG, Davidson, BPCG).
    /// @note For PW calculations, each act() call may involve FFT transforms (Veff),
    ///       BLAS3 gemm operations (Nonlocal), or element-wise vector ops (Ekinetic).
    /// @see Operator::act(), HamiltPW, DiagoCG::diag()
    virtual hpsi_info hPsi(hpsi_info& input) const;

    virtual void init(const int ik_in);

    virtual void add(Operator* next);

    virtual int get_ik() const
    {
        return this->ik;
    }

    /// do operation : |hpsi_choosed> = V|psi_choosed>
    /// V is the target operator act on choosed psi, the consequence should be added to choosed hpsi
    ///  interface type 1: pointer-only (default)
    ///  @note PW: nbasis = max_npw * npol, nbands = nband * npol, npol = npol. Strange but PAY ATTENTION!!!
    virtual void act(const int nbands,
                     const int nbasis,
                     const int npol,
                     const T* tmpsi_in,
                     T* tmhpsi,
                     const int ngk_ik = 0,
                     const bool is_first_node = false) const {};

    /// developer-friendly interfaces for act() function
    /// interface type 2: input and change the Psi-type HPsi
	// virtual void act(const psi::Psi<T, Device>& psi_in, psi::Psi<T, Device>& psi_out) const {};
	virtual void act(const psi::Psi<T, Device>& psi_in, 
			psi::Psi<T, Device>& psi_out, 
			const int nbands) const {};

    /// interface type 3: return a Psi-type HPsi
    // virtual psi::Psi<T> act(const psi::Psi<T,Device>& psi_in) const { return psi_in; };

    Operator* next_op = nullptr;

    /// type 1 (default): pointer-only
    ///         act(const T* psi_in, T* psi_out)
    /// type 2: use the `Psi`class
    ///         act(const Psi& psi_in, Psi& psi_out)
    int get_act_type() const
    {
        return this->act_type;
    }

    calculation_type get_cal_type() const
    {
        return this->cal_type;
    }

  protected:
    int ik = 0;
    int act_type = 1; ///< determine which act() interface would be called in hPsi()

    mutable bool in_place = false;

    // calculation type, only different type can be in main chain table
    enum calculation_type cal_type;
    Operator* next_sub_op = nullptr;
    bool is_first_node = true;

    // if this Operator is first node in chain table, hpsi would not be empty
    mutable psi::Psi<T, Device>* hpsi = nullptr;

    /*This function would analyze hpsi_info and choose how to arrange hpsi storage
    In hpsi_info, if the third parameter hpsi_pointer is set, which indicates memory of hpsi is arranged by developer;
    if hpsi_pointer is not set(nullptr), which indicates memory of hpsi is arranged by Operator, this case is rare.
    two cases would occurred:
    1. hpsi_pointer != nullptr && psi_pointer == hpsi_pointer , psi would be replaced by hpsi, hpsi need a temporary
    memory
    2. hpsi_pointer != nullptr && psi_pointer != hpsi_pointer , this is the commonly case
    */
    T* get_hpsi(const hpsi_info& info) const;

    Device* ctx = {};
    using set_memory_op = base_device::memory::set_memory_op<T, Device>;
};

} // end namespace hamilt

#endif
