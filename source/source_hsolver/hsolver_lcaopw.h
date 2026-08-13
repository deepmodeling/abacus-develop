#ifndef HSOLVERLIP_H
#define HSOLVERLIP_H

#include "source_estate/elecstate.h"
#include "source_hamilt/hamilt.h"
#include "source_base/macros.h"

/// General_Exx_Info forward declaration, full definition in general_exx_info.h
struct General_Exx_Info;

namespace hsolver
{

// LCAO-in-PW does not support GPU now.
template <typename T>
class HSolverLIP
{
  private:
    // Note GetTypeReal<T>::type will
    // return T if T is real type(float, double),
    // otherwise return the real type of T(complex<float>, std::complex<double>)
    using Real = typename GetTypeReal<T>::type;

  public:
    HSolverLIP(ModulePW::PW_Basis_K* wfc_basis_in,
               const bool use_uspp_in,
               const std::string basis_type_in,
               const std::string calculation_in)
        : wfc_basis(wfc_basis_in), use_uspp(use_uspp_in), basis_type(basis_type_in), calculation(calculation_in) {};

    /// @brief solve function for lcao_in_pw
    /// @param pHamilt interface to hamilt
    /// @param psi reference to psi
    /// @param pes interface to elecstate
    /// @param transform transformation matrix between lcao and pw
    /// @param skip_charge
    void solve(hamilt::Hamilt<T>* pHamilt,
               psi::Psi<T>& psi,
               elecstate::ElecState* pes,
               psi::Psi<T>& transform,
               const bool skip_charge,
               const double tpiba,
               const int nat,
               const General_Exx_Info& exx_info);

  private:
    ModulePW::PW_Basis_K* wfc_basis = nullptr;

    const bool use_uspp; // true if ultrasoft pseudopotentials are in use

    const std::string basis_type;  // "lcao_in_pw" for this solver
    const std::string calculation; // "scf", "nscf", "md", "relax", ...
};

} // namespace hsolver

#endif