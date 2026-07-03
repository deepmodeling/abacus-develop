#ifndef STERNHEIMER_ABACUS_DFT_ZERO_ORDER_H
#define STERNHEIMER_ABACUS_DFT_ZERO_ORDER_H

#include "source_lcao/module_ri/sternheimer_abacus_fd_adapter.h"
#include "source_lcao/module_ri/sternheimer_dft_zero_order.h"

#include <string>
#include <vector>

class UnitCell;

namespace elecstate
{
class ElecState;
class Potential;
}

namespace ModulePW
{
class PW_Basis;
}

namespace ModuleRI
{

struct SternheimerABACUSDFTZeroOrderResult
{
    SternheimerABACUSFDGridData grid_data;
    SternheimerFDZeroOrderStates fd_states;
    std::vector<double> dft_eigenvalues;
    std::vector<double> dft_occupations;
    SternheimerDFTZeroOrderComparison comparison;
    std::string hamiltonian_mode;
};

SternheimerABACUSDFTZeroOrderResult compare_sternheimer_abacus_fd_zero_order_to_dft(
    const elecstate::Potential& potential,
    const ModulePW::PW_Basis& pw_basis,
    const elecstate::ElecState& elec_state,
    int spin,
    int k_index,
    int num_bands,
    double eigenvalue_tolerance,
    int max_dense_size = 4096,
    double kinetic_prefactor = 1.0);

SternheimerABACUSDFTZeroOrderResult compare_sternheimer_abacus_fd_zero_order_to_dft(
    const elecstate::Potential& potential,
    const ModulePW::PW_Basis& pw_basis,
    const UnitCell& ucell,
    const elecstate::ElecState& elec_state,
    int spin,
    int k_index,
    int num_bands,
    double eigenvalue_tolerance,
    int max_dense_size = 4096,
    double kinetic_prefactor = 1.0);

} // namespace ModuleRI

#endif
