#include "source_lcao/module_ri/sternheimer_abacus_dft_zero_order.h"

#include "source_lcao/module_ri/sternheimer_abacus_dft_reference.h"

namespace ModuleRI
{

SternheimerABACUSDFTZeroOrderResult compare_sternheimer_abacus_fd_zero_order_to_dft(
    const elecstate::Potential& potential,
    const ModulePW::PW_Basis& pw_basis,
    const elecstate::ElecState& elec_state,
    const int spin,
    const int k_index,
    const int num_bands,
    const double eigenvalue_tolerance,
    const int max_dense_size,
    const double kinetic_prefactor)
{
    SternheimerABACUSDFTZeroOrderResult result;
    result.grid_data = make_sternheimer_fd_grid(pw_basis);
    const std::vector<double> local_potential = copy_sternheimer_local_potential(potential, pw_basis, spin);
    const SternheimerFDHamiltonian hamiltonian(result.grid_data.grid, local_potential, kinetic_prefactor);

    result.fd_states
        = solve_sternheimer_fd_zero_order_dense(hamiltonian, num_bands, result.grid_data.volume_element, max_dense_size);
    result.dft_eigenvalues = copy_sternheimer_dft_eigenvalues(elec_state, k_index, num_bands);
    result.dft_occupations = copy_sternheimer_dft_occupations(elec_state, k_index, num_bands);
    result.comparison
        = compare_sternheimer_fd_zero_order_to_dft(result.fd_states, result.dft_eigenvalues, eigenvalue_tolerance);

    return result;
}

SternheimerABACUSDFTZeroOrderResult compare_sternheimer_abacus_fd_zero_order_to_dft(
    const elecstate::Potential& potential,
    const ModulePW::PW_Basis& pw_basis,
    const UnitCell& ucell,
    const elecstate::ElecState& elec_state,
    const int spin,
    const int k_index,
    const int num_bands,
    const double eigenvalue_tolerance,
    const int max_dense_size,
    const double kinetic_prefactor)
{
    SternheimerABACUSDFTZeroOrderResult result;
    result.grid_data = make_sternheimer_fd_grid(pw_basis);
    const SternheimerFDHamiltonian hamiltonian
        = make_sternheimer_fd_hamiltonian(potential, pw_basis, ucell, spin, kinetic_prefactor);

    result.fd_states
        = solve_sternheimer_fd_zero_order_dense(hamiltonian, num_bands, result.grid_data.volume_element, max_dense_size);
    result.dft_eigenvalues = copy_sternheimer_dft_eigenvalues(elec_state, k_index, num_bands);
    result.dft_occupations = copy_sternheimer_dft_occupations(elec_state, k_index, num_bands);
    result.comparison
        = compare_sternheimer_fd_zero_order_to_dft(result.fd_states, result.dft_eigenvalues, eigenvalue_tolerance);

    return result;
}

} // namespace ModuleRI
