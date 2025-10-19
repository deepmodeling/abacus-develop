#ifndef WRITE_DOS_LCAO_H
#define WRITE_DOS_LCAO_H

#include "source_base/matrix.h" // use matrix
#include "source_cell/klist.h"  // use K_Vectors
#include "source_psi/psi.h"     // use psi::Psi<T>
#include "source_hamilt/hamilt.h" // use hamilt::Hamilt<T>
#include "source_basis/module_ao/parallel_orbitals.h" // use Parallel_Orbitals
#include "source_estate/fp_energy.h" // use elecstate::Efermi


namespace ModuleIO
{
	/// @brief calculate density of states(DOS), 
    /// partial density of states(PDOS),
    ///  and mulliken charge for LCAO base
    template <typename T>
    void write_dos_lcao(
        const psi::Psi<T>* psi,      // LCAO wave functions
		hamilt::Hamilt<T>* p_ham,    // Hamiltonian
        const Parallel_Orbitals &pv, // Parallel scheme for LCAO wave functions
        const UnitCell& ucell,       // Unit cell information
		const K_Vectors& kv,         // k-point information in Brillouin zone
		const int nbands,            // Number of bands
		const elecstate::Efermi &energy_fermi,  // Fermi energy
        const ModuleBase::matrix& ekb,          // Eigenvalues per k-point and band
        const ModuleBase::matrix& wg,           // Weights of eigenvalues
        const double& dos_edelta_ev,            // Delta energy
        const double& dos_scale,                
        const double& bcoeff,
        const bool out_app_flag,
        const int istep,
        std::ofstream &ofs_running);
}
#endif
