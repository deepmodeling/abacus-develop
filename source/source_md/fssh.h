#ifndef FSSH_MD_H_
#define FSSH_MD_H_

#include "md_base.h"
#include "source_lcao/module_operator_lcao/fssh_driver.h"
#include "source_io/module_parameter/parameter.h"

/// @brief Fewest Switches Surface Hopping (FSSH) Molecular Dynamics Class.
/// Inherits from standard MD_base, utilizing Verlet integration for nuclei
/// combined with stochastic electronic quantum jumps.
class FsshMD : public MD_base {
public:
    FsshMD(const Parameter& param_in, UnitCell& unit_in);
    virtual ~FsshMD();

    /// @brief Standard Verlet-style first half integration
    void first_half(std::ofstream& ofs) override;

    /// @brief Standard Verlet-style second half integration
    void second_half() override;

    /// @brief The core Quantum-Classical interface hooking into ESolver
    void execute_hopping(ModuleESolver::ESolver* p_esolver, const Parameter& param_in);

private:
    FsshDriver fssh_engine_;                 // The encapsulated FSSH physical engine
    ModuleBase::ComplexMatrix coef_old_;     // Cache for previous electronic wavefunctions
    bool fssh_initialized_;                  // Flag for initial setup

    // Velocity synchronization between MD_base flat arrays and UnitCell structured arrays
    void sync_vel_to_ucell();
    void sync_vel_from_ucell();

    /// @brief Hijacks the force calculation by calling the external abacus-fd tool.
    /// @details Dumps current STRU, runs abacus-fd via system call, parses excited_forces.txt,
    ///          and overwrites MD_base::force with the active state's numerical force.
    void update_force_active_state_fd(const Parameter& param_in);
};

#endif // FSSH_MD_H_