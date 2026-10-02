#ifndef TD_INFO_H
#define TD_INFO_H
#include "source_base/timer.h"
#include "source_basis/module_nao/two_center_integrator.h"
#include "source_hamilt/module_hcontainer/hcontainer.h"
#include "source_io/module_hs/pos_op_mat.h"
#include "source_lcao/module_ri/abfs_vector3_order.h"

#include <map>
// Class to store TDDFT infos, mainly for periodic system.
class TD_info
{
  public:
    TD_info(const UnitCell* ucell_in, const Parallel_Orbitals& pv, const LCAO_Orbitals& orb, const int restart_step);
    ~TD_info();

    /// @brief switch to control the output of HR
    static bool out_mat_R;

    /// @brief pointer to the only TD_info object itself
    static TD_info* td_vel_op;

    /// @brief switch to control the output of current
    static int out_current;

    /// @brief switch to control the format of the output current, in total or in each k-point
    static bool out_current_k;

    /// @brief if need to calculate more than once
    static bool evolve_once;

    /// @brief Restart step
    static int estep_shift;

    /** @brief Propagation vector potential in Hartree atomic units, published by the ESolver. */
    static ModuleBase::Vector3<double> A_prop_ha;

    /** @brief Bind the manager-selected propagation A and refresh hybrid phases at an explicit step. */
    void set_A_prop(const int step, const ModuleBase::Vector3<double>& A_ha);

    // allocate memory for current term.
    void initialize_current_term(const hamilt::HContainer<std::complex<double>>* HR, const Parallel_Orbitals* paraV);

    hamilt::HContainer<std::complex<double>>* get_current_term_pointer(const int& i) const
    {
        return this->current_term[i];
    }
    // allocate memory for phase_hybrid.
    template <typename TR>
    void initialize_phase_hybrid(const UnitCell& ucell, const hamilt::HContainer<TR>* hR);

    const std::map<ModuleBase::Vector3<int>, std::complex<double>>& get_phase_hybrid() const
    {
        return this->phase_hybrid;
    }

    void calculate_grad_overlap(const Parallel_Orbitals& paraV,
                                const UnitCell& ucell,
                                const Grid_Driver& GridD,
                                const std::vector<double>& orb_cutoff,
                                const TwoCenterIntegrator* intor);
    std::vector<hamilt::HContainer<double>*> get_grad_overlap() const
    {
        return this->grad_overlap;
    }
    // set velocity HR.
    void set_velocity_HR(hamilt::HContainer<std::complex<double>>* HR)
    {
        this->velocity_HR = HR;
    }
    hamilt::HContainer<std::complex<double>>* get_velocity_HR_pointer() const
    {
        return this->velocity_HR;
    }

    int get_istep()
    {
        return istep;
    }
    // For TDDFT velocity gauge, to fix the output of HR
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>> HR_sparse_td_vel[2];

    //r_calculator
    Position_op r_calculator;

  private:
    /// @brief lattice vectors, used to calculate the extra phase for hybrid gauge
    ModuleBase::Vector3<double> a1, a2, a3;
    double lat0;

    /// @brief store time-dependent phase for hybrid gauge
    std::map<ModuleBase::Vector3<int>, std::complex<double>> phase_hybrid;

    /// @brief store isteps now
    static int istep;

    /// @brief store the dS/dD matrix
    std::vector<hamilt::HContainer<double>*> grad_overlap = {nullptr, nullptr, nullptr};

    /// @brief destory HSR data stored
    void destroy_HS_R_td_sparse();

    /// @brief part of Momentum operator, -i∇ - i[r,Vnl]. Used to calculate current.
    std::vector<hamilt::HContainer<std::complex<double>>*> current_term = {nullptr, nullptr, nullptr};

    /// @brief store kinetic hamilton
    hamilt::HContainer<std::complex<double>>* velocity_HR = nullptr;
};

#endif
