#ifndef TD_INFO_H
#define TD_INFO_H
#include "source_base/timer.h"
#include "source_basis/module_nao/two_center_integrator.h"
#include "source_estate/module_pot/td_field_manager.h"
#include "source_hamilt/module_hcontainer/hcontainer.h"
#include "source_io/module_hs/pos_op_mat.h"
#include "source_lcao/module_ri/abfs_vector3_order.h"

#include <map>
// Class to store TDDFT infos, mainly for periodic system.
class TD_info
{
  public:
    TD_info(const UnitCell* ucell_in,
            const Parallel_Orbitals& pv,
            const LCAO_Orbitals& orb,
            const elecstate::TDFieldManager& field_manager);
    ~TD_info();

    /// @brief switch to control the output of HR
    static bool out_mat_R;

    /// @brief switch to control the output of current
    static int out_current;

    /// @brief switch to control the format of the output current, in total or in each k-point
    static bool out_current_k;

    /// @brief if need to calculate more than once
    static bool evolve_once;

    /** @brief Read the solver-owned field state without copying or publishing it. */
    const elecstate::TDFieldManager& field() const { return field_manager_; }

    /** @brief Read the current endpoint directly from the shared field manager. */
    const ModuleBase::Vector3<double>& endpoint_A() const;

    /** @brief Refresh hybrid phases using the manager's current endpoint. */
    void update_field();

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

    // For TDDFT velocity gauge, to fix the output of HR
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>> HR_sparse_td_vel[2];

    // r_calculator
    Position_op r_calculator;

  private:
    const elecstate::TDFieldManager& field_manager_;
    /// @brief lattice vectors, used to calculate the extra phase for hybrid gauge
    ModuleBase::Vector3<double> a1, a2, a3;
    double lat0;

    /// @brief store time-dependent phase for hybrid gauge
    std::map<ModuleBase::Vector3<int>, std::complex<double>> phase_hybrid;

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
