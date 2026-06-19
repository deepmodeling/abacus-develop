#ifndef TD_INFO_H
#define TD_INFO_H
#include "source_lcao/module_ri/abfs-vector3_order.h"
#include "source_base/timer.h"
#include "source_lcao/module_hcontainer/hcontainer.h"
#include "source_io/module_hs/cal_r_overlap_R.h"
#include "source_basis/module_nao/two_center_integrator.h"

#include <map>
// Class to store TDDFT infos, mainly for periodic system.
class TD_info
{
  public:
    TD_info(const UnitCell* ucell_in,const Parallel_Orbitals& pv, const LCAO_Orbitals& orb);
    ~TD_info();

    /// @brief switch to control the output of HR
    static bool out_mat_R;

    /// @brief pointer to the only TD_info object itself
    static TD_info* td_vel_op;

    /// @brief switch to control the output of At
    static bool out_vecpot;

    /// @brief switch to control the output of current
    static int out_current;

    /// @brief switch to control the format of the output current, in total or in each k-point
    static bool out_current_k;

    /// @brief switch to control the source of At
    static bool init_vecpot_file;

    /// @brief if need to calculate more than once
    static bool evolve_once;

    /// @brief Restart step
    static int estep_shift;

    /// @brief Store the vector potential for tddft calculation
    static ModuleBase::Vector3<double> cart_At;

    /// @brief calculate the At in cartesian coordinate
    void cal_cart_At(const ModuleBase::Vector3<double>& At);

    /// @brief output RT-TDDFT info for restart
    void out_restart_info(const int nstep, 
                          const ModuleBase::Vector3<double>& At_current, 
                          const ModuleBase::Vector3<double>& At_laststep);

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
    cal_r_overlap_R r_calculator;

  private:
    /// @brief lattice vectors, used to calculate the extra phase for hybrid gauge
    ModuleBase::Vector3<double>a1, a2, a3;
    double lat0;

    /// @brief store time-dependent phase for hybrid gauge
    std::map<ModuleBase::Vector3<int>, std::complex<double>> phase_hybrid;

    /// @brief read At from output file
    void read_cart_At();

    /// @brief output cart_At to output file
    void output_cart_At(const std::string& out_dir);

    /// @brief store isteps now
    static int istep;

    /// @brief total steps of read in At
    static int max_istep;

    /// @brief store the read in At_data
    static std::vector<ModuleBase::Vector3<double>> At_from_file;

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
