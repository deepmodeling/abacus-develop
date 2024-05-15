#ifndef TD_VELOCITY_H
#define TD_VELOCITY_H
#include <map>

#include "module_base/timer.h"
#include "module_base/abfs-vector3_order.h"
//Class to store TDDFT velocity gague infos.
class TD_Velocity
{
  public:
    TD_Velocity();
    ~TD_Velocity();

    void init();

    /// @brief Judge if in tddft calculation or not 
    static bool tddft_velocity;

    /// @brief switch to control the output of HR
    static bool out_mat_R;

    /// @brief pointer to the only TD_Velocity object itself
    static TD_Velocity* td_vel_op;

    /// @brief switch to control the output of At
    static bool out_vecpot;

    /// @brief switch to control the source of At
    static bool init_vecpot_file;

    /// @brief store the read in At_data
    static std::vector<ModuleBase::Vector3<double>> At_from_file;

    /// @brief total steps of read in At
    static int max_istep;

    //For TDDFT velocity gague, to fix the output of HR
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>> HR_sparse_td_vel[2];

    /// @brief Store the vector potential for tddft calculation
    ModuleBase::Vector3<double> cart_At;

    /// @brief calculate the At in cartesian coordinate
    void cal_cart_At(ModuleBase::Vector3<double> a0, 
                     ModuleBase::Vector3<double> a1, 
                     ModuleBase::Vector3<double> a2,
                     ModuleBase::Vector3<double> At);

    /// @brief read At from output file
    void read_cart_At(void);


  private:

    /// @brief destory HSR data stored
    void destroy_HS_R_td_sparse(void);

    /// @brief output cart_At to output file
    void output_cart_At(const std::string &out_dir);

    /// @brief store isteps now
    static int istep;

};

#endif