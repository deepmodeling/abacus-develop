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

    /// @brief Judge if in tddft calculation or not 
    static bool tddft_velocity;

    /// @brief switch to control the output of HR
    static bool out_mat_R;

    /// @brief pointer to the only TD_Velocity object itself
    static TD_Velocity* td_vel_op;

    //For TDDFT velocity gague, to fix the output of HR
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>> HR_sparse_td_vel[2];


  private:

    /// @brief destory HSR data stored
    void destroy_HS_R_td_sparse(void);

};

#endif