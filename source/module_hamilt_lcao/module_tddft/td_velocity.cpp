#include "module_base/timer.h"
#include "td_velocity.h"


bool TD_Velocity::tddft_velocity = false;
bool TD_Velocity::out_mat_R = false;
TD_Velocity* TD_Velocity::td_vel_op = nullptr;

TD_Velocity::TD_Velocity()
{
  return;
}
TD_Velocity::~TD_Velocity()
{
  this->destroy_HS_R_td_sparse();
}

void TD_Velocity::destroy_HS_R_td_sparse(void)
{
  std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>> empty_HR_sparse_td_vel_up;
  std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>> empty_HR_sparse_td_vel_down;
  HR_sparse_td_vel[0].swap(empty_HR_sparse_td_vel_up);
  HR_sparse_td_vel[1].swap(empty_HR_sparse_td_vel_down);
}