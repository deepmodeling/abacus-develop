#include <complex>

#include "source_pw/module_pwdft/onsite_proj.h"
#include "source_pw/module_pwdft/dftu_base.h"

template <typename T, typename Device>
void projectors::OnsiteProjector<T, Device>::cal_force_onsite_dftu(int ik, int npm, T* force,
                                                        const Plus_U_Base& dftu, int nks,
                                                        const double* wg_ik) const
{
    const int isk_val = this->isk_ ? this->isk_[ik] : 0;
    const std::complex<double>* pot_onsite_ptr = dftu.get_pot_uterm_pw_spin(isk_val);
    const int pot_onsite_size = dftu.get_size_pot_uterm_pw_spin();
    this->fs_tools->cal_force_dftu(ik, npm, force,
        dftu.get_orbital_corr_vec().data(), pot_onsite_ptr, pot_onsite_size, wg_ik);
}

template <typename T, typename Device>
double projectors::OnsiteProjector<T, Device>::cal_stress_onsite_dftu(int ik, int npm,
                                                           const Plus_U_Base& dftu, int nks,
                                                           const double* wg_ik) const
{
    const int isk_val = this->isk_ ? this->isk_[ik] : 0;
    const std::complex<double>* pot_onsite_ptr = dftu.get_pot_uterm_pw_spin(isk_val);
    const int pot_onsite_size = dftu.get_size_pot_uterm_pw_spin();
    return this->fs_tools->cal_stress_dftu(ik, npm,
        dftu.get_orbital_corr_vec().data(), pot_onsite_ptr, pot_onsite_size, wg_ik);
}

template <typename T, typename Device>
void projectors::OnsiteProjector<T, Device>::cal_force_onsite_dspin(int ik, int npm, T* force,
                                                         const ModuleBase::Vector3<double>* lambda,
                                                         const double* wg_ik) const
{
    this->fs_tools->cal_force_dspin(ik, npm, force, lambda, wg_ik);
}

template <typename T, typename Device>
double projectors::OnsiteProjector<T, Device>::cal_stress_onsite_dspin(int ik, int npm,
                                                            const ModuleBase::Vector3<double>* lambda,
                                                            const double* wg_ik) const
{
    return this->fs_tools->cal_stress_dspin(ik, npm, lambda, wg_ik);
}

// explicit method instantiation
template
void projectors::OnsiteProjector<double, base_device::DEVICE_CPU>::cal_force_onsite_dftu(
    int, int, double*, const Plus_U_Base&, int, const double*) const;

template
double projectors::OnsiteProjector<double, base_device::DEVICE_CPU>::cal_stress_onsite_dftu(
    int, int, const Plus_U_Base&, int, const double*) const;

template
void projectors::OnsiteProjector<double, base_device::DEVICE_CPU>::cal_force_onsite_dspin(
    int, int, double*, const ModuleBase::Vector3<double>*, const double*) const;

template
double projectors::OnsiteProjector<double, base_device::DEVICE_CPU>::cal_stress_onsite_dspin(
    int, int, const ModuleBase::Vector3<double>*, const double*) const;

#if ((defined __CUDA) || (defined __ROCM))
template
void projectors::OnsiteProjector<double, base_device::DEVICE_GPU>::cal_force_onsite_dftu(
    int, int, double*, const Plus_U_Base&, int, const double*) const;

template
double projectors::OnsiteProjector<double, base_device::DEVICE_GPU>::cal_stress_onsite_dftu(
    int, int, const Plus_U_Base&, int, const double*) const;

template
void projectors::OnsiteProjector<double, base_device::DEVICE_GPU>::cal_force_onsite_dspin(
    int, int, double*, const ModuleBase::Vector3<double>*, const double*) const;

template
double projectors::OnsiteProjector<double, base_device::DEVICE_GPU>::cal_stress_onsite_dspin(
    int, int, const ModuleBase::Vector3<double>*, const double*) const;
#endif
