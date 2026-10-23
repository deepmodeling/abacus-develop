#ifndef TD_TEST_UTILS_H
#define TD_TEST_UTILS_H

#include "source_pw/module_rt/td_spin_pw.h"
#include "source_pw/module_rt/td_vel_pw.h"

#include <gtest/gtest.h>

/** @brief Restore the global input used by legacy PW initialization. */
class TestParameters
{
  public:
    TestParameters();
    ~TestParameters();

  private:
    Input_para input_;
    System_para system_;
};

namespace td_test
{
using Complex = std::complex<double>;
using Device = base_device::DEVICE_CPU;
using Vector = std::vector<Complex>;

/** @brief Small real PW basis and relativistic projectors, without an SCF calculation. */
class Observables : public testing::Test
{
  protected:
    TestParameters parameters_;
    UnitCell cell_;
    ModulePW::PW_Basis_K basis_;
    K_Vectors kv_;
    Structure_Factor sf_;
    pseudopot_cell_vnl pp_;
    Charge charge_;
    elecstate::ElecState state_;
    std::unique_ptr<psi::Psi<Complex>> psi_;

    void SetUp() override;
    void TearDown() override;
    Vector projectors(const ModuleBase::Vector3<double>& a);
    Vector nonlocal(const Vector& beta, const Vector& x, int atom) const;
    Vector hamiltonian(const ModuleBase::Vector3<double>& a, const Vector& x);
};

/** @brief Apply a Pauli matrix to a packed two-component vector. */
Vector pauli(const Vector& x, int axis);
Complex inner(const Vector& x, const Vector& y);
} // namespace td_test
#endif
