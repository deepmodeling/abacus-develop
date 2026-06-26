#ifndef TD_MOVING_GAUGE_H
#define TD_MOVING_GAUGE_H

#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_basis/module_nao/two_center_integrator.h"
#include "source_cell/unitcell.h"
#include "source_lcao/module_hcontainer/hcontainer.h"
#include "source_lcao/module_hcontainer/hcontainer_funcs.h"

#include <complex>
#include <vector>

namespace module_rt
{

class TD_MovingGauge
{
  public:
    TD_MovingGauge() = default;
    ~TD_MovingGauge();

    // Initialize the R-space derivative matrices D_R (x, y, z)
    // using the provided sR_template for consistent sparse atomic pair topology
    // D_{K,\mu\nu}(R) = <\phi_{\mu 0}|∂\phi_{\nu R}/∂\tau_K> where tau_K is the position of atom K
    template <typename T_sR>
    void init_DR(const hamilt::HContainer<T_sR>* sR_template,
                 const UnitCell* ucell,
                 const Parallel_Orbitals* paraV,
                 TwoCenterIntegrator* intor);

    // Update the R-space matrix D_R (x, y, z)
    template <typename T_sR>
    void update_DR(const hamilt::HContainer<T_sR>* sR_template,
                   const UnitCell* ucell,
                   const Parallel_Orbitals* paraV,
                   TwoCenterIntegrator* intor);

    // Fourier transform D(R) to D(k)
    // Note: folding_HR performs an accumulation (+=) operation, need to ensure Dk matrices are zeroed before calling
    // D_{K,\mu\nu}(k) = \sum_R e^{ikR} D_{K,\mu\nu}(R)
    template <typename TK>
    void get_D_k(int K, const ModuleBase::Vector3<double>& kvec_d, TK* Dk_x, TK* Dk_y, TK* Dk_z, int hk_ld) const;

    // Calculate the moving spatial gauge matrix P_k and accumulate it to the input P_k matrix
    // Note: The unit is converted to Rydberg atomic units, and multiplied by 2 internally
    // P_{\mu\nu}(k) = -i \sum_K vel_K \cdot D_{K,\mu\nu}(k) where vel_K is the velocity of atom K
    template <typename TK>
    void get_P_k(const UnitCell* ucell, const ModuleBase::Vector3<double>& kvec_d, TK* P_k, int matrix_size, int hk_ld)
        const;

  private:
    int nat_ = 0;

    std::vector<hamilt::HContainer<double>*> DR_x_;
    std::vector<hamilt::HContainer<double>*> DR_y_;
    std::vector<hamilt::HContainer<double>*> DR_z_;
};

} // namespace module_rt

#endif // TD_MOVING_GAUGE_H
