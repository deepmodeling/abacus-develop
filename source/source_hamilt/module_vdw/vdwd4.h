#ifndef VDWD4_H
#define VDWD4_H

#include "vdw.h"

#include <array>
#include <string>
#include <vector>

namespace vdw
{

class Vdwd4 : public Vdw
{
  public:
    Vdwd4(const UnitCell& unit_in, const std::string& xc_name, const Input_para& input);

    ~Vdwd4() override = default;

  private:
    std::string xc_name_;
    std::string model_name_;
    double cutoff_disp2_ = 0.0; // Bohr, two-body dispersion cutoff
    double cutoff_disp3_ = 0.0; // Bohr, three-body ATM cutoff
    double cutoff_cn_ = 0.0;    // Bohr, coordination-number cutoff

    bool has_force_cache_ = false;
    bool has_stress_cache_ = false;

    void set_force_from_gradient(const std::vector<double>& gradient_ha_bohr);
    void set_stress_from_sigma(const std::array<double, 9>& sigma_ha);

    void cal_energy() override;
    void cal_force() override;
    void cal_stress() override;

    void build_structure(std::vector<int>& numbers,
                         std::vector<double>& positions,
                         std::vector<double>& lattice,
                         std::array<bool, 3>& periodic) const;

    // Optional output buffers are caller-owned and non-owning here.
    // Vdwd4 writes to them during compute() and never stores their pointers.
    void compute(double& energy_ha,
                 std::vector<double>* gradient_ha_bohr,
                 std::array<double, 9>* sigma_ha);
};

} // namespace vdw

#endif // VDWD4_H
