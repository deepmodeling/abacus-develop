#ifndef VDW_H
#define VDW_H

#include <memory>
#include <vector>
#include "module_cell/unitcell.h"
#include "vdw_parameters.h"
#include "vdwd2_parameters.h"
#include "vdwd3_parameters.h"

namespace vdw
{

template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args &&... args) {
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

class Vdw
{
  public:
    Vdw(const UnitCell &unit_in) : ucell_(unit_in) {};

    virtual ~Vdw() = default;

    inline double get_energy(bool cal=true) {
        if (cal) { cal_energy(); }
        return energy_;
    }
    inline const std::vector<ModuleBase::Vector3<double>> &get_force(bool cal=true) {
        if (cal) { cal_force(); }
        return force_;
    }
    inline const ModuleBase::Matrix3 &get_stress(bool cal=true) {
        if (cal) { cal_stress(); }
        return stress_;
    }

  protected:
    const UnitCell &ucell_;

    double energy_ = 0;
    std::vector<ModuleBase::Vector3<double>> force_;
    ModuleBase::Matrix3 stress_;

    virtual void cal_energy() = 0;
    virtual void cal_force() = 0;
    virtual void cal_stress() = 0;
};

/**
 * @brief make vdw correction object
 * 
 * @param ucell UnitCell instance
 * @param input Parameter instance
 * @param plog optional, for logging the parameter setting process
 * @return std::unique_ptr<Vdw> 
 */
std::unique_ptr<Vdw> make_vdw(const UnitCell &ucell, 
                              const Input_para &input,
                              std::ofstream* plog = nullptr);

} // namespace vdw

#endif // VDW_H
