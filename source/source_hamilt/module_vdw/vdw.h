#ifndef VDW_H
#define VDW_H

#include <memory>
#include <utility>
#include <vector>

#include "source_cell/unitcell.h"
#include "vdw_parameters.h"
#include "vdwd2_parameters.h"
#include "vdwd3_parameters.h"

namespace vdw
{

template <typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args)
{
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

struct VdwRequest
{
    VdwRequest(const bool force_in, const bool stress_in) : force(force_in), stress(stress_in) {}

    bool force;
    bool stress;
};

struct VdwResult
{
    VdwResult() : energy(0.0), has_force(false), has_stress(false)
    {
        stress.Zero();
    }

    double energy;
    std::vector<ModuleBase::Vector3<double>> force;
    ModuleBase::Matrix3 stress;
    bool has_force;
    bool has_stress;
};

class Vdw
{
  public:
    Vdw(const UnitCell& unit_in) : ucell_(unit_in) {}

    virtual ~Vdw() = default;

    VdwResult evaluate(const VdwRequest& request)
    {
        VdwResult result;
        evaluate_impl(request, result);
        return result;
    }

  protected:
    const UnitCell& ucell_;

    virtual void evaluate_impl(const VdwRequest& request, VdwResult& result) = 0;
};

/**
 * @brief make vdw correction object
 *
 * @param ucell UnitCell instance
 * @param input Parameter instance
 * @param plog optional, for logging the parameter setting process
 * @return std::unique_ptr<Vdw>
 */
std::unique_ptr<Vdw> make_vdw(const UnitCell& ucell,
                              const Input_para& input,
                              std::ofstream* plog = nullptr);

} // namespace vdw

#endif // VDW_H
