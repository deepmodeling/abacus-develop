//==========================================================
// AUTHOR : Peize Lin
// DATE : 2014-04-25
// UPDATE : 2019-04-26
//==========================================================

#ifndef VDWD2_PARAMETERS_H
#define VDWD2_PARAMETERS_H

#include <map>
#include "module_io/input.h"
#include "module_cell/unitcell.h"
#include "vdw_parameters.h"

namespace vdw
{

class Vdwd2Parameters : public VdwParameters
{

  public:
    Vdwd2Parameters() : VdwParameters()
    {
        C6_ = C6_default_;
        R0_ = R0_default_;
    }

    ~Vdwd2Parameters() = default;

    void C6_input(const std::string &file, const std::string &unit);
    void R0_input(const std::string &file, const std::string &unit);

    void initset(const UnitCell &ucell); // init sets of vdwd2 once this correction is called
    void initial_parameters(const Input &input); // initial parameters of Vdwd2 with INPUT file

    inline const std::map<std::string, double> C6() const { return C6_; }
    inline const std::map<std::string, double> R0() const { return R0_; }
    inline double damping() const { return damping_; }
    inline double scaling() const { return scaling_; }

  private:
    double scaling_;
    double damping_;
    double radius_;
    std::map<std::string, double> C6_;
    std::map<std::string, double> R0_;
    static const std::map<std::string, double> C6_default_;
    static const std::map<std::string, double> R0_default_;
};

} // namespace vdw

#endif // VDWD2_PARAMETERS_H
