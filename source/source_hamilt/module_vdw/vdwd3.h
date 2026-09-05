#ifndef VDWD3_H
#define VDWD3_H

#include "vdwd3_types.h"
#include "vdw.h"

#include <fstream>
#include <string>

namespace vdw
{

class Vdwd3 : public Vdw
{
  public:
    Vdwd3(const UnitCell& unit_in,
          const std::string& xc_name,
          const Input_para& input,
          std::ofstream* plog = nullptr);

    ~Vdwd3() override = default;

  private:
    d3::Parameters parameters_;
    d3::Cutoffs cutoffs_;
    std::string canonical_method_;

    void evaluate_impl(const VdwRequest& request, VdwResult& result) override;
    d3::Structure build_structure() const;
    void write_parameters(std::ofstream* plog) const;
};

} // namespace vdw

#endif // VDWD3_H
