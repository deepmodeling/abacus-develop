#ifndef VSEP_IN_PW
#define VSEP_IN_PW

#include "source_base/matrix.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/sep_cell.h"

#include <vector>

class VSep
{
  public:
    VSep() noexcept;
    ~VSep() noexcept;

    void init_vsep(const ModulePW::PW_Basis& rho_basis, const Sep_Cell& sep_cell);
    void generate_vsep_r(const ModulePW::PW_Basis& rho_basis, const ModuleBase::ComplexMatrix& sf_in, const Sep_Cell& sep_cell);

    ModuleBase::matrix vsep_form;
    std::vector<double> vsep_r;

  private:
    int nrxx = 0;
};
//
// namespace GlobalC
// {
// extern VSep vsep_cell;
// }

#endif /* ifndef VSEP_IN_PW */
