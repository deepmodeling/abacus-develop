#include "meta_lcao.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

namespace hamilt
{

template class Meta<OperatorLCAO<double, double>>;

template class Meta<OperatorLCAO<std::complex<double>, double>>;

template class Meta<OperatorLCAO<std::complex<double>, std::complex<double>>>;

}