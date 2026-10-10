#include "xc_kxc.h"

#ifdef __LIBXC
#include <xc.h>
#include <stdexcept>
#include <string>

void LR::require_xc_kxc(const xc_func_type& functional, const bool need_kxc)
{
    if (!need_kxc) { return; }
    // Only mixed functionals evaluate these children as separate XC contributions.
    // Other auxiliary objects may merely provide internal parameters or lower derivatives.
    if (functional.mix_coef != nullptr)
    {
        for (int component = 0; component < functional.n_func_aux; ++component)
        {
            require_xc_kxc(*functional.func_aux[component], need_kxc);
        }
    }
    const xc_func_info_type* info = xc_func_get_info(&functional);
    const int flags = xc_func_info_get_flags(info);
    if ((flags & XC_FLAGS_HAVE_KXC) == 0)
    {
        const std::string name = xc_func_info_get_name(info);
        const std::string version = xc_version_string();
        const std::string message = "LR-TDDFT gradients require LibXC KXC, but functional '" + name
            + "' in the loaded LibXC " + version + " does not provide it. "
            "If KXC was disabled at build time, rebuild LibXC with -DDISABLE_KXC=OFF. "
            "If the functional itself lacks third derivatives, use a supported functional or LibXC version.";
        throw std::runtime_error(message);
    }
}
#endif
