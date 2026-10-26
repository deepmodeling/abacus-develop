#include "nonlocal_xc.h"
#ifdef __LIBVDWXC
#include "nonlocal_vdwxc.h"
#endif

#include "source_base/tool_quit.h"

namespace vdw
{

std::unique_ptr<NonlocalXC> make_nonlocal_xc(const ModulePW::PW_Basis* pw, const std::string& method)
{
    if (method == "rvv10")
    {
#ifdef __LIBVDWXC
        if (pw == nullptr)
        {
            ModuleBase::WARNING_QUIT("NonlocalXC", "libvdwxc rVV10 requires a density PW basis");
        }
        return std::unique_ptr<NonlocalXC>(new LibvdwxcNonlocalXC(*pw, RVV10_B_RPW86_PBE, RVV10_C_RPW86_PBE));
#else
        ModuleBase::WARNING_QUIT("NonlocalXC", "xc_nonlocal=rvv10 requires a CMake build with ENABLE_LIBVDWXC=ON");
#endif
    }

    const std::string message = "Unsupported nonlocal XC method: " + method;
    ModuleBase::WARNING_QUIT("NonlocalXC", message);
    return nullptr;
}

} // namespace vdw
