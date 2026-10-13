#include "nonlocal_xc.h"
#ifdef __LIBVDWXC
#include "nonlocal_xc_libvdwxc.h"
#endif

#include "source_base/tool_quit.h"

namespace vdw
{

std::unique_ptr<NonlocalXC> make_nonlocal_xc(const ModulePW::PW_Basis* pw,
                                             const std::string& method,
                                             double b,
                                             double c)
{
    if (method == "rvv10")
    {
#ifdef __LIBVDWXC
        if (pw == nullptr)
        {
            ModuleBase::WARNING_QUIT("NonlocalXC", "libvdwxc rVV10 requires a density PW basis");
        }
        return std::unique_ptr<NonlocalXC>(new LibvdwxcNonlocalXC(*pw, b, c));
#else
        ModuleBase::WARNING_QUIT("NonlocalXC", "xc_nonlocal=rvv10 requires a CMake build with ENABLE_LIBVDWXC=ON");
#endif
    }

    ModuleBase::WARNING_QUIT("NonlocalXC", "Unsupported nonlocal XC method: " + method);
    return nullptr;
}

} // namespace vdw
