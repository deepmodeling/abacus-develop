#ifndef LIBXC_HOST_H
#define LIBXC_HOST_H

#include <xc.h>

namespace XC_Functional_Libxc
{
// A CUDA-enabled Libxc defaults to device pointers. Keep host callers explicit.
inline int init_func_host(xc_func_type* func, const int id, const int spin)
{
#ifdef XC_FLAGS_ON_HOST
    return xc_func_init_flags(func, id, spin, XC_FLAGS_ON_HOST);
#else
    return xc_func_init(func, id, spin);
#endif
}
} // namespace XC_Functional_Libxc

#endif
