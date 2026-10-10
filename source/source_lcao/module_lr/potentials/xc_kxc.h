#ifndef ABACUS_LR_XC_KXC_H
#define ABACUS_LR_XC_KXC_H

#ifdef __LIBXC
struct xc_func_type;
namespace LR
{
/// Check initialized runtime functionals and their mixed components before requesting KXC.
void require_xc_kxc(const xc_func_type& functional, bool need_kxc);
}
#endif

#endif
