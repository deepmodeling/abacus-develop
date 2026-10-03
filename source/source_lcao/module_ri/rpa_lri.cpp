// The implementation is kept in small internal fragments so the public header
// stays declaration-only and no single source fragment exceeds the repository
// size guideline. These fragments are included into one translation unit so the
// class-template definitions remain visible for explicit instantiation.
#include "rpa_lri_impl_00.h"
#include "rpa_lri_impl_01.h"
#include "rpa_lri_impl_02.h"
#include "rpa_lri_impl_03.h"
#include "rpa_lri_impl_04.h"
#include "rpa_lri_impl_05.h"
#include "rpa_lri_impl_06.h"
#include "rpa_lri_impl_07.h"

#ifdef __EXX
template class RPA_LRI<double, double>;
template class RPA_LRI<std::complex<double>, double>;
#endif
