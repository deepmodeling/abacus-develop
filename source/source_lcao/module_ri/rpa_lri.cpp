// The implementation is kept in small internal fragments so the public header
// stays declaration-only and no single source fragment exceeds the repository
// size guideline. These fragments are included into one translation unit so the
// class-template definitions remain visible for explicit instantiation. A few
// large functions continue across adjacent fragments, so include order matters.
#include "rpa_lri_detail.h"
#include "rpa_lri_eigenvector.h"
#include "rpa_lri_coulomb.h"
#include "rpa_lri_overlap.h"
#include "rpa_lri_overlap_io.h"
#include "rpa_lri_eigenvector_io.h"
#include "rpa_lri_output.h"
#include "rpa_lri_output_v1.h"

template class RPA_LRI<double, double>;
template class RPA_LRI<std::complex<double>, double>;
