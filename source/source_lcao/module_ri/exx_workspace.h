#ifndef EXX_WORKSPACE_H
#define EXX_WORKSPACE_H

#include <cstddef>

namespace RI
{
template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
class Exx;
}

// The destination must be fresh, with its own parallel context for the same geometry.
// Tensor storage is shared read-only; map containers and electronic state are independent.
template<typename Tdata>
void share_exx_geometry(const RI::Exx<int, int, 3, Tdata>& source,
                       RI::Exx<int, int, 3, Tdata>& destination);

#endif
