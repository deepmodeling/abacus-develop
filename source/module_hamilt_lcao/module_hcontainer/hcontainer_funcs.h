#pragma once

#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"

namespace hamilt
{
/**
 * @brief calculate the Hk matrix with specific k vector
 * @param hR the HContainer of <I,J,R> atom pairs
 * @param hk the data pointer of Hk matrix, the size of hk would be nrow * ncol
 * @param kvec_d_in the k vector in Direct coordinate
 * @param hk_ld the leading dimension number of hk, ncol for row-major, nrow for column-major
 * @param hk_type the data-type of hk, 0 is row-major, 1 is column-major
*/
template<typename TR>
void folding_HR(const hamilt::HContainer<TR>& hR,
                std::complex<double>* hk,
                const ModuleBase::Vector3<double>& kvec_d_in,
                const int ncol,
                const int hk_type);

void folding_HR(const hamilt::HContainer<double>& hR,
                double* hk,
                const ModuleBase::Vector3<double>& kvec_d_in,
                const int ncol,
                const int hk_type);

} // namespace hamilt