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

#ifdef __MPI
/**
 * @brief transfer the HContainer from serial object to parallel object
 * @param hR_s the HContainer of <I,J,R> atom pairs in serial object
 * @param hR_p the HContainer of <I,J,R> atom pairs in parallel object
 * @param my_rank the rank of current process
*/
template<typename TR>
void transferSerial2Parallels(const hamilt::HContainer<TR>& hR_s,
                             hamilt::HContainer<TR>* hR_p,
                             const int serial_rank);
/**
 * @brief transfer the HContainer from parallel object to serial object
 * @param hR_p the HContainer of <I,J,R> atom pairs in parallel object
 * @param hR_s the HContainer of <I,J,R> atom pairs in serial object
 * @param my_rank the rank of current process
*/
template<typename TR>
void transferParallels2Serial(const hamilt::HContainer<TR>& hR_p,
                             hamilt::HContainer<TR>* hR_s,
                             const int serial_rank);

/**
 * @brief transfer the HContainer from all serial objects to all parallel objects
 * @param hR_s the HContainer of <I,J,R> atom pairs in serial object
 * @param hR_p the HContainer of <I,J,R> atom pairs in parallel object
*/
template<typename TR>
void transferSerials2Parallels(const hamilt::HContainer<TR>& hR_s,
                             hamilt::HContainer<TR>* hR_p);

/**
 * @brief transfer the HContainer from all serial objects to all parallel objects
 * @param hR_p the HContainer of <I,J,R> atom pairs in parallel object
 * @param hR_s the HContainer of <I,J,R> atom pairs in serial object
*/
template<typename TR>
void transferParallels2Serials(const hamilt::HContainer<TR>& hR_p,
                             hamilt::HContainer<TR>* hR_s);

/**
 * @brief gather the HContainer from all parallel objects to target serial object
 * the serial object should be empty before gather
 * @param hR_p the HContainer of <I,J,R> atom pairs in parallel object
 * @param hR_s the empty HContainer of <I,J,R> atom pairs in serial object
 * @param serial_rank the rank of target serial object
*/
template<typename TR>
void gatherParallels(const hamilt::HContainer<TR>& hR_p,
                     hamilt::HContainer<TR>* hR_s,
                     const int serial_rank);

#endif // __MPI

} // namespace hamilt