#include "phi_operator_gpu.h"

namespace ModuleGint
{

PhiOperatorGpu::PhiOperatorGpu(cudaStream_t stream)
: stream_(stream)
{
    atoms_num_info_h_.reserve(BatchBigGrid::get_max_batch_size());
    atoms_num_info_h_.reserve(BatchBigGrid::get_max_batch_size());
    atoms_iat_h_.reserve(BatchBigGrid::get_max_atoms_num());
    atoms_iat_d_.reserve(BatchBigGrid::get_max_atoms_num());
    atoms_bgrids_rcoords_h_.reserve(BatchBigGrid::get_max_atoms_num());
    atoms_bgrids_rcoords_d_.reserve(BatchBigGrid::get_max_atoms_num());
}
  


}