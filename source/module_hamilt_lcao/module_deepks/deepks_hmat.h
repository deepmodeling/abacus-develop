#ifndef DEEPKS_HMAT_H 
#define DEEPKS_HMAT_H 

#ifdef __DEEPKS

#include "module_base/matrix.h"
#include "module_base/complexmatrix.h"
#include "module_base/timer.h"
#include "module_basis/module_ao/parallel_orbitals.h"

#include <torch/script.h>
#include <torch/torch.h>
#include <unordered_map>

namespace DeePKS_domain
{
    //Collect data in h_in to matrix h_out. Note that left lower trianger in h_out is filled
    template <typename TK, typename TH>
    void collect_h_mat(
        const Parallel_Orbitals &pv,
        const std::vector<std::vector<TK>>& h_in,
        std::vector<TH> &h_out,
        const int nlocal,
        const int nks);

    // write h_mat to file h_file for checking // not used in the code now
    template <typename TH>
    void check_h_mat(
		const std::vector<TH> &H,
		const std::string &h_file,
		const int nlocal,
        const int nks);
}

#endif
#endif
