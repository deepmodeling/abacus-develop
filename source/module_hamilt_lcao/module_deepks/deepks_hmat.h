#ifndef DEEPKS_HMAT_H 
#define DEEPKS_HMAT_H 

#ifdef __DEEPKS

#include "module_base/matrix.h"
#include "module_base/timer.h"
#include "module_basis/module_ao/parallel_orbitals.h"

#include <torch/script.h>
#include <torch/torch.h>
#include <unordered_map>

namespace DeePKS_domain
{
	void save_h_mat(
			const double *h_mat_in,
			const int nloc);

	void save_h_mat(
			const std::complex<double> *h_mat_in,
			const int nloc);

    //Collect data in h_in to matrix h_out. Note that left lower trianger in h_out is filled
    void collect_h_mat(
        const Parallel_Orbitals &pv,
        const std::vector<double>& h_in,
        ModuleBase::matrix &h_out,
        const int nlocal);

    void check_h_mat(
        const ModuleBase::matrix &H,
        const std::string &h_file,
        const int nlocal);
}

#endif
#endif
