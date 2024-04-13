#ifndef SPARSE_FORMAT_HSR_H 
#define SPARSE_FORMAT_HSR_H

#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_matrix.h"

namespace sparse_format
{

	void cal_HSR(
            const Parallel_Orbitals &pv,
			LCAO_Matrix &lm,
			Grid_Driver &grid,
			const int &current_spin, 
			const double &sparse_thr, 
			const int (&nmp)[3], 
			hamilt::Hamilt<std::complex<double>>* p_ham);

	void cal_HContainer_d(
            const Parallel_Orbitals &pv,
			const int &current_spin, 
			const double &sparse_threshold, 
			const hamilt::HContainer<double>& hR, 
			std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>>& target);

	void cal_HContainer_cd(
            const Parallel_Orbitals &pv,
			const int &current_spin, 
			const double &sparse_threshold, 
			const hamilt::HContainer<std::complex<double>>& hR, 
			std::map<Abfs::Vector3_Order<int>, 
			std::map<size_t, std::map<size_t, std::complex<double>>>>& target);

	void clear_zero_elements(
			LCAO_Matrix &lm,
			const int &current_spin, 
			const double &sparse_thr);

}

#endif
