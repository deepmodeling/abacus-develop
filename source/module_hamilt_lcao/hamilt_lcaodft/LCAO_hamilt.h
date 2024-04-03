#ifndef LCAO_HAMILT_H
#define LCAO_HAMILT_H

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "LCAO_gen_fixedH.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "module_hamilt_general/hamilt.h"
#include "module_hamilt_lcao/module_gint/grid_technique.h"

#include "module_hamilt_lcao/module_gint/gint_gamma.h"
#include "module_hamilt_lcao/module_gint/gint_k.h"


#ifdef __EXX
#include <RI/global/Tensor.h>
#endif

class LCAO_Hamilt
{
    public:

    LCAO_Hamilt();

    ~LCAO_Hamilt();

	void grid_prepare(
			const Grid_Technique& gt, 
			Gint_Gamma &gint_gamma,
			Gint_k &gint_k,
			const ModulePW::PW_Basis& rhopw, 
			const ModulePW::PW_Basis_Big& bigpw);

    // jingan add 2021-6-4
    void set_R_range_sparse(LCAO_Matrix &lm);

    void cal_HContainer_sparse_d(const int &current_spin, 
        const double &sparse_threshold, 
        const hamilt::HContainer<double>& hR, 
        std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>>& target);

    void cal_HContainer_sparse_cd(const int &current_spin, 
        const double &sparse_threshold, 
        const hamilt::HContainer<std::complex<double>>& hR, 
        std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>>& target);

    void cal_dSTN_R_sparse(const int &current_spin, const double &sparse_threshold);

    // mohan update 2024-04-01
	void cal_STN_R_sparse_for_S(
			const Parallel_Orbitals pv,
			std::vector<double> &slocR,
			std::vector<std::complex<double>> &slocR_soc,
			std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> &SR_sparse,
			std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>> &SR_soc_sparse,
			const double &sparse_threshold);

    void cal_STN_R_sparse_for_T(const double &sparse_threshold);

    void cal_HR_dftu_sparse(const int &current_spin, const double &sparse_threshold);

    void cal_HR_dftu_soc_sparse(const int &current_spin, const double &sparse_threshold);

#ifdef __EXX
    template<typename Tdata> void cal_HR_exx_sparse(
            const int &current_spin,
            const double &sparse_threshold,
            const int (&nmp)[3],
            const std::vector< std::map <int, std::map < std::pair<int, std::array<int,3>>, RI::Tensor<Tdata> > >>& Hexxs);
#endif

	void cal_HSR_sparse(
			const int &current_spin, 
			const double &sparse_threshold, 
			const int (&nmp)[3], 
			hamilt::Hamilt<std::complex<double>>* p_ham);

    void cal_SR_sparse(const double &sparse_threshold, hamilt::Hamilt<std::complex<double>>* p_ham);

    void clear_zero_elements(const int &current_spin, const double &sparse_threshold);

    void destroy_all_HSR_sparse(void);

    void cal_TR_sparse(const double &sparse_threshold);

    void destroy_TR_sparse(void);

	void cal_dH_sparse(
			const int &current_spin, 
			const double &sparse_threshold,
			Gint_k &gint_k); // mohan add 2024-04-01

    void destroy_dH_R_sparse(void);

    // use overlap matrix to generate fixed Hamiltonian
    LCAO_gen_fixedH genH;

    LCAO_Matrix* LM;
};

#endif
