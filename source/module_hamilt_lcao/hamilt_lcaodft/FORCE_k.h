#ifndef FORCE_LCAO_K_H
#define FORCE_LCAO_K_H

#include "FORCE_gamma.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/matrix.h"
#include "module_elecstate/module_dm/density_matrix.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_hamilt.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_matrix.h"
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_charge.h"

class Force_LCAO_k : public Force_LCAO_gamma
{
    public:
    template <typename T>
    friend class Force_Stress_LCAO;

    Force_LCAO_k();
    ~Force_LCAO_k();

    private:
    LCAO_Hamilt* UHM;

    // orthonormal force + contribution from T and VNL
    void ftable_k(const bool isforce,
                  const bool isstress,
                  Record_adj& ra,
                  const psi::Psi<std::complex<double>>* psi,
                  Local_Orbital_Charge& loc,
                  const elecstate::ElecState* pelec,
                  ModuleBase::matrix& foverlap,
                  ModuleBase::matrix& ftvnl_dphi,
                  ModuleBase::matrix& fvnl_dbeta,
                  ModuleBase::matrix& fvl_dphi,
                  ModuleBase::matrix& soverlap,
                  ModuleBase::matrix& stvnl_dphi,
                  ModuleBase::matrix& svnl_dbeta,
#ifdef __DEEPKS
                  ModuleBase::matrix& svl_dphi,
                  ModuleBase::matrix& svnl_dalpha,
#else
                  ModuleBase::matrix& svl_dphi,
#endif
				  LCAO_Hamilt &uhm,
				  Parallel_Orbitals &pv,
				  LCAO_Matrix &lm,
                  const K_Vectors &kv);

    // get the ds, dt, dvnl.
    void allocate_k(const Parallel_Orbitals& pv,
                    LCAO_Matrix &lm,
                    const int& nks,
                    const std::vector<ModuleBase::Vector3<double>>& kvec_d);

    void finish_k(LCAO_Matrix &lm);

    // calculate the force due to < dphi | beta > < beta | phi >
	void cal_ftvnl_dphi_k(const elecstate::DensityMatrix<std::complex<double>, double>* DM,
			const Parallel_Orbitals &pv,
			LCAO_Matrix &lm,
			const bool isforce,
			const bool isstress,
			Record_adj& ra,
			ModuleBase::matrix& ftvnl_dphi,
			ModuleBase::matrix& stvnl_dphi);

    // calculate the overlap force
    void cal_foverlap_k(const bool isforce,
                        const bool isstress,
                        Record_adj &ra,
                        const psi::Psi<std::complex<double>> *psi,
						Local_Orbital_Charge &loc,
						Parallel_Orbitals &pv,
						LCAO_Matrix &lm,
						const elecstate::DensityMatrix<std::complex<double>, double> *DM,
                        ModuleBase::matrix &foverlap,
                        ModuleBase::matrix &soverlap,
                        const elecstate::ElecState *pelec,
                        const int &nks,
                        const K_Vectors &kv);

    // calculate the force due to < phi | Vlocal | dphi >
	void cal_fvl_dphi_k(const bool isforce,
		  	const bool isstress,
			LCAO_Matrix &lm,
			const elecstate::Potential* pot_in,
			ModuleBase::matrix& fvl_dphi,
			ModuleBase::matrix& svl_dphi,
			double** DM_R);

    // new method to calculate the force due to < phi | dbeta > < beta | phi > , developed by wenfei-li
    void cal_fvnl_dbeta_k(const elecstate::DensityMatrix<std::complex<double>, double>* DM,
                          const bool isforce,
						  const bool isstress, 
						  const Parallel_Orbitals &pv,
						  ModuleBase::matrix& fvnl_dbeta,
                          ModuleBase::matrix& svnl_dbeta);

	void test(
			Parallel_Orbitals &pv,
			double* mm, 
			const std::string& name);
};

#endif
