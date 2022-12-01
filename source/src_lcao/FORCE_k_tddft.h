#ifndef FORCE_LCAO_K_TDDFT_H
#define FORCE_LCAO_K_TDDFT_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/matrix.h"
#include "FORCE_gamma.h"
#include "FORCE_k.h"
#include "LCAO_matrix.h"
#include "src_lcao/LCAO_hamilt.h"
#include "src_lcao/local_orbital_charge.h"

class Force_LCAO_k_tddft : public Force_LCAO_k
{
  public:
    friend class Force_Stress_LCAO;
    friend class Force_Stress_LCAO_TDDFT;

    Force_LCAO_k_tddft();
    ~Force_LCAO_k_tddft();

  protected:
    //LCAO_Hamilt* UHM;
    //const Parallel_Orbitals* ParaV;

    // orthonormal force + contribution from T and VNL
    void ftable_k(const bool isforce,
                  const bool isstress,
                  Record_adj& ra,
                  const psi::Psi<std::complex<double>>* psi,
                  Local_Orbital_Charge& loc,
                  ModuleBase::matrix& foverlap,
                  ModuleBase::matrix& ftddft,
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
                  LCAO_Hamilt& uhm,
                  ModuleBase::Vector3<double>* vel);

    // get the ds, dt, dvnl.
    void allocate_k(const Parallel_Orbitals& pv);
	// calculate the force due to < dphi | beta > < beta | phi >
   
   // void cal_ftvnl_dphi_k(double** dm2d, const bool isforce, const bool isstress, Record_adj& ra,
    //    ModuleBase::matrix& ftvnl_dphi, ModuleBase::matrix& stvnl_dphi);

	// calculate the overlap force
    //void cal_foverlap_k(const bool isforce, const bool isstress, Record_adj &ra, const psi::Psi<std::complex<double>>* psi,
      //  Local_Orbital_Charge& loc, ModuleBase::matrix& foverlap, ModuleBase::matrix& soverlap);

	// calculate the force due to < phi | Vlocal | dphi >
	//void cal_fvl_dphi_k(const bool isforce, const bool isstress, ModuleBase::matrix& fvl_dphi, ModuleBase::matrix& svl_dphi, double **DM_R);

  // old method to calculate the force due to < phi | dbeta > < beta | phi >
	//void cal_fvnl_dbeta_k(double** dm2d, const bool isforce, const bool isstress, ModuleBase::matrix& fvnl_dbeta, ModuleBase::matrix& svnl_dbeta);
	// new method to calculate the force due to < phi | dbeta > < beta | phi > , developed by wenfei-li
	//void cal_fvnl_dbeta_k_new(double** dm2d, const bool isforce, const bool isstress, ModuleBase::matrix& fvnl_dbeta, ModuleBase::matrix& svnl_dbeta);


    // calculate the overlap force
    void cal_ftddft_k(const bool isforce,
                      Record_adj& ra,
                      const psi::Psi<std::complex<double>>* psi,
                      Local_Orbital_Charge& loc,
                      ModuleBase::matrix& ftddft,
                      ModuleBase::Vector3<double>* vel);

    // calculate the force due to < phi | dbeta > < beta | phi >
	/*void calFvnlDbeta(
		double** dm2d, 
		const bool &isforce, 
		const bool &isstress, 
		ModuleBase::matrix& fvnl_dbeta, 
		ModuleBase::matrix& svnl_dbeta,
		const int &vnl_method);*/

    void finish_k(void);
};
#endif
