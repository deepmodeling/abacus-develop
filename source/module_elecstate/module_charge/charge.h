#ifndef CHARGE_H
#define CHARGE_H

#include "module_base/complexmatrix.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_pw/pw_basis.h"
#include "src_parallel/parallel_global.h"

//==========================================================
// Electron Charge Density
//==========================================================
class Charge
{

  public:
    Charge();
    ~Charge();

    //==========================================================
    // MEMBER VARIABLES :
    // init_chg : "atomic" or "file"
    // NAME : total number of electrons
    // NAME : rho (nspin,ncxyz), the charge density in real space
    // NAME : rho_save (nspin,ncxyz), for charge mixing
    // NAME : rhog, charge density in G space
    // NAME : rhog_save, chage density in G space
    // NAME : rho_core [nrxx], the core charge in real space
    // NAME : rhog_core [ngm], the core charge in reciprocal space
    //==========================================================

    double **rho;
    double **rho_save;

    std::complex<double> **rhog;
    std::complex<double> **rhog_save;

    double **kin_r; // kinetic energy density in real space, for meta-GGA
    double **kin_r_save; // kinetic energy density in real space, for meta-GGA
                         // wenfei 2021-07-28

    double *rho_core;
    std::complex<double> *rhog_core;

    double *start_mag_type;
    double *start_mag_atom;
    int prenspin = 1;

    void init_rho();
    // mohan update 2021-02-20
    void allocate(const int &nspin_in, const int &nrxx_in, const int &ngmc_in);

    void atomic_rho(const int spin_number_need, double **rho_in, ModulePW::PW_Basis *rho_basis) const;

    void set_rho_core(const ModuleBase::ComplexMatrix &structure_factor);

    void cal_nelec(); // calculate total number of electrons  Yu Liu add 2021-07-03

    void renormalize_rho(void);

    void save_rho_before_sum_band(void);

	// for non-linear core correction
    void non_linear_core_correction
    (
        const bool &numeric,
        const int mesh,
        const double *r,
        const double *rab,
        const double *rhoc,
        double *rhocg,
        ModulePW::PW_Basis* rho_basis
    ) const;

	double check_ne(const double *rho_in) const;

    void init_final_scf(); //LiuXh add 20180619

	public:

    void rho_mpi(void);

    // mohan add 2021-02-20
    int nrxx; // number of r vectors in this processor
    int ngmc; // number of g vectors in this processor
    int nspin; // number of spins
  private:
    double sum_rho(void) const;

    bool allocate_rho;

    bool allocate_rho_final_scf; // LiuXh add 20180606
};

#endif // charge
