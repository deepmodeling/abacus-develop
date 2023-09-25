/************************************
//LiaoChen Modify on 2010-3-22
***********************************/
#ifndef LCAO_gen_fixedH_H
#define LCAO_gen_fixedH_H

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_basis/module_ao/ORB_gen_tables.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_matrix.h"
#include "module_basis/module_ao/ORB_gen_tables.h"

class LCAO_gen_fixedH
{
    friend class Force_LCAO_gamma;
    friend class Force_LCAO_k;
    friend class LCAO_Hamilt;

  public:
    LCAO_Matrix* LM;

    LCAO_gen_fixedH();
    ~LCAO_gen_fixedH();

    void calculate_NL_no(double* HlocR);
    // void calculate_NL_no(std::complex<double>* HlocR);
    void calculate_T_no(double* HlocR);
    // void calculate_T_no(std::complex<double>* HlocR);
    void calculate_S_no(double* SlocR);
    // void calculate_S_no(std::complex<double>* SlocR);
    void build_ST_new(const char& dtype,
                      const bool& cal_deri,
                      const UnitCell& ucell,
                      double* SHlocR,
                      bool cal_syns = false,
                      double dmax = 0.0);
	// cal_syns : calculate asynchronous overlap matrix for Hefei-NAMD

  private:
    // can used in gamma algorithm.
    void build_Nonlocal_beta_new(double* Hloc);

    void build_Nonlocal_mu_new(double* HlocR, const bool& calc_deri);
};

#endif
