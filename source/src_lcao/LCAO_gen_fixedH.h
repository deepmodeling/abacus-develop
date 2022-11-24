/************************************
//LiaoChen Modify on 2010-3-22
***********************************/
#ifndef LCAO_gen_fixedH_H
#define LCAO_gen_fixedH_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_orbital/ORB_gen_tables.h"
#include "../module_neighbor/sltk_grid_driver.h"
#include "src_lcao/LCAO_matrix.h"

class LCAO_gen_fixedH
{
	friend class Force_LCAO_gamma;
	friend class Force_LCAO_k;
	friend class energy;//qifeng 2019/9/10
	friend class Mulliken_Charge;//qifeng  2019/9/10

public:

    LCAO_Matrix* LM;

	LCAO_gen_fixedH();
	~LCAO_gen_fixedH();

	void calculate_NL_no(double* HlocR);
	//void calculate_NL_no(std::complex<double>* HlocR);
	void calculate_T_no(double* HlocR);
	//void calculate_T_no(std::complex<double>* HlocR);
	void calculate_S_no(double* SlocR);
	//void calculate_S_no(std::complex<double>* SlocR);

	private:

	void build_ST_new(const char& dtype, const bool& cal_deri, const UnitCell &ucell, double* SHlocR);
	//void build_ST_new(const char& dtype, const bool& cal_deri, const UnitCell &ucell, std::complex<double>* SHlocR);	
	
	// can used in gamma algorithm.
	void build_Nonlocal_beta_new (double* Hloc);

	void build_Nonlocal_mu_new (double* HlocR, const bool &calc_deri); 

};

#endif
