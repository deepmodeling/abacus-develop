/************************************
//LiaoChen Modify on 2010-3-22
***********************************/
#ifndef LCAO_gen_fixedH_H
#define LCAO_gen_fixedH_H

#include "../src_pw/tools.h"
#include "../module_ORB/ORB_gen_tables.h"
#include "../src_global/sltk_grid_driver.h"

class LCAO_gen_fixedH
{
	friend class Force_LCAO_gamma;
	friend class Force_LCAO_k;
	friend class energy;//qifeng 2019/9/10
	friend class Mulliken_Charge;//qifeng  2019/9/10

	public:

	LCAO_gen_fixedH();
	~LCAO_gen_fixedH();

	void calculate_NL_no(void);
	void calculate_T_no(void);
	void calculate_S_no(void);

	private:

	void build_ST_new(const char& dtype, const bool& cal_deri);	
	
	// can used in gamma algorithm.
	void build_Nonlocal_beta (const bool& calc_deri);

	// used if k point is used,
	// also can used in gamma algorithm.
	void build_Nonlocal_mu (const bool &calc_deri); 
	void test_Nonlocal(void);

};

#endif
