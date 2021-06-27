#ifndef MULLIKEN_CHARGE_H
#define MULLIKEN_CHARGE_H

#include "src_pw/tools.h"
#include "module_ORB/ORB_gen_tables.h"
#include "module_neighbor/sltk_grid_driver.h"
#include "src_lcao/LCAO_matrix.h"
#include "src_lcao/global_fp.h"
#include "src_lcao/wfc_dm_2d.h"
#include "src_global/lapack_connector.h"
#include "src_global/scalapack_connector.h"
#include "src_global/matrix.h"
#include "src_global/complexmatrix.h"
#include "src_pdiag/pdiag_double.h"
#include "src_pdiag/GenELPA.h"
#include <vector>

// by qifeng
class Mulliken_Charge
{
	public:

	Mulliken_Charge();
	~Mulliken_Charge();

	double**  DecMulP ;
	double**  MecMulP ;
	double***  ADecMulP ;
	Wfc_Dm_2d M;

	complex<double> *mug;

	void cal_mulliken(void);

	void stdout_mulliken(void);

	private:

};
#endif
