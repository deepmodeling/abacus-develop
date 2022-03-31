#ifndef MULLIKEN_CHARGE_H
#define MULLIKEN_CHARGE_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#ifdef __LCAO
#include "../module_orbital/ORB_gen_tables.h"
#include "module_orbital/ORB_control.h"
#include "../module_neighbor/sltk_grid_driver.h"
#include "../src_lcao/LCAO_matrix.h"
#include "../src_lcao/global_fp.h"
#include "../src_pdiag/pdiag_double.h"
#include "../src_pdiag/GenELPA.h"
#endif
#include "../module_base/lapack_connector.h"
#include "../module_base/scalapack_connector.h"
#include "../module_base/matrix.h"
#include "../module_base/complexmatrix.h"
#include <vector>

// by qifeng
class Mulliken_Charge
{
	public:

	Mulliken_Charge(std::vector<ModuleBase::matrix> *wfc_gamma_in, 
        std::vector<ModuleBase::ComplexMatrix> *wfc_k_in);
	~Mulliken_Charge();

	double**  DecMulP ;
	double**  MecMulP ;
	double***  ADecMulP ;
    std::vector<ModuleBase::matrix> *wfc_gamma;
    std::vector<ModuleBase::ComplexMatrix> *wfc_k;

	std::complex<double> *mug;

	void cal_mulliken(LCAO_Hamilt &uhm);

	void stdout_mulliken(LCAO_Hamilt &uhm );

	private:

};
#endif
