#ifndef LOOP_CELL_H
#define LOOP_CELL_H
#include<vector>
#include "module_base/matrix.h"
#include "module_base/complexmatrix.h"
#include "module_orbital/ORB_control.h"
#include "src_lcao/LCAO_matrix.h"
#include "module_esolver/esolver.h"

class LOOP_cell
{
	public:

	LOOP_cell(Parallel_Orbitals &pv);
	~LOOP_cell();

	void opt_cell(ORB_control &orb_con, ModuleESolver::ESolver *p_esolver);

private:
    LCAO_Matrix LM;
};

#endif
