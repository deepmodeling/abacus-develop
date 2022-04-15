#ifndef IONS_MOVE_BFGS_H
#define IONS_MOVE_BFGS_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/matrix.h"
//xiaohui modified 2013-04-25
//#include "../src_algorithms/myfunc.h"
#include "bfgs_basic.h"
class Ions_Move_BFGS : public BFGS_Basic
{
public:
	Ions_Move_BFGS();
	~Ions_Move_BFGS();
	
	void allocate(void);
	void start(const ModuleBase::matrix& force,  const double &energy_in);

	private:

	bool init_done;
	void bfgs_routine(void);
	void terminate_bfgs(void);
	void restart_bfgs(void);

};

#endif
