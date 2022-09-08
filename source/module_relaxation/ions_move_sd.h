#ifndef IONS_MOVE_SD_H
#define IONS_MOVE_SD_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/matrix.h"
class Ions_Move_SD
{
public:
    Ions_Move_SD();
    ~Ions_Move_SD();

    void allocate(void);
    void start(const ModuleBase::matrix& force, const double &etot);

private:

	double energy_saved;
    double* pos_saved;
    double* grad_saved;

	void cal_tradius_sd(void)const;
};

#endif
