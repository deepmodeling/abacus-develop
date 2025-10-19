#ifndef SETUP_PSI_LCAO_H
#define SETUP_PSI_LCAO_H

#include "source_psi/psi_init.h"
#include "source_cell/klist.h"
#include "source_io/module_parameter/input_parameter.h"

template <typename T>
class Setup_Psi_LCAO
{
    public:

    Setup_Psi_LCAO();
    ~Setup_Psi_LCAO();

    static void before_runner(
		const K_Vectors &kv,
		const Input_para &inp);

};


#endif
