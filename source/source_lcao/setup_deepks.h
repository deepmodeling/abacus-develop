#ifndef SETUP_DEEPKS_H
#define SETUP_DEEPKS_H

#include "source_cell/unitcell.h" // use unitcell
#include "source_io/module_parameter/input_parameter.h" // Input_para
#include "source_basis/module_ao/parallel_orbitals.h" // parallel orbitals
#include "source_basis/module_ao/ORB_read.h" // orb

#ifdef __MLALGO
#include "source_lcao/module_deepks/LCAO_deepks.h" // deepks
#endif


template <typename TK>
class Setup_DeePKS
{
    public:

    Setup_DeePKS();
    ~Setup_DeePKS();

#ifdef __MLALGO
    LCAO_Deepks<TK> ld;
#endif

	void before_runner(
			const UnitCell& ucell, // unitcell
			const int nks, // k points
            const LCAO_Orbitals &orb, // orbital info
			Parallel_Orbitals &pv, // parallel orbitals
			const Input_para &inp);

};


#endif
