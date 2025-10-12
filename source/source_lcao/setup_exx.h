#ifndef SETUP_EXX_NAO_H
#define SETUP_EXX_NAO_H

#include "source_cell/unitcell.h" // use unitcell
#include "source_cell/klist.h" // k points
#include "source_io/module_parameter/input_parameter.h" // Input_para
#include "source_basis/module_ao/parallel_orbitals.h" // parallel orbitals
#include "source_basis/module_ao/ORB_read.h" // orb

// for EXX
#ifdef __EXX
#include "source_lcao/module_ri/Exx_LRI_interface.h"
#include "source_lcao/module_ri/Mix_DMk_2D.h"
#endif

template <typename TK>
class Exx_NAO
{
    public:

    Exx_NAO();
    ~Exx_NAO();

#ifdef __EXX
    std::shared_ptr<Exx_LRI_Interface<TK, double>> exd = nullptr;
    std::shared_ptr<Exx_LRI_Interface<TK, std::complex<double>>> exc = nullptr;
#endif

    void init();

	void before_runner(
			UnitCell& ucell, // unitcell
			K_Vectors &kv, // k points
            const LCAO_Orbitals &orb, // orbital info
			const Parallel_Orbitals &pv, // parallel orbitals
			const Input_para& inp);

};


#endif
