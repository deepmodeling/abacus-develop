#ifndef SETUP_PSI_H
#define SETUP_PSI_H

#include "source_psi/psi_init.h"
#include "source_cell/klist.h"
#include "source_io/module_parameter/input_parameter.h"
#include "source_basis/module_ao/parallel_orbitals.h" // use para_orb

template <typename T>
class Setup_Psi
{
    public:

    Setup_Psi();
    ~Setup_Psi();

	static void allocate_psi(
		psi::Psi<T>* &psi,
		const K_Vectors &kv,
        const Parallel_Orbitals &para_orb,
		const Input_para &inp);

    static void deallocate_psi(psi::Psi<T>* &psi);

};


#endif
