#ifndef SETUP_EXX_NAO_H
#define SETUP_EXX_NAO_H

/*
#include "source_cell/unitcell.h" // use unitcell
#include "source_estate/elecstate_lcao.h"// use ElecStateLCAO
#include "source_psi/psi.h" // use electronic wave functions
#include "source_estate/module_charge/charge.h" // use charge
*/

// for EXX
#ifdef __EXX
#include "source_lcao/module_ri/Exx_LRI_interface.h"
#include "source_lcao/module_ri/Mix_DMk_2D.h"
#endif

class Exx_NAO
{
    public:

#ifdef __EXX
    std::shared_ptr<Exx_LRI_Interface<TK, double>> exd = nullptr;
    std::shared_ptr<Exx_LRI_Interface<TK, std::complex<double>>> exc = nullptr;
#endif

    void init();

    void before_runner();



};


#endif
