//=======================
#ifndef RUN_MD_PW_H
#define RUN_MD_PW_H

#include "../src_pw/electrons.h"
#include "../src_pw/tools.h"
#include "../src_pw/charge_extra.h"
#include "../src_pw/sto_elec.h" //mohan added 2021-01-28
#include "../src_ions/ions_move_methods.h"
#include "../src_ions/lattice_change_methods.h"
#include "src_pw/global.h"

class Run_MD_PW
{
public:
    Run_MD_PW(){};
    ~Run_MD_PW(){};

    void md_ions_pw();
    void md_cells_pw();

private:
    Electrons elec;
    Stochastic_Elec elec_sto;
    int istep;
    Ions_Move_Methods IMM;
    Charge_Extra CE;
    Lattice_Change_Methods LCM;
};

#endif