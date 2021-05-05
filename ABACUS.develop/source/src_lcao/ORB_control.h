//==========================================================
// AUTHOR : mohan, ywcui
// Last Update: 2021-02-10
//==========================================================
#ifndef ORB_CONTROL_H 
#define ORB_CONTROL_H 

#include "../src_pw/tools.h"

#include "ORB_gen_tables.h"
#include "../src_pdiag/pdiag_double.h"

class ORB_control 
{

	public:

    ORB_control();
    ~ORB_control();

    // Generate the S(overlap),T,NL matrix.
    void set_orb_tables();
    void clear_after_ions();

};
#endif
