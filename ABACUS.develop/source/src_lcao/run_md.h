#ifndef RUN_MD_H
#define RUN_MD_H 

#include "LOOP_elec.h"
#include "../src_ions/ions_move_methods.h"
#include "../src_pw/charge_extra.h"
#include "../src_pw/md.h"
#include "../src_pw/mdNVT.h"
#include "../src_pw/mdNVE.h"
#include "../src_ions/lattice_change_methods.h"

class Run_MD 
{

	public:

	Run_MD();
	~Run_MD();

	LOOP_elec LOE;

	void opt_cell(void);
	void opt_ions(void);

	//2014-06-06, xiaohui
	mdnvt MDNVT ;
	mdNVE MDNVE ;

	private:

	Ions_Move_Methods IMM;

	//bool force_stress(void);
	Lattice_Change_Methods LCM;

	bool force_stress(const int &istep, int &force_step, int &stress_step);

	int istep;

	// electron charge density extropolation method
	Charge_Extra CE;

	//choose md ensemble, zheng daye
	int mdtype;

	void final_scf(void);

};

#endif
