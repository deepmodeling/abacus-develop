#include "module_base/global_variable.h"
#include "module_base/global_function.h"

#include "module_basis/module_ao/ORB_control.h"	
#include "module_basis/module_ao/ORB_read.h"

#include "module_cell/unitcell.h"

#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_cell/module_neighbor/sltk_atom_arrange.h"
#include "klist.h"
#include "parallel_orbitals.h"

#include "../LCAO_deepks.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace Test_Deepks
{
    extern Grid_Driver GridD;
}

class test_deepks
{

public:
    test_deepks();
    ~test_deepks();

	LCAO_Orbitals ORB;
	ORB_gen_tables OGT;
	ORB_control ooo;

	UnitCell ucell;

    Parallel_Orbitals ParaO;
	Test_Deepks::K_Vectors kv;
	LCAO_Deepks ld;

	int failed_check = 0;
	int total_check = 0;

	int my_rank = 0;

	double lcao_ecut = 0; // (Ry)
	double lcao_dk = 0.01;
	double lcao_dr = 0.01;
	double lcao_rmax = 30; // (a.u.)

	int out_mat_r = 0;

	int lmax=2;
	int ntype=0;
	int nnr;

	std::vector<ModuleBase::matrix> dm;
	std::vector<ModuleBase::ComplexMatrix> dm_k;

//preparation
	void preparation();
	void set_parameters();//set some global variables
	void setup_cell();

	void count_ntype(); //from STRU, count types of elements
	void set_ekcut();	//from LCAO files, read and set ekcut

	void prep_neighbour();
	void setup_kpt();
	void set_orbs(const double &lat0_in);

	void cal_nnr();
	void folding_nnr(const Test_Deepks::K_Vectors &kv);

//checking
	void check_dstable(void);
	void check_psialpha(void);

	void read_dm(void);
	void read_dm_k(const int nks);

	void check_pdm(void);
	void check_gdmx(void);

	void check_descriptor(void);
	void check_gvx(void);

	void check_edelta(void);

	void check_v_delta(void);

	void check_e_deltabands(void);
	void check_f_delta(void);

	//compares numbers stored in two files
	void compare_with_ref(
		const std::string f1,
		const std::string f2);
};
