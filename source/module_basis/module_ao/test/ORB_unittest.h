#ifndef _ORBUNITTEST_
#define _ORBUNITTEST_

#include "gtest/gtest.h"
#include "module_basis/module_ao/ORB_control.h"
#include "module_base/global_function.h"	
#include "module_hamilt_lcao/hamilt_lcaodft/center2_orb-orb11.h"
//#include "mock_center2.h"
#include <fstream>
#include <iomanip>
#include <iostream>

#include <vector>
#include <map>
#include <set>

class test_orb : public testing::Test
{
protected:
	void SetUp() override;
	void TearDown() override;
public:

	LCAO_Orbitals ORB;
	ORB_gen_tables OGT;
	ORB_gaunt_table Center2_MGT;	//gaunt table used in center2orb
	ORB_control ooo;
	std::ofstream ofs_running;

	std::map < size_t,
        std::map<size_t,
            std::map<size_t,
                std::map<size_t,
                    std::map<size_t,
                        std::map<size_t, 
							std::unique_ptr<Center2_Orb::Orb11>>>>>>> test_center2_orb11;
/*
	std::map < size_t,
        std::map<size_t,
            std::map<size_t,
                std::map<size_t,
                    std::map<size_t,
                        std::map<size_t, 
							std::unique_ptr<MockCenter2Orb11>>>>>>> mock_center2_orb11;
*/
	void count_ntype(); //from STRU, count types of elements
	void set_files();   //from STRU, read names of LCAO files
	void set_ekcut();	//from LCAO files, read and set ekcut
	void set_orbs();	//interface to Read_PAO
	void set_center2orbs();	//interface to Center2orb	
	template <class c2o>
	void set_single_c2o(int TA, int TB, int LA, int NA, int LB, int NB);
	double randr(double Rmax);
	void gen_table_center2();
	

	bool force_flag = 0;
	int my_rank = 0;
	int ntype_read;

	double lcao_ecut = 0; // (Ry)
	double lcao_dk = 0.01;
	double lcao_dr = 0.01;
	double lcao_rmax = 30; // (a.u.)

	int out_descriptor = 0;
	int out_mat_r = 0;

	int lmax=1;
	double lat0 = 1.0;
	string case_dir = "./GaAs/";
    string *orbital_fn;
    string descriptor_file;
};
#endif
