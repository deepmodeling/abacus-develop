#ifndef CAL_R_OVERLAP_R_H
#define CAL_R_OVERLAP_R_H

#include<vector>
using std::vector;
#include<map>
using std::map;
#include<set>
using std::set;

#include "src_lcao/center2_orb.h"
#include "src_lcao/center2_orb-orb11.h"
#include "src_lcao/center2_orb-orb21.h"

#include "module_ORB/ORB_table_phi.h"
#include "module_ORB/ORB_gaunt_table.h"
#include "module_ORB/ORB_atomic_lm.h"
#include "module_ORB/ORB_read.h"
#include "module_base/vector3.h"
#include "src_global/ylm.h"
#include "src_lcao/global_fp.h"

#include "src_pw/global.h"


// output r_R matrix, added by Jingan
class cal_r_overlap_R
{

	public:

	cal_r_overlap_R();
	~cal_r_overlap_R();
	

	ORB_table_phi MOT;

	ORB_gaunt_table MGT;

	Numerical_Orbital_Lm orb_r;

	vector<vector<vector<Numerical_Orbital_Lm>>> orbital_phi;
	
	int R_x_num;
    int R_y_num;
    int R_z_num;
	int R_minX;
	int R_minY;
	int R_minZ;
	
	Vector3<double> ****psi_r_psi;
	bool allocate_psi_r_psi = false;

	map<size_t,
		map<size_t,
			map<size_t,
				map<size_t,
					map<size_t,
						map<size_t,
						Center2_Orb::Orb11>>>>>> center2_orb11;
	
	map<size_t,
		map<size_t,
			map<size_t,
				map<size_t,
					map<size_t,
						map<size_t,
						Center2_Orb::Orb21>>>>>> center2_orb21_r;
						
	void init();
	void out_r_overlap_R(const int nspin);
	
	int iw2it(int iw);
	int iw2ia(int iw);
	int iw2iL(int iw);
	int iw2iN(int iw);
	int iw2im(int iw);					
	
};
#endif

