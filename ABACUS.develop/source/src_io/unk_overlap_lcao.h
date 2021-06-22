#ifndef UNKOVERLAP_LCAO
#define UNKOVERLAP_LCAO

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
#include "src_global/vector3.h"
#include "src_global/ylm.h"

#include "src_pw/global.h"

class unkOverlap_lcao
{
public:

	ORB_table_phi MOT;
	ORB_gaunt_table MGT;
	Numerical_Orbital_Lm orb_r;  // The newly created r vector exists in the form of atomic orbitals and is expanded by solid spherical harmonic functions
	
	vector<vector<vector<Vector3<double>>>> orb1_orb2_R;
	vector<vector<vector<double>>> psi_psi;
	vector<vector<vector<Vector3<double>>>> psi_r_psi;
	bool allocate_flag; // Used to initialize the array
	complex<double>***lcao_wfc_global; // Wave function coefficients under the global lcao basis set
	int** cal_tag; // For parallel scenarios
	
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


	unkOverlap_lcao();
	~unkOverlap_lcao();
	
	
	void init();
	int iw2it(int iw);
	int iw2ia(int iw);
	int iw2iL(int iw);
	int iw2iN(int iw);
	int iw2im(int iw);
	complex<double> unkdotp_LCAO(const int ik_L, const int ik_R, const int iband_L, const int iband_R, const Vector3<double> dk);
	void cal_R_number();
	void cal_orb_overlap();
	void get_lcao_wfc_global_ik(complex<double> **ctot, complex<double> **cc);
	void prepare_midmatrix_pblas(const int ik_L, const int ik_R, const Vector3<double> dk, complex<double> *&midmatrix);
	complex<double> det_berryphase(const int ik_L, const int ik_R, const Vector3<double> dk, const int occ_bands);
	
	void test();
	
	
};






#endif
