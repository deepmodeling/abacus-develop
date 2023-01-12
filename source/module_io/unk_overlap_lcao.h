#ifndef UNKOVERLAP_LCAO
#define UNKOVERLAP_LCAO

#include <vector>
#include <map>
#include <set>

#include "../src_lcao/center2_orb.h"
#include "../src_lcao/center2_orb-orb11.h"
#include "../src_lcao/center2_orb-orb21.h"
#include "../src_lcao/local_orbital_wfc.h"

#include "module_orbital/ORB_table_phi.h"
#include "module_orbital/ORB_gaunt_table.h"
#include "module_orbital/ORB_atomic_lm.h"
#include "module_orbital/ORB_read.h"
#include "module_orbital/parallel_orbitals.h"
#include "../module_base/vector3.h"
#include "../module_base/ylm.h"


class unkOverlap_lcao
{
public:

	ORB_table_phi MOT;
	ORB_gaunt_table MGT;
	Numerical_Orbital_Lm orb_r;  // 新建的r矢量,以原子轨道形式存在,以solid球谐函数展开
	
	std::vector<std::vector<std::vector<ModuleBase::Vector3<double>>>> orb1_orb2_R;
	std::vector<std::vector<std::vector<double>>> psi_psi;
	std::vector<std::vector<std::vector<ModuleBase::Vector3<double>>>> psi_r_psi;
	bool allocate_flag; // 用于初始化数组的
	std::complex<double>***lcao_wfc_global; // 全局的lcao基组下的波函数系数
	int** cal_tag; // 用于并行方案
	
	std::map<size_t,
		std::map<size_t,
			std::map<size_t,
				std::map<size_t,
					std::map<size_t,
						std::map<size_t,
						Center2_Orb::Orb11>>>>>> center2_orb11;
						
	std::map<size_t,
		std::map<size_t,
			std::map<size_t,
				std::map<size_t,
					std::map<size_t,
						std::map<size_t,
						Center2_Orb::Orb21>>>>>> center2_orb21_r;


	unkOverlap_lcao();
	~unkOverlap_lcao();
	
	
	void init(std::complex<double>*** wfc_k_grid);
	int iw2it(int iw);
	int iw2ia(int iw);
	int iw2iL(int iw);
	int iw2iN(int iw);
	int iw2im(int iw);
	std::complex<double> unkdotp_LCAO(const int ik_L, const int ik_R, const int iband_L, const int iband_R, const ModuleBase::Vector3<double> dk);
	void cal_R_number();
	void cal_orb_overlap();
	void get_lcao_wfc_global_ik(std::complex<double> **ctot, std::complex<double> **cc);
	void prepare_midmatrix_pblas(const int ik_L, const int ik_R, const ModuleBase::Vector3<double> dk, std::complex<double> *&midmatrix, const Parallel_Orbitals &pv);
    std::complex<double> det_berryphase(const int ik_L, const int ik_R,
        const ModuleBase::Vector3<double> dk, const int occ_bands,
        Local_Orbital_wfc &lowf,
		const psi::Psi<std::complex<double>>* psi_in);

	void test(std::complex<double>*** wfc_k_grid);
	
	
};






#endif
