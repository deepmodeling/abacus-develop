#ifndef EXX_ABFS_H
#define EXX_ABFS_H

#include <vector>
using std::vector;
#include <map>
using std::map;
#include <string>

#include "../module_orbital/ORB_atomic_lm.h"
#include "../module_base/element_basis_index.h"
#include "../module_base/matrix.h"
#include "../module_base/vector3.h"

class Exx_Abfs
{
public:

	void test_all() const;	
	void generate_matrix() const;
	void test_abfs1() const;
	void test_abfs2() const;
	void cal_exx() const;
	
	class Util;
	class Abfs_Index;
	class Jle;
	class Inverse_Matrix_Double;
	class IO;
	class Construct_Orbs;
	class PCA;
	class DM;
	class Parallel;
	class Screen
	{public:
		class Schwarz;
		class Cauchy;
	};
	
	class Matrix_Orbs11;
	class Matrix_Orbs21;
	class Matrix_Orbs22;
	
	class Matrix_Lcaoslcaos_Lcaoslcaos;
//	template< typename type_Orb22 > class Matrix_Lcaoslcaos_Lcaoslcaos2;	// Peize Lin test
	
	int rmesh_times = 5;				// Peize Lin test
	int kmesh_times = 1;				// Peize Lin test
//	static int Lmax;					// Peize Lin test
	
public:
	
	static std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<std::vector<ModuleBase::matrix>>>>>> cal_I( 
		const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,ModuleBase::matrix>>>> &ms, 
		const ModuleBase::Element_Basis_Index::IndexLNM &index );
	
	static std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<ModuleBase::matrix>>>>> cal_C( 
		const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<ModuleBase::matrix>>>>> &matrix_A, 
		const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<std::vector<ModuleBase::matrix>>>>>> &matrix_L );
	
	void cal_CVC( 
		const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<ModuleBase::matrix>>>>> &matrix_C, 
		const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,ModuleBase::matrix>>>> &matrix_V ) const;
		
	static std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,ModuleBase::matrix>>>> cal_lcaos2_lcaos2_proj_asa( 
		const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<ModuleBase::matrix>>>>> &ms_lcaos2_asa, 
		const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<std::vector<ModuleBase::matrix>>>>>> &ms_asa_asa_I,
		const ModuleBase::Element_Basis_Index::Range &range,
		const ModuleBase::Element_Basis_Index::IndexLNM &index);
		
	static std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<ModuleBase::matrix>>>>> cal_lcaos2_jys_proj_asa(
		const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<ModuleBase::matrix>>>>> &ms_lcaos2_asa,
		const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<std::vector<ModuleBase::matrix>>>>>> &ms_asa_asa_I,
		const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,ModuleBase::matrix>>>> &ms_asa_jys);	
};

#endif
