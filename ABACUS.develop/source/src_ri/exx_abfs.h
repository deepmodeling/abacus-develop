#ifndef EXX_ABFS_H
#define EXX_ABFS_H

#include<vector>
using std::vector;
#include<map>
using std::map;
#include<string>

#include "module_ORB/ORB_atomic_lm.h"
#include "src_global/element_basis_index.h"
#include "src_global/matrix.h"
#include "src_global/vector3.h"

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
	static int Lmax;					// Peize Lin test
	
public:
	
	static map<size_t,map<size_t,map<size_t,map<size_t,vector<vector<matrix>>>>>> cal_I( 
		const map<size_t,map<size_t,map<size_t,map<size_t,matrix>>>> &ms, 
		const Element_Basis_Index::IndexLNM &index );
	
	static map<size_t,map<size_t,map<size_t,map<size_t,vector<matrix>>>>> cal_C( 
		const map<size_t,map<size_t,map<size_t,map<size_t,vector<matrix>>>>> &matrix_A, 
		const map<size_t,map<size_t,map<size_t,map<size_t,vector<vector<matrix>>>>>> &matrix_L );
	
	void cal_CVC( 
		const map<size_t,map<size_t,map<size_t,map<size_t,vector<matrix>>>>> &matrix_C, 
		const map<size_t,map<size_t,map<size_t,map<size_t,matrix>>>> &matrix_V ) const;
		
	static map<size_t,map<size_t,map<size_t,map<size_t,matrix>>>> cal_lcaos2_lcaos2_proj_asa( 
		const map<size_t,map<size_t,map<size_t,map<size_t,vector<matrix>>>>> &ms_lcaos2_asa, 
		const map<size_t,map<size_t,map<size_t,map<size_t,vector<vector<matrix>>>>>> &ms_asa_asa_I,
		const Element_Basis_Index::Range &range,
		const Element_Basis_Index::IndexLNM &index);
		
	static map<size_t,map<size_t,map<size_t,map<size_t,vector<matrix>>>>> cal_lcaos2_jys_proj_asa(
		const map<size_t,map<size_t,map<size_t,map<size_t,vector<matrix>>>>> &ms_lcaos2_asa,
		const map<size_t,map<size_t,map<size_t,map<size_t,vector<vector<matrix>>>>>> &ms_asa_asa_I,
		const map<size_t,map<size_t,map<size_t,map<size_t,matrix>>>> &ms_asa_jys);	
};

#endif
