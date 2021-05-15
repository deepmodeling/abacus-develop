#ifndef EXX_ABFS_IO_H
#define EXX_ABFS_IO_H

#include "exx_abfs.h"

#include <map>
#include <vector>
#include "module_ORB/ORB_atomic_lm.h"
#include "src_global/matrix.h"
#include "src_global/element_basis_index.h"

class LCAO_Orbitals;

class Exx_Abfs::IO
{
public:
	static void print_matrix( 
		const string &file_name_prefix, 
		const map<size_t,map<size_t,map<size_t,map<size_t,vector<matrix>>>>> &matrixes_Q, 
		const map<size_t,map<size_t,map<size_t,map<size_t,matrix>>>> &matrixes_S,
		const map<size_t,map<size_t,map<size_t,map<size_t,matrix>>>> &matrixes_V,
		const Element_Basis_Index::Range &range_jles, 
		const Element_Basis_Index::IndexLNM &index_jles, 
		const Element_Basis_Index::Range &range_lcaos,
		const Element_Basis_Index::IndexLNM &index_lcaos );
		
	static vector<vector<vector<Numerical_Orbital_Lm>>> construct_abfs( 
		const LCAO_Orbitals &orbs,
		const vector<string> &files_abfs,
		const double kmesh_times=1 );				// close dK, keep Kcut	
		
	static vector<vector<vector<Numerical_Orbital_Lm>>> construct_abfs( 
		const vector<vector<vector<Numerical_Orbital_Lm>>> &abfs_pre, 	
		const LCAO_Orbitals &orbs,
		const vector<string> &files_abfs,
		const double kmesh_times=1 );				// close dK, keep Kcut
		
	template<typename T>
	static void output_binary( const T &data, const string &file_name );
	template<typename T>
	static T input_binary( const string &file_name );
	template<typename T>
	static void output_text( const T &data, const string &file_name );
	template<typename T>
	static T input_text( const string &file_name );
	template<typename T>
	static void bcast( T &data, const int rank_src, MPI_Comm mpi_comm );
	
private:
	static vector<vector<Numerical_Orbital_Lm>> construct_abfs_T(
		const string & file_name,
		const int &T,
		const int &nk,
		const double &dk,
		const double &dr_uniform);
};

#endif	// EXX_ABFS_IO_H