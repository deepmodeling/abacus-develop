#ifndef EXX_ABFS_IO_H
#define EXX_ABFS_IO_H

#include "exx_abfs.h"

#include <map>
#include <vector>
#include "../module_orbital/ORB_atomic_lm.h"
#include "../module_base/matrix.h"
#include "../module_base/element_basis_index.h"
#ifdef __MPI
#include "mpi.h"
#endif

class LCAO_Orbitals;

class Exx_Abfs::IO
{
public:
	static void print_matrix( 
		const std::string &file_name_prefix, 
		const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<ModuleBase::matrix>>>>> &matrixes_Q, 
		const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,ModuleBase::matrix>>>> &matrixes_S,
		const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,ModuleBase::matrix>>>> &matrixes_V,
		const ModuleBase::Element_Basis_Index::Range &range_jles, 
		const ModuleBase::Element_Basis_Index::IndexLNM &index_jles, 
		const ModuleBase::Element_Basis_Index::Range &range_lcaos,
		const ModuleBase::Element_Basis_Index::IndexLNM &index_lcaos );
		
	static std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> construct_abfs( 
		const LCAO_Orbitals &orbs,
		const std::vector<std::string> &files_abfs,
		const double kmesh_times=1 );				// close dK, keep Kcut	
		
	static std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> construct_abfs( 
		const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &abfs_pre, 	
		const LCAO_Orbitals &orbs,
		const std::vector<std::string> &files_abfs,
		const double kmesh_times=1 );				// close dK, keep Kcut
		
	template<typename T>
	static void output_binary( const T &data, const std::string &file_name );
	template<typename T>
	static T input_binary( const std::string &file_name );
	template<typename T>
	static void output_text( const T &data, const std::string &file_name );
	template<typename T>
	static T input_text( const std::string &file_name );
#ifdef __MPI
	template<typename T>
	static void bcast( T &data, const int rank_src, MPI_Comm mpi_comm );
#endif
	
private:
	static std::vector<std::vector<Numerical_Orbital_Lm>> construct_abfs_T(
		const std::string & file_name,
		const int &T,
		const int &nk,
		const double &dk,
		const double &dr_uniform);
};

#endif	// EXX_ABFS_IO_H
