#ifndef EXX_ABFS_CONSTRUCT_ORBS_H
#define EXX_ABFS_CONSTRUCT_ORBS_H

#include "exx_abfs.h"

#include <limits>
#include "module_ORB/ORB_atomic_lm.h"

class LCAO_Orbitals;

class Exx_Abfs::Construct_Orbs
{
public:
	static vector<vector<vector<Numerical_Orbital_Lm>>> change_orbs( 
		const LCAO_Orbitals &orb_in,
		const double kmesh_times );
	static vector<vector<vector<Numerical_Orbital_Lm>>> change_orbs( 
		const vector<vector<vector<Numerical_Orbital_Lm>>> &orb_in,
		const double kmesh_times );

	static vector<vector<vector<Numerical_Orbital_Lm>>> abfs_same_atom( 
		const vector<vector<vector<Numerical_Orbital_Lm>>> &lcaos,
		const double kmesh_times_mot,
		const double times_threshold=0);
		
	static vector<vector<vector<Numerical_Orbital_Lm>>> orth_orbs( 
		const vector<vector<vector<Numerical_Orbital_Lm>>> &orbs,
		const double norm_threshold=std::numeric_limits<double>::min() );
		
private:
	static vector<vector<vector<vector<double>>>> psi_mult_psi( 
		const vector<vector<vector<Numerical_Orbital_Lm>>> &lcaos );
		
	static vector<vector<vector<vector<double>>>> psir_mult_psir( 
		const vector<vector<vector<Numerical_Orbital_Lm>>> &lcaos );
		
	static vector<vector<vector<vector<double>>>> orth( 
		const vector<vector<vector<vector<double>>>> &psis,
		const vector<vector<vector<Numerical_Orbital_Lm>>> &lcaos,
		const double norm_threshold = std::numeric_limits<double>::min() );

	static vector<vector<vector<vector<double>>>> pca(
		const vector<vector<vector<Numerical_Orbital_Lm>>> &abfs,
		const vector<vector<vector<Numerical_Orbital_Lm>>> &orbs,
		const double kmesh_times_mot,
		const double times_threshold );
		
	static vector<vector<vector<vector<double>>>> div_r( 
		const vector<vector<vector<vector<double>>>> &psirs,
		const vector<double> &r_radial );		
		
	static vector<vector<vector<Numerical_Orbital_Lm>>> orbital(
		const vector<vector<vector<vector<double>>>> &psis,
		const vector<vector<vector<Numerical_Orbital_Lm>>> &orbs_info,
		const double kmesh_times);		
		
	static vector<vector<vector<vector<double>>>> get_psi(
		const vector<vector<vector<Numerical_Orbital_Lm>>> &orbs );
};

#endif	// EXX_ABFS_IO_ASA_H
