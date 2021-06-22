#ifndef LCAO_ORBITALS_H
#define LCAO_ORBITALS_H

#include "ORB_atomic.h"
#include "ORB_atomic_lm.h"
#include "ORB_nonlocal.h"

////////////////////////////////////////////////////////////
/// advices for reconstructions:
/// -------------------------------
/// each set of orbitals should have: lmax, dr, dk, rmax, lmax, etc.
///
/// the orbitals include : NAO, non-local projectors, descriptors, etc.
///
/// mohan note 2021-02-13
///////////////////////////////////////////////////////////

class LCAO_Orbitals
{
	public:

	LCAO_Orbitals();
	~LCAO_Orbitals();

	void Read_Orbitals(
		ofstream &ofs_in, // mohan add 2021-05-07
		const int &ntype_in,
		const int &lmax_in,
		const int &out_descriptor, //  mohan add 2021-04-25
		const int &out_r_matrix, // mohan add 2021-04-26
		const bool &force_flag, // mohan add 2021-05-07
		const int &my_rank); // mohan add 2021-04-26

	void Read_PAO(
		ofstream &ofs_in,
		const int& it,
		const bool &force_flag, // mohan add 2021-05-07
		const int& my_rank); // mohan add 2021-04-26

#ifdef __NORMAL

#else
	/// in order to get rid of the .NONLOCAL file.
	void Set_NonLocal(const int &it, int &n_projectors);

	/// read in the NONLOCAL projector from file.
	void Read_NonLocal(const int& it, int &n_projectors, const int &my_rank);
#endif


	void Read_Descriptor(
		ofstream &ofs_in,
		const bool &force_flag, // mohan add 2021-05-07
		const int &my_rank);	//caoyu add 2020-3-16

#ifdef __MPI
	void bcast_files(const int &ntype_in, const int &my_rank);
#endif

	const double& get_ecutwfc(void) const {return ecutwfc;}
	const int& get_kmesh(void) const{return kmesh;}
	const double& get_dk(void) const {return dk;}
	const double& get_dR(void) const {return dR;}
	const double& get_Rmax(void) const {return Rmax;}
	const int& get_lmax(void) const {return lmax;}
	const int& get_lmax_d(void) const { return lmax_d; }		///<lmax of descriptor basis
	const int& get_nchimax(void) const {return nchimax;}
	const int& get_nchimax_d(void) const { return nchimax_d; }	///<nchimax of descriptor basis
	const int& get_ntype(void) const {return ntype;}
	const double& get_dr_uniform(void) const { return dr_uniform; }

	//caoyu add 2021-05-24
	const double& get_rcutmax_Phi(void) const { return rcutmax_Phi; }
	const double& get_rcutmax_Beta(void) const { return rcutmax_Beta; }

	/// numerical atomic orbitals
	Numerical_Orbital* Phi;
	
	/// nonlocal projectors (1-dimension array)
	Numerical_Nonlocal* Beta;
	
	//caoyu add 2021-3-10
	/// descriptor bases, saved as one-type atom orbital
	Numerical_Orbital* Alpha;

	// initialized in input.cpp
	double ecutwfc;
	double dk;
	double dR;
	double Rmax;
	int *nproj; //mohan add 2010-12-19
	int nprojmax; // mohan add 2010-03-07
	
	double dr_uniform;

	// initalized in unitcell_pseudo
	// assume ntype < 20.
	bool read_in_flag;
	std::vector<string> orbital_file;
	std::vector<string> nonlocal_file;
	string descriptor_file;	//caoyu add 2020-3-16

private:

	int ntype; // number of elements
	int kmesh; // number of points on kmesh

	int lmax;
	int nchimax;

	int lmax_d;	//caoyu add 2021-03-17
	int nchimax_d;	//caoyu add 2021-03-17

	double rcutmax_Phi;	//caoyu add 2021-05-24
	double rcutmax_Beta;	//caoyu add 2021-05-24

	void read_orb_file(
		ofstream &ofs_in,
		ifstream &ifs, 
		const int &it, 
		int &lmax, 
		int &nchimax, 
		Numerical_Orbital* ao,
		const bool &force_flag, // mohan add 2021-05-07
		const int &my_rank);	//caoyu add 2021-04-26

};

/// PLEASE avoid using 'ORB' as global variable 
///
///mohan note 2021 - 03 - 23
extern LCAO_Orbitals ORB;
#endif
