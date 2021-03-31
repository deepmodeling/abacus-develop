#ifndef LCAO_ORBITALS_H
#define LCAO_ORBITALS_H

#include "../src_pw/tools.h"
#include "ORB_atomic.h"
#include "ORB_atomic_lm.h"
#include "ORB_nonlocal.h"

//---------------------------------------------------------------------
// advices for reconstructions:
// each set of orbitals should have: lmax, dr, dk, rmax, lmax, etc.
// the orbitals include : NAO, non-local projectors, descriptors, etc.
// mohan note 2021-02-13
//---------------------------------------------------------------------

class LCAO_Orbitals
{
	public:

	LCAO_Orbitals();
	~LCAO_Orbitals();

	void Read_Orbitals(void);

	void Read_PAO(const int& it);

	// in order to get rid of the .NONLOCAL file.
	void Set_NonLocal(const int &it, int &n_projectors);

	// read in the NONLOCAL projector from file.
	void Read_NonLocal(const int& it, int &n_projectors);

	void set_nl_index(void);

	void Read_Descriptor(void);		//caoyu add 2020-3-16

#ifdef __MPI
	void bcast_files(void);
#endif

	const double& get_ecutwfc(void) const {return ecutwfc;}
	const int& get_kmesh(void) const{return kmesh;}
	const double& get_dk(void) const {return dk;}
	const double& get_dR(void) const {return dR;}
	const double& get_Rmax(void) const {return Rmax;}
	const int& get_lmax(void) const {return lmax;}
	const int& get_lmax_d(void) const { return lmax_d; }		//lmax of descriptor basis		//caoyu add 2021-03-17
	const int& get_nchimax(void) const {return nchimax;}
	const int& get_nchimax_d(void) const { return nchimax_d; }	//nchimax of descriptor basis		//caoyu add 2021-03-17
	const int& get_ntype(void) const {return ntype;}
	const double& get_dr_uniform(void) const {return dr_uniform;}

	// numerical atomic orbitals
	Numerical_Orbital* Phi;
	
	// nonlocal projectors (1-dimension array)
	Numerical_Nonlocal* Beta;
	
	//caoyu add 2021-3-10
	// descriptor bases, saved as one-type atom orbital
	Numerical_Orbital* Alpha;

	// initialized in input.cpp
	double ecutwfc;
	double dk;
	double dR;
	double Rmax;
	int *nproj; //mohan add 2010-12-19
	int nprojmax; // mohan add 2010-03-07
	int nkb; // total number of projectors.
	IntArray itiaib2ib_all;
	IntArray ib2_ylm;
	
	double dr_uniform;

	// initalized in unitcell_pseudo
	// assume ntype < 20.
	bool read_in_flag;
	std::vector<string> orbital_file;
	std::vector<string> nonlocal_file;
	string descriptor_file;	//caoyu add 2020-3-16

	private:

	int kmesh;
	int lmax;
	int nchimax;
	int lmax_d;	//caoyu add 2021-03-17
	int nchimax_d;	//caoyu add 2021-03-17
	int ntype;

};

// PLEASE avoid using 'ORB' as global variable 
// mohan note 2021-03-23
extern LCAO_Orbitals ORB;
#endif
