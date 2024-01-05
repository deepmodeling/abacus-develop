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

	// static function to get global instance
	static const LCAO_Orbitals& get_const_instance();

	void Read_Orbitals(
		std::ofstream &ofs_in, // mohan add 2021-05-07
		const int &ntype_in,
		const int &lmax_in,
		const bool &deepks_setorb, //  mohan add 2021-04-25
		const int &out_mat_r, // mohan add 2021-04-26
		const bool &force_flag, // mohan add 2021-05-07
		const int &my_rank); // mohan add 2021-04-26

	void Read_PAO(
		std::ofstream &ofs_in,
		const int& it,
		const bool &force_flag, // mohan add 2021-05-07
		const int& my_rank); // mohan add 2021-04-26



	void Read_Descriptor(
		std::ofstream &ofs_in,
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

	/// numerical atomic orbitals
	Numerical_Orbital* Phi;
	
	
	//caoyu add 2021-3-10
	/// descriptor bases, saved as one-type atom orbital
	Numerical_Orbital* Alpha;

	// initialized in input.cpp
	double ecutwfc;
	double dk;
	double dR;
	double Rmax;
	
	double dr_uniform;

	// initalized in UnitCell
	// assume ntype < 20.
	bool read_in_flag;
	std::vector<std::string> orbital_file;
	std::vector<std::string> nonlocal_file;
	std::string descriptor_file;	//caoyu add 2020-3-16

private:

	int ntype; // number of elements
	int kmesh; // number of points on kmesh

	int lmax;
	int nchimax;

	int lmax_d;	//max l of descriptor orbitals
	int nchimax_d;	//max number of descriptor orbitals per l

	double rcutmax_Phi;	//caoyu add 2021-05-24

	void read_orb_file(
		std::ofstream &ofs_in,
		std::ifstream &ifs, 
		const int &it, 
		int &lmax, 
		int &nchimax, 
		Numerical_Orbital* ao,
		const bool &force_flag, // mohan add 2021-05-07
		const int &my_rank);	//caoyu add 2021-04-26

    friend class TwoCenterBundle; // for the sake of TwoCenterBundle::to_LCAO_Orbitals
};

/// PLEASE avoid using 'ORB' as global variable 
///
///mohan note 2021 - 03 - 23
namespace GlobalC
{
extern LCAO_Orbitals ORB;
}
#endif
