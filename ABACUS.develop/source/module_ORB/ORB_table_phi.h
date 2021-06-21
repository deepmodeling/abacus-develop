#ifndef ORB_TABLE_PHI_H 
#define ORB_TABLE_PHI_H 

#include "ORB_read.h"
#include "ORB_atomic_lm.h"
#include "../src_global/sph_bessel_recursive.h"
#include "../src_global/intarray.h"
#include <set>

class ORB_table_phi
{
	public:

	ORB_table_phi();
	~ORB_table_phi();

	void allocate (
		const int &ntype, ///< number of atom types
		const int &lmax_in,///< max L used to calculate overlap
		const int &kmesh_in,///< kpoints, for integration in k space
		const double &Rmax_in,///< max value of radial table
		const double &dR_in,///< delta R, for making radial table
		const double& dk_in///< delta k, for integration in k space
	);

	void init_Table(
		const int &job,
		LCAO_Orbitals &orb);

	void Destroy_Table(LCAO_Orbitals &orb);

	/**
	 * Five dimension:
	 *-----------------------------
	 * (1) 0: normal (S(R)) ; 1: derivative( dS/dR )
	 *
	 * (2) pairs type number.
	 *
	 * (3) pairs chi
	 *
	 * (4) Max angular momentum: L.
	 *
	 * (5) Distance between atoms: R.
	 */
	double***** Table_SR;
	double***** Table_TR;

	bool destroy_sr;
	bool destroy_tr;
	/**
	 * \brief make table of Spherical bessel
	 *
	 * Sph_Bes : jlx[kmesh][Rmesh][L],
	 * L should be 2*Lmax, which is max L of all type
	 */
	// Peize Lin update 2016-01-26
	void init_Lmax(
		const int orb_num, 
		const int mode, 
		int &Lmax_used, 
		int &Lmax,
		const int &Lmax_exx) const;

	void init_Table_Spherical_Bessel(
		const int orb_num, 
		const int mode, 
		int &Lmax_used, 
		int &Lmax,
		const int &Lmax_exx);

	// Peize Lin add 2017-04-24, and change all jlx in this class
	Sph_Bessel_Recursive::D2* pSB = nullptr;

	///
	/// make the index, in order to get the element from Table_SR and Table_TR quickly.
	///

	///
	/// OV stands for 'overlap'
	///
	/// T stands for atom type.
	///
	void init_OV_Tpair(LCAO_Orbitals& orb);
	///
	/// O stands for orbitals.
	///
	void init_OV_Opair(LCAO_Orbitals& orb);

	int OV_nTpairs;
    IntArray OV_Tpair;
    IntArray OV_Opair;
    IntArray OV_L2plus1;

	static int get_rmesh( const double &R1, const double &R2);

	static double dr;
	int Rmesh;


	void cal_ST_Phi12_R(
		const int &job,
		const int &l,
		const Numerical_Orbital_Lm &n1,
		const Numerical_Orbital_Lm &n2,
		const int &rmesh,
		double *rs,
		double *drs) const;

	// Peize Lin add 2017-10-13
	void cal_ST_Phi12_R(
		const int &job,
		const int &l,
		const Numerical_Orbital_Lm &n1,
		const Numerical_Orbital_Lm &n2,
		const set<size_t> &radials,				// only calculate ir in radials
		double *rs,
		double *drs) const;

	private:

	// variables
    int ntype;

	int lmax;

	double Rmax;

	double dk;

	int nlm;

	int kmesh;

	double *kpoint;

	double *r;

	double *rab;

	double *kab;	
};
#endif
