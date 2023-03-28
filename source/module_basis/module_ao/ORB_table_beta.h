#ifndef ORB_TABLE_BETA_H 
#define ORB_TABLE_BETA_H 

#include "ORB_read.h" // use LCAO_Orbitals
#include "ORB_atomic_lm.h" // use Numerical_Orbital_Lm
#include "module_base/sph_bessel_recursive.h" // use ModuleBase::Sph_Bessel_Recursive
#include "module_base/intarray.h"

class ORB_table_beta
{
	public:

	ORB_table_beta();
	~ORB_table_beta();

	void allocate (
		const int &ntype,
		const int &lmax_in,
		const int &kmesh_in,
		const double &Rmax_in,
		const double &dR_in,
		const double &dk_in);

	double***** Table_NR;
	bool destroy_nr;
	
	///
	/// NL stands for 'nonlocal', 
	/// T stands for atom type.
	///O stands for orbitals.
	///
	void init_NL_Tpair(
		const Numerical_Orbital* phi_,
		const Numerical_Nonlocal* beta_
	);

    void init_NL_Opair(
		LCAO_Orbitals &orb,
		const int nprojmax,
		const int* nproj);

	int NL_nTpairs;
	ModuleBase::IntArray NL_Tpair;
	ModuleBase::IntArray NL_Opair;
	ModuleBase::IntArray NL_L2plus1;

	void init_Table_Beta(
		ModuleBase::Sph_Bessel_Recursive::D2 *pSB,
		const Numerical_Orbital* phi_,
		const Numerical_Nonlocal* beta_,
		const int* nproj_);

	void Destroy_Table_Beta(const int& ntype,
	const Numerical_Orbital* phi_,
	const int* nproj_);

	static int get_rmesh( const double &R1, const double &R2);

	static double dr;

	int Rmesh;

	private:

	void cal_VNL_PhiBeta_R(
		ModuleBase::Sph_Bessel_Recursive::D2 *pSB, // mohan add 2021-03-06
        const int &l,
        const Numerical_Orbital_Lm &n1,
        const Numerical_Nonlocal_Lm &n2,
        const int &rmesh,
        double *rs,
		double *drs);

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
