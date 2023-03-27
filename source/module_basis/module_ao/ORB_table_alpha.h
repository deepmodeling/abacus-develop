#ifndef ORB_TABLE_ALPHA_H
#define ORB_TABLE_ALPHA_H

#include "ORB_atomic_lm.h"
#include "ORB_read.h" 
#include "module_base/sph_bessel_recursive.h"
#include "module_base/intarray.h"

//caoyu add 2021-03-17

class ORB_table_alpha
{
public:
	ORB_table_alpha();
	~ORB_table_alpha();

	void allocate(
		const int &ntype,
		const int &lmax_in,
		const int &kmesh_in,
		const double &Rmax_in,
		const double &dR_in,
		const double &dk_in);

	/// overlap between lcao basis phi and descriptor basis alpha
	double *****Table_DSR;

	bool table_allocated;

	
	/// O stands for orbitals.
	void init_DS_Opair(LCAO_Orbitals &orb);

	void init_DS_2Lplus1(LCAO_Orbitals &orb);

	ModuleBase::IntArray DS_Opair;

	int *DS_2Lplus1;

	void init_Table_Alpha(ModuleBase::Sph_Bessel_Recursive::D2 *pSB, LCAO_Orbitals &orb);

	void Destroy_Table_Alpha(LCAO_Orbitals &orb);

	int get_rmesh(const double &R1, const double &R2) const;

	double dr;

	int Rmesh;

	int ntype;

	int lmax;

	void print_Table_DSR(LCAO_Orbitals &orb);		//caoyu add 2021-03-20

private:

	void cal_S_PhiAlpha_R(
		ModuleBase::Sph_Bessel_Recursive::D2 *pSB, // mohan add 2021-03-06
		const int &l,
		const Numerical_Orbital_Lm &n1,
		const Numerical_Orbital_Lm &n2,
		const int &rmesh,
		double *rs,
		double *drs);

	// variables
	double Rmax;
	double dk;
	int nlm;
	int kmesh;
	double *kpoint;
	//double *r;
	//double *rab;
	double *kab;

	//------------------------------
	// sizes of table
	//int ntype_; // table is created with this->ntype
	int lmax_d_;
	std::vector<int> lmax_;
	std::vector<int> nchi_pairs_;

	// automatically deallocate Table_DSR using lmax_d_, lmax_ & nchi_pairs_
	// called by destructor
	void _destroy_table();

};
#endif
