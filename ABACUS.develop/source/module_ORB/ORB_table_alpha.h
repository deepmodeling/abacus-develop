#ifndef ORB_TABLE_ALPHA_H
#define ORB_TABLE_ALPHA_H

#include "ORB_atomic_lm.h"
#include "../src_global/sph_bessel_recursive.h"

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

	bool destroy_nr;

	
	/// O stands for orbitals.
	void init_DS_Opair(void);

	void init_DS_2Lplus1(void);

	IntArray DS_Opair;

	int *DS_2Lplus1;

	void init_Table_Alpha(Sph_Bessel_Recursive::D2 *pSB);

	void Destroy_Table_Alpha(void);

	static int get_rmesh(const double &R1, const double &R2);

	static double dr;

	int Rmesh;

	int ntype;

	int lmax;

	//void print_Table_DSR(void);		//caoyu add 2021-03-20

private:

	void cal_S_PhiAlpha_R(
		Sph_Bessel_Recursive::D2 *pSB, // mohan add 2021-03-06
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
	double *r;
	double *rab;
	double *kab;
};
#endif
