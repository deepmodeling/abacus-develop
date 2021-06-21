#ifndef ORB_GEN_TABLES_H
#define ORB_GEN_TABLES_H

#include "ORB_gaunt_table.h"
#include "ORB_table_beta.h"
#include "ORB_table_phi.h"
#include "ORB_table_alpha.h"		//caoyu add 2020-3-18
#include "ORB_read.h"
#include "../src_global/vector3.h"
#include "../src_global/matrix.h"

/// used to be 'Use_Overlap_Table',
/// now the name is 'ORB_gen_tables'
class ORB_gen_tables
{
	public:

	friend class ORB_control;
	
	ORB_gen_tables();
	~ORB_gen_tables();

	void gen_tables( 
		ofstream &ofs_in, // mohan add 2021-05-07
		const int &job0, 
		LCAO_Orbitals &orb,
		const int &Lmax_exx,
		const int& out_descriptor///<[in] whether to generate descriptors
	);
	void set_unit(const double& v) { lat0 = v; }
	
	void snap_psipsi(
		double olm[],
		const int &job, ///<[in]0 for matrix element of either S or T, 1 for its derivatives
	    const char &dtype, ///<[in] derivative type, 'S' for overlap, 'T' for kinetic energy, 'D' for descriptor in deepks
		const Vector3<double> &R1,
    	const int &I1,
    	const int &l1,
    	const int &m1,
    	const int &n1,
    	const Vector3<double> &R2,
    	const int &I2,
    	const int &l2,
    	const int &m2,
		const int &n2,
		const int &nspin,
		complex<double> *olm1=NULL)const;
		

	void snap_psibeta(
		double nlm[],
		const int& job/**<[in]	job = 0 for vnl matrix elements, job = 1 for its derivatives*/,
		const Vector3<double> &R1,
		const int &I1,
		const int &l1,
		const int &m1,
		const int &n1,
		const Vector3<double> &R2,
		const int &I2,
		const int &l2,
		const int &m2,
		const int &n2,
		const Vector3<double> &Rnl,
		const int &type,
		const matrix &dion, // mohan add 2021-04-25
		const int &nspin, // mohan add 2021-05-07
		const ComplexArray &d_so, // mohan add 2021-04-25
		const int &count_soc, // mohan add 2021-05-07
		int* index1_soc, // mohan add 2021-05-07
		int* index2_soc, // mohan add 2021-05-07
		const int &nproj_in, // mohan add 2021-05-07
		complex<double> *nlm1=NULL,
		const int is=0)const;

	/// set as public because in hamilt_linear, 
	/// we need to destroy the tables: SR,TR,NR
	/// after ionic optimization is done.
	ORB_table_phi MOT;
	ORB_table_beta tbeta;

	/// if we want to add table for descriptors,
	/// we should consider here -- mohan 2021-02-09
	ORB_table_alpha talpha;		//caoyu add 2021-03-17

	private:

	ORB_gaunt_table MGT;

	double get_distance(const Vector3<double> &R1, const Vector3<double> &R2)const;

	double lat0;

};

/// PLEASE try to get rid of UOT, which is a global variable
/// mohan add 2021-03-30
extern ORB_gen_tables UOT;

#endif
