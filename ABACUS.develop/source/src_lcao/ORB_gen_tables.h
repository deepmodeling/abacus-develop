//=========================================================
//AUTHOR : Mohan 
//DATE : 2009-04-22
//=========================================================
#ifndef USE_OVERLAP_TABLE_H
#define USE_OVERLAP_TABLE_H

#include "src_pw/tools.h"
#include "ORB_gaunt_table.h"
#include "src_global/ylm.h"
#include "ORB_table_beta.h"
#include "ORB_table_phi.h"
#include "ORB_table_alpha.h"		//caoyu add 2020-3-18

//------------------------------------
// used to be 'Use_Overlap_Table',
// now the name is 'ORB_gen_tables'
//------------------------------------
class ORB_gen_tables
{
	public:

	friend class ORB_control;
	
	ORB_gen_tables();
	~ORB_gen_tables();

	void gen_tables( const int &job0 );
	void set_unit( const double &v ){lat0=v;}
	
	void snap_psipsi(
		double olm[],
		const int &job, // 0 for matrix element of either S or T, 1 for its derivatives
	    const char &dtype, // derivative type, 'S' for overlap, 'T' for kinetic energy
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
		complex<double> *olm1=NULL)const;
		
	//job = 0 for vnl matrix elements
	//job = 1 for its derivatives
	void snap_psibeta(
		double nlm[],
		const int& job,
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
		complex<double> *nlm1=NULL,
		const int is=0)const;

	//caoyu add 2021-03-17
	//job = 0 for vnl matrix elements
	//job = 1 for its derivatives
	void snap_psialpha(
		double olm[],
		const int& job,
		const Vector3<double>& R1,
		const int& I1,
		const int& l1,
		const int& m1,
		const int& n1,
		const Vector3<double>& R2,
		const int& I2,
		const int& l2,
		const int& m2,
		const int& n2,
		complex<double>* olm1 = NULL,
		const int is = 0)const;

	// set as public because in hamilt_linear, 
	// we need to destroy the tables: SR,TR,NR
	// after ionic optimization is done.
	ORB_table_phi MOT;
	ORB_table_beta tbeta;

	// if we want to add table for descriptors,
	// we should consider here -- mohan 2021-02-09
	ORB_table_alpha talpha;		//caoyu add 2021-03-17

	private:

	ORB_gaunt_table MGT;

	double get_distance(const Vector3<double> &R1, const Vector3<double> &R2)const;

	double lat0;

};

// PLEASE try to get rid of UOT, which is a global variable
// mohan add 2021-03-30
extern ORB_gen_tables UOT;

#endif
