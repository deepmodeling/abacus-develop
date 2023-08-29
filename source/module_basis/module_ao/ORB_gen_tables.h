#ifndef ORB_GEN_TABLES_H
#define ORB_GEN_TABLES_H

// temporary switch
// if USE_NEW_TWO_CENTER is defined, the new two-center method is used
// to replace the old UOT.snap_xxx
// otherwise the new method is calculated but merely used to compare with the old method

// if __EXX is defined, the new two-center method is off anyway.
#ifndef __EXX
#define USE_NEW_TWO_CENTER
#endif

#include <memory>
#include "ORB_gaunt_table.h"
#include "ORB_table_beta.h"
#include "ORB_table_phi.h"
#include "ORB_table_alpha.h"		//caoyu add 2020-3-18
#include "ORB_read.h"
#include "module_base/vector3.h"
#include "module_base/matrix.h"
#include "module_cell/setup_nonlocal.h"
#include "module_base/intarray.h"
#include "module_base/complexarray.h"
#include "module_basis/module_nao/two_center_bundle.h"

/// used to be 'Use_Overlap_Table',
/// now the name is 'ORB_gen_tables'
class ORB_gen_tables
{
	public:

	friend class ORB_control;

	ORB_gen_tables();
	~ORB_gen_tables();

	// static function to get global instance
	static const ORB_gen_tables& get_const_instance();

	void gen_tables(
		std::ofstream &ofs_in, // mohan add 2021-05-07
		LCAO_Orbitals &orb,
		const int &Lmax_exx,
		const bool& deepks_setorb,///<[in] whether to generate descriptors
		const int &nprojmax,
		const int* nproj,
		const Numerical_Nonlocal* beta_
	);
	void set_unit(const double& v) { lat0 = v; }

	void snap_psipsi(
		const LCAO_Orbitals &orb,
		double olm[],
		const int &job, ///<[in]0 for matrix element of either S or T, 1 for its derivatives
	    const char &dtype, ///<[in] derivative type, 'S' for overlap, 'T' for kinetic energy, 'D' for descriptor in deepks
		const ModuleBase::Vector3<double> &R1,
    	const int &I1,
    	const int &l1,
    	const int &m1,
    	const int &n1,
    	const ModuleBase::Vector3<double> &R2,
    	const int &I2,
    	const int &l2,
    	const int &m2,
		const int &n2,
		bool cal_syns = false,
		double dmax = 0.0)const;

	void snap_psibeta_half(
		const LCAO_Orbitals &orb,
		const InfoNonlocal &infoNL_,
		std::vector<std::vector<double>> &nlm,
		const ModuleBase::Vector3<double> &R1,
		const int &T1,
		const int &L1,
		const int &m1,
		const int &N1,
		const ModuleBase::Vector3<double> &R0, // The projector.
		const int &T0,
		const bool &calc_deri)const; // mohan add 2021-04-25);
	/// set as public because in hamilt_linear,
#ifdef __DEEPKS
    void snap_psialpha_half(
			const LCAO_Orbitals& orb,
			std::vector<std::vector<double>> &nlm,
			const int& job,
			const ModuleBase::Vector3<double>& R1,
			const int& T1,
			const int& L1,
			const int& m1,
			const int& N1,
			const ModuleBase::Vector3<double>& R0, // The projector.
			const int& T0,
			const int& I0
        ) const;
#endif
    /// set as public because in hamilt_linear,
	/// we need to destroy the tables: SR,TR,NR
	/// after ionic optimization is done.
	ORB_table_phi MOT;
	ORB_table_beta tbeta;

	/// if we want to add table for descriptors,
	/// we should consider here -- mohan 2021-02-09
	ORB_table_alpha talpha;		//caoyu add 2021-03-17

    /// a temporary solution in refactoring
    std::unique_ptr<TwoCenterBundle> two_center_bundle;

	private:

	ORB_gaunt_table MGT;

	double get_distance(const ModuleBase::Vector3<double> &R1, const ModuleBase::Vector3<double> &R2)const;

	double lat0;

};

/// PLEASE try to get rid of GlobalC::UOT, which is a global variable
/// mohan add 2021-03-30
namespace GlobalC
{
extern ORB_gen_tables UOT;
}

#endif
