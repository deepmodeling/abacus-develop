#ifndef INFONONLOCAL_H
#define INFONONLOCAL_H

#include "atom_spec.h"
#include "../src_pw/tools.h"
#include "../module_orbital/ORB_nonlocal.h"
class InfoNonlocal
{
		public:
		InfoNonlocal();
		~InfoNonlocal();
		///
		///NON-LOCAL part for LCAO
		///
		Numerical_Nonlocal* Beta;/// nonlocal projectors (1-dimension array)
		int *nproj; //mohan add 2010-12-19
		int nprojmax; // mohan add 2010-03-07
		double rcutmax_Beta;	//caoyu add 2021-05-24
		const double& get_rcutmax_Beta(void) const { return rcutmax_Beta; }
		/// in order to get rid of the .NONLOCAL file.
		void Set_NonLocal(
			const int &it, 
			Atom* atom, 
			int &n_projectors,
			const int& kmesh,
			const double& dk,
			const double& dr_uniform);
		/// read in the NONLOCAL projector from file.
		void Read_NonLocal(
			const int &it, 
			Atom* atom, 
			int &n_projectors, 
			const int &my_rank,
			const int& kmesh,
			const double& dk,
			const double& dr_uniform,
			const string& nonlocalFile);
};

#endif