#ifndef ATOM_PSEUDO_H
#define ATOM_PSEUDO_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/vector3.h"
#include "../module_io/output.h"
#include "../module_base/complexarray.h"
#include "../module_base/complexmatrix.h"
#include "pseudo_nc.h"


class Atom_pseudo : public pseudo_nc
{
public:

	Atom_pseudo();
	~Atom_pseudo();

	// mohan add 2021-05-07
	ModuleBase::ComplexArray d_so; //(:,:,:), spin-orbit case
	ModuleBase::matrix d_real; //(:,:), non-spin-orbit case
	int nproj;
	int nproj_soc; // dimension of D_ij^so
	int non_zero_count_soc[4];
	int *index1_soc[4];
	int *index2_soc[4];

	void set_d_so( // mohan add 2021-05-07
		ModuleBase::ComplexMatrix &d_so_in,
		const int &nproj_in,
		const int &nproj_in_so,
		const bool has_so);


	inline void get_d(const int& is, const int& p1, const int& p2, const std::complex<double>*& tmp_d)
	{
		tmp_d = &this->d_so(is, p1, p2);
		return;
	}
	inline void get_d(const int& is, const int& p1, const int& p2, const double*& tmp_d)
	{
		tmp_d = &this->d_real(p1, p2);
		return;
	}
	

#ifdef __MPI
	void bcast_atom_pseudo(void); // for upf201
#endif

};

#endif
