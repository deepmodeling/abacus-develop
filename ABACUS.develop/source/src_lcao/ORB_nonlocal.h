#ifndef NUMERICAL_NONLOCAL_H
#define NUMERICAL_NONLOCAL_H

#include "../src_global/complexarray.h"
#include "../src_global/complexmatrix.h"
#include "ORB_nonlocal_lm.h"
//=========================================================
//CLASS  Numerical_Nonlocal
//Note : contain nonlocal ps(:pseudopotential) information 
//about atoms
//Feature : set and store information about ps infomation 
//			related to atoms
//AUTHOR : liaochen
//DATE : 2008-03-04
//=========================================================

class Numerical_Nonlocal
{
public:

	Numerical_Nonlocal();
	~Numerical_Nonlocal();

	// mohan add 2011-03-07
	const int& getL_Beta(const int &nb){ return this->LfromBeta[nb]; }

	const int& getLmax() const { return this->lmax; }

   	const int& getType() const { return this->type; }

	const string& getLabel() const { return this->label; }

	const string& getType_ps() const { return this->type_ps; }

	const complex<double>& getCoefficient_D_so(const int& is, const int& L1, const int& L2) 
						const { return this->Coefficient_D_so(is, L1, L2); }

	const int& get_index1_soc(const int& is, const int& no)
						const { return this->index1_soc[is][no]; }

	const int& get_index2_soc(const int& is, const int& no)
						const { return this->index2_soc[is][no]; }

	const int& get_count_soc(const int& is)
						const { return this->non_zero_count_soc[is]; }

	void set_type_info(
			const int& type_in,
			const string& label_in,
			const string& type_ps_in,
			const int& lmax_in,
			ComplexMatrix& Coefficient_D_in_so,
			const int& nproj_in,
			const int& nproj_in_so,
			int* LfromBeta_in,
			const Numerical_Nonlocal_Lm* ps_orbital_in,
			const bool has_so);

	Numerical_Nonlocal_Lm* Proj; // length: nproj(only store radial function )

	const double& get_rcut_max(void) const { return rcut_max; }

	const int& get_nproj_soc(void) const {return nproj_soc;}

	private:
	
	string label; // atom type

	int type; // number of atom types

	int lmax; // max value of L angular momentum

	double rcut_max;

	string type_ps; //local or nonlocal

	int nproj;

	int *LfromBeta;

	// spin-orbital related
	int nproj_soc; //dimension of D_ij^so

	ComplexArray Coefficient_D_so;   //(:,:,:),  spin-orbit case,  added by zhengdy-soc

	int non_zero_count_soc[4];

	int *index1_soc[4];
 
	int *index2_soc[4];
};

#endif
