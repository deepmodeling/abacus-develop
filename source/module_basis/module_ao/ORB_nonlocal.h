#ifndef NUMERICAL_NONLOCAL_H
#define NUMERICAL_NONLOCAL_H

#include "module_base/complexarray.h"
#include "module_base/complexmatrix.h"
#include "ORB_nonlocal_lm.h"
/**
 * \class Numerical_Nonlocal
 *CLASS  Numerical_Nonlocal
 *----------------------------
 * Note : contain nonlocal ps(:pseudopotential) information 
 * about atoms
 *
 * Feature : set and store information about ps infomation
 * related to atoms
 *
 * AUTHOR : liaochen
 *
 * DATE : 2008-03-04
 */
class Numerical_Nonlocal
{
public:

	Numerical_Nonlocal();
	~Numerical_Nonlocal();

	const int& getLmax() const { return this->lmax; }

   	const int& getType() const { return this->type; }

	const std::string& getLabel() const { return this->label; }

	const std::string& getType_ps() const { return this->type_ps; }


	void set_type_info(
			const int& type_in,
			const std::string& label_in,
			const std::string& type_ps_in,
			const int& lmax_in,
			const int& nproj_in,
			const Numerical_Nonlocal_Lm* ps_orbital_in);

	Numerical_Nonlocal_Lm* Proj; ///< length: nproj(only store radial function )

	const double& get_rcut_max(void) const { return rcut_max; }
    const int& get_nproj() const { return nproj; }

	private:
	
	std::string label; /// <element label

	int type; ///< element index

	int lmax; ///< max value of L angular momentum

	double rcut_max;

	std::string type_ps; ///<local or nonlocal

	int nproj;

};

#endif
