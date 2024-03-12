//=========================================================
//AUTHOR : liaochen
//DATE : 2008-11-12
//UPDATE : Peize Lin change all pointer to std::vector 2016-05-14
//=========================================================
#ifndef NUMERICAL_ORBITAL_LM_H
#define NUMERICAL_ORBITAL_LM_H

#include <vector>
using std::vector;
#include "module_base/global_function.h"
#include "module_base/global_variable.h"

/**
 * CLASS Num_orbital_lm
 * -------------------------------
 * Note : contain information about each orbital : psi(l,m)
 *  all features of orbital's shape
 */

class Numerical_Orbital_Lm
{
	friend class Numerical_Orbital;

	public:

	std::vector<double> psi_uniform;// mohan add 2009-5-10
	std::vector<double> dpsi_uniform; //liaochen add 2010/5/11
	std::vector<double> ddpsi_uniform; //wenfei add 2022/7/13
	
	int nr_uniform;// mohan add 2009-5-10
	double dr_uniform;// mohan add 2009-5-10
	double zty; ///< the valus of psi at 0.
	
	Numerical_Orbital_Lm();
	~Numerical_Orbital_Lm();

	// Peize Lin add 2017-12-12
	enum class Psi_Type{ Psi, Psif, Psik, Psik2 };

	/// EXPLAIN : set information about Numerical_Orbital_Lm
	void set_orbital_info
	(
 		const std::string &label_in,
	 	const int &index_atom_type_in,
		const int &angular_momentum_l_in,
		const int &index_chi_in,
	    const int &nr_in,
		const double *rab_in,
		const double *r_radial_in,
		const Psi_Type &psi_type,				// Peize Lin add 2017-12-12
		const double *psi_in,
		const int &nk_in,
		const double &dk_in,
		// Peize Lin delete lat0 2016-02-03
		const double &dr_uniform,
		bool flag_plot,						// Peize Lin add flag_plot 2016-08-31
		bool flag_sbpool,					// Peize Lin add flag_sbpool 2017-10-02
		const bool &force_flag // mohan add 2021-05-07
	);

private:

	void copy_parameter(
		const std::string &label_in,
		const int &index_atom_type_in,
		const int &angular_momentum_l_in,
		const int &index_chi_in,
		const int &nr_in,
		const double *rab_in,
		const double *r_radial_in,
		const int &nk_in,
		const double &dk_in,
		const double &dr_uniform_in);

	void cal_kradial(void);
	void cal_kradial_sbpool(void);
	void cal_rradial_sbpool(void);
	//void norm_test()const;
	void plot()const;
	//void use_uniform(const double &dr_uniform_in);
	void extra_uniform(const double &dr_uniform_in, const bool &force_flag);

	std::string label;
	int index_atom_type;
	int angular_momentum_l;
	int index_chi;

	int nr;
	int nk;

	double rcut;
	double kcut;
	double dk;

	std::vector<double> r_radial; ///<points of r
	std::vector<double> k_radial;

	std::vector<double> rab;

	std::vector<double> psi;   ///< psi(r)
	std::vector<double> psir;  ///< psi(r) * r
	std::vector<double> psif;  ///< psi(k)
	std::vector<double> psik;  ///< psi(k) * k
	std::vector<double> psik2; ///< psi(k) * k^2

public:

	const std::string& getLabel() const { return label; }
	const int& getType() const     { return index_atom_type; }
	const int& getL() const        { return angular_momentum_l; }
	const int& getChi() const      { return index_chi; }

	const double* getPsiuniform() const  { return ModuleBase::GlobalFunc::VECTOR_TO_PTR(psi_uniform); }
	const double* getDpsiuniform() const { return ModuleBase::GlobalFunc::VECTOR_TO_PTR(dpsi_uniform); }
	const int& getNruniform() const { return nr_uniform; }
	const double& getDruniform() const { return dr_uniform; }

	const int& getNr() const { return nr; }
	const int& getNk() const { return nk; }

	const double& getRcut() const { return rcut; }
	const double& getKcut() const { return kcut; }
	
	const double* getRadial() const             { return ModuleBase::GlobalFunc::VECTOR_TO_PTR(r_radial); }
	const std::vector<double>& get_r_radial() const  { return r_radial; }
	const double& getRadial(const int ir) const { return r_radial[ir]; }

	const double* getRab() const { return ModuleBase::GlobalFunc::VECTOR_TO_PTR(rab); }
	const std::vector<double>& get_rab() const { return rab; }
	const double& getRab(const int ir) const { return rab[ir]; }
	
	const double& getDk()const { return dk; }
	const double* getKpoint() const { return ModuleBase::GlobalFunc::VECTOR_TO_PTR(k_radial); }
	const double& getKpoint(const int ik) const { return k_radial[ik]; }
	const std::vector<double>& get_k_radial() const { return k_radial; }

	const double* getPsi() const { return ModuleBase::GlobalFunc::VECTOR_TO_PTR(psi);}
	const double& getPsi(const int ir) const { return psi[ir];}	
	const std::vector<double>& get_psi() const { return psi; }
	const double* getPsi_r() const { return ModuleBase::GlobalFunc::VECTOR_TO_PTR(psir); }
	const double& getPsi_r(const int ir) const { return psir[ir]; }

	const double* getPsif() const               { return ModuleBase::GlobalFunc::VECTOR_TO_PTR(psif); }
	const double& getPsif(const int ik) const   { return psif[ik]; }
	const std::vector<double>& get_psif() const      { return psif; }
	const double* getPsi_k() const              { return ModuleBase::GlobalFunc::VECTOR_TO_PTR(psik); }
	const double& getPsi_k(const int ik) const  { return psik[ik]; }
	const std::vector<double>& get_psi_k() const      { return psik; }
	const double* getPsi_k2() const             { return ModuleBase::GlobalFunc::VECTOR_TO_PTR(psik2); }
	const double& getPsi_k2(const int ik) const { return psik2[ik]; }
	const std::vector<double>& get_psi_k2() const      { return psik2; }
};

#endif

