#ifndef NUMERICAL_NONLOCAL_LM
#define NUMERICAL_NONLOCAL_LM

#include <string>

/**
 * \class Numerical_Nonlocal_Lm
 * CLASS Numerical_Nonlocal_Lm
 * -------------------------------
 * Note : contain information about each projector
 * all features of projector's shape
 *
 * AUTHOR : liaochen
 *
 * DATE : 2008-03-04
 */

class Numerical_Nonlocal_Lm
{

	public:

	double* beta_uniform;
	double* dbeta_uniform;
	int nr_uniform;
	double dr_uniform;

	public:

	Numerical_Nonlocal_Lm();
	~Numerical_Nonlocal_Lm();	

	const int& getL() const { return this->angular_momentum_l; }
	const int& getType() const { return this->index_atom_type; }
	const double& getRcut() const { return this->rcut; }

    const int& getNr() const { return this->nr; }
	const double* getRadial() const { return this->r_radial; }
	const double& getRadial(const int &ir) const { return this->r_radial[ir]; }
	const double* getBeta_r() const { return this->beta_r; }
	const double& getBeta_r(const int &ir) const { return this->beta_r[ir]; }

	const double& getDk()const { return this->dk; }
	const double* getKpoint()const { return this->k_radial; }
	const double& getKpoint(const int &ik) const { return this->k_radial[ik]; }
	const double* getBeta_k() const { return this->beta_k; }
	const double& getBeta_k(const int &ik) const { return this->beta_k[ik]; }
	
    // enables deep copy
	Numerical_Nonlocal_Lm& operator= (const Numerical_Nonlocal_Lm& nol );

	void set_NL_proj(
 		const std::string &label,
    	const int &index_atom_type_in,
    	const int &angular_momentum_l_in,
    	const int &nr_in,
    	const double *rab_in,
    	const double *r_radial_in,
    	const double *beta_r_in,
    	const int &nk_in,
    	const double &dk_in,
		const double &dr_uniform_in);

	void plot(const int &my_rank)const;

	private:

	void freemem(void);
	void renew(void);
	//void extra_uniform(const double &dr_uniform);
	void get_kradial(void);

	std::string label;
	int index_atom_type;
	int angular_momentum_l;
	int index_proj;

	int nr;
	int nk;

	double rcut;
	double kcut;
	double dk;

	double* r_radial; //points of r
	double* k_radial;

	double* rab;
	double* beta_r; // |beta(r) * r>
	double* beta_k;
};

#endif
