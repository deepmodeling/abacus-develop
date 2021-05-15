#ifndef LCAO_DESCRIPTOR_H
#define LCAO_DESCRIPTOR_H

#include "../src_global/intarray.h"
#include "../src_global/complexmatrix.h"

//caoyu add 2021-03-29
class LCAO_Descriptor
{
public:

    LCAO_Descriptor(int lm, int inlm);
    ~LCAO_Descriptor();

	// cal S_alpha_mu: overlap between lcao basis Phi and descriptor basis Al
    void build_S_descriptor(const bool &calc_deri); 

	// cal pdm: S_alpha_mu * inv(Sloc) * DM * inv(Sloc) * S_nu_beta
    void cal_projected_DM(void);

	// cal d: EIGENVALUE of pdm in block of I_n_l
    void cal_descriptor(void);
    void print_descriptor(void);

private:

	// overlap between lcao and descriptor basis
	double** S_mu_alpha;	//[Alpha.totnchi][NLOCAL*(2*lmax+1)]	caoyu modified 2021-05-07

	//d(S) for f_delta
	double** DS_mu_alpha_x;
	double** DS_mu_alpha_y;
	double** DS_mu_alpha_z;
	
	// projected density matrix
	double** pdm;	//[2l+1][2l+1]	caoyu modified 2021-05-07
	
	//gdmx: dD/dX		\sum_{mu,nu} 4*c_mu*c_nu * <dpsi_mu/dx|alpha_m><alpha_m'|psi_nu>
	double** gdmx;	//[2l+1][2l+1]	
	double** gdmy;
	double** gdmz;

	// descriptors
    double *d;

    int n_descriptor;

	// \sum_L{Nchi(L)*(2L+1)}
	int des_per_atom;

	const int lmaxd;

	const int inlmax;

	IntArray* mu_index;
	IntArray* inl_index;	//caoyu add 2021-05-07
	int* inl_ld;	//inl_ld[inl_index] = l_descriptor

    void init_mu_index(void);
    
	void set_S_mu_alpha(
		const int &iw1_all, 
		const int& inl,
		const int& im,
		const double& v);

    void print_projected_DM(
		ofstream &ofs, 
		ComplexMatrix &des, 
		const int &it, 
		const int &ia, 
		const int &l, 
		const int& n);
	
	void set_DS_mu_alpha(
	const int& iw1_all,
	const int& inl,
	const int& im,
	const double& vx,
	const double& vy,
	const double& vz);

	void cal_gdmx();	//dD/dX
	void init_gdmx();
	void del_gdmx();
	
	void getdm(double* dm);

	void cal_v_delta();	//pytorch term remaining!
	void cal_f_delta();	//pytorch term remaining!
	
};


#endif
