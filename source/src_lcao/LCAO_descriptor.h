#ifndef LCAO_DESCRIPTOR_H
#define LCAO_DESCRIPTOR_H

#include "../module_base/intarray.h"
#include "../module_base/complexmatrix.h"
#include <torch/script.h>
#include "../src_pw/global.h"

//caoyu add 2021-03-29
class LCAO_Descriptor
{
public:

    explicit LCAO_Descriptor();
    ~LCAO_Descriptor();

	void init(int lm, int nm, int tot_inl);	// index of descriptor in all atoms
	// cal S_alpha_mu: overlap between lcao basis Phi and descriptor basis Al
    void build_S_descriptor(const bool &calc_deri);

	// cal pdm: S_alpha_mu * inv(Sloc) * DM * inv(Sloc) * S_nu_beta
    void cal_projected_DM(void);

	// cal d: EIGENVALUE of pdm in block of I_n_l
    void cal_descriptor(void);
	void print_descriptor(void);


	void cal_v_delta(const string& model_file);//<psi|V_delta|psi>
	void cal_f_delta(matrix& dm);	//pytorch term remaining!
	void print_H_V_delta();
	void print_F_delta();

	/*These 3 func save the [dm_eig], [e_base], [f_base]
	of current configuration as .npy file, when deepks_scf = 1.
	After a full group of consfigurations are calculated,
    we need a python script to 'load' and 'torch.cat' these .npy files,
    and get l_e_delta and l_f_delta corresponding to the exact e,f data.*/
	void save_npy_d();
	void save_npy_e(double& ebase);	//Ry
	void save_npy_f(matrix& fbase);//Ry

	//deepks E_delta(Ry)
	double E_delta = 0.0;
	//deepks V_delta, to be added to Hamiltonian matrix
	double* H_V_delta;
	//deepks F_delta(Ry/Bohr), to be added to atom force
	matrix	F_delta;

private:
	torch::jit::script::Module module;

	// overlap between lcao and descriptor basis
	double** S_mu_alpha;	//[tot_Inl][NLOCAL][2l+1]	caoyu modified 2021-05-07

	//d(S) for f_delta:	<\psi_mu|d\alpha^I_nlm> , [tot_Inl][NLOCAL][2l+1]
	double** DS_mu_alpha_x;
	double** DS_mu_alpha_y;
	double** DS_mu_alpha_z;

	// projected density matrix
	double** pdm;	//[tot_Inl][2l+1][2l+1]	caoyu modified 2021-05-07
	std::vector<torch::Tensor> pdm_tensor;

	// descriptors
    double *d;
	vector<torch::Tensor> d_tensor;

	//gedm:dE/dD, [tot_Inl][2l+1][2l+1]	(E: Hartree)
	std::vector<torch::Tensor> gedm_tensor;

	//gdmx: dD/dX		\sum_{mu,nu} 4*c_mu*c_nu * <dpsi_mu/dx|alpha_m><alpha_m'|psi_nu>
	double*** gdmx;	//[natom][tot_Inl][2l+1][2l+1]
	double*** gdmy;
	double*** gdmz;

	//dE/dD, autograd from loaded model(E: Ry)
	double** gedm;	//[tot_Inl][2l+1][2l+1]


	int n_descriptor;

	// \sum_L{Nchi(L)*(2L+1)}
	int des_per_atom;

	int lmaxd = 0;
	int nmaxd = 0;
	int inlmax = 0;

	IntArray* alpha_index;
	IntArray* inl_index;	//caoyu add 2021-05-07
	int* inl_l;	//inl_l[inl_index] = l of descriptor with inl_index

	void init_index(void);	// index of descriptor in all atoms

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

	void init_gdmx();
	void load_model(const string& model_file);
	void cal_gedm();	//need to load model in this step
	void cal_gdmx(matrix& dm);	//dD/dX
	void del_gdmx();

	void getdm(double* dm);

	void cal_descriptor_tensor();

};
extern LCAO_Descriptor ld;

#endif
