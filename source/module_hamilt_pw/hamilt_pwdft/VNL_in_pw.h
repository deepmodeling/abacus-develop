#ifndef VNL_IN_PW_H
#define VNL_IN_PW_H

#include "VL_in_pw.h"
#include "module_base/complexarray.h"
#include "module_base/complexmatrix.h"
#include "module_base/intarray.h"
#include "module_base/realarray.h"
#include "module_cell/unitcell.h"
#include "module_hamilt_pw/hamilt_pwdft/soc.h"
#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"
#include "module_psi/psi.h"
#ifdef __LCAO
#include "module_basis/module_ao/ORB_gen_tables.h"
#endif

//==========================================================
// Calculate the non-local pseudopotential in reciprocal
// space using plane wave as basis set.
//==========================================================
class pseudopot_cell_vnl: public pseudopot_cell_vl
{

public:

	pseudopot_cell_vnl();
    ~pseudopot_cell_vnl();
    void init(const int ntype,
              Structure_Factor* psf_in,
              const ModulePW::PW_Basis_K* wfc_basis = nullptr,
              const bool allocate_vkb = 1);

    double cell_factor; //LiuXh add 20180619

	int nkb; // total number of beta functions considering all atoms

	int lmaxkb; // max angular momentum for non-local projectors

    void init_vnl(UnitCell& cell, const ModulePW::PW_Basis* rho_basis);

    template <typename FPTYPE, typename Device>
    void getvnl(Device * ctx, const int &ik, std::complex<FPTYPE>* vkb_in)const;

    void getvnl(const int &ik, ModuleBase::ComplexMatrix& vkb_in)const;

	// void getvnl_alpha(const int &ik);

	void init_vnl_alpha(void);

	void initgradq_vnl(const UnitCell &cell);

	void getgradq_vnl(const int ik);

    //===============================================================
    // MEMBER VARIABLES :
    // NAME : nqx(number of interpolation points)
    // NAME : nqxq(size of interpolation table)
    // NAME : nhm(max number of different beta functions per atom)
    // NAME : lmaxq
    // NAME : dq(space between points in the pseudopotential tab)
    //===============================================================
    // private:

    int nhm;
    int nbetam; // max number of beta functions

    int lmaxq;

	ModuleBase::matrix indv;		// indes linking  atomic beta's to beta's in the solid
	ModuleBase::matrix nhtol;      	// correspondence n <-> angular momentum l
	ModuleBase::matrix nhtolm;     	// correspondence n <-> combined lm index for (l,m)
	ModuleBase::matrix nhtoj;		// new added

	ModuleBase::realArray dvan;		//(:,:,:),  the D functions of the solid
	ModuleBase::ComplexArray dvan_so;	//(:,:,:),  spin-orbit case,  added by zhengdy-soc

	ModuleBase::realArray tab;		//(:,:,:), interpolation table for PPs: 4pi/sqrt(V) * \int betar(r)jl(qr)rdr
	ModuleBase::realArray tab_alpha;
	ModuleBase::realArray tab_at;	//(:,:,:), interpolation table for atomic wfc
	ModuleBase::realArray tab_dq;   //4pi/sqrt(V) * \int betar(r)*djl(qr)/d(qr)*r^2 dr

	ModuleBase::realArray deeq;		//(:,:,:,:), the integral of V_eff and Q_{nm}
	bool multi_proj = false;
    float *s_deeq = nullptr;
    double *d_deeq = nullptr;
	ModuleBase::ComplexArray deeq_nc;	//(:,:,:,:), the spin-orbit case
    std::complex<float> *c_deeq_nc = nullptr; // GPU array of deeq_nc
    std::complex<double>* z_deeq_nc = nullptr; // GPU array of deeq_nc

    // liuyu add 2023-10-03
    // uspp
    int* indv_ijkb0 = nullptr;      // first beta (index in the solid) for each atom
    ModuleBase::IntArray ijtoh;     // correspondence beta indexes ih,jh -> composite index ijh
    ModuleBase::realArray qq_at;    // the integral of q functions in the solid (ONE PER ATOM)
    ModuleBase::realArray qq_nt;    // the integral of q functions in the solid (ONE PER NTYP) used to be the qq array
    ModuleBase::ComplexArray qq_so; // Q_{nm} for spin-orbit case
    ModuleBase::realArray ap;       // the expansion coefficients
    ModuleBase::IntArray lpx;       // for each input limi,ljmj is the number of LM in the sum
    ModuleBase::IntArray lpl;       // for each input limi,ljmj points to the allowed LM
    ModuleBase::realArray qrad;     // radial FT of Q functions

    float* s_qq_nt = nullptr;
    double* d_qq_nt = nullptr;
    std::complex<float>* c_qq_so = nullptr;  // GPU array of qq_so
    std::complex<double>* z_qq_so = nullptr; // GPU array of qq_so

    mutable ModuleBase::ComplexMatrix vkb;    // all beta functions in reciprocal space
    mutable ModuleBase::ComplexArray gradvkb; // gradient of beta functions
    std::complex<double>*** vkb1_alpha;
    std::complex<double>*** vkb_alpha;
    Structure_Factor* psf = nullptr;

    // other variables
    std::complex<double> Cal_C(int alpha, int lu, int mu, int L, int M);

    double CG(int l1, int m1, int l2, int m2, int L, int M);

    void print_vnl(std::ofstream& ofs);

    /**
     * @brief Compute the radial Fourier transform of the Q functions
     *
     * The interpolation table for the radial Fourier transform is stored in qrad.
     *
     * The formula implemented here is:
     *   \[ q(g,i,j) = \sum_\text{lm} (-i)^l \text{ap}(\text{lm},i,j)
     *   \text{yr}_\text{lm}(g) \text{qrad}(g,l,i,j) \]
     *
     * @param ng [in] the number of G vectors
     * @param ih [in] the first index of Q
     * @param jh [in] the second index of Q
     * @param itype [in] the atomic type
     * @param qnorm [in] the norm of q+g vectors
     * @param ylm [in] the real spherical harmonics
     * @param qg [out] the Fourier transform of interest
     */
    void radial_fft_q(const int ng,
                      const int ih,
                      const int jh,
                      const int itype,
                      const double* qnorm,
                      const ModuleBase::matrix ylm,
                      std::complex<double>* qg) const;
    template <typename FPTYPE, typename Device>
    void radial_fft_q(Device* ctx,
                      const int ng,
                      const int ih,
                      const int jh,
                      const int itype,
                      const FPTYPE* qnorm,
                      const FPTYPE* ylm,
                      std::complex<FPTYPE>* qg) const;

    /**
     * @brief calculate the effective coefficient matrix for non-local pseudopotential projectors
     *
     * @param veff effective potential
     * @param rho_basis potential FFT grids
     * @param cell UnitCell
     */
    void cal_effective_D(const ModuleBase::matrix& veff, const ModulePW::PW_Basis* rho_basis, UnitCell& cell);
#ifdef __LCAO
    ORB_gaunt_table MGT;
#endif

    template <typename FPTYPE>
    FPTYPE* get_nhtol_data() const;
    template <typename FPTYPE>
    FPTYPE* get_nhtolm_data() const;
    template <typename FPTYPE>
    FPTYPE* get_indv_data() const;
    template <typename FPTYPE>
    FPTYPE* get_tab_data() const;
    template <typename FPTYPE>
    FPTYPE* get_deeq_data() const;
    template <typename FPTYPE>
    FPTYPE* get_qq_nt_data() const;
    template <typename FPTYPE>
    std::complex<FPTYPE>* get_qq_so_data() const;
    template <typename FPTYPE>
    std::complex<FPTYPE>* get_vkb_data() const;
    template <typename FPTYPE>
    std::complex<FPTYPE>* get_deeq_nc_data() const;

  private:
    float *s_nhtol = nullptr, *s_nhtolm = nullptr, *s_indv = nullptr, *s_tab = nullptr;
    std::complex<float>* c_vkb = nullptr;

    double *d_nhtol = nullptr, *d_nhtolm = nullptr, *d_indv = nullptr, *d_tab = nullptr;
    std::complex<double>* z_vkb = nullptr;

    const ModulePW::PW_Basis_K* wfcpw = nullptr;

    Soc soc;

    /**
     * @brief Compute interpolation table qrad
     *
     * Compute interpolation table qrad(i,nm,l,nt) = Q^{(L)}_{nm,nt}(q_i)
     * of angular momentum L, for atom of type nt, on grid q_i, where
     * nm = combined index for n,m=1,nh(nt)
     *
     * @param cell UnitCell
     */
    void compute_qrad(UnitCell& cell);

    /**
     * @brief computes the integral of the effective potential with the Q function
     *
     * @param veff effective potential
     * @param rho_basis potential FFT grids
     * @param cell UnitCell
     */
    void newq(const ModuleBase::matrix& veff, const ModulePW::PW_Basis* rho_basis, UnitCell& cell);

    /**
     * @brief calculate D functions in the soc case when tvanp is true
     *
     * @param iat the index of atom
     * @param cell UnitCell
     */
    void newd_so(const int& iat, UnitCell& cell);

    /**
     * @brief calculate D functions in the noncolin case when tvanp is true
     *
     * @param iat the index of atom
     * @param cell UnitCell
     */
    void newd_nc(const int& iat, UnitCell& cell);
};

#endif // VNL_IN_PW
