#ifndef LCAO_MATRIX_H
#define LCAO_MATRIX_H

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/vector3.h"
#include "module_base/complexmatrix.h"
#include "module_basis/module_ao/parallel_orbitals.h"

// add by jingan for map<> in 2021-12-2, will be deleted in the future
#include "module_base/abfs-vector3_order.h"
#ifdef __EXX
#include <RI/global/Tensor.h>
#endif

class LCAO_Matrix
{
    public:

    LCAO_Matrix();
    ~LCAO_Matrix();

    void divide_HS_in_frag(const bool isGamma, Parallel_Orbitals& pv, const int& nks);
    
    // folding the fixed Hamiltonian (T+Vnl) if
	// k-point algorithm is used.
	void folding_fixedH(const int &ik, 
                        const std::vector<ModuleBase::Vector3<double>>& kvec_d, 
                        bool cal_syns = false);

    Parallel_Orbitals *ParaV;
    
#ifdef __EXX
    using TAC = std::pair<int, std::array<int, 3>>;
    std::vector< std::map<int, std::map<TAC, RI::Tensor<double>>>> *Hexxd;
    std::vector< std::map<int, std::map<TAC, RI::Tensor<std::complex<double>>>>>* Hexxc;
    /// @brief Hexxk for all k-points, only for the 1st scf loop ofrestart load
    std::vector<std::vector<double>> Hexxd_k_load;
    std::vector<std::vector<std::complex<double>>> Hexxc_k_load;
#endif

    void allocate_HS_k(const long &nloc);

private:

    void allocate_HS_gamma(const long &nloc);


    public:
    //------------------------------
    // H, S, Hfixed
    // used in gamma only algorithm.
    // thse matrix are used to
    // diagonalize.
    //------------------------------
    std::vector<double> Hloc;
    std::vector<double> Sloc;
    std::vector<double> Hloc_fixed;

    //------------------------------
    // 1. Hamiltonian(vl),
    // 2. overlap matrix Sloc2
    // 3. fixed (vna+T+Vnl) matrix.
    // used in kpoint algorithm.
    // these matrix are used to
    // diagonalize.
    //------------------------------
    std::vector<std::complex<double>> Hloc2;
    std::vector<std::complex<double>> Sloc2;
    std::vector<std::complex<double>> Hloc_fixed2;
    //with soc, zhengdy-soc
/*	ModuleBase::ComplexMatrix Hloc2_soc;
    ModuleBase::ComplexMatrix Sloc2_soc;
    ModuleBase::ComplexMatrix Hloc_fixed2_soc;*/


    //------------------------------
    // Store H(mu,nu')
    // nu' : nu in near unitcell R.
    // used in kpoint algorithm.
    // these matrixed are used
    // for 'folding_matrix' in lcao_nnr,
    // HlocR -> Hloc2,
    // SlocR -> Sloc2,
    //------------------------------
    std::vector<double> SlocR;
    std::vector<double> Hloc_fixedR;

    //with soc, zhengdy-soc
    std::vector<std::complex<double>> SlocR_soc;
    std::vector<std::complex<double>> Hloc_fixedR_soc;

    //LiuXh add 2019-07-15
    double ****Hloc_fixedR_tr;
    double ****SlocR_tr;
    double ****HR_tr;


    std::complex<double> ****Hloc_fixedR_tr_soc;
    std::complex<double> ****SlocR_tr_soc;
    std::complex<double> ****HR_tr_soc;

    // jingan add 2021-6-4, modify 2021-12-2
    // Sparse form of HR and SR, the format is [R_direct_coor][orbit_row][orbit_col]
    
    // For HR_sparse[2], when nspin=1, only 0 is valid, when nspin=2, 0 means spin up, 1 means spin down
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> HR_sparse[2];
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> SR_sparse;
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> TR_sparse;

    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> dHRx_sparse[2];
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> dHRy_sparse[2];
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> dHRz_sparse[2];

    // For nspin = 4
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>> HR_soc_sparse;
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>> SR_soc_sparse;
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>> TR_soc_sparse;

    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>> dHRx_soc_sparse;
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>> dHRy_soc_sparse;
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>> dHRz_soc_sparse;

    // Record all R direct coordinate information, even if HR or SR is a zero matrix
    std::set<Abfs::Vector3_Order<int>> all_R_coor;

    // Records the R direct coordinates of HR and SR output, This variable will be filled with data when HR and SR files are output.
    std::set<Abfs::Vector3_Order<int>> output_R_coor;


    //========================================
    // FORCE
    //========================================

    //-----------------------------------------
    // force in LCAO
    // used in gamma only algorithm.
    //-----------------------------------------
    double* DSloc_x;
    double* DSloc_y;
    double* DSloc_z;

    //-----------------------------------------
    // force in LCAO
    // used in k-points algorithm.
    //-----------------------------------------
    double* DSloc_Rx;
    double* DSloc_Ry;
    double* DSloc_Rz;

    //-----------------------------------------
    // dT + part of dVNL
    // used in gamma only algorithm.
    //-----------------------------------------
    double* DHloc_fixed_x;
    double* DHloc_fixed_y;
    double* DHloc_fixed_z;

    //-----------------------------------------
    // dT + part of dVNL
    // used in kpoint algorithm.
    //-----------------------------------------
    double* DHloc_fixedR_x;
    double* DHloc_fixedR_y;
    double* DHloc_fixedR_z;

    //----------------------------------------
    // r_mu - r_nu
    //----------------------------------------

    double* DH_r;//zhengdy added 2017-07


    double* stvnl11;
    double* stvnl12;
    double* stvnl13;
    double* stvnl22;
    double* stvnl23;
    double* stvnl33;

    double* DSloc_11;
    double* DSloc_12;
    double* DSloc_13;
    double* DSloc_22;
    double* DSloc_23;
    double* DSloc_33;

    double* DHloc_fixed_11;
    double* DHloc_fixed_12;
    double* DHloc_fixed_13;
    double* DHloc_fixed_22;
    double* DHloc_fixed_23;
    double* DHloc_fixed_33;

	template <typename T>
		static void set_mat2d(
				const int& global_ir, 
				const int& global_ic, 
				const T& v, 
				const Parallel_Orbitals& pv, 
				T* mat);

	void set_HSgamma(
			const int& iw1_all, 
			const int& iw2_all, 
			const double& v, 
			double* HSloc);

	void set_HSk(
			const int &iw1_all, 
			const int &iw2_all, 
			const std::complex<double> &v, 
			const char &dtype, 
			const int spin = 0);

	void set_force (
			const int& iw1_all, 
			const int& iw2_all, 
			const double& vx, 
			const double& vy, 
			const double& vz, 
			const char &dtype);

	void set_stress (
			const int& iw1_all, 
			const int& iw2_all, 
			const double& vx, 
			const double& vy,
			const double& vz, 
			const char &dtype, 
			const ModuleBase::Vector3<double> &dtau);

	void set_HR_tr(
			const int &Rx, 
			const int &Ry, 
			const int &Rz, 
			const int &iw1_all, 
			const int &iw2_all, 
			const double &v);

	void set_HR_tr_soc(
			const int &Rx, 
			const int &Ry, 
			const int &Rz, 
			const int &iw1_all, 
			const int &iw2_all, 
			const std::complex<double> &v); //LiuXh add 2019-07-16

    void zeros_HSgamma(const char &mtype);

    void zeros_HSk(const char &mtype);

    void zeros_HSR(const char &mtype);

    void print_HSgamma(const char &mtype, std::ostream &os=std::cout);

	void print_HSk(
			const char &mtype, 
			const char &vtype = 'C', 
			const double &accuracy = 1.0e-5, 
			std::ostream &os=std::cout);

    void update_Hloc(void);

    void update_Hloc2(const int &ik);

    void allocate_HS_R(const int &nnr);

    void output_HSk(const char &mtype, std::string &fn);

    //LiuXh add 2019-07-15
    void allocate_Hloc_fixedR_tr(void);

    void allocate_HR_tr(void);

    void allocate_SlocR_tr(void);

    void destroy_Hloc_fixedR_tr(void);

    // jingan add 2021-6-4, modify 2021-12-2
    void destroy_HS_R_sparse(void);

    void destroy_T_R_sparse(void);

    void destroy_dH_R_sparse(void);

};

#include "LCAO_matrix.hpp"

#endif
