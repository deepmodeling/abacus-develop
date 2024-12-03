#ifndef LCAO_DEEPKS_IO_H 
#define LCAO_DEEPKS_IO_H

#ifdef __DEEPKS

#include <string>
#include <vector>
#include "module_base/tool_title.h"
#include "module_base/matrix.h"
#include "module_base/complexmatrix.h"

#include <torch/script.h>
#include <torch/torch.h>

namespace LCAO_deepks_io
{

    /// This file contains subroutines that contains interface with libnpy
    /// since many arrays must be saved in numpy format
    /// It also contains subroutines for printing density matrices
    /// which is used in unit tests

    /// There are 2 subroutines for printing density matrices:
    /// 1. print_dm : for gamma only
    /// 2. print_dm_k : for multi-k

    /// others print quantities in .npy format

    /// 3. save_npy_d : descriptor ->dm_eig.npy
    /// 4. save_npy_gvx : gvx ->grad_vx.npy
    /// 5. save_npy_e : energy
    /// 6. save_npy_f : force
    /// 7. save_npy_s : stress
    /// 8. save_npy_o: orbital
    /// 9. save_npy_orbital_precalc: orbital_precalc -> orbital_precalc.npy
    /// 10. save_npy_h : Hamiltonian
    /// 11. save_npy_v_delta_precalc : v_delta_precalc
    /// 12. save_npy_psialpha : psialpha
    /// 13. save_npy_gevdm : grav_evdm , can use psialpha and gevdm to calculate v_delta_precalc

/// print density matrices
void print_dm(const std::vector<double> &dm,
		const int nlocal,
		const int nrow);

void print_dm_k(const int nks,
		const int nlocal,
		const int nrow,
		const std::vector<std::vector<std::complex<double>>>& dm);

void load_npy_gedm(const int nat,
		const int des_per_atom,
		double** gedm,
		double& e_delta,
		const int rank);

    ///----------------------------------------------------------------------
    /// The following 4 functions save the `[dm_eig], [e_base], [f_base], [grad_vx]`
    /// of current configuration as `.npy` file, when `deepks_scf = 1`.
    /// After a full group of consfigurations are calculated,
    /// we need a python script to `load` and `torch.cat` these `.npy` files,
    /// and get `l_e_delta,npy` and `l_f_delta.npy` corresponding to the exact E, F data.
    ///
    /// Unit of energy: Ry
    ///
    /// Unit of force: Ry/Bohr
    ///----------------------------------------------------------------------

void save_npy_d(const int nat,
		const int des_per_atom,
		const int inlmax,
        const int* inl_l,
		const bool deepks_equiv,
		const std::vector<torch::Tensor> &d_tensor,
		const std::string& out_dir,
		const int rank);

void save_npy_gvx(const int nat,
		const int des_per_atom,
		const torch::Tensor &gvx_tensor,
        const std::string& out_dir,
		const int rank);

void save_npy_gvepsl(const int nat,
		const int des_per_atom,
		const torch::Tensor &gvepsl_tensor,
		const std::string& out_dir,
		const int rank);

void save_npy_e(const double &e,  /**<[in] \f$E_{base}\f$ or \f$E_{tot}\f$, in Ry*/
		const std::string &e_file,
		const int rank);

void save_npy_f(const ModuleBase::matrix &f, /**<[in] \f$F_{base}\f$ or \f$F_{tot}\f$, in Ry/Bohr*/
		const std::string &f_file,
		const int nat,
		const int rank);

void save_npy_s(const ModuleBase::matrix &stress, /**<[in] \f$S_{base}\f$ or \f$S_{tot}\f$, in Ry/Bohr^3*/
		const std::string &s_file,
		const double &omega,
		const int rank);

/// QO added on 2021-12-15
void save_npy_o(const ModuleBase::matrix &bandgap, /**<[in] \f$E_{base}\f$ or \f$E_{tot}\f$, in Ry*/
		const std::string &o_file,
		const int nks,
		const int rank);

void save_npy_orbital_precalc(const int nat, 
		const int nks,
		const int des_per_atom,
		const torch::Tensor& orbital_precalc_tensor,
        const std::string& out_dir,
		const int rank);

/// xinyuan added on 2023-2-20
/// for gamma only
void save_npy_h(const ModuleBase::matrix &hamilt,
		const std::string &h_file,
		const int nlocal,
		const int rank);

/// for multi-k
void save_npy_h(const std::vector<ModuleBase::ComplexMatrix> &hamilt,
        const std::string &h_file,
        const int nlocal,
        const int nks,
        const int rank);

void save_npy_v_delta_precalc(const int nat,
		const int nks,
		const int nlocal,
		const int des_per_atom,
		const torch::Tensor& v_delta_precalc_tensor,
		const std::string& out_dir,
		const int rank);

void save_npy_psialpha(const int nat,
		const int nks,
		const int nlocal,
		const int inlmax,
		const int lmaxd,
		const torch::Tensor &psialpha_tensor,
		const std::string& out_dir,
		const int rank);

void save_npy_gevdm(const int nat,
		const int inlmax,
		const int lmaxd,
		const torch::Tensor& gevdm_tensor,
		const std::string& out_dir,
		const int rank);
};

#endif
#endif
