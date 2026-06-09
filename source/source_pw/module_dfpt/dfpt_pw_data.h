// ============================================================
// This code is added by Mohan Chen on 2026-05-18.
// This code is currently in the design phase and has not been
// put into production yet. It may change in the future.
// Please use this code with caution. Only developers who know
// what they are doing should use this code.
// ============================================================

#ifndef DFPT_PW_DATA_H
#define DFPT_PW_DATA_H

#include "source_base/matrix.h"
#include "source_base/vector3.h"
#include "source_psi/psi.h"
#include "source_cell/qlist.h"
#include <vector>
#include <complex>

namespace ModuleDFPT {

class DFPT_PW_Data {
public:
    DFPT_PW_Data();
    ~DFPT_PW_Data();
    
    void init(ModuleCell::QList* qlist, int nk, int nbands, int npw_max, 
              int nrxx, int nspin, int nat);
    
    void clean();
    
    int get_nq() const;
    ModuleBase::Vector3<double> get_qvec(int q_idx) const;
    int get_nirr(int q_idx) const;
    std::vector<int> get_irrep_modes(int q_idx, int irrep) const;
    
    void set_dpsi(int q_idx, int k_idx, int band_idx, 
                  const std::vector<std::complex<double>>& psi);
    std::vector<std::complex<double>> get_dpsi(int q_idx, int k_idx, int band_idx) const;
    psi::Psi<std::complex<double>>& get_dpsi_obj(int q_idx);
    
    void set_drho_r(int q_idx, int spin, const std::vector<double>& rho);
    std::vector<double> get_drho_r(int q_idx, int spin) const;
    void set_drho_g(int q_idx, int spin, const std::vector<std::complex<double>>& rho);
    std::vector<std::complex<double>> get_drho_g(int q_idx, int spin) const;
    
    void set_dv_r(int q_idx, int spin, const std::vector<double>& v);
    std::vector<double> get_dv_r(int q_idx, int spin) const;
    
    void set_dynmat(int q_idx, const ModuleBase::matrix& dm);
    ModuleBase::matrix get_dynmat(int q_idx) const;
    void set_phon_freq(int q_idx, const std::vector<double>& freq);
    std::vector<double> get_phon_freq(int q_idx) const;
    
    void set_dielectric(const ModuleBase::matrix& eps);
    ModuleBase::matrix get_dielectric() const;
    void set_born(int atom_idx, const ModuleBase::matrix& z);
    ModuleBase::matrix get_born(int atom_idx) const;
    
    void set_compute_q0(bool flag) { compute_q0_ = flag; }
    bool get_compute_q0() const { return compute_q0_; }
    void set_loto(bool flag) { loto_ = flag; }
    bool get_loto() const { return loto_; }
    
    void set_is_metal(bool flag) { is_metal_ = flag; }
    bool get_is_metal() const { return is_metal_; }
    void set_dmu(double dmu) { dmu_ = dmu; }
    double get_dmu() const { return dmu_; }
    
    void set_max_iter(int iter) { max_iter_ = iter; }
    int get_max_iter() const { return max_iter_; }
    void set_conv_thr(double thr) { conv_thr_ = thr; }
    double get_conv_thr() const { return conv_thr_; }
    void set_current_iter(int iter) { current_iter_ = iter; }
    int get_current_iter() const { return current_iter_; }
    void set_converged(bool flag) { converged_ = flag; }
    bool get_converged() const { return converged_; }
    
    void add_residual(double r) { residuals_.push_back(r); }
    std::vector<double> get_residuals() const { return residuals_; }
    
private:
    ModuleCell::QList* qlist_ = nullptr;
    
    int nk_ = 0;
    int nbands_ = 0;
    int npw_max_ = 0;
    int nrxx_ = 0;
    int nspin_ = 1;
    int nat_ = 0;
    
    std::vector<psi::Psi<std::complex<double>>> dpsi_;
    
    std::vector<std::vector<std::vector<double>>> drho_r_;
    std::vector<std::vector<std::vector<std::complex<double>>>> drho_g_;
    
    std::vector<std::vector<std::vector<double>>> dv_r_;
    
    std::vector<ModuleBase::matrix> dynmat_;
    std::vector<std::vector<double>> phon_freq_;
    
    bool compute_q0_ = false;
    bool loto_ = false;
    ModuleBase::matrix dielectric_;
    std::vector<ModuleBase::matrix> born_;
    
    bool is_metal_ = false;
    double dmu_ = 0.0;
    
    int max_iter_ = 100;
    double conv_thr_ = 1e-8;
    int current_iter_ = 0;
    bool converged_ = false;
    std::vector<double> residuals_;
    
    bool is_initialized_ = false;
    
    void allocate_memory();
    void deallocate_memory();
};

} // namespace ModuleDFPT

#endif // DFPT_PW_DATA_H