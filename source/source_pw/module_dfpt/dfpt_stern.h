// ============================================================
// This code is added by Mohan Chen on 2026-05-18.
// This code is currently in the design phase and has not been
// put into production yet. It may change in the future.
// Please use this code with caution. Only developers who know
// what they are doing should use this code.
// ============================================================

#ifndef DFPT_STERN_H
#define DFPT_STERN_H

#include "dfpt_pw_data.h"
#include "source_psi/psi.h"

namespace ModuleDFPT {

class DFPT_Stern {
public:
    DFPT_Stern();
    ~DFPT_Stern();
    
    void init(int nk, int nbands, int npw_max, const ModuleBase::matrix& eig, 
              const ModuleBase::matrix& wg, double alpha);
    
    void solve(const psi::Psi<std::complex<double>>& psi, 
               const std::vector<std::complex<double>>& dv_psi,
               int q_idx, int k_idx, int band_idx, double omega, 
               DFPT_PW_Data& data);
    
private:
    int nk_ = 0;
    int nbands_ = 0;
    int npw_max_ = 0;
    double alpha_ = 1.0;
    
    ModuleBase::matrix eig_;
    ModuleBase::matrix wg_;
    
    void apply_pv(const std::vector<std::complex<double>>& x, 
                  std::vector<std::complex<double>>& px);
    
    void apply_op(const std::vector<std::complex<double>>& x,
                  std::vector<std::complex<double>>& y,
                  int k_idx, int band_idx);
    
    void cg_solve(const std::vector<std::complex<double>>& b,
                  std::vector<std::complex<double>>& x,
                  int k_idx, int band_idx, double& residual);
};

} // namespace ModuleDFPT

#endif // DFPT_STERN_H