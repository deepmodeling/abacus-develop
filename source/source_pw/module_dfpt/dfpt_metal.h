// ============================================================
// This code is added by Mohan Chen on 2026-05-18.
// This code is currently in the design phase and has not been
// put into production yet. It may change in the future.
// Please use this code with caution. Only developers who know
// what they are doing should use this code.
// ============================================================

#ifndef DFPT_METAL_H
#define DFPT_METAL_H

#include "dfpt_pw_data.h"
#include "source_psi/psi.h"

namespace ModuleDFPT {

class DFPT_Metal {
public:
    DFPT_Metal();
    ~DFPT_Metal();
    
    void init(double sigma, const std::string& smearing_type);
    
    void dfdeps(const ModuleBase::matrix& eig, double efermi, 
                ModuleBase::matrix& dfdeps);
    
    void compute_dmu(int q_idx, const psi::Psi<std::complex<double>>& psi,
                     const ModuleBase::matrix& wg, const ModuleBase::matrix& dfdeps,
                     DFPT_PW_Data& data);
    
    void compute_drho_metal(int q_idx, const psi::Psi<std::complex<double>>& psi,
                            const ModuleBase::matrix& wg, const ModuleBase::matrix& dfdeps,
                            double dmu, DFPT_PW_Data& data);

private:
    double sigma_ = 0.0;
    std::string smearing_type_ = "gaussian";
    
    double fd_dfdeps(double e, double efermi);
    
    double gauss_dfdeps(double e, double efermi);
};

} // namespace ModuleDFPT

#endif // DFPT_METAL_H