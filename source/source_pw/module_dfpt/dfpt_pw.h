// ============================================================
// This code is added by Mohan Chen on 2026-05-18.
// This code is currently in the design phase and has not been
// put into production yet. It may change in the future.
// Please use this code with caution. Only developers who know
// what they are doing should use this code.
// ============================================================

#ifndef DFPT_PW_H
#define DFPT_PW_H

#include <string>
#include <vector>
#include "source_cell/unitcell.h"
#include "source_psi/psi.h"

class Plus_U;

namespace ModuleDFPT {

class DFPT_PW {
public:
    DFPT_PW();
    ~DFPT_PW();
    
    void init(UnitCell& ucell, const psi::Psi<std::complex<double>>& psi, 
              double nelec, double ecutwfc, const Plus_U* dftu);
    
    void run();
    
    /// DFT+U reservation accessors (U0): with_u() reports whether a DFT+U
    /// provider is wired (dft_plus_u enabled upstream); u_active() further
    /// requires the provider to be usable (locale initialized, i.e. the LCAO
    /// orbital files are present).
    bool get_with_u() const;
    bool get_u_active() const;
    
    std::vector<double> get_phonon_freq(int q_idx) const;
    
    ModuleBase::matrix get_dielectric_tensor() const;
    
    ModuleBase::matrix get_born_charges(int atom_idx) const;
    
    void set_parameters(const std::string& param_file);
    
    void set_qmesh(int nqx, int nqy, int nqz);
    
    void set_conv_thr(double thr);
    
    void set_max_iter(int max_iter);
    
private:
    class Impl;
    Impl* pimpl_;
};

} // namespace ModuleDFPT

#endif // DFPT_PW_H