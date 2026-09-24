#ifndef CHG_MIX_CFG_H
#define CHG_MIX_CFG_H

#include <string>

/// Configuration for charge mixing, aggregating the INPUT mixing parameters
/// together with the runtime globals (nspin, scf_thr_type, double_grid,
/// gamma_only_pw, domag, domag_z) that the mixing logic needs, so that
/// Charge_Mixing does not read PARAM/GlobalV directly. Callers fill this
/// from the parsed input once per run.
struct MixingConfig
{
    std::string mixing_mode;   ///< mixing mode: "plain", "broyden", "pulay"
    double mixing_beta;        ///< mixing beta for density
    int mixing_ndim;           ///< mixing ndim for broyden and pulay
    double mixing_gg0;         ///< mixing gg0 for Kerker screen
    bool mixing_tau;           ///< whether to use tau mixing
    double mixing_beta_mag;    ///< mixing beta for magnetism
    double mixing_gg0_mag;     ///< mixing gg0 for Kerker screen for magnetism
    double mixing_gg0_min;     ///< minimum kerker coefficient
    double mixing_angle;       ///< mixing angle for nspin=4
    bool mixing_dmr;           ///< whether to mix real space density matrix
    int nspin;                 ///< number of spins
    int scf_thr_type;          ///< 1: reciprocal, 2: real space threshold
    bool double_grid;          ///< whether double grid is used
    bool gamma_only_pw;        ///< whether gamma-only plane wave is used
    bool domag;                ///< whether magnetism (non-collinear) is considered
    bool domag_z;              ///< whether only the z-component magnetism is considered
    int scf_nmax;              ///< max SCF iterations, sizes the drho oscillation history (PARAM.inp.scf_nmax)
};

#endif // CHG_MIX_CFG_H
