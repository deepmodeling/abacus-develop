#ifndef DFTU_TYPES_H
#define DFTU_TYPES_H

namespace DFTU_BASE
{

/// Frequency configuration for numbered occupation-matrix output.
///
/// Each output ionic step (out_freq_ion) owns one file occ_matg{#}.txt.
/// Within that file, one section is appended at every electronic step
/// selected by out_freq_elec, scf_nmax or convergence.
struct OccmatOutputCfg
{
    int out_freq_ion;  ///< ionic-step interval; 0 disables numbered files
    int out_freq_elec; ///< electronic-iteration interval recorded inside a numbered file
    int scf_nmax;      ///< maximum number of electronic iterations
    bool out_occ_mat;  ///< master switch of occupation-matrix output (out_occ_mat INPUT parameter)
    /// DFT+U calculation mode (INPUT parameter dft_plus_u). The occupation-
    /// matrix IO is meaningful only when dft_plus_u > 0; otherwise the
    /// occupation matrix is never computed and the writers must early-out
    /// to avoid touching an empty l_channel vector. This honours the
    /// documented contract that out_occ_mat only takes effect for DFT+U.
    int dft_plus_u;
};

/// Text format used by write_occup_m().
enum OccmatTextFormat
{
    OCMAT_FMT_LEGACY,   ///< whitespace-token layout parsed by read_occup_m()
    OCMAT_FMT_READABLE  ///< compact human-readable layout for snapshot files
};

/// Storage layout of the real (2l+1)*npol matrix when nspin == 4.
enum OccmatSocLayout
{
    /// 4 contiguous Pauli-component blocks [charge, sigma_x, sigma_y, sigma_z],
    /// each (2l+1)x(2l+1); used by the PW path
    SOC_LAYOUT_PAULI,
    /// real symmetric matrix in spin basis (m + spin-polarization index);
    /// imaginary parts are not stored; used by the LCAO path
    SOC_LAYOUT_SPIN_BASIS_REAL
};

} // namespace DFTU_BASE

#endif // DFTU_TYPES_H
