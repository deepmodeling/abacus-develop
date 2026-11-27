#include "density_matrix.h"

#include "source_base/libm/libm.h"
#include "source_base/tool_title.h"
#include "source_base/tool_quit.h"
#include "source_base/constants.h"
#include "source_base/timer.h"
#include "source_cell/klist.h"

namespace module_dm
{

template <>
void exp_mul_dmk<double>(
    const std::complex<double> kphase,
    const std::vector<std::complex<double>>& dmk_row,
    double* dmr_mat)
{
    const std::size_t mat_size = dmk_row.size();
    for (std::size_t i = 0; i < mat_size; i++)
    {
        dmr_mat[i] += kphase.real() * dmk_row[i].real() - kphase.imag() * dmk_row[i].imag();
    }
}

template <>
void exp_mul_dmk<std::complex<double>>(
    const std::complex<double> kphase,
    const std::vector<std::complex<double>>& dmk_row,
    std::complex<double>* dmr_mat)
{
    BlasConnector::axpy(dmk_row.size(), kphase, dmk_row.data(), 1, dmr_mat, 1);
}

template <>
void xyz_to_updown<double>(
    const std::complex<double> spin_block[4],
    const int icol,
    const int spin_stride[4],
    double* dmr_mat)
{
    dmr_mat[icol + spin_stride[0]] = spin_block[0].real() + spin_block[3].real();  // rho_0 = (rho_upup + rho_downdown).real()
    dmr_mat[icol + spin_stride[1]] = spin_block[1].real() + spin_block[2].real();  // rho_x = (rho_updown + rho_downup).real()
    // rho_y: the stored DM block is the complex conjugate of the physical 1-RDM P (dm_from_psi builds
    // DM_{ab}=sum conj(c_a) c_b = conj(P), so spin_block[1]=DM_{ud}=conj(P_{ud})). Extracting m_y from the
    // CONJUGATED block therefore carries the opposite sign of the bare-textbook formula; m_x/m_z read
    // Re() and are conjugation-invariant. Using the bare formula (PR #7664) sign-flips m_y and quenches
    // in-plane non-collinear moments (e.g. Mn3Sn 120-deg AFM); see issue #7831.
    dmr_mat[icol + spin_stride[2]] = spin_block[1].imag() - spin_block[2].imag();  // rho_y = Im(P_updown) - Im(P_downup)
    dmr_mat[icol + spin_stride[3]] = spin_block[0].real() - spin_block[3].real();  // rho_z = (rho_upup - rho_downdown).real()
}

template <>
void xyz_to_updown<std::complex<double>>(
    const std::complex<double> spin_block[4],
    const int icol,
    const int spin_stride[4],
    std::complex<double>* dmr_mat)
{
    dmr_mat[icol + spin_stride[0]] = spin_block[0] + spin_block[3];  // rho_0 = (rho_upup + rho_downdown)
    dmr_mat[icol + spin_stride[1]] = spin_block[1] + spin_block[2];  // rho_x = (rho_updown + rho_downup)
    // rho_y sign accounts for the conjugated stored DM block (conj(P)); see the <double> specialization above.
    dmr_mat[icol + spin_stride[2]] = -ModuleBase::IMAG_UNIT * (spin_block[1] - spin_block[2]);  // rho_y = -i*(rho_updown - rho_downup)
    dmr_mat[icol + spin_stride[3]] = spin_block[0] - spin_block[3];  // rho_z = (rho_upup - rho_downdown)
}

DmrBlock get_dmr_block(
    const Parallel_Orbitals* pv,
    const int iat1,
    const int iat2)
{
    DmrBlock block;
    block.row0 = pv->atom_begin_row[iat1];
    block.col0 = pv->atom_begin_col[iat2];
    block.nrows = pv->get_nrow_atom(iat1);
    block.ncols = pv->get_ncol_atom(iat2);
    assert(block.row0 != -1 && block.col0 != -1 && "Atom-pair not belong this process");
    return block;
}

template <typename TK, typename TR>
void build_kphase(
    hamilt::AtomPair<TR>& atom_pair,
    const std::vector<ModuleBase::Vector3<double>>& kvec_d,
    const int nk,
    const std::map<ModuleBase::Vector3<int>, std::complex<double>>& phase_hybrid,
    std::vector<std::vector<TK>>& kphase_vec,
    std::vector<TR*>& dmr_mats)
{
    const int R_size = atom_pair.get_R_size();
    kphase_vec.assign(nk, std::vector<TK>(R_size));
    dmr_mats.assign(R_size, nullptr);
    for (int iR = 0; iR < R_size; ++iR)
    {
        const ModuleBase::Vector3<int> R_index = atom_pair.get_R_index(iR);
        hamilt::BaseMatrix<TR>* const dmr_R = atom_pair.find_matrix(R_index);
#ifdef __DEBUG
        if (dmr_R == nullptr)
        {
            std::cout << "dmr_R is nullptr" << std::endl;
            continue;
        }
#endif
        dmr_mats[iR] = dmr_R->get_pointer();
        for (int ik = 0; ik < nk; ++ik)
        {
            const ModuleBase::Vector3<double> dR(R_index[0], R_index[1], R_index[2]);
            const double arg = (kvec_d[ik] * dR) * ModuleBase::TWO_PI;
            double sinp;
            double cosp;
            ModuleBase::libm::sincos(arg, &sinp, &cosp);
            kphase_vec[ik][iR] = TK(cosp, sinp);
            if (!phase_hybrid.empty())
            {
                kphase_vec[ik][iR] *= phase_hybrid.at(R_index);
            }
        }
    }
}

template <typename TK>
void transpose_dmk_block(
    const TK* dmk_col_major,
    const int ld_hk,
    const DmrBlock& block,
    TK* dmk_row)
{
    for (int icol = 0; icol < block.ncols; ++icol)
    {
        for (int irow = 0; irow < block.nrows; ++irow)
        {
            dmk_row[irow * block.ncols + icol] = dmk_col_major[icol * ld_hk + irow];
        }
    }
}

template <typename TK, typename TR>
void add_dmr_real(
    const DensityMatrix<TK, TR>& dm,
    const DmrBlock& block,
    const int ik_begin,
    const std::vector<std::vector<TK>>& kphase_vec,
    const int ld_hk,
    const int ik_in,
    std::vector<TR*>& dmr_mats)
{
    const int R_size = dmr_mats.size();
    std::vector<TK> dmk_row(block.size());
    if (ik_in >= 0)
    {
        // single k-point
        const TK* const dmk_col = dm.dmk[ik_in + ik_begin].data() + block.col0 * ld_hk + block.row0;
        transpose_dmk_block(dmk_col, ld_hk, block, dmk_row.data());
        for (int iR = 0; iR < R_size; ++iR)
        {
            exp_mul_dmk(kphase_vec[ik_in][iR], dmk_row, dmr_mats[iR]);
        }
    }
    else
    {
        // all k-points
        for (int ik = 0; ik < dm._nk; ++ik)
        {
            const TK* const dmk_col = dm.dmk[ik + ik_begin].data() + block.col0 * ld_hk + block.row0;
            transpose_dmk_block(dmk_col, ld_hk, block, dmk_row.data());
            for (int iR = 0; iR < R_size; ++iR)
            {
                exp_mul_dmk(kphase_vec[ik][iR], dmk_row, dmr_mats[iR]);
            }
        }
    }
}

template <typename TK, typename TR>
void add_dmr_soc(
    const DensityMatrix<TK, TR>& dm,
    const DmrBlock& block,
    const int ik_begin,
    const std::vector<std::vector<TK>>& kphase_vec,
    const int ld_hk,
    const int ik_in,
    const int col_stride,
    std::vector<TR*>& dmr_mats)
{
    const int mat_size = block.size();
    const int R_size = dmr_mats.size();
    std::vector<TK> soc_dmr_R(mat_size * R_size, TK(0.0, 0.0));

    // transpose DMK block to row-major and axpy into the per-R buffer
    std::vector<TK> dmk_row(mat_size);
    if (ik_in >= 0)
    {
        // single k-point
        const TK* const dmk_col = dm.dmk[ik_in + ik_begin].data() + block.col0 * ld_hk + block.row0;
        transpose_dmk_block(dmk_col, ld_hk, block, dmk_row.data());
        for (int iR = 0; iR < R_size; ++iR)
        {
            BlasConnector::axpy(mat_size,
                                kphase_vec[ik_in][iR],
                                dmk_row.data(),
                                1,
                                &soc_dmr_R[iR * mat_size],
                                1);
        }
    }
    else
    {
        // all k-points
        for (int ik = 0; ik < dm._nk; ++ik)
        {
            const TK* const dmk_col = dm.dmk[ik + ik_begin].data() + block.col0 * ld_hk + block.row0;
            transpose_dmk_block(dmk_col, ld_hk, block, dmk_row.data());
            for (int iR = 0; iR < R_size; ++iR)
            {
                BlasConnector::axpy(mat_size,
                                    kphase_vec[ik][iR],
                                    dmk_row.data(),
                                    1,
                                    &soc_dmr_R[iR * mat_size],
                                    1);
            }
        }
    }

    // spin-block column offsets for the 2x2 (upup, updown, downup, downdown) components
    int spin_stride[4]{};
    constexpr int npol = 2;
    for (int is = 0; is < npol; ++is)
    {
        for (int is2 = 0; is2 < npol; ++is2)
        {
            spin_stride[is * npol + is2] = col_stride * is + is2;
        }
    }

    // transform each 2x2 spin block to Pauli components and write back
    TK spin_block[4]{};
    for (int iR = 0; iR < R_size; ++iR)
    {
        const TK* soc_mat = &soc_dmr_R[iR * mat_size];
        TR* dmr_mat = dmr_mats[iR];
        for (int irow = 0; irow < block.nrows; irow += 2)
        {
            for (int icol = 0; icol < block.ncols; icol += 2)
            {
                spin_block[0] = soc_mat[icol + spin_stride[0]];
                spin_block[1] = soc_mat[icol + spin_stride[1]];
                spin_block[2] = soc_mat[icol + spin_stride[2]];
                spin_block[3] = soc_mat[icol + spin_stride[3]];
                xyz_to_updown(spin_block, icol, spin_stride, dmr_mat);
            }
            soc_mat += block.ncols * 2;
            dmr_mat += block.ncols * 2;
        }
    }
}

// explicit instantiations for build_kphase
template void build_kphase<std::complex<double>, double>(
    hamilt::AtomPair<double>&,
    const std::vector<ModuleBase::Vector3<double>>&,
    const int,
    const std::map<ModuleBase::Vector3<int>, std::complex<double>>&,
    std::vector<std::vector<std::complex<double>>>&,
    std::vector<double*>&);

template void build_kphase<std::complex<double>, std::complex<double>>(
    hamilt::AtomPair<std::complex<double>>&,
    const std::vector<ModuleBase::Vector3<double>>&,
    const int,
    const std::map<ModuleBase::Vector3<int>, std::complex<double>>&,
    std::vector<std::vector<std::complex<double>>>&,
    std::vector<std::complex<double>*>&);

// explicit instantiations for transpose_dmk_block
template void transpose_dmk_block<std::complex<double>>(
    const std::complex<double>*,
    const int,
    const DmrBlock&,
    std::complex<double>*);

// explicit instantiations for add_dmr_real
template void add_dmr_real<std::complex<double>, double>(
    const DensityMatrix<std::complex<double>, double>&,
    const DmrBlock&,
    const int,
    const std::vector<std::vector<std::complex<double>>>&,
    const int,
    const int,
    std::vector<double*>&);

template void add_dmr_real<std::complex<double>, std::complex<double>>(
    const DensityMatrix<std::complex<double>, std::complex<double>>&,
    const DmrBlock&,
    const int,
    const std::vector<std::vector<std::complex<double>>>&,
    const int,
    const int,
    std::vector<std::complex<double>*>&);

// explicit instantiations for add_dmr_soc
template void add_dmr_soc<std::complex<double>, double>(
    const DensityMatrix<std::complex<double>, double>&,
    const DmrBlock&,
    const int,
    const std::vector<std::vector<std::complex<double>>>&,
    const int,
    const int,
    const int,
    std::vector<double*>&);

template void add_dmr_soc<std::complex<double>, std::complex<double>>(
    const DensityMatrix<std::complex<double>, std::complex<double>>&,
    const DmrBlock&,
    const int,
    const std::vector<std::vector<std::complex<double>>>&,
    const int,
    const int,
    const int,
    std::vector<std::complex<double>*>&);

} // namespace module_dm
