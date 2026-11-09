#include <cassert>
#include <cmath>
#include <vector>
#include <map>
#include <set>
#include <tuple>
#include <complex>
#include <fstream>
#include <memory>
#include "source_cell/unitcell.h"
#include "source_base/sph_bessel_tf.h"
#include "source_basis/module_nao/two_center_integrator.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_cell/module_neighbor/sltk_atom_arrange.h"
#include "source_io/module_hs/angmom_op_mat.h"
#include "source_io/module_hs/hs_sparse_io.h"
#include "source_base/parallel_common.h"
#include "source_base/timer.h"
/**
 * 
 * FIXME: the following part will be transfered to TwoCenterIntegrator soon
 * 
 * Notation
 * --------
 * ylm: complex spherical harmonics
 * slm: solid (real) spherical harmonics
 * 
 * Changelog
 * ---------
 * Switch to support the solid spherical harmonics to keep consistent with
 * implementations of other parts.
 * Formulation (in Chinese):
 * https://my.feishu.cn/wiki/D0enwcUKfiJgtSkJ5scc9Dagntc
 */

namespace
{
// L+ylm = sqrt((l-m)(l+m+1))ylm+1, return the sqrt((l-m)(l+m+1))
double _lambda_plus(const int l, const int m)
{
    return std::sqrt((l - m) * (l + m + 1)); // NOTE: complex spherical harmonics
}

// L-ylm = sqrt((l+m)(l-m+1))ylm-1, return the sqrt((l+m)(l-m+1))
double _lambda_minus(const int l, const int m)
{
    return std::sqrt((l + m) * (l - m + 1)); // NOTE: complex spherical harmonics
}

const std::complex<double> kImag = {0., 1.};
const double kInvSqrt2 = std::sqrt(2) * 0.5;
/// @brief Threshold for judging whether the analytical coefficient
///        sqrt(l(l+1)-m(m+1)) is exactly zero. Tighter than the output
///        sparse threshold (1e-10) because this guards an analytical
///        zero (m at extremal value), not a numerical matrix element.
constexpr double kCoeffZeroThreshold = 1e-12;

/// @brief Elements with a magnitude below this value are omitted from the
///        CSR output, matching the threshold used by the other sparse
///        matrix writers.
constexpr double kLMatrixThreshold = 1e-10;

/**
 * @brief Map a signed magnetic quantum number to its position inside one
 *        angular-momentum block.
 * @details ABACUS orders the real spherical harmonics as
 *          0, +1, -1, +2, -2, ... (see Ylm::get_ylm_real and
 *          TwoCenterIntegrator::snap), so the position is 0 for m = 0,
 *          2m - 1 for m > 0 and -2m for m < 0.
 * @param m Signed magnetic quantum number.
 * @return Zero-based position inside the (2l+1)-dimensional block.
 */
int signed_m_to_position(const int m)
{
    if (m > 0)
    {
        return 2 * m - 1;
    }
    if (m < 0)
    {
        return -2 * m;
    }
    return 0;
}

/**
 * @brief First local orbital index of every angular-momentum block of one atom.
 * @details For each l the atom stores l_nchi[l] radial functions, each carrying
 *          2l+1 real spherical harmonics. The returned vector holds the offset
 *          of block l, so that the local orbital index of a given (l, iz, m) is
 *          offset[l] + iz * (2l+1) + signed_m_to_position(m).
 * @param atom Atom whose orbital layout is inspected.
 * @return Offsets indexed by l, with atom.nwl + 1 entries.
 */
std::vector<int> orbital_block_offsets(const Atom& atom)
{
    std::vector<int> offsets(atom.nwl + 1, 0);
    for (int l = 0; l < atom.nwl; ++l)
    {
        const int block_size = atom.l_nchi[l] * (2 * l + 1);
        offsets[l + 1] = offsets[l] + block_size;
    }
    return offsets;
}

/**
 * @brief Assemble one component of L(R) into a real sparse R matrix.
 * @details The two-center integrals are evaluated on the fly for every
 *          neighboring atom pair; the resulting matrix elements are stored
 *          under the global NAO orbital indices of the two centers. Only
 *          entries above kLMatrixThreshold are kept.
 *
 *          Since the NAO basis is real, L = -i r x grad gives matrix elements
 *          that are purely imaginary, so only the imaginary part is recorded.
 * @param calculator Two-center integrator.
 * @param neighbor_searcher Neighbor list driver.
 * @param matrix Output sparse matrix keyed by R coordinate; entries hold the
 *        imaginary part of <phi_i|L|phi_j>.
 * @param all_R_coor Output set of R coordinates encountered during the scan.
 * @param ucell Unit cell.
 * @param dir Angular-momentum component, 'x', 'y' or 'z'.
 */
void assemble_angular_momentum(
    const std::unique_ptr<TwoCenterIntegrator>& calculator,
    Grid_Driver& neighbor_searcher,
    ModuleIO::SparseRMatrix<double>& matrix,
    std::set<ModuleIO::RCoordinate>& all_R_coor,
    const UnitCell& ucell,
    const char dir)
{
    const bool valid_dir = (dir == 'x' || dir == 'y' || dir == 'z');
    assert(valid_dir);

    std::complex<double> val = 0;
    ModuleBase::Vector3<double> taui;
    ModuleBase::Vector3<double> dtau;
    AdjacentAtomInfo adjinfo;
    for (int it = 0; it < ucell.ntype; it++)
    {
        const Atom& atyp_i = ucell.atoms[it];
        const std::vector<int> offsets_i = orbital_block_offsets(atyp_i);
        for (int ia = 0; ia < atyp_i.na; ia++)
        {
            const int iat_i = ucell.itia2iat(it, ia);
            taui = ucell.get_tau(iat_i);
            neighbor_searcher.Find_atom(ucell, taui, it, ia, &adjinfo);
            for (int ia_adj = 0; ia_adj < adjinfo.adj_num + 1; ia_adj++)
            {
                const int jt = adjinfo.ntype[ia_adj];
                const int ja = adjinfo.natom[ia_adj];
                const Atom& atyp_j = ucell.atoms[jt];
                const ModuleBase::Vector3<int> iR = adjinfo.box[ia_adj];
                const int iat_j = ucell.itia2iat(jt, ja);
                dtau = ucell.cal_dtau(iat_i, iat_j, iR) * ucell.lat0;

                const ModuleIO::RCoordinate R_coor(iR.x, iR.y, iR.z);
                all_R_coor.insert(R_coor);
                ModuleIO::SparseRBlock<double>& block = matrix[R_coor];

                const std::vector<int> offsets_j = orbital_block_offsets(atyp_j);
                for (int li = 0; li < atyp_i.nwl + 1; li++)
                {
                    const int block_size_i = 2 * li + 1;
                    for (int iz = 0; iz < atyp_i.l_nchi[li]; iz++)
                    {
                        for (int mi = -li; mi <= li; mi++)
                        {
                            const int pos_i = signed_m_to_position(mi);
                            const int iw_i = offsets_i[li] + iz * block_size_i + pos_i;
                            const int row = ucell.itiaiw2iwt(it, ia, iw_i);
                            for (int lj = 0; lj < atyp_j.nwl + 1; lj++)
                            {
                                const int block_size_j = 2 * lj + 1;
                                for (int jz = 0; jz < atyp_j.l_nchi[lj]; jz++)
                                {
                                    for (int mj = -lj; mj <= lj; mj++)
                                    {
                                        if (dir == 'x')
                                        {
                                            val = ModuleIO::cal_LxijR(calculator,
                                                it, ia, li, iz, mi, jt, ja, lj, jz, mj, dtau);
                                        }
                                        else if (dir == 'y')
                                        {
                                            val = ModuleIO::cal_LyijR(calculator,
                                                it, ia, li, iz, mi, jt, ja, lj, jz, mj, dtau);
                                        }
                                        else
                                        {
                                            val = ModuleIO::cal_LzijR(calculator,
                                                it, ia, li, iz, mi, jt, ja, lj, jz, mj, dtau);
                                        }

                                        if (std::fabs(val.imag()) > kLMatrixThreshold)
                                        {
                                            const int pos_j = signed_m_to_position(mj);
                                            const int iw_j = offsets_j[lj] + jz * block_size_j + pos_j;
                                            const int col = ucell.itiaiw2iwt(jt, ja, iw_j);
                                            block[row][col] = val.imag();
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
} // namespace

std::complex<double> ModuleIO::cal_LzijR(
    const std::unique_ptr<TwoCenterIntegrator>& calculator,
    const int it, const int ia, const int il, const int iz, const int mi,
    const int jt, const int ja, const int jl, const int jz, const int mj,
    const ModuleBase::Vector3<double>& vR)
{
    if(mj == 0)
    {
        return std::complex<double>(0.);
    }
    double val_ = 0;
    calculator->calculate(it, il, iz, mi, jt, jl, jz, -mj, vR, &val_);
    return kImag * static_cast<double>(mj) * val_;
}

std::complex<double> ModuleIO::cal_LxijR(
    const std::unique_ptr<TwoCenterIntegrator>& calculator,
    const int it, const int ia, const int il, const int iz, const int im,
    const int jt, const int ja, const int jl, const int jz, const int jm,
    const ModuleBase::Vector3<double>& vR)
{
    const double lmbdp = _lambda_plus(jl, jm);
    const double lmbdm = _lambda_minus(jl, jm);
    // two-center-integral placeholders
    double valp = 0.;
    double valm = 0.;
    if (jm > 1)
    {
        if (std::fabs(lmbdp) > kCoeffZeroThreshold)
        {
            calculator->calculate(it, il, iz, im, jt, jl, jz, -(jm+1), vR, &valp);
        }
        if (std::fabs(lmbdm) > kCoeffZeroThreshold)
        {
            calculator->calculate(it, il, iz, im, jt, jl, jz, -(jm-1), vR, &valm);
        }
        return kImag * 0.5 * (lmbdp * valp + lmbdm * valm);
    }
    if (jm == 1)
    {
        if (std::fabs(lmbdp) > kCoeffZeroThreshold)
        {
            calculator->calculate(it, il, iz, im, jt, jl, jz, -2, vR, &valp);
        }
        return kImag * 0.5 * lmbdp * valp;
    }
    if (jm == 0)
    {
        const double lmbd = _lambda_plus(jl, 0); // std::sqrt(jl*(jl+1))
        if (std::fabs(lmbd) > kCoeffZeroThreshold)
        {
            calculator->calculate(it, il, iz, im, jt, jl, jz, -1, vR, &valp);
        }
        return kImag * kInvSqrt2 * lmbd * valp;
    }
    if (jm == -1)
    {
        if (std::fabs(lmbdp) > kCoeffZeroThreshold)
        {
            calculator->calculate(it, il, iz, im, jt, jl, jz, 0, vR, &valp);
        }
        if (std::fabs(lmbdm) > kCoeffZeroThreshold)
        {
            calculator->calculate(it, il, iz, im, jt, jl, jz, 2, vR, &valm);
        }
        return -kImag * 0.5 * (std::sqrt(2) * lmbdp * valp + lmbdm * valm);
    }
    else
    {
        assert(jm < -1); // defensive check
        if (std::fabs(lmbdp) > kCoeffZeroThreshold)
        {
            calculator->calculate(it, il, iz, im, jt, jl, jz, -(jm+1), vR, &valp);
        }
        if (std::fabs(lmbdm) > kCoeffZeroThreshold)
        {
            calculator->calculate(it, il, iz, im, jt, jl, jz, -(jm-1), vR, &valm);
        }
        return -kImag * 0.5 * (lmbdp * valp + lmbdm * valm);
    }
}

std::complex<double> ModuleIO::cal_LyijR(
    const std::unique_ptr<TwoCenterIntegrator>& calculator,
    const int it, const int ia, const int il, const int iz, const int im,
    const int jt, const int ja, const int jl, const int jz, const int jm,
    const ModuleBase::Vector3<double>& vR)
{   
    const double lmbdp = _lambda_plus(jl, jm);
    const double lmbdm = _lambda_minus(jl, jm);
    // two-center-integral placeholders
    double valp = 0.;
    double valm = 0.;
    if (jm > 1)
    {
        if (std::fabs(lmbdp) > kCoeffZeroThreshold)
        {
            calculator->calculate(it, il, iz, im, jt, jl, jz, jm+1, vR, &valp);
        }
        if (std::fabs(lmbdm) > kCoeffZeroThreshold)
        {
            calculator->calculate(it, il, iz, im, jt, jl, jz, jm-1, vR, &valm);
        }
        return -kImag * 0.5 * (lmbdp * valp - lmbdm * valm);
    }
    if (jm == 1)
    {
        if (std::fabs(lmbdp) > kCoeffZeroThreshold)
        {
            calculator->calculate(it, il, iz, im, jt, jl, jz, 2, vR, &valp);
        }
        if (std::fabs(lmbdm) > kCoeffZeroThreshold)
        {
            calculator->calculate(it, il, iz, im, jt, jl, jz, 0, vR, &valm);
        }
        return -kImag * 0.5 * (lmbdp * valp - std::sqrt(2) * lmbdm * valm);
    }
    if (jm == 0)
    {
        const double lmbd = _lambda_plus(jl, 0); // std::sqrt(l*(l+1))
        if (std::fabs(lmbd) > kCoeffZeroThreshold)
        {
            calculator->calculate(it, il, iz, im, jt, jl, jz, 1, vR, &valp);
        }
        return -kImag * kInvSqrt2 * lmbd * valp;
    }
    if (jm == -1)
    {
        if (std::fabs(lmbdm) > kCoeffZeroThreshold)
        {
            calculator->calculate(it, il, iz, im, jt, jl, jz, -2, vR, &valm);
        }
        return -kImag * 0.5 * lmbdm * valm;
    }
    else
    {
        assert(jm < -1); // defensive check
        if (std::fabs(lmbdp) > kCoeffZeroThreshold)
        {
            calculator->calculate(it, il, iz, im, jt, jl, jz, jm+1, vR, &valp);
        }
        if (std::fabs(lmbdm) > kCoeffZeroThreshold)
        {
            calculator->calculate(it, il, iz, im, jt, jl, jz, jm-1, vR, &valm);
        }
        return kImag * 0.5 * (lmbdp * valp - lmbdm * valm);
    }
}

ModuleIO::Angmom_op::Angmom_op(
    const std::string& orbital_dir,
    const UnitCell& ucell,
    const double& search_radius,
    const int tdestructor,
    const int tgrid,
    const int tatom,
    const bool searchpbc,
    const std::string& out_level,
    const bool gamma_only,
    std::ofstream* ptr_log,
    const int rank)
{
    
    this->ofs_ = ptr_log;
    if (this->ofs_ == nullptr)
    {
        this->fallback_ofs_.open("/dev/null");
        this->ofs_ = &this->fallback_ofs_;
    }

    *ofs_ << " Calculate angular momentum expectation values Lx, Ly, Lz (<a|L|b>) in NAO basis." << std::endl;

    int ntype_ = ucell.ntype;
    Parallel_Common::bcast_int(ntype_);
    std::vector<std::string> forb(ntype_);
    if (rank == 0)
    {
        for (int i = 0; i < ucell.ntype; ++i)
        {
            forb[i] = orbital_dir + ucell.orbital_fn[i];
        }
    }
    Parallel_Common::bcast_string(forb.data(), ntype_);
    
    this->orb_ = std::unique_ptr<RadialCollection>(new RadialCollection);
    this->orb_->build(ucell.ntype, forb.data(), 'o');
    
    ModuleBase::SphericalBesselTransformer sbt(true);
    this->orb_->set_transformer(sbt);
    
    const double rcut_max = orb_->rcut_max();
    const int ngrid = int(rcut_max / 0.01) + 1;
    const double cutoff = 2.0 * rcut_max;
    this->orb_->set_uniform_grid(true, ngrid, cutoff, 'i', true);
    
    this->calculator_ = std::unique_ptr<TwoCenterIntegrator>(new TwoCenterIntegrator);
    this->calculator_->tabulate(*orb_, *orb_, 'S', ngrid, cutoff);
    
    // Initialize Ylm coefficients
    ModuleBase::Ylm::set_coefficients();
    
    // for neighbor list search
    double temp = -1.0;
    if (search_radius < rcut_max)
    {
        *ofs_ << "Find the `search_radius` from the input file being smaller than the \n"
                 "`rcut_max` of the orbitals.\n"
              << "Reset the `search_radius` (" << search_radius << ") "
              << "to `rcut_max` ("<< rcut_max << ")." 
              << std::endl;
        // we don't really set, but use std::max to mask :)
    }
    temp = atom_arrange::set_sr_NL(*ofs_,
                                   out_level,
                                   std::max(search_radius, rcut_max),
                                   ucell.infoNL->get_rcutmax_Beta(),
                                   gamma_only);
    temp = std::max(temp, std::max(search_radius, rcut_max));
    this->neighbor_searcher_ = std::unique_ptr<Grid_Driver>(new Grid_Driver(tdestructor, tgrid));
    atom_arrange::search(searchpbc,
                         *ofs_,
                         *neighbor_searcher_,
                         ucell,
                         temp,
                         tatom);
}

void ModuleIO::Angmom_op::calculate(
    const std::string& prefix,
    const std::string& outdir,
    const UnitCell& ucell,
    const Parallel_Orbitals& pv,
    const int precision,
    const int rank,
    const int istep,
    const int out_app_flag)
{
    ModuleBase::TITLE("Angmom_op", "calculate");
    ModuleBase::timer::start("Angmom_op", "calculate");

    // The CSR writer below covers the scalar (npol = 1) basis only. For
    // nspin = 4 the orbital dimension is doubled by npol, so reject the
    // request explicitly instead of writing a matrix with mismatched sizes.
    const int npol = ucell.get_npol();
    if (npol != 1)
    {
        ModuleBase::WARNING_QUIT("Angmom_op::calculate",
                                 "out_mat_l is not supported for nspin = 4 (npol = 2).");
    }

    if (rank != 0)
    {
        ModuleBase::timer::end("Angmom_op", "calculate");
        return;
    }

    /// @brief Filename suffix: modes 0 and 1 reuse the fixed name (mode 0
    ///        overwrites, mode 1 appends each written step); only mode 2
    ///        carries "g{istep+1}" (e.g., "g1") following the out_freq_ion
    ///        convention.
    const std::string step_suffix = (istep >= 0 && out_app_flag == 2)
        ? "g" + std::to_string(istep + 1) : "";
    const std::string dir = "xyz";

    for (char d : dir)
    {
        ModuleIO::SparseRMatrix<double> matrix;
        std::set<ModuleIO::RCoordinate> all_R_coor;
        assemble_angular_momentum(this->calculator_,
                                  *this->neighbor_searcher_,
                                  matrix,
                                  all_R_coor,
                                  ucell,
                                  d);

        std::string fn = outdir + "l" + d + step_suffix + "_nao.csr";
        ModuleIO::SparseWriteOptions options;
        options.filename = fn;
        options.label = std::string("L") + d;
        options.threshold = kLMatrixThreshold;
        options.binary = false;
        options.precision = precision;
        options.istep = istep;
        // Rank 0 owns the fully assembled matrix, so no MPI reduction is
        // needed and every rank-local entry is already global.
        options.reduce = false;
        options.out_app_flag = out_app_flag;
        options.csr_format = true;
        options.ucell = &ucell;

        ModuleIO::save_sparse(matrix, all_R_coor, pv, options);
        *ofs_ << " Write L(R) (" << d << " component) matrix in NAO basis to file: " << fn << std::endl;
    }

    ModuleBase::timer::end("Angmom_op", "calculate");
}
