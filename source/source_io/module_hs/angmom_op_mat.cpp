#include <cassert>
#include <cmath>
#include <vector>
#include <map>
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
#include "source_base/formatter.h"
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

void ModuleIO::Angmom_op::kernel(
    std::ofstream* ofs,
    const UnitCell& ucell,
    const char dir,
    const int precision)
{
    if (ofs == nullptr || !ofs->is_open())
    {
        return;
    }
    // an easy sanity check
    assert(dir == 'x' || dir == 'y' || dir == 'z');

    // it, ia, il, iz, im, iRx, iRy, iRz, jt, ja, jl, jz, jm
    // the iRx, iRy, iRz are the indices of the supercell in which the two-center-integral
    // it and jt are indexes of atomtypes,
    // ia and ja are indexes of atoms within the atomtypes,
    // il and jl are indexes of the angular momentum,
    // iz and jz are indexes of the zeta functions
    // im and jm are indexes of the magnetic quantum numbers.
    std::string fmtstr = "%4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d";
    fmtstr += " %" + std::to_string(precision*2) + "." + std::to_string(precision) + "e";
    fmtstr += " %" + std::to_string(precision*2) + "." + std::to_string(precision) + "e\n";
    FmtCore fmt(fmtstr);

    // placeholders
    std::complex<double> val = 0;
    ModuleBase::Vector3<double> taui; // the origin position
    ModuleBase::Vector3<double> dtau; // the displacement
    AdjacentAtomInfo adjinfo; // adjacent atom information carrier
    for (int it = 0; it < ucell.ntype; it++)
    {
        const Atom& atyp_i = ucell.atoms[it];
        for (int ia = 0; ia < atyp_i.na; ia++)
        {
            taui = ucell.get_tau(ucell.itia2iat(it, ia));
            neighbor_searcher_->Find_atom(ucell, taui, it, ia, &adjinfo);
            for (int ia_adj = 0; ia_adj < adjinfo.adj_num + 1; ia_adj++) // "+1" is to include itself
            {
                int jt = adjinfo.ntype[ia_adj]; // ityp
                int ja = adjinfo.natom[ia_adj]; // iat with in atomtype
                const Atom& atyp_j = ucell.atoms[jt];
                const ModuleBase::Vector3<int> iR = adjinfo.box[ia_adj];
                dtau = ucell.cal_dtau(ucell.itia2iat(it, ia), 
                                      ucell.itia2iat(jt, ja), 
                                      iR) * ucell.lat0; // convert to unit of Bohr

                // nested loop: calculate the two-center-integral
                for (int li = 0; li < atyp_i.nwl + 1; li++)
                {
                    for (int iz = 0; iz < atyp_i.l_nchi[li]; iz++)
                    {
                        for (int mi = -li; mi <= li; mi++)
                        {
                            for (int lj = 0; lj < atyp_j.nwl + 1; lj++)
                            {
                                for (int jz = 0; jz < atyp_j.l_nchi[lj]; jz++)
                                {
                                    for (int mj = -lj; mj <= lj; mj++)
                                    {
                                        if (dir == 'x')
                                        {
                                            val = cal_LxijR(calculator_, 
                                                it, ia, li, iz, mi, jt, ja, lj, jz, mj, dtau);
                                        }
                                        else if (dir == 'y')
                                        {
                                            val = cal_LyijR(calculator_, 
                                                it, ia, li, iz, mi, jt, ja, lj, jz, mj, dtau);
                                        }
                                        else if (dir == 'z')
                                        {
                                            val = cal_LzijR(calculator_, 
                                                it, ia, li, iz, mi, jt, ja, lj, jz, mj, dtau);
                                        }

                                        *ofs << fmt.format(
                                            it, ia, li, iz, mi,
                                            iR.x, iR.y, iR.z,
                                            jt, ja, lj, jz, mj,
                                            val.real(), val.imag());
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

void ModuleIO::Angmom_op::calculate(
    const std::string& prefix,
    const std::string& outdir,
    const UnitCell& ucell,
    const int precision,
    const int rank,
    const int istep)
{
    ModuleBase::TITLE("Angmom_op", "calculate");
    ModuleBase::timer::start("Angmom_op", "calculate");

    if (rank != 0)
    {
        ModuleBase::timer::end("Angmom_op", "calculate");
        return;
    }
    std::ofstream ofout;
    const std::string dir = "xyz";
    /// @brief Filename suffix: if istep >= 0, append "g{istep+1}" (e.g., "g1")
    ///        to follow the out_freq_ion convention; otherwise no suffix.
    const std::string step_suffix = (istep >= 0) ? "g" + std::to_string(istep + 1) : "";
    const std::string title = "# it ia il iz im iRx iRy iRz jt ja jl jz jm Re[<a|L|b>] Im[<a|L|b>]\n"
                              "# it: atomtype index of the first atom\n"
                              "# ia: atomic index of the first atom within the atomtype\n"
                              "# il: angular momentum index of the first atom\n"
                              "# iz: zeta function index of the first atom\n"
                              "# im: magnetic quantum number of the first atom\n"
                              "# iRx, iRy, iRz: the indices of the supercell\n"
                              "# jt: atomtype index of the second atom\n"
                              "# ja: atomic index of the second atom within the atomtype\n"
                              "# jl: angular momentum index of the second atom\n"
                              "# jz: zeta function index of the second atom\n"
                              "# jm: magnetic quantum number of the second atom\n"
                              "# Re[<a|L|b>], Im[<a|L|b>]: the real and imaginary parts "
                              "of the value of the matrix element\n";
    
    for (char d : dir)
    {
        std::string fn = outdir + "l" + d + step_suffix + "_nao.txt";
        ofout.open(fn, std::ios::out);
        ofout << title;
        this->kernel(&ofout, ucell, d, precision);
        ofout.close();
        *ofs_ << " Write L(R) (" << d << " component) matrix in NAO basis to file: " << fn << std::endl;
    }
    ModuleBase::timer::end("Angmom_op", "calculate");
}
