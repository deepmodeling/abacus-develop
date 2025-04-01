#include <cmath>
#include <vector>
#include <map>
#include <tuple>
#include <complex>
#include <fstream>
#include <memory>
#include "module_cell/unitcell.h"
#include "module_base/spherical_bessel_transformer.h"
#include "module_basis/module_nao/two_center_integrator.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_cell/module_neighbor/sltk_atom_arrange.h"
#include "module_parameter/parameter.h"
#include "module_io/cal_pLpR.h"

/**
 * 
 * FIXME: the following part will be transfered to TwoCenterIntegrator soon
 * 
 */

// L+|l, m> = sqrt((l-m)(l+m+1))|l, m+1>, return the sqrt((l-m)(l+m+1))
double _lplus_on_ylm(const int l, const int m)
{
    return std::sqrt((l - m) * (l + m + 1));
}

// L-|l, m> = sqrt((l+m)(l-m+1))|l, m-1>, return the sqrt((l+m)(l-m+1))
double _lminus_on_ylm(const int l, const int m)
{
    return std::sqrt((l + m) * (l - m + 1));
}

std::complex<double> ModuleIO::cal_LzijR(
    const std::unique_ptr<TwoCenterIntegrator>& calculator,
    const int it, const int ia, const int il, const int iz, const int mi,
    const int jt, const int ja, const int jl, const int jz, const int mj,
    const ModuleBase::Vector3<double>& vR)
{
    double val_ = 0;
    calculator->calculate(it, il, iz, mi, jt, jl, jz, mj, vR, &val_);
    return std::complex<double>(mi) * val_;
}

std::complex<double> ModuleIO::cal_LyijR(
    const std::unique_ptr<TwoCenterIntegrator>& calculator,
    const int it, const int ia, const int il, const int iz, const int im,
    const int jt, const int ja, const int jl, const int jz, const int jm,
    const ModuleBase::Vector3<double>& vR)
{
    // Ly = -i/2 * (L+ - L-)
    const double plus_ = _lplus_on_ylm(jl, jm);
    const double minus_ = _lminus_on_ylm(jl, jm);
    double val_plus = 0, val_minus = 0;
    if (plus_ != 0)
    {
        calculator->calculate(it, il, iz, im, jt, jl, jz, jm + 1, vR, &val_plus);
        val_plus *= plus_;
    }
    if (minus_ != 0)
    {
        calculator->calculate(it, il, iz, im, jt, jl, jz, jm - 1, vR, &val_minus);
        val_minus *= minus_;
    }
    return std::complex<double>(0, -0.5) * (val_plus - val_minus);
}

std::complex<double> ModuleIO::cal_LxijR(
    const std::unique_ptr<TwoCenterIntegrator>& calculator,
    const int it, const int ia, const int il, const int iz, const int im,
    const int jt, const int ja, const int jl, const int jz, const int jm,
    const ModuleBase::Vector3<double>& vR)
{   
    // Lx = 1/2 * (L+ + L-)
    const double plus_ = _lplus_on_ylm(jl, jm);
    const double minus_ = _lminus_on_ylm(jl, jm);
    double val_plus = 0, val_minus = 0;
    if (plus_ != 0)
    {
        calculator->calculate(it, il, iz, im, jt, jl, jz, jm + 1, vR, &val_plus);
        val_plus *= plus_;
    }
    if (minus_ != 0)
    {
        calculator->calculate(it, il, iz, im, jt, jl, jz, jm - 1, vR, &val_minus);
        val_minus *= minus_;
    }
    return std::complex<double>(0.5) * (val_plus + val_minus);
}

ModuleIO::AngularMomentumExpectationCalculator::AngularMomentumExpectationCalculator(
    const std::string& orbital_dir,
    const UnitCell& ucell,
    const double& search_radius,
    const int tdestructor,
    const int tgrid,
    const int tatom,
    const bool searchpbc,
    std::ofstream* ptr_log)
{
    std::vector<std::string> forb(ucell.ntype);
    for (int i = 0; i < ucell.ntype; ++i)
    {
        forb[i] = orbital_dir + ucell.orbital_fn[i];
    }
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

    // ofs_running
    this->ofs_ = ptr_log;

    // for neighbor list search
    this->neighbor_searcher_ = std::unique_ptr<Grid_Driver>(new Grid_Driver(tdestructor, tgrid));
    atom_arrange::search(
        searchpbc,
        *ofs_,
        *neighbor_searcher_,
        ucell,
        search_radius,
        tatom);
}

void ModuleIO::AngularMomentumExpectationCalculator::calculate(
    std::ofstream* ofs,
    const UnitCell& ucell,
    const char dir)
{
    if (ofs == nullptr)
    {
        return;
    }
    ModuleBase::Vector3<double> ri, rj, dr;
    for (int it = 0; it < ucell.ntype; it++)
    {
        const Atom& atyp_i = ucell.atoms[it];
        for (int ia = 0; ia < atyp_i.na; ia++)
        {
            ri = atyp_i.tau[ia];
            neighbor_searcher_->Find_atom(ucell, ri, it, ia);
            for (int ia_adj = 0; ia_adj < neighbor_searcher_->getAdjacentNum(); ia_adj++)
            {
                rj = neighbor_searcher_->getAdjacentTau(ia_adj);
                int jt = neighbor_searcher_->getType(ia_adj);
                dr = (ri - rj) * ucell.lat0;
                const ModuleBase::Vector3<int> iR = neighbor_searcher_->getBox(ia_adj);
            }
        }
    }
}