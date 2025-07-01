#include <cmath>
#include <vector>
#include <map>
#include <tuple>
#include <complex>
#include <fstream>
#include <memory>
#include "source_cell/unitcell.h"
#include "source_base/spherical_bessel_transformer.h"
#include "source_basis/module_nao/two_center_integrator.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_cell/module_neighbor/sltk_atom_arrange.h"
#include "module_parameter/parameter.h"
#include "source_io/cal_pLpR.h"
#include "source_base/formatter.h"
#include "source_base/parallel_common.h"
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

ModuleIO::AngularMomentumCalculator::AngularMomentumCalculator(
    const std::string& orbital_dir,
    const UnitCell& ucell,
    const double& search_radius,
    const int tdestructor,
    const int tgrid,
    const int tatom,
    const bool searchpbc,
    std::ofstream* ptr_log,
    const int rank)
{
    
    // ofs_running
    this->ofs_ = ptr_log;
    *ofs_ << "\n\n\n\n";
    *ofs_ << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    *ofs_ << " |                                                                    |" << std::endl;
    *ofs_ << " |  Angular momentum expectation value calculation:                   |" << std::endl;
    *ofs_ << " |  This is a post-processing step. The expectation value of operator |" << std::endl;
    *ofs_ << " |  Lx, Ly, Lz (<a|L|b>, in which a and b are ABACUS numerical atomic |" << std::endl;
    *ofs_ << " |  orbitals) will be calculated.                                     |" << std::endl;
    *ofs_ << " |  The result will be printed to file with name ${suffix}_Lx/y/z.dat |" << std::endl;
    *ofs_ << " |                                                                    |" << std::endl;
    *ofs_ << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
    *ofs_ << "\n\n\n\n";

    int ntype_ = ucell.ntype;
#ifdef __MPI
    Parallel_Common::bcast_int(ntype_);
#endif
    std::vector<std::string> forb(ntype_);
    if (rank == 0)
    {
        for (int i = 0; i < ucell.ntype; ++i)
        {
            forb[i] = orbital_dir + ucell.orbital_fn[i];
        }
    }
#ifdef __MPI
    Parallel_Common::bcast_string(forb.data(), ntype_);
#endif
    
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
    temp = atom_arrange::set_sr_NL(*ofs_,
                                   PARAM.inp.out_level,
                                   search_radius,
                                   ucell.infoNL.get_rcutmax_Beta(),
                                   PARAM.globalv.gamma_only_local);
    temp = std::max(temp, search_radius);
    this->neighbor_searcher_ = std::unique_ptr<Grid_Driver>(new Grid_Driver(tdestructor, tgrid));
    atom_arrange::search(searchpbc,
                         *ofs_,
                         *neighbor_searcher_,
                         ucell,
                         temp,
                         tatom);
}

void ModuleIO::AngularMomentumCalculator::kernel(
    std::ofstream* ofs,
    const UnitCell& ucell,
    const char dir,
    const int precision)
{
    if (!ofs->is_open())
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
    std::string fmtstr = "%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d";
    fmtstr += "%" + std::to_string(precision*2) + "." + std::to_string(precision) + "e";
    fmtstr += "%" + std::to_string(precision*2) + "." + std::to_string(precision) + "e\n";
    FmtCore fmt(fmtstr);

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
                const Atom& atyp_j = ucell.atoms[jt];
                int ja = neighbor_searcher_->getNatom(ia_adj);
                dr = (ri - rj) * ucell.lat0;
                const ModuleBase::Vector3<int> iR = neighbor_searcher_->getBox(ia_adj);
                // the two-center-integral

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
                                        std::complex<double> val = 0;
                                        if (dir == 'x')
                                        {
                                            val = cal_LxijR(calculator_, 
                                                it, ia, li, iz, mi, jt, ja, lj, jz, mj, dr);
                                        }
                                        else if (dir == 'y')
                                        {
                                            val = cal_LyijR(calculator_, 
                                                it, ia, li, iz, mi, jt, ja, lj, jz, mj, dr);
                                        }
                                        else if (dir == 'z')
                                        {
                                            val = cal_LzijR(calculator_, 
                                                it, ia, li, iz, mi, jt, ja, lj, jz, mj, dr);
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

void ModuleIO::AngularMomentumCalculator::calculate(
    const std::string& prefix,
    const std::string& outdir,
    const UnitCell& ucell,
    const int precision,
    const int rank)
{
    if (rank != 0)
    {
        return;
    }
    std::ofstream ofout;
    const std::string dir = "xyz";
    const std::string title = "# it ia il iz im iRx iRy iRz jt ja jl jz jm <a|L|b>\n"
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
                              "# <a|L|b>: the value of the matrix element\n";
    
    for (char d : dir)
    {
        std::string fn = outdir + prefix + "_L" + d + ".dat";
        ofout.open(fn, std::ios::out);
        ofout << title;
        this->kernel(&ofout, ucell, d, precision);
        ofout.close();
    }
}