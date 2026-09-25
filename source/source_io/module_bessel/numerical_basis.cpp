#include "numerical_basis.h"

#include "source_io/module_parameter/parameter.h"
#include "source_base/constants.h"
#include "source_base/global_variable.h"
#include "source_base/intarray.h"
#include "source_base/math_ylmreal.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_base/vector3.h"
#include "source_cell/module_symmetry/symmetry.h"
#include "numerical_basis_jyjy.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <functional>
#include <memory>
#include <stdexcept>
#include <vector>
Numerical_Basis::Numerical_Basis()
{
}
Numerical_Basis::~Numerical_Basis()
{
}

int Numerical_Basis::siab_abacus_m(const int conventional_m)
{
    if (conventional_m == 0)
    {
        return 0;
    }
    return conventional_m > 0 ? 2 * conventional_m - 1 : -2 * conventional_m;
}

Numerical_Basis::SIABPrimitiveParameters Numerical_Basis::siab_parameters_from_input(const int rcut_index,
                                                                                      const int lmax)
{
    if (rcut_index < 0 || rcut_index >= static_cast<int>(PARAM.inp.bessel_nao_rcuts.size()))
    {
        throw std::out_of_range("SIAB primitive rcut index is out of range");
    }
    const double primitive_ecut_ry = PARAM.inp.bessel_nao_ecut == "default"
                                         ? PARAM.inp.ecutwfc
                                         : std::stod(PARAM.inp.bessel_nao_ecut);
    if (lmax < -1)
    {
        throw std::invalid_argument("SIAB primitive lmax must be -1 or non-negative");
    }
    return SIABPrimitiveParameters{
        primitive_ecut_ry,
        PARAM.inp.bessel_nao_rcuts[rcut_index],
        PARAM.inp.bessel_nao_smooth,
        PARAM.inp.bessel_nao_sigma,
        PARAM.inp.bessel_nao_tolerence,
        lmax,
    };
}

void Numerical_Basis::initialize_siab_basis(const UnitCell& ucell,
                                            const SIABPrimitiveParameters& parameters)
{
    if (!std::isfinite(parameters.ecut_ry) || !std::isfinite(parameters.rcut_bohr)
        || !std::isfinite(parameters.sigma) || !std::isfinite(parameters.tolerance)
        || !(parameters.ecut_ry > 0.0) || !(parameters.rcut_bohr > 0.0)
        || (parameters.smooth && !(parameters.sigma > 0.0)) || !(parameters.tolerance > 0.0))
    {
        throw std::invalid_argument("SIAB primitive parameters are invalid");
    }

    const int target_lmax = parameters.lmax < 0 ? ucell.lmax : parameters.lmax;
    const std::vector<int> basis_shape = Numerical_Basis::siab_basis_shape_signature(ucell);
    const bool same_parameters = this->initialized_siab_parameters.ecut_ry == parameters.ecut_ry
                                 && this->initialized_siab_parameters.rcut_bohr == parameters.rcut_bohr
                                 && this->initialized_siab_parameters.smooth == parameters.smooth
                                 && this->initialized_siab_parameters.sigma == parameters.sigma
                                 && this->initialized_siab_parameters.tolerance == parameters.tolerance
                                 && this->initialized_siab_parameters.lmax == parameters.lmax;
    if (this->init_label && same_parameters && this->initialized_siab_ntype == ucell.ntype
        && this->initialized_siab_lmax == target_lmax && this->initialized_siab_nmax == ucell.nmax
        && this->initialized_siab_basis_shape == basis_shape)
    {
        return;
    }

    this->bessel_basis.init(false,
                            parameters.ecut_ry,
                            ucell.ntype,
                            target_lmax,
                            parameters.smooth,
                            parameters.sigma,
                            parameters.rcut_bohr,
                            parameters.tolerance,
                            ucell);
    this->mu_index = this->init_mu_index(ucell);
    this->initialized_siab_parameters = parameters;
    this->initialized_siab_ntype = ucell.ntype;
    this->initialized_siab_lmax = target_lmax;
    this->initialized_siab_nmax = ucell.nmax;
    this->initialized_siab_basis_shape = basis_shape;
    this->init_label = true;
}

std::vector<int> Numerical_Basis::siab_basis_shape_signature(const UnitCell& ucell)
{
    std::vector<int> signature;
    for (int type = 0; type < ucell.ntype; ++type)
    {
        signature.push_back(ucell.atoms[type].na);
        signature.push_back(ucell.atoms[type].nwl);
        signature.push_back(static_cast<int>(ucell.atoms[type].l_nchi.size()));
        signature.insert(signature.end(), ucell.atoms[type].l_nchi.begin(), ucell.atoms[type].l_nchi.end());
    }
    return signature;
}

void Numerical_Basis::fill_reciprocal_primitive(
    const int l,
    const int abacus_m,
    const int primitive_index,
    const std::complex<double>* structure_factor,
    const ModuleBase::realArray& flq,
    const ModuleBase::matrix& ylm,
    const std::vector<double>& gpow,
    const double normalization,
    std::vector<std::complex<double>>& values)
{
    const int npw = static_cast<int>(gpow.size());
    if (values.size() < static_cast<std::size_t>(npw))
    {
        throw std::invalid_argument("SIAB reciprocal primitive storage is too small");
    }

    const int lm = l * l + abacus_m;
    const std::complex<double> lphase = normalization * std::pow(ModuleBase::IMAG_UNIT, -l);
    for (int ig = 0; ig < npw; ++ig)
    {
        values[ig] = lphase * structure_factor[ig] * ylm(lm, ig) * flq(l, primitive_index, ig) * gpow[ig];
    }
}

std::vector<Numerical_Basis::SIABPrimitiveReciprocalBlock> Numerical_Basis::siab_primitive_reciprocal_values(
    const int ik,
    const ModulePW::PW_Basis_K* wfcpw,
    const Structure_Factor& sf,
    const UnitCell& ucell,
    const SIABPrimitiveParameters& parameters)
{
    if (wfcpw == nullptr)
    {
        throw std::invalid_argument("SIAB reciprocal primitives require a PW basis");
    }
    if (ik < 0 || ik >= wfcpw->nks)
    {
        throw std::out_of_range("SIAB primitive k-point index is out of range");
    }
    if (!(ucell.omega > 0.0))
    {
        throw std::invalid_argument("SIAB reciprocal primitives require a positive cell volume");
    }

    this->initialize_siab_basis(ucell, parameters);
    const int target_lmax = parameters.lmax < 0 ? ucell.lmax : parameters.lmax;

    const int npw = wfcpw->npwk[ik];
    std::unique_ptr<ModuleBase::realArray> flq;
    std::unique_ptr<ModuleBase::matrix> ylm;
    std::vector<double> gpow;
    if (npw > 0)
    {
        std::vector<ModuleBase::Vector3<double>> gk(npw);
        for (int ig = 0; ig < npw; ++ig)
        {
            gk[ig] = wfcpw->getgpluskcar(ik, ig) * ucell.tpiba;
        }
        flq.reset(new ModuleBase::realArray(this->cal_flq(gk, target_lmax)));
        ylm.reset(new ModuleBase::matrix(Numerical_Basis::cal_ylm(gk, target_lmax)));
        gpow.assign(npw, 1.0);
    }
    const double normalization = 4.0 * ModuleBase::PI / std::sqrt(ucell.omega);
    const int nprimitive = this->bessel_basis.get_ecut_number();

    std::vector<SIABPrimitiveReciprocalBlock> blocks;
    int offset = 0;
    int global_atom = 0;
    for (int type = 0; type < ucell.ntype; ++type)
    {
        for (int atom = 0; atom < ucell.atoms[type].na; ++atom, ++global_atom)
        {
            std::unique_ptr<std::complex<double>[]> sk(sf.get_sk(ik, type, atom, wfcpw));
            for (int l = 0; l <= target_lmax; ++l)
            {
                for (int conventional_m = -l; conventional_m <= l; ++conventional_m)
                {
                    SIABPrimitiveReciprocalBlock block;
                    block.type_index = type;
                    block.element = ucell.atoms[type].label;
                    block.atom_index = global_atom;
                    block.l = l;
                    block.m = conventional_m;
                    block.n_primitive = nprimitive;
                    block.offset = offset;
                    block.values.reserve(nprimitive);

                    const int abacus_m = Numerical_Basis::siab_abacus_m(conventional_m);
                    for (int ie = 0; ie < nprimitive; ++ie)
                    {
                        std::vector<std::complex<double>> reciprocal(static_cast<std::size_t>(npw), ModuleBase::ZERO);
                        if (npw > 0)
                        {
                            Numerical_Basis::fill_reciprocal_primitive(l,
                                                                      abacus_m,
                                                                      ie,
                                                                      sk.get(),
                                                                      *flq,
                                                                      *ylm,
                                                                      gpow,
                                                                      normalization,
                                                                      reciprocal);
                        }

                        block.values.push_back(std::move(reciprocal));
                    }
                    blocks.push_back(std::move(block));
                    offset += nprimitive;
                }
            }
        }
    }
    return blocks;
}

std::vector<Numerical_Basis::SIABPrimitiveGridBlock> Numerical_Basis::siab_primitive_grid_values(
    const int ik,
    const ModulePW::PW_Basis_K* wfcpw,
    const Structure_Factor& sf,
    const UnitCell& ucell,
    const SIABPrimitiveParameters& parameters)
{
    if (wfcpw == nullptr)
    {
        throw std::invalid_argument("SIAB primitive grid requires a PW basis");
    }
    const auto reciprocal_blocks
        = this->siab_primitive_reciprocal_values(ik, wfcpw, sf, ucell, parameters);
    const double real_space_scale = 1.0 / std::sqrt(ucell.omega);
    std::vector<SIABPrimitiveGridBlock> blocks;
    blocks.reserve(reciprocal_blocks.size());
    for (const SIABPrimitiveReciprocalBlock& reciprocal_block: reciprocal_blocks)
    {
        SIABPrimitiveGridBlock block;
        block.type_index = reciprocal_block.type_index;
        block.element = reciprocal_block.element;
        block.atom_index = reciprocal_block.atom_index;
        block.l = reciprocal_block.l;
        block.m = reciprocal_block.m;
        block.n_primitive = reciprocal_block.n_primitive;
        block.offset = reciprocal_block.offset;
        block.values.reserve(reciprocal_block.values.size());
        for (const std::vector<std::complex<double>>& reciprocal: reciprocal_block.values)
        {
            std::vector<std::complex<double>> local_grid(static_cast<std::size_t>(wfcpw->nrxx), ModuleBase::ZERO);
            if (wfcpw->gamma_only)
            {
                std::vector<double> local_real(static_cast<std::size_t>(std::max(1, wfcpw->nrxx)), 0.0);
                wfcpw->recip2real(reciprocal.data(), local_real.data(), ik);
                for (int ir = 0; ir < wfcpw->nrxx; ++ir)
                {
                    local_grid[static_cast<std::size_t>(ir)]
                        = std::complex<double>(local_real[static_cast<std::size_t>(ir)] * real_space_scale, 0.0);
                }
            }
            else
            {
                std::vector<std::complex<double>> local_complex(
                    static_cast<std::size_t>(std::max(1, wfcpw->nrxx)), ModuleBase::ZERO);
                wfcpw->recip2real(reciprocal.data(), local_complex.data(), ik);
                for (int ir = 0; ir < wfcpw->nrxx; ++ir)
                {
                    local_grid[static_cast<std::size_t>(ir)]
                        = local_complex[static_cast<std::size_t>(ir)] * real_space_scale;
                }
            }
            block.values.push_back(std::move(local_grid));
        }
        blocks.push_back(std::move(block));
    }
    return blocks;
}

//============================================================
// MEMBER FUNCTION :
// NAME : init
// DESCRIPTION : Two main functions:
// (1) start_from_file = true;
// Firstly, use check(1) to call bessel_basis.init
// to generate TableOne.
// Secondly readin C4 from file.
// Thirdly generate 3D atomic wfc in G space, put the
// results in psi.
//
// (2) If output overlap Q, start_from_file = false;
// Firstly, use check(0) to call bessel_basis,init
// to generate TableOne
// Secondly output overlap, use psi(evc) and jlq3d.
//============================================================
// void Numerical_Basis::start_from_file_k(const int& ik, ModuleBase::ComplexMatrix& psi, const Structure_Factor& sf,
//                                         const ModulePW::PW_Basis_K* wfcpw, const UnitCell& ucell)
// {
//     ModuleBase::TITLE("Numerical_Basis", "start_from_file_k");

//     if (!this->init_label)
//     {
//         // true stands for : start_from_file
//         this->bessel_basis.init(true, std::stod(PARAM.inp.bessel_nao_ecut), ucell.ntype, ucell.lmax,
//                                 PARAM.inp.bessel_nao_smooth, PARAM.inp.bessel_nao_sigma, PARAM.globalv.bessel_nao_rcut,
//                                 PARAM.inp.bessel_nao_tolerence, ucell);
//         this->mu_index = this->init_mu_index(ucell);
//         this->init_label = true;
//     }
//     this->numerical_atomic_wfc(ik, wfcpw, psi, sf, ucell);
// }

// The function is called in run_fp.cpp.
std::vector<ModuleBase::IntArray> Numerical_Basis::init_mu_index(const UnitCell& ucell)
{
    GlobalV::ofs_running << " Initialize the mu index" << std::endl;
    std::vector<ModuleBase::IntArray> mu_index_(ucell.ntype);

    int mu = 0;
    for (int it = 0; it < ucell.ntype; it++)
    {
        mu_index_[it].create(ucell.atoms[it].na, ucell.atoms[it].nwl + 1, ucell.nmax,
                             2 * (ucell.atoms[it].nwl + 1) + 1);
        mu_index_[it].zero_out();

        GlobalV::ofs_running << "Type " << it + 1 << " number_of_atoms " << ucell.atoms[it].na << " number_of_L "
                             << ucell.atoms[it].nwl + 1 << " number_of_n " << ucell.nmax << " number_of_m "
                             << 2 * (ucell.atoms[it].nwl + 1) + 1 << std::endl;

        for (int ia = 0; ia < ucell.atoms[it].na; ia++)
        {
            for (int l = 0; l < ucell.atoms[it].nwl + 1; l++)
            {
                for (int n = 0; n < ucell.atoms[it].l_nchi[l]; n++)
                {
                    for (int m = 0; m < 2 * l + 1; m++)
                    {
                        mu_index_[it](ia, l, n, m) = mu;
                        mu++;
                    }
                }
            }
        }
    }
    return mu_index_;
}

void Numerical_Basis::numerical_atomic_wfc(const int& ik, const ModulePW::PW_Basis_K* wfcpw,
                                           ModuleBase::ComplexMatrix& psi, const Structure_Factor& sf,
                                           const UnitCell& ucell)
{
    ModuleBase::TITLE("Numerical_Basis", "numerical_atomic_wfc");
    const int np = wfcpw->npwk[ik];
    std::vector<ModuleBase::Vector3<double>> gk(np);
    for (int ig = 0; ig < np; ig++) {
        gk[ig] = wfcpw->getgpluskcar(ik, ig);
}

    const int total_lm = (ucell.lmax + 1) * (ucell.lmax + 1);
    ModuleBase::matrix ylm(total_lm, np);
    ModuleBase::YlmReal::Ylm_Real(total_lm, np, gk.data(), ylm);

    std::vector<double> flq(np);
    for (int it = 0; it < ucell.ntype; it++)
    {
        // OUT("it",it);
        for (int ia = 0; ia < ucell.atoms[it].na; ia++)
        {
            // OUT("ia",ia);
            std::complex<double>* sk = sf.get_sk(ik, it, ia, wfcpw);
            for (int l = 0; l < ucell.atoms[it].nwl + 1; l++)
            {
                // OUT("l",l);
                std::complex<double> lphase = pow(ModuleBase::IMAG_UNIT, l);
                for (int ic = 0; ic < ucell.atoms[it].l_nchi[l]; ic++)
                {
                    // OUT("ic",ic);
                    for (int ig = 0; ig < np; ig++)
                    {
                        flq[ig] = this->bessel_basis.Polynomial_Interpolation(it, l, ic, gk[ig].norm() * ucell.tpiba);
                    }

                    for (int m = 0; m < 2 * l + 1; m++)
                    {
                        // OUT("m",m);
                        const int lm = l * l + m;
                        for (int ig = 0; ig < np; ig++)
                        {
                            psi(this->mu_index[it](ia, l, ic, m), ig) = lphase * sk[ig] * ylm(lm, ig) * flq[ig];
                        }
                    }
                }
            }
            delete[] sk;
            sk = nullptr;
        }
    }
}
