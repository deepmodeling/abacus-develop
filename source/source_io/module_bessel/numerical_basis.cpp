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
void Numerical_Basis::output_overlap(const psi::Psi<std::complex<double>>& psi,
                                     const Structure_Factor& sf,
                                     const K_Vectors& kv,
                                     const ModulePW::PW_Basis_K* wfcpw,
                                     const UnitCell& ucell,
                                     const int& index)
{
    ModuleBase::TITLE("Numerical_Basis", "output_overlap");
    ModuleBase::GlobalFunc::NEW_PART("Overlap Data For Spillage Minimization");
    const SIABPrimitiveParameters parameters = Numerical_Basis::siab_parameters_from_input(index);
    const double bessel_nao_rcut = parameters.rcut_bohr;

    //---------------------------------------------------------
    // if the numerical_basis hasn't been initialized yet,
    // then we initial here.
    //---------------------------------------------------------
    this->initialize_siab_basis(ucell, parameters);
    ModuleBase::GlobalFunc::MAKE_DIR(PARAM.inp.spillage_outdir);
    for (int derivative_order = 0; derivative_order <= 1; ++derivative_order) // Peize Lin add 2020.04.23
    {
        std::ofstream ofs;
        std::stringstream ss;
        ss << PARAM.inp.spillage_outdir << "/";

        if (PARAM.inp.bessel_nao_rcuts.size() > 1)
        {
            ss << "orb_matrix_rcut" << bessel_nao_rcut << "deriv";
        }
        else
        {
            ss << "orb_matrix.";
        } // to make it compatible with old version of orbital generation
        ss << derivative_order << ".dat";

        if (GlobalV::MY_RANK == 0)
        {
            ofs.open(ss.str().c_str());
        }

        // ALLOCATE MEMORY FOR THE OVERLAP MATRIX
        // OVERLAP : < J_mu | Psi >
        std::vector<ModuleBase::ComplexArray> overlap_Q(kv.get_nks());
        // OVERLAP : < J_mu | J_nu >
        std::vector<ModuleBase::ComplexArray> overlap_Sq(kv.get_nks());

        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "number of k points", kv.get_nks());
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "number of bands", PARAM.inp.nbands);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "number of local orbitals", PARAM.globalv.nlocal);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "number of eigenvalues of Jl(x)",
                                    this->bessel_basis.get_ecut_number());

        // CALCULATE THE OVERLAP MATRIX
        // nks now is the reduced k-points.
        for (int ik = 0; ik < kv.get_nks(); ik++)
        {
            const int npw = kv.ngk[ik];
            GlobalV::ofs_running << " --------------------------------------------------------" << std::endl;
            GlobalV::ofs_running << " Print the overlap matrixs Q and S for this kpoint" << std::endl;
            GlobalV::ofs_running << std::setw(8) << "ik" << std::setw(8) << "npw" << std::endl;
            GlobalV::ofs_running << std::setw(8) << ik + 1 << std::setw(8) << npw << std::endl;
            GlobalV::ofs_running << " --------------------------------------------------------" << std::endl;

            // search for all k-points.
            psi.fix_k(ik);
            overlap_Q[ik] = this->cal_overlap_Q(ik, npw, wfcpw, psi, static_cast<double>(derivative_order), sf, ucell);
            ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "cal_overlap_Q");

            // (2) generate Sq matrix if necessary.
            if (PARAM.inp.out_spillage == 2)
            {
#ifndef __LCAO
                // compute <jY|jY> in plane-wave basis
                overlap_Sq[ik] = this->cal_overlap_Sq(ik, npw, static_cast<double>(derivative_order), sf, wfcpw, ucell);
#else
                // compute <jY|jY> with two-center integration
                assert(derivative_order == 0 || derivative_order == 1);
                char type = (derivative_order == 0) ? 'S' : 'T';
                std::vector<int> natom;
                std::vector<int> lmax;
                std::vector<std::vector<ModuleBase::Vector3<double>>> tau_cart;
                for (int it = 0; it < ucell.ntype; ++it)
                {
                    natom.push_back(ucell.atoms[it].na);
                    lmax.push_back(ucell.atoms[it].nwl);
                    tau_cart.emplace_back();

                    for (int ia = 0; ia < ucell.atoms[it].na; ++ia)
                    {
                        tau_cart[it].push_back(ucell.atoms[it].tau[ia] * ucell.lat0);
                    }
                }

                overlap_Sq[ik] = NumericalBasis::cal_overlap_Sq(
                    type, ucell.lmaxmax, this->bessel_basis.get_ecut_number(), bessel_nao_rcut, tau_cart,
                    ucell.lat0 * ucell.latvec, NumericalBasis::indexgen(natom, lmax));
#endif
                ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "cal_overlap_Sq");
            }
        }

        const ModuleBase::matrix overlap_V
            = this->cal_overlap_V(wfcpw, psi, static_cast<double>(derivative_order), kv, ucell.tpiba);

        // ALTHOUGH THIS FUNCTION NAMES output_overlap, IT ACTUALLY OUTPUTS THE OVERLAP MATRIX HERE
#ifdef __MPI
        for (int ik = 0; ik < kv.get_nks(); ik++)
        {
            Parallel_Reduce::reduce_pool(overlap_Q[ik].ptr, overlap_Q[ik].getSize());
            // Parallel_Reduce::reduce_pool(overlap_Sq[ik].ptr, overlap_Sq[ik].getSize());
        }
        Parallel_Reduce::reduce_pool(overlap_V.c, overlap_V.nr * overlap_V.nc); // Peize Lin add 2020.04.23
#endif
        // exception handling following, for FileNotOpenFailure
        if (ofs.good()) {
            this->output_info(ofs, bessel_basis, kv, ucell); // header of orb_matrix* file
        } else {
            ModuleBase::WARNING_QUIT("Numerical_Basis", "Failed to open file for writing the overlap matrix.");
}
        // because one stage of file io complete, re-check the file status.
        if (ofs.good()) {
            this->output_k(ofs, kv); // <WEIGHTS_OF_KPOINTS>...</WEIGHTS_OF_KPOINTS>
        } else {
            ModuleBase::WARNING_QUIT("Numerical_Basis", "Failed to write k-points to file.");
}
        // because one stage of file io complete, re-check the file status.
        if (ofs.good()) {
            this->output_overlap_Q(ofs, overlap_Q, kv); // <OVERLAP_Q>...</OVERLAP_Q>
        } else {
            ModuleBase::WARNING_QUIT("Numerical_Basis", "Failed to write overlap Q to file.");
}
        // because one stage of file io complete, re-check the file status.
        if (PARAM.inp.out_spillage == 2)
        {
            // caution: this is the largest matrix to be output, always flush
            if (ofs.good()) {
                this->output_overlap_Sq(ss.str(), ofs, overlap_Sq, kv); // <OVERLAP_Sq>...</OVERLAP_Sq>
            } else {
                ModuleBase::WARNING_QUIT("Numerical_Basis", "Failed to write overlap S to file.");
}
        }
        // because one stage of file io complete, re-check the file status.
        if (ofs.good()) {
            this->output_overlap_V(ofs, overlap_V); // <OVERLAP_V>...</OVERLAP_V>
                                                    // Peize Lin add 2020.04.23
        } else {
            ModuleBase::WARNING_QUIT("Numerical_Basis", "Failed to write overlap V to file.");
}
        if (GlobalV::MY_RANK == 0) {
            ofs.close();
}
    }
    return;
}

ModuleBase::ComplexArray Numerical_Basis::cal_overlap_Q(const int& ik, const int& np, const ModulePW::PW_Basis_K* wfcpw,
                                                        const psi::Psi<std::complex<double>>& psi,
                                                        const double derivative_order, const Structure_Factor& sf,
                                                        const UnitCell& ucell) const
{
    ModuleBase::TITLE("Numerical_Basis", "cal_overlap_Q");
    ModuleBase::timer::tick("Numerical_Basis", "cal_overlap_Q");

    GlobalV::ofs_running << " OUTPUT THE OVERLAP BETWEEN SPHERICAL BESSEL FUNCTIONS AND BLOCH WAVE FUNCTIONS"
                         << std::endl;
    GlobalV::ofs_running << " Q = < J_mu, q | Psi_n, k > " << std::endl;

    ModuleBase::ComplexArray overlap_Q(PARAM.inp.nbands, PARAM.globalv.nlocal, this->bessel_basis.get_ecut_number());
    overlap_Q.zero_out();

    const double normalization = (4 * ModuleBase::PI) / sqrt(ucell.omega); // Peize Lin add normalization 2015-12-29

    std::vector<ModuleBase::Vector3<double>> gk(np);
    for (int ig = 0; ig < np; ig++)
    {
        gk[ig] = wfcpw->getgpluskcar(ik, ig) * ucell.tpiba;
    }

    const std::vector<double> gpow = Numerical_Basis::cal_gpow(gk, derivative_order);

    const ModuleBase::realArray flq = this->cal_flq(gk, ucell.lmax);

    const ModuleBase::matrix ylm = Numerical_Basis::cal_ylm(gk, ucell.lmax);

    GlobalV::ofs_running << "\n " << std::setw(5) << "ik" << std::setw(8) << "Type1" << std::setw(8) << "Atom1"
                         << std::setw(8) << "L" << std::endl;

    for (int T = 0; T < ucell.ntype; T++)
    {
        // OUT("T",T);
        for (int I = 0; I < ucell.atoms[T].na; I++)
        {
            // OUT("I",I);
            std::complex<double>* sk = sf.get_sk(ik, T, I, wfcpw);
            for (int L = 0; L < ucell.atoms[T].nwl + 1; L++)
            {
                GlobalV::ofs_running << " " << std::setw(5) << ik + 1 << std::setw(8) << ucell.atoms[T].label
                                     << std::setw(8) << I + 1 << std::setw(8) << L << std::endl;
                // OUT("l",l);
                for (int ie = 0; ie < this->bessel_basis.get_ecut_number(); ie++)
                {
                    const int N = 0;
                    assert(ucell.nmax == 1);
                    for (int m = 0; m < 2 * L + 1; m++)
                    {
                        std::vector<std::complex<double>> primitive(np);
                        Numerical_Basis::fill_reciprocal_primitive(L,
                                                                  m,
                                                                  ie,
                                                                  sk,
                                                                  flq,
                                                                  ylm,
                                                                  gpow,
                                                                  normalization,
                                                                  primitive);
                        for (int ib = 0; ib < PARAM.inp.nbands; ib++)
                        {
                            std::complex<double> overlap_tmp = ModuleBase::ZERO;
                            for (int ig = 0; ig < np; ig++)
                            {
                                overlap_tmp += conj(primitive[ig]) * psi(ib, ig); // psi is bloch orbitals
                            }
                            overlap_Q(ib, this->mu_index[T](I, L, N, m), ie) = overlap_tmp;
                        }
                    }
                } // end ie
            }     // end l
            delete[] sk;
            sk = nullptr;
        }
    }

    ModuleBase::timer::tick("Numerical_Basis", "cal_overlap_Q");
    return overlap_Q;
}

ModuleBase::ComplexArray Numerical_Basis::cal_overlap_Sq(const int& ik, const int& np, const double derivative_order,
                                                         const Structure_Factor& sf, const ModulePW::PW_Basis_K* wfcpw,
                                                         const UnitCell& ucell) const
{
    ModuleBase::TITLE("Numerical_Basis", "cal_overlap_Sq");
    ModuleBase::timer::tick("Numerical_Basis", "cal_overlap_Sq");

    GlobalV::ofs_running << " OUTPUT THE OVERLAP BETWEEN SPHERICAL BESSEL FUNCTIONS" << std::endl;
    GlobalV::ofs_running << " S = < J_mu,q1 | J_nu,q2 >" << std::endl;

    const int enumber = this->bessel_basis.get_ecut_number();
    ModuleBase::ComplexArray overlap_Sq(PARAM.globalv.nlocal, PARAM.globalv.nlocal, enumber, enumber);
    overlap_Sq.zero_out();

    const double normalization
        = (4 * ModuleBase::PI) * (4 * ModuleBase::PI) / ucell.omega; // Peize Lin add normalization 2015-12-29

    std::vector<ModuleBase::Vector3<double>> gk(np);
    for (int ig = 0; ig < np; ig++) {
        gk[ig] = wfcpw->getgpluskcar(ik, ig) * ucell.tpiba;
}

    const std::vector<double> gpow = Numerical_Basis::cal_gpow(gk, derivative_order);

    const ModuleBase::realArray flq = this->cal_flq(gk, ucell.lmax);

    const ModuleBase::matrix ylm = Numerical_Basis::cal_ylm(gk, ucell.lmax);

    GlobalV::ofs_running << "\n " << std::setw(5) << "ik" << std::setw(8) << "Type1" << std::setw(8) << "Atom1"
                         << std::setw(8) << "L1" << std::setw(8) << "Type2" << std::setw(8) << "Atom2" << std::setw(8)
                         << "L2" << std::endl;

    for (int T1 = 0; T1 < ucell.ntype; T1++) // 1.1
    {
        for (int I1 = 0; I1 < ucell.atoms[T1].na; I1++) // 1.2
        {
            std::complex<double>* sk1 = sf.get_sk(ik, T1, I1, wfcpw);
            for (int T2 = 0; T2 < ucell.ntype; T2++) // 2.1
            {
                for (int I2 = 0; I2 < ucell.atoms[T2].na; I2++) // 2.2
                {
                    std::complex<double>* sk2 = sf.get_sk(ik, T2, I2, wfcpw);
                    for (int l1 = 0; l1 < ucell.atoms[T1].nwl + 1; l1++) // 1.3
                    {
                        const std::complex<double> lphase1
                            = normalization * pow(ModuleBase::IMAG_UNIT, l1); // Peize Lin add normalization 2015-12-29
                        for (int l2 = 0; l2 < ucell.atoms[T2].nwl + 1; l2++)  // 2.3
                        {
                            GlobalV::ofs_running << " " << std::setw(5) << ik + 1 << std::setw(8)
                                                 << ucell.atoms[T1].label << std::setw(8) << I1 + 1 << std::setw(8)
                                                 << l1 << std::setw(8) << ucell.atoms[T2].label << std::setw(8)
                                                 << I2 + 1 << std::setw(8) << l2 << std::setw(8) << std::endl;

                            const std::complex<double> lphase2 = pow(ModuleBase::IMAG_UNIT, l2);
                            for (int ic1 = 0; ic1 < ucell.nmax; ic1++) // 1.5
                            {
                                for (int ic2 = 0; ic2 < ucell.nmax; ic2++) // 2.5
                                {
                                    for (int m1 = 0; m1 < 2 * l1 + 1; m1++) // 1.6
                                    {
                                        const int lm1 = l1 * l1 + m1;
                                        const int iwt1 = this->mu_index[T1](I1, l1, ic1, m1);

                                        std::vector<std::complex<double>> about_ig1(np, std::complex<double>(0.0, 0.0));
                                        for (int ig = 0; ig < np; ig++) {
                                            about_ig1[ig] = conj(lphase1 * sk1[ig] * ylm(lm1, ig))
                                                            * gpow[ig]; // Peize Lin add for dpsi 2020.04.23
}

                                        for (int m2 = 0; m2 < 2 * l2 + 1; m2++) // 2.6
                                        {
                                            const int lm2 = l2 * l2 + m2;
                                            const int iwt2 = this->mu_index[T2](I2, l2, ic2, m2);

                                            std::vector<std::complex<double>> about_ig2(np,
                                                                                        std::complex<double>(0.0, 0.0));
                                            for (int ig = 0; ig < np; ++ig) {
                                                about_ig2[ig] = lphase2 * sk2[ig] * ylm(lm2, ig) * about_ig1[ig];
}

                                            /* same as:
                                            for (int ig=0; ig<np; ig++)
                                                for (int ie1=0; ie1 < enumber; ie1++)
                                                    for (int ie2=0; ie2 < enumber; ie2++)
                                                        overlap_Sq( iwt1, iwt2, ie1, ie2) +=
                                                            about_ig2[ig] * flq(l1,ie1,ig) * flq(l2,ie2,ig);
                                            */

                                            ModuleBase::ComplexMatrix about_ig3_1(enumber, np);
                                            std::copy(&flq(l1, 0, 0), &flq(l1, 0, 0) + enumber * np, about_ig3_1.c);

                                            ModuleBase::ComplexMatrix about_ig3_2(enumber, np);
                                            for (int ie2 = 0; ie2 < enumber; ++ie2) {
                                                std::transform(&flq(l2, ie2, 0), &flq(l2, ie2, 0) + np,
                                                               about_ig2.data(), about_ig3_2.c + ie2 * np,
                                                               std::multiplies<std::complex<double>>());
}

                                            BlasConnector::gemm('N', 'T', enumber, enumber, np, 1.0, about_ig3_1.c, np,
                                                                about_ig3_2.c, np, 1.0, &overlap_Sq(iwt1, iwt2, 0, 0),
                                                                enumber);
                                        }
                                    }
                                }
                            }
                        }
                    }
                    delete[] sk2;
                    sk2 = nullptr;
                }
            }
            delete[] sk1;
            sk1 = nullptr;
        }
    }

    ModuleBase::timer::tick("Numerical_Basis", "cal_overlap_Sq");
    return overlap_Sq;
}

// Peize Lin add for dpsi 2020.04.23
ModuleBase::matrix Numerical_Basis::cal_overlap_V(const ModulePW::PW_Basis_K* wfcpw,
                                                  const psi::Psi<std::complex<double>>& psi,
                                                  const double derivative_order, const K_Vectors& kv,
                                                  const double tpiba)
{
    ModuleBase::matrix overlap_V(kv.get_nks(), PARAM.inp.nbands);
    for (int ik = 0; ik < kv.get_nks(); ++ik)
    {
        std::vector<ModuleBase::Vector3<double>> gk(kv.ngk[ik]);
        for (int ig = 0; ig < gk.size(); ig++) {
            gk[ig] = wfcpw->getgpluskcar(ik, ig) * tpiba;
}

        const std::vector<double> gpow = Numerical_Basis::cal_gpow(gk, derivative_order);

        for (int ib = 0; ib < PARAM.inp.nbands; ++ib) {
            for (int ig = 0; ig < kv.ngk[ik]; ++ig) {
                overlap_V(ik, ib) += norm(psi(ik, ib, ig)) * gpow[ig];
}
}
    }
    return overlap_V;
}

ModuleBase::realArray Numerical_Basis::cal_flq(const std::vector<ModuleBase::Vector3<double>>& gk,
                                               const int ucell_lmax) const
{
    const int np = gk.size();
    const int enumber = this->bessel_basis.get_ecut_number();

    // get flq(G) = \int f(r)jl(G*r) from interpolation table.
    ModuleBase::realArray flq(ucell_lmax + 1, enumber, np);
    for (int il = 0; il < ucell_lmax + 1; il++)
    {
        for (int ie = 0; ie < enumber; ie++)
        {
            for (int ig = 0; ig < np; ig++)
            {
                flq(il, ie, ig) = this->bessel_basis.Polynomial_Interpolation2(il, ie, gk[ig].norm());
            }
        }
    }
    return flq;
}

ModuleBase::matrix Numerical_Basis::cal_ylm(const std::vector<ModuleBase::Vector3<double>>& gk, const int ucell_lmax)
{
    const int total_lm = (ucell_lmax + 1) * (ucell_lmax + 1);
    ModuleBase::matrix ylm(total_lm, gk.size());
    ModuleBase::YlmReal::Ylm_Real(total_lm, gk.size(), gk.data(), ylm);
    return ylm;
}

std::vector<double> Numerical_Basis::cal_gpow(const std::vector<ModuleBase::Vector3<double>>& gk,
                                              const double derivative_order)
{
    constexpr double thr = 1E-12;
    std::vector<double> gpow(gk.size(), 0.0);
    for (int ig = 0; ig < gpow.size(); ++ig)
    {
        if (derivative_order >= 0)
        {
            gpow[ig] = std::pow(gk[ig].norm2(), derivative_order);
        }
        else
        {
            if (gk[ig].norm2() >= thr) {
                gpow[ig] = std::pow(gk[ig].norm2(), derivative_order);
}
        }
    }
    return gpow;
}

std::vector<ModuleBase::IntArray> Numerical_Basis::init_mu_index(const UnitCell& ucell)
{
    GlobalV::ofs_running << " Initialize the mu index" << std::endl;
    std::vector<ModuleBase::IntArray> mu_index_(ucell.ntype);

    int mu = 0;
    for (int it = 0; it < ucell.ntype; it++)
    {
        mu_index_[it].create(ucell.atoms[it].na, ucell.atoms[it].nwl + 1, ucell.nmax,
                             2 * (ucell.atoms[it].nwl + 1) + 1); // m ==> 2*l+1

        mu_index_[it].zero_out();

        // mohan added 2021-01-03
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
