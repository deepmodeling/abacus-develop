/**
 * @file py_esolver_lcao.cpp
 * @brief Python bindings for ESolver_KS_LCAO
 *
 * This file provides pybind11 bindings for the LCAO ESolver,
 * enabling Python-controlled SCF workflows with breakpoint support.
 */

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>

#include "py_esolver_lcao.hpp"

// ABACUS headers for actual implementation
#include "source_estate/module_charge/charge.h"
#include "source_estate/fp_energy.h"
#include "source_estate/elecstate.h"
#include "source_estate/module_dm/density_matrix.h"
#include "source_lcao/hamilt_lcao.h"
#include "source_lcao/module_hcontainer/hcontainer.h"
#include "source_basis/module_ao/parallel_orbitals.h"

#include <complex>
#include <stdexcept>
#include <iostream>

namespace py = pybind11;
using namespace pybind11::literals;

namespace py_esolver
{

// ============================================================================
// PyChargeAccessor Implementation
// ============================================================================

void PyChargeAccessor::set_from_charge(const Charge* chr)
{
    if (chr == nullptr)
    {
        chr_ptr_ = nullptr;
        rho_ptr_ = nullptr;
        nspin_ = 0;
        nrxx_ = 0;
        ngmc_ = 0;
        return;
    }

    chr_ptr_ = chr;
    rho_ptr_ = nullptr;  // Use chr_ptr_ instead
    nspin_ = chr->nspin;
    nrxx_ = chr->nrxx;
    ngmc_ = chr->ngmc;
}

void PyChargeAccessor::set_data(const double* rho_ptr, int nspin, int nrxx)
{
    chr_ptr_ = nullptr;  // Not using Charge object
    rho_ptr_ = rho_ptr;
    nspin_ = nspin;
    nrxx_ = nrxx;
    ngmc_ = 0;
}

py::array_t<double> PyChargeAccessor::get_rho() const
{
    if (!is_valid())
    {
        throw std::runtime_error("Charge data not available. Run SCF first.");
    }

    // Create numpy array with shape (nspin, nrxx)
    std::vector<ssize_t> shape = {static_cast<ssize_t>(nspin_), static_cast<ssize_t>(nrxx_)};

    auto result = py::array_t<double>(shape);
    auto buf = result.request();
    double* ptr = static_cast<double*>(buf.ptr);

    // Copy data from either chr_ptr_ or rho_ptr_
    if (chr_ptr_ != nullptr && chr_ptr_->rho != nullptr)
    {
        // Copy from Charge object (rho is double** with shape [nspin][nrxx])
        for (int is = 0; is < nspin_; ++is)
        {
            if (chr_ptr_->rho[is] != nullptr)
            {
                std::copy(chr_ptr_->rho[is], chr_ptr_->rho[is] + nrxx_, ptr + is * nrxx_);
            }
        }
    }
    else if (rho_ptr_ != nullptr)
    {
        // Copy from flat array (legacy mode)
        std::copy(rho_ptr_, rho_ptr_ + nspin_ * nrxx_, ptr);
    }
    else
    {
        throw std::runtime_error("No valid charge data source.");
    }

    return result;
}

py::array_t<std::complex<double>> PyChargeAccessor::get_rhog() const
{
    if (chr_ptr_ == nullptr || chr_ptr_->rhog == nullptr)
    {
        throw std::runtime_error("Reciprocal-space charge density not available.");
    }

    // Create numpy array with shape (nspin, ngmc)
    std::vector<ssize_t> shape = {static_cast<ssize_t>(nspin_), static_cast<ssize_t>(ngmc_)};

    auto result = py::array_t<std::complex<double>>(shape);
    auto buf = result.request();
    std::complex<double>* ptr = static_cast<std::complex<double>*>(buf.ptr);

    // Copy from Charge object (rhog is complex** with shape [nspin][ngmc])
    for (int is = 0; is < nspin_; ++is)
    {
        if (chr_ptr_->rhog[is] != nullptr)
        {
            std::copy(chr_ptr_->rhog[is], chr_ptr_->rhog[is] + ngmc_, ptr + is * ngmc_);
        }
    }

    return result;
}

py::array_t<double> PyChargeAccessor::get_rho_core() const
{
    if (chr_ptr_ == nullptr || chr_ptr_->rho_core == nullptr)
    {
        throw std::runtime_error("Core charge density not available.");
    }

    // Create numpy array with shape (nrxx,)
    std::vector<ssize_t> shape = {static_cast<ssize_t>(nrxx_)};

    auto result = py::array_t<double>(shape);
    auto buf = result.request();
    double* ptr = static_cast<double*>(buf.ptr);

    std::copy(chr_ptr_->rho_core, chr_ptr_->rho_core + nrxx_, ptr);

    return result;
}

// ============================================================================
// PyEnergyAccessor Implementation
// ============================================================================

void PyEnergyAccessor::set_from_fenergy(const elecstate::fenergy* f_en)
{
    if (f_en == nullptr)
    {
        etot_ = 0.0;
        eband_ = 0.0;
        hartree_energy_ = 0.0;
        etxc_ = 0.0;
        ewald_energy_ = 0.0;
        demet_ = 0.0;
        exx_ = 0.0;
        evdw_ = 0.0;
        return;
    }

    etot_ = f_en->etot;
    eband_ = f_en->eband;
    hartree_energy_ = f_en->hartree_energy;
    etxc_ = f_en->etxc;
    ewald_energy_ = f_en->ewald_energy;
    demet_ = f_en->demet;
    exx_ = f_en->exx;
    evdw_ = f_en->evdw;
}

void PyEnergyAccessor::set_energies(double etot, double eband, double hartree,
                                     double etxc, double ewald, double demet,
                                     double exx, double evdw)
{
    etot_ = etot;
    eband_ = eband;
    hartree_energy_ = hartree;
    etxc_ = etxc;
    ewald_energy_ = ewald;
    demet_ = demet;
    exx_ = exx;
    evdw_ = evdw;
}

py::dict PyEnergyAccessor::get_all_energies() const
{
    py::dict result;
    result["etot"] = etot_;
    result["eband"] = eband_;
    result["hartree_energy"] = hartree_energy_;
    result["etxc"] = etxc_;
    result["ewald_energy"] = ewald_energy_;
    result["demet"] = demet_;
    result["exx"] = exx_;
    result["evdw"] = evdw_;
    return result;
}

// ============================================================================
// PyHamiltonianAccessor Implementation (template)
// ============================================================================

template <typename TK, typename TR>
void PyHamiltonianAccessor<TK, TR>::set_from_hamilt(hamilt::HamiltLCAO<TK, TR>* hamilt_lcao, int nks, const Parallel_Orbitals* pv)
{
    hamilt_ptr_ = hamilt_lcao;
    pv_ = pv;
    nks_ = nks;

    if (hamilt_lcao == nullptr)
    {
        nbasis_ = 0;
        nloc_ = 0;
        nrow_ = 0;
        ncol_ = 0;
        return;
    }

    // Get dimensions from Parallel_Orbitals if available
    if (pv != nullptr)
    {
        nrow_ = pv->get_row_size();
        ncol_ = pv->get_col_size();
        nloc_ = nrow_ * ncol_;
        nbasis_ = pv->get_global_row_size();
    }

    // Initialize pointer arrays for compatibility mode
    hk_ptrs_.resize(nks, nullptr);
    sk_ptrs_.resize(nks, nullptr);
    matrix_dims_.resize(nks, {nrow_, ncol_});
}

template <typename TK, typename TR>
void PyHamiltonianAccessor<TK, TR>::set_dimensions(int nbasis, int nks)
{
    nbasis_ = nbasis;
    nks_ = nks;
    hk_ptrs_.resize(nks, nullptr);
    sk_ptrs_.resize(nks, nullptr);
    matrix_dims_.resize(nks, {0, 0});
}

template <typename TK, typename TR>
void PyHamiltonianAccessor<TK, TR>::set_Hk_data(int ik, const TK* data, int nrow, int ncol)
{
    if (ik >= 0 && ik < nks_)
    {
        hk_ptrs_[ik] = data;
        matrix_dims_[ik] = {nrow, ncol};
    }
}

template <typename TK, typename TR>
void PyHamiltonianAccessor<TK, TR>::set_Sk_data(int ik, const TK* data, int nrow, int ncol)
{
    if (ik >= 0 && ik < nks_)
    {
        sk_ptrs_[ik] = data;
        matrix_dims_[ik] = {nrow, ncol};
    }
}

template <typename TK, typename TR>
py::array_t<TK> PyHamiltonianAccessor<TK, TR>::get_Hk(int ik) const
{
    if (!is_valid() || ik < 0 || ik >= nks_)
    {
        throw std::runtime_error("Invalid k-point index or Hamiltonian not available.");
    }

    if (hk_ptrs_[ik] == nullptr)
    {
        throw std::runtime_error("H(k) data not set for this k-point.");
    }

    auto [nrow, ncol] = matrix_dims_[ik];
    std::vector<ssize_t> shape = {nrow, ncol};

    auto result = py::array_t<TK>(shape);
    auto buf = result.request();
    TK* ptr = static_cast<TK*>(buf.ptr);

    std::copy(hk_ptrs_[ik], hk_ptrs_[ik] + nrow * ncol, ptr);

    return result;
}

template <typename TK, typename TR>
py::array_t<TK> PyHamiltonianAccessor<TK, TR>::get_Sk(int ik) const
{
    if (!is_valid() || ik < 0 || ik >= nks_)
    {
        throw std::runtime_error("Invalid k-point index or overlap matrix not available.");
    }

    if (sk_ptrs_[ik] == nullptr)
    {
        throw std::runtime_error("S(k) data not set for this k-point.");
    }

    auto [nrow, ncol] = matrix_dims_[ik];
    std::vector<ssize_t> shape = {nrow, ncol};

    auto result = py::array_t<TK>(shape);
    auto buf = result.request();
    TK* ptr = static_cast<TK*>(buf.ptr);

    std::copy(sk_ptrs_[ik], sk_ptrs_[ik] + nrow * ncol, ptr);

    return result;
}

template <typename TK, typename TR>
py::dict PyHamiltonianAccessor<TK, TR>::get_HR() const
{
    // Placeholder: will be implemented when full ABACUS integration is available
    py::dict result;
    return result;
}

template <typename TK, typename TR>
py::dict PyHamiltonianAccessor<TK, TR>::get_SR() const
{
    // Placeholder: will be implemented when full ABACUS integration is available
    py::dict result;
    return result;
}

// Explicit template instantiations
template class PyHamiltonianAccessor<double, double>;
template class PyHamiltonianAccessor<std::complex<double>, double>;

// ============================================================================
// PyDensityMatrixAccessor Implementation (template)
// ============================================================================

template <typename TK, typename TR>
void PyDensityMatrixAccessor<TK, TR>::set_from_dm(elecstate::DensityMatrix<TK, TR>* dm)
{
    dm_ptr_ = dm;

    if (dm == nullptr)
    {
        nks_ = 0;
        nrow_ = 0;
        ncol_ = 0;
        return;
    }

    nks_ = dm->get_DMK_nks();
    nrow_ = dm->get_DMK_nrow();
    ncol_ = dm->get_DMK_ncol();

    // Initialize pointer arrays for compatibility mode
    dmk_ptrs_.resize(nks_, nullptr);
}

template <typename TK, typename TR>
void PyDensityMatrixAccessor<TK, TR>::set_dimensions(int nks, int nrow, int ncol)
{
    nks_ = nks;
    nrow_ = nrow;
    ncol_ = ncol;
    dmk_ptrs_.resize(nks, nullptr);
}

template <typename TK, typename TR>
void PyDensityMatrixAccessor<TK, TR>::set_DMK_data(int ik, const TK* data)
{
    if (ik >= 0 && ik < nks_)
    {
        dmk_ptrs_[ik] = data;
    }
}

template <typename TK, typename TR>
py::array_t<TK> PyDensityMatrixAccessor<TK, TR>::get_DMK(int ik) const
{
    if (!is_valid() || ik < 0 || ik >= nks_)
    {
        throw std::runtime_error("Invalid k-point index or density matrix not available.");
    }

    if (dmk_ptrs_[ik] == nullptr)
    {
        throw std::runtime_error("DM(k) data not set for this k-point.");
    }

    std::vector<ssize_t> shape = {static_cast<ssize_t>(nrow_), static_cast<ssize_t>(ncol_)};

    auto result = py::array_t<TK>(shape);
    auto buf = result.request();
    TK* ptr = static_cast<TK*>(buf.ptr);

    std::copy(dmk_ptrs_[ik], dmk_ptrs_[ik] + nrow_ * ncol_, ptr);

    return result;
}

template <typename TK, typename TR>
std::vector<py::array_t<TK>> PyDensityMatrixAccessor<TK, TR>::get_DMK_all() const
{
    std::vector<py::array_t<TK>> result;
    for (int ik = 0; ik < nks_; ++ik)
    {
        result.push_back(get_DMK(ik));
    }
    return result;
}

template <typename TK, typename TR>
py::dict PyDensityMatrixAccessor<TK, TR>::get_DMR() const
{
    // Placeholder: will be implemented when full ABACUS integration is available
    py::dict result;
    return result;
}

// Explicit template instantiations
template class PyDensityMatrixAccessor<double, double>;
template class PyDensityMatrixAccessor<std::complex<double>, double>;

// ============================================================================
// PyESolverLCAO Implementation (template)
// ============================================================================

template <typename TK, typename TR>
PyESolverLCAO<TK, TR>::PyESolverLCAO()
{
    // Constructor - initialization deferred to initialize()
}

template <typename TK, typename TR>
PyESolverLCAO<TK, TR>::~PyESolverLCAO()
{
    // Destructor - cleanup will be implemented in Phase 3
}

template <typename TK, typename TR>
void PyESolverLCAO<TK, TR>::initialize(const std::string& input_dir)
{
    // Placeholder: will be implemented in Phase 3
    // This will:
    // 1. Read INPUT file from input_dir
    // 2. Initialize UnitCell
    // 3. Create ESolver_KS_LCAO instance
    initialized_ = true;
    std::cout << "[PyESolverLCAO] Initialized with input directory: " << input_dir << std::endl;
}

template <typename TK, typename TR>
void PyESolverLCAO<TK, TR>::before_all_runners()
{
    if (!initialized_)
    {
        throw std::runtime_error("ESolver not initialized. Call initialize() first.");
    }
    // Placeholder: will call esolver_->before_all_runners() in Phase 3
    std::cout << "[PyESolverLCAO] before_all_runners called" << std::endl;
}

template <typename TK, typename TR>
void PyESolverLCAO<TK, TR>::before_scf(int istep)
{
    if (!initialized_)
    {
        throw std::runtime_error("ESolver not initialized. Call initialize() first.");
    }
    istep_ = istep;
    scf_started_ = true;
    conv_esolver_ = false;
    niter_ = 0;
    // Placeholder: will call esolver_->before_scf() in Phase 3
    std::cout << "[PyESolverLCAO] before_scf called for step " << istep << std::endl;
}

template <typename TK, typename TR>
void PyESolverLCAO<TK, TR>::run_scf_iteration(int iter)
{
    if (!scf_started_)
    {
        throw std::runtime_error("SCF not started. Call before_scf() first.");
    }
    niter_ = iter;
    // Placeholder: will implement actual SCF iteration in Phase 3
    // 1. iter_init()
    // 2. hamilt2rho()
    // 3. iter_finish()
    std::cout << "[PyESolverLCAO] SCF iteration " << iter << std::endl;
}

template <typename TK, typename TR>
void PyESolverLCAO<TK, TR>::run_scf(int max_iter)
{
    before_scf(istep_);

    for (int iter = 1; iter <= max_iter; ++iter)
    {
        run_scf_iteration(iter);
        if (conv_esolver_)
        {
            break;
        }
    }
}

template <typename TK, typename TR>
void PyESolverLCAO<TK, TR>::after_scf(int istep)
{
    if (!scf_started_)
    {
        throw std::runtime_error("SCF not started. Call before_scf() first.");
    }
    // Placeholder: will call esolver_->after_scf() in Phase 3
    std::cout << "[PyESolverLCAO] after_scf called for step " << istep << std::endl;
    scf_started_ = false;
}

template <typename TK, typename TR>
PyChargeAccessor PyESolverLCAO<TK, TR>::get_charge() const
{
    PyChargeAccessor accessor;
    // Note: esolver_ connection will be implemented when full ABACUS integration is available
    // For now, return empty accessor
    return accessor;
}

template <typename TK, typename TR>
PyEnergyAccessor PyESolverLCAO<TK, TR>::get_energy() const
{
    PyEnergyAccessor accessor;
    // Note: esolver_ connection will be implemented when full ABACUS integration is available
    // For now, return empty accessor
    return accessor;
}

template <typename TK, typename TR>
PyHamiltonianAccessor<TK, TR> PyESolverLCAO<TK, TR>::get_hamiltonian() const
{
    PyHamiltonianAccessor<TK, TR> accessor;
    // Note: esolver_ connection will be implemented when full ABACUS integration is available
    // For now, return empty accessor
    return accessor;
}

template <typename TK, typename TR>
PyDensityMatrixAccessor<TK, TR> PyESolverLCAO<TK, TR>::get_density_matrix() const
{
    PyDensityMatrixAccessor<TK, TR> accessor;
    // Note: esolver_ connection will be implemented when full ABACUS integration is available
    // For now, return empty accessor
    return accessor;
}

template <typename TK, typename TR>
py::array_t<TK> PyESolverLCAO<TK, TR>::get_psi(int ik) const
{
    // Note: Will return wave function coefficients when full ABACUS integration is available
    return py::array_t<TK>();
}

template <typename TK, typename TR>
py::array_t<double> PyESolverLCAO<TK, TR>::get_eigenvalues(int ik) const
{
    // Note: Will return eigenvalues when full ABACUS integration is available
    return py::array_t<double>();
}

template <typename TK, typename TR>
py::array_t<double> PyESolverLCAO<TK, TR>::get_occupations(int ik) const
{
    // Note: Will return occupation numbers when full ABACUS integration is available
    return py::array_t<double>();
}

template <typename TK, typename TR>
int PyESolverLCAO<TK, TR>::get_nks() const
{
    // Note: Will return actual nks when full ABACUS integration is available
    return 0;
}

template <typename TK, typename TR>
py::array_t<double> PyESolverLCAO<TK, TR>::get_kvec_d(int ik) const
{
    std::vector<ssize_t> shape = {3};
    auto result = py::array_t<double>(shape);
    auto buf = result.request();
    double* ptr = static_cast<double*>(buf.ptr);
    ptr[0] = ptr[1] = ptr[2] = 0.0;
    return result;
}

template <typename TK, typename TR>
py::array_t<double> PyESolverLCAO<TK, TR>::get_wk() const
{
    // Note: Will return k-point weights when full ABACUS integration is available
    return py::array_t<double>();
}

template <typename TK, typename TR>
int PyESolverLCAO<TK, TR>::get_nbasis() const
{
    // Note: Will return actual nbasis when full ABACUS integration is available
    return 0;
}

template <typename TK, typename TR>
int PyESolverLCAO<TK, TR>::get_nbands() const
{
    // Note: Will return actual nbands when full ABACUS integration is available
    return 0;
}

template <typename TK, typename TR>
int PyESolverLCAO<TK, TR>::get_nspin() const
{
    // Note: Will return actual nspin when full ABACUS integration is available
    return 1;
}

template <typename TK, typename TR>
int PyESolverLCAO<TK, TR>::get_nat() const
{
    // Note: Will return actual nat when full ABACUS integration is available
    return 0;
}

// Explicit template instantiations
template class PyESolverLCAO<double, double>;
template class PyESolverLCAO<std::complex<double>, double>;

} // namespace py_esolver

// ============================================================================
// Pybind11 Module Definition
// ============================================================================

void bind_charge_accessor(py::module& m)
{
    py::class_<py_esolver::PyChargeAccessor>(m, "ChargeAccessor",
        R"pbdoc(
        Accessor for charge density data.

        Provides access to real-space charge density (rho) and related quantities.
        )pbdoc")
        .def(py::init<>())
        .def("get_rho", &py_esolver::PyChargeAccessor::get_rho,
            R"pbdoc(
            Get real-space charge density as numpy array.

            Returns
            -------
            numpy.ndarray
                Charge density with shape (nspin, nrxx)
            )pbdoc")
        .def("get_rhog", &py_esolver::PyChargeAccessor::get_rhog,
            R"pbdoc(
            Get reciprocal-space charge density as numpy array.

            Returns
            -------
            numpy.ndarray
                Charge density in G-space with shape (nspin, ngmc)
            )pbdoc")
        .def("get_rho_core", &py_esolver::PyChargeAccessor::get_rho_core,
            R"pbdoc(
            Get core charge density as numpy array.

            Returns
            -------
            numpy.ndarray
                Core charge density with shape (nrxx,)
            )pbdoc")
        .def_property_readonly("nspin", &py_esolver::PyChargeAccessor::get_nspin,
            "Number of spin channels")
        .def_property_readonly("nrxx", &py_esolver::PyChargeAccessor::get_nrxx,
            "Number of real-space grid points")
        .def_property_readonly("ngmc", &py_esolver::PyChargeAccessor::get_ngmc,
            "Number of G-vectors for charge density")
        .def("is_valid", &py_esolver::PyChargeAccessor::is_valid,
            "Check if charge data is available");
}

void bind_energy_accessor(py::module& m)
{
    py::class_<py_esolver::PyEnergyAccessor>(m, "EnergyAccessor",
        R"pbdoc(
        Accessor for energy data.

        Provides access to various energy components from the calculation.
        All energies are in Rydberg units.
        )pbdoc")
        .def(py::init<>())
        .def_property_readonly("etot", &py_esolver::PyEnergyAccessor::get_etot,
            "Total energy (Ry)")
        .def_property_readonly("eband", &py_esolver::PyEnergyAccessor::get_eband,
            "Band energy (Ry)")
        .def_property_readonly("hartree_energy", &py_esolver::PyEnergyAccessor::get_hartree_energy,
            "Hartree energy (Ry)")
        .def_property_readonly("etxc", &py_esolver::PyEnergyAccessor::get_etxc,
            "Exchange-correlation energy (Ry)")
        .def_property_readonly("ewald_energy", &py_esolver::PyEnergyAccessor::get_ewald_energy,
            "Ewald energy (Ry)")
        .def_property_readonly("demet", &py_esolver::PyEnergyAccessor::get_demet,
            "-TS term for metals (Ry)")
        .def_property_readonly("exx", &py_esolver::PyEnergyAccessor::get_exx,
            "Exact exchange energy (Ry)")
        .def_property_readonly("evdw", &py_esolver::PyEnergyAccessor::get_evdw,
            "van der Waals energy (Ry)")
        .def("get_all_energies", &py_esolver::PyEnergyAccessor::get_all_energies,
            "Get all energies as a dictionary");
}

template <typename TK>
void bind_hamiltonian_accessor(py::module& m, const std::string& suffix)
{
    using HamiltAccessor = py_esolver::PyHamiltonianAccessor<TK>;

    std::string class_name = "HamiltonianAccessor" + suffix;

    py::class_<HamiltAccessor>(m, class_name.c_str(),
        R"pbdoc(
        Accessor for Hamiltonian matrix data.

        Provides access to H(k), S(k), H(R), and S(R) matrices.
        )pbdoc")
        .def(py::init<>())
        .def_property_readonly("nbasis", &HamiltAccessor::get_nbasis,
            "Number of basis functions")
        .def_property_readonly("nks", &HamiltAccessor::get_nks,
            "Number of k-points")
        .def("get_Hk", &HamiltAccessor::get_Hk,
            R"pbdoc(
            Get H(k) matrix for specific k-point.

            Parameters
            ----------
            ik : int
                K-point index

            Returns
            -------
            numpy.ndarray
                Hamiltonian matrix at k-point ik
            )pbdoc", "ik"_a)
        .def("get_Sk", &HamiltAccessor::get_Sk,
            R"pbdoc(
            Get S(k) overlap matrix for specific k-point.

            Parameters
            ----------
            ik : int
                K-point index

            Returns
            -------
            numpy.ndarray
                Overlap matrix at k-point ik
            )pbdoc", "ik"_a)
        .def("get_HR", &HamiltAccessor::get_HR,
            "Get H(R) in sparse format")
        .def("get_SR", &HamiltAccessor::get_SR,
            "Get S(R) in sparse format")
        .def("is_valid", &HamiltAccessor::is_valid,
            "Check if Hamiltonian data is available");
}

template <typename TK>
void bind_density_matrix_accessor(py::module& m, const std::string& suffix)
{
    using DMAccessor = py_esolver::PyDensityMatrixAccessor<TK>;

    std::string class_name = "DensityMatrixAccessor" + suffix;

    py::class_<DMAccessor>(m, class_name.c_str(),
        R"pbdoc(
        Accessor for density matrix data.

        Provides access to DM(k) and DM(R) matrices.
        )pbdoc")
        .def(py::init<>())
        .def_property_readonly("nks", &DMAccessor::get_nks,
            "Number of k-points")
        .def_property_readonly("nrow", &DMAccessor::get_nrow,
            "Number of rows in density matrix")
        .def_property_readonly("ncol", &DMAccessor::get_ncol,
            "Number of columns in density matrix")
        .def("get_DMK", &DMAccessor::get_DMK,
            R"pbdoc(
            Get DM(k) for specific k-point.

            Parameters
            ----------
            ik : int
                K-point index

            Returns
            -------
            numpy.ndarray
                Density matrix at k-point ik
            )pbdoc", "ik"_a)
        .def("get_DMK_all", &DMAccessor::get_DMK_all,
            "Get all DM(k) matrices as a list")
        .def("get_DMR", &DMAccessor::get_DMR,
            "Get DM(R) in sparse format")
        .def("is_valid", &DMAccessor::is_valid,
            "Check if density matrix data is available");
}

template <typename TK, typename TR>
void bind_esolver_lcao(py::module& m, const std::string& suffix)
{
    using ESolver = py_esolver::PyESolverLCAO<TK, TR>;

    std::string class_name = "ESolverLCAO" + suffix;

    py::class_<ESolver>(m, class_name.c_str(),
        R"pbdoc(
        Python wrapper for ESolver_KS_LCAO.

        This class provides a Python interface for LCAO calculations
        with support for breakpoints and state inspection during SCF.

        Example
        -------
        >>> esolver = ESolverLCAO_gamma()
        >>> esolver.initialize("./")
        >>> esolver.before_all_runners()
        >>> esolver.before_scf(0)
        >>> for iter in range(1, 101):
        ...     esolver.run_scf_iteration(iter)
        ...     energy = esolver.get_energy()
        ...     print(f"Iter {iter}: E = {energy.etot}")
        ...     if esolver.is_converged():
        ...         break
        >>> # Breakpoint before after_scf - inspect state here
        >>> charge = esolver.get_charge()
        >>> hamiltonian = esolver.get_hamiltonian()
        >>> esolver.after_scf(0)
        )pbdoc")
        .def(py::init<>())

        // Initialization
        .def("initialize", &ESolver::initialize,
            R"pbdoc(
            Initialize ESolver from INPUT file.

            Parameters
            ----------
            input_dir : str
                Directory containing INPUT, STRU, and other input files
            )pbdoc", "input_dir"_a)
        .def("before_all_runners", &ESolver::before_all_runners,
            "Initialize calculation environment")

        // SCF Control
        .def("before_scf", &ESolver::before_scf,
            R"pbdoc(
            Prepare for SCF calculation.

            Parameters
            ----------
            istep : int, optional
                Ion step index (default: 0)
            )pbdoc", "istep"_a = 0)
        .def("run_scf_iteration", &ESolver::run_scf_iteration,
            R"pbdoc(
            Run a single SCF iteration.

            Parameters
            ----------
            iter : int
                Iteration number (1-based)
            )pbdoc", "iter"_a)
        .def("run_scf", &ESolver::run_scf,
            R"pbdoc(
            Run complete SCF loop.

            Parameters
            ----------
            max_iter : int, optional
                Maximum number of iterations (default: 100)
            )pbdoc", "max_iter"_a = 100)
        .def("after_scf", &ESolver::after_scf,
            R"pbdoc(
            Finalize SCF calculation.

            Parameters
            ----------
            istep : int, optional
                Ion step index (default: 0)
            )pbdoc", "istep"_a = 0)

        // Status
        .def("is_converged", &ESolver::is_converged,
            "Check if SCF is converged")
        .def_property_readonly("niter", &ESolver::get_niter,
            "Current iteration number")
        .def_property_readonly("drho", &ESolver::get_drho,
            "Charge density difference")
        .def_property_readonly("istep", &ESolver::get_istep,
            "Current ion step")

        // Data Accessors
        .def("get_charge", &ESolver::get_charge,
            "Get charge density accessor")
        .def("get_energy", &ESolver::get_energy,
            "Get energy accessor")
        .def("get_hamiltonian", &ESolver::get_hamiltonian,
            "Get Hamiltonian accessor")
        .def("get_density_matrix", &ESolver::get_density_matrix,
            "Get density matrix accessor")

        // Wave functions
        .def("get_psi", &ESolver::get_psi,
            "Get wave function coefficients for k-point ik", "ik"_a)
        .def("get_eigenvalues", &ESolver::get_eigenvalues,
            "Get eigenvalues for k-point ik", "ik"_a)
        .def("get_occupations", &ESolver::get_occupations,
            "Get occupation numbers for k-point ik", "ik"_a)

        // K-points
        .def_property_readonly("nks", &ESolver::get_nks,
            "Number of k-points")
        .def("get_kvec_d", &ESolver::get_kvec_d,
            "Get k-vector in direct coordinates", "ik"_a)
        .def("get_wk", &ESolver::get_wk,
            "Get k-point weights")

        // System info
        .def_property_readonly("nbasis", &ESolver::get_nbasis,
            "Number of basis functions")
        .def_property_readonly("nbands", &ESolver::get_nbands,
            "Number of bands")
        .def_property_readonly("nspin", &ESolver::get_nspin,
            "Number of spin channels")
        .def_property_readonly("nat", &ESolver::get_nat,
            "Number of atoms");
}

PYBIND11_MODULE(_esolver_pack, m)
{
    m.doc() = R"pbdoc(
        PyABACUS ESolver Module
        -----------------------

        This module provides Python bindings for ABACUS ESolver_KS_LCAO,
        enabling Python-controlled SCF workflows with breakpoint support.

        Main Classes
        ------------
        ESolverLCAO_gamma : ESolver for gamma-only calculations
        ESolverLCAO_multi_k : ESolver for multi-k calculations

        Accessor Classes
        ----------------
        ChargeAccessor : Access charge density data
        EnergyAccessor : Access energy components
        HamiltonianAccessor_gamma/multi_k : Access Hamiltonian matrices
        DensityMatrixAccessor_gamma/multi_k : Access density matrices

        Example
        -------
        >>> from pyabacus.esolver import ESolverLCAO_gamma
        >>> esolver = ESolverLCAO_gamma()
        >>> esolver.initialize("./")
        >>> esolver.before_all_runners()
        >>> esolver.before_scf(0)
        >>> # Run SCF with breakpoint support
        >>> for iter in range(1, 101):
        ...     esolver.run_scf_iteration(iter)
        ...     if esolver.is_converged():
        ...         break
        >>> # Inspect state before after_scf
        >>> charge = esolver.get_charge()
        >>> energy = esolver.get_energy()
        >>> esolver.after_scf(0)
    )pbdoc";

    // Bind accessor classes
    bind_charge_accessor(m);
    bind_energy_accessor(m);
    bind_hamiltonian_accessor<double>(m, "_gamma");
    bind_hamiltonian_accessor<std::complex<double>>(m, "_multi_k");
    bind_density_matrix_accessor<double>(m, "_gamma");
    bind_density_matrix_accessor<std::complex<double>>(m, "_multi_k");

    // Bind ESolver classes
    bind_esolver_lcao<double, double>(m, "_gamma");
    bind_esolver_lcao<std::complex<double>, double>(m, "_multi_k");
}
