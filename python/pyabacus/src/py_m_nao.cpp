#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "module_base/vector3.h"
#include "module_basis/module_nao/radial_collection.h"
#include "module_basis/module_nao/two_center_integrator.h"

namespace py = pybind11;
using namespace pybind11::literals;
template <typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;

void bind_m_nao(py::module& m)
{
    // Create the submodule for Module NAO
    py::module m_nao = m.def_submodule("ModuleNAO");
    m_nao.doc() = "Module for Numerical Atomic Orbitals (NAO) in ABACUS";

    // Bind the RadialCollection class
    py::class_<RadialCollection>(m_nao, "RadialCollection")
        .def(py::init<>(), R"pbdoc(
            A class that holds all numerical radial functions of the same kind. 
            
            An instance of this class could be the collection of all radial functions 
            of numerical atomic orbitals, or all Kleinman-Bylander beta functions from 
            all elements involved in a calculation.
            )pbdoc")
        .def(
            "build",
            [](RadialCollection& self, int nfile, const py::list& file_list, char ftype) {
                std::vector<std::string> files;
                files.reserve(nfile);
                for (auto file: file_list)
                {
                    files.push_back(file.cast<std::string>());
                }
                self.build(nfile, files.data(), ftype);
            },
            "Builds the collection from (orbital) files",
            "nfile"_a,
            "file_list"_a,
            "ftype"_a = '\0')
        .def("set_transformer",
             &RadialCollection::set_transformer,
             "Sets a spherical Bessel transformers for all RadialSet objects.",
             "sbt"_a,
             "update"_a = 0)
        .def("set_uniform_grid",
             &RadialCollection::set_uniform_grid,
             "Sets a common uniform grid for all RadialSet objects",
             "for_r_space"_a,
             "ngrid"_a,
             "cutoff"_a,
             "mode"_a = 'i',
             "enable_fft"_a = false)
        .def(
            "set_grid",
            [](RadialCollection& self,
               const bool for_r_space,
               const int ngrid,
               py::array_t<double> grid,
               const char mode = 'i') {
                py::buffer_info grid_info = grid.request();
                if (grid_info.size != ngrid)
                {
                    throw std::runtime_error("grid array must be of size ngrid");
                }
                self.set_grid(for_r_space, ngrid, static_cast<double*>(grid_info.ptr), mode);
            },
            "Sets a common grid for all RadialSet objects",
            "for_r_space"_a,
            "ngrid"_a,
            "grid"_a,
            "mode"_a = 'i')
        // Getters
        .def("symbol", &RadialCollection::symbol, "itype"_a)
        .def_property_readonly("ntype", &RadialCollection::ntype)
        .def("lmax", overload_cast_<const int>()(&RadialCollection::lmax, py::const_), "itype"_a)
        .def_property_readonly("lmax", overload_cast_<>()(&RadialCollection::lmax, py::const_))
        .def("rcut_max", overload_cast_<const int>()(&RadialCollection::rcut_max, py::const_), "itype"_a)
        .def_property_readonly("rcut_max", overload_cast_<>()(&RadialCollection::rcut_max, py::const_))
        .def("nzeta", &RadialCollection::nzeta, "itype"_a, "l"_a)
        .def("nzeta_max", overload_cast_<const int>()(&RadialCollection::nzeta_max, py::const_), "itype"_a)
        .def_property_readonly("nzeta_max", overload_cast_<>()(&RadialCollection::nzeta_max, py::const_))
        .def("nchi", overload_cast_<const int>()(&RadialCollection::nchi, py::const_), "itype"_a)
        .def_property_readonly("nchi", overload_cast_<>()(&RadialCollection::nchi, py::const_));
    // Bind the TwoCenterIntegrator class
    py::class_<TwoCenterIntegrator>(m_nao, "TwoCenterIntegrator")
        .def(py::init<>(), R"pbdoc(
    A class to compute two-center integrals.

    This class computes two-center integrals of the form:

                        /    
                 I(R) = | dr phi1(r) (op) phi2(r - R)
                        /               

    as well as their gradients, where op is 1 (overlap) or minus Laplacian (kinetic), and phi1, 
    phi2 are "atomic-orbital-like" functions of the form:

                 phi(r) = chi(|r|) * Ylm(r/|r|)

    where chi is some numerical radial function and Ylm is some real spherical harmonics.

    This class is designed to efficiently compute the two-center integrals
    between two "collections" of the above functions with various R, e.g., the
    overlap integrals between all numerical atomic orbitals and all
    Kleinman-Bylander nonlocal projectors, the overlap & kinetic integrals between all numerical atomic orbitals, etc.
    This is done by tabulating the radial part of the integrals on an r-space grid and the real Gaunt coefficients in advance.
    )pbdoc")
        .def("tabulate",
             &TwoCenterIntegrator::tabulate,
             R"pbdoc(
        Tabulates the radial part of a two-center integral.

        Parameters:
        bra (RadialFunctions): The radial functions of the first collection.
        ket (RadialFunctions): The radial functions of the second collection.
        op (char): Operator, could be 'S' (overlap) or 'T' (kinetic).
        nr (int): Number of r-space grid points.
        cutoff (float): r-space cutoff radius.
        )pbdoc",
             "bra"_a,
             "ket"_a,
             "op"_a,
             "nr"_a,
             "cutoff"_a)
        .def(
            "calculate",
            [](TwoCenterIntegrator& self,
               const int itype1,
               const int l1,
               const int izeta1,
               const int m1,
               const int itype2,
               const int l2,
               const int izeta2,
               const int m2,
               py::array_t<double> pvR,
               bool cal_grad = false) {
                py::buffer_info pvR_info = pvR.request();
                if (pvR_info.size != 3)
                {
                    throw std::runtime_error("Radial part must have 3 elements");
                }
                double* cvR = static_cast<double*>(pvR_info.ptr);
                ModuleBase::Vector3<double> vR(cvR[0], cvR[1], cvR[2]);
                double out[1] = {0.0};
                double* grad_out = nullptr;
                if (cal_grad)
                {
                    grad_out = new double[3];
                }
                self.calculate(itype1, l1, izeta1, m1, itype2, l2, izeta2, m2, vR, out, grad_out);
                py::array_t<double> out_array({1}, out);
                if (cal_grad)
                {
                    py::array_t<double> grad_out_array({3}, grad_out);
                    return py::make_tuple(out_array, grad_out_array);
                }
                else
                {
                    py::array_t<double> grad_out_array({0});
                    return py::make_tuple(out_array, grad_out_array);
                }
            },
            R"pbdoc(
    Compute the two-center integrals.

    This function calculates the two-center integral

                        /    
                 I(R) = | dr phi1(r) (op_) phi2(r - R)
                        /               

    or its gradient by using the tabulated radial part and real Gaunt coefficients.

    Parameters
    ----------
    itype1 : int
        Element index of orbital 1.
    l1 : int
        Angular momentum of orbital 1.
    izeta1 : int
        Zeta number of orbital 1.
    m1 : int
        Magnetic quantum number of orbital 1.
    itype2 : int
        Element index of orbital 2.
    l2 : int
        Angular momentum of orbital 2.
    izeta2 : int
        Zeta number of orbital 2.
    m2 : int
        Magnetic quantum number of orbital 2.
    pvR : array_like
        R2 - R1, the displacement vector between the two centers.
    cal_grad : bool, optional
        The gradient will not be computed if cal_grad is false.
    
    Returns
    -------
    out_array : array_like
        The two-center integral.
    grad_out_array : array_like
        Gradient of the two-center integral.
        )pbdoc",
            "itype1"_a,
            "l1"_a,
            "izeta1"_a,
            "m1"_a,
            "itype2"_a,
            "l2"_a,
            "izeta2"_a,
            "m2"_a,
            "pvR"_a,
            "cal_grad"_a = false)
        .def(
            "snap",
            [](TwoCenterIntegrator& self,
               const int itype1,
               const int l1,
               const int izeta1,
               const int m1,
               const int itype2,
               py::array_t<double> pvR,
               const bool deriv) {
                py::buffer_info pvR_info = pvR.request();
                if (pvR_info.size != 3)
                {
                    throw std::runtime_error("Radial part must have 3 elements");
                }
                double* cvR = static_cast<double*>(pvR_info.ptr);
                ModuleBase::Vector3<double> vR(cvR[0], cvR[1], cvR[2]);
                // TODO: check deriv & out memory allocation
                std::vector<std::vector<double>> out;
                self.snap(itype1, l1, izeta1, m1, itype2, vR, deriv, out);
                return out;
            },
            R"pbdoc(
    Compute a batch of two-center integrals.

    This function calculates the two-center integrals (and optionally their gradients)
    between one orbital and all orbitals of a certain type from the other collection.
        )pbdoc",
            "itype1"_a,
            "l1"_a,
            "izeta1"_a,
            "m1"_a,
            "itype2"_a,
            "pvR"_a,
            "deriv"_a = false);
    // Bind the NumericalRadial class
    py::class_<NumericalRadial>(m_nao, "NumericalRadial")
        .def(py::init<>())
        .def(
            "build",
            [](NumericalRadial& self,
               const int l,
               const bool for_r_space,
               const int ngrid,
               py::array_t<double> grid,
               py::array_t<double> value,
               const int p = 0,
               const int izeta = 0,
               const std::string symbol = "",
               const int itype = 0,
               const bool init_sbt = true) {
                py::buffer_info grid_info = grid.request();
                py::buffer_info value_info = value.request();
                if (grid_info.size != ngrid)
                {
                    throw std::runtime_error("grid array must be of size ngrid");
                }
                self.build(l,
                           for_r_space,
                           ngrid,
                           static_cast<double*>(grid_info.ptr),
                           static_cast<double*>(value_info.ptr),
                           p,
                           izeta,
                           symbol,
                           itype,
                           init_sbt);
            },
            "l"_a,
            "for_r_space"_a,
            "ngrid"_a,
            "grid"_a,
            "value"_a,
            "p"_a = 0,
            "izeta"_a = 0,
            "symbol"_a = "",
            "itype"_a = 0,
            "init_sbt"_a = true)
        .def("set_transformer", &NumericalRadial::set_transformer, "sbt"_a, "update"_a = 0)
        .def(
            "set_grid",
            [](NumericalRadial& self,
               const bool for_r_space,
               const int ngrid,
               py::array_t<double> grid,
               const char mode = 'i') {
                py::buffer_info grid_info = grid.request();
                if (grid_info.size != ngrid)
                {
                    throw std::runtime_error("grid array must be of size ngrid");
                }
                self.set_grid(for_r_space, ngrid, static_cast<double*>(grid_info.ptr), mode);
            },
            "for_r_space"_a,
            "ngrid"_a,
            "grid"_a,
            "mode"_a = 'i')
        .def("set_uniform_grid",
             &NumericalRadial::set_uniform_grid,
             "for_r_space"_a,
             "ngrid"_a,
             "cutoff"_a,
             "mode"_a = 'i',
             "enable_fft"_a = false)
        .def(
            "set_value",
            [](&NumericalRadial self, const bool for_r_space, py::array_t<double> value, const int p) {
                py::buffer_info value_info = value.request();
                self.set_value(for_r_space, static_cast<double*>(value_info.ptr), p);
            },
            "for_r_space"_a,
            "value"_a,
            "p"_a)
        .def("wipe", &NumericalRadial::wipe, "r_space"_a = true, "k_space"_a = true)
        .def("normalize", &NumericalRadial::normalize, "for_r_space"_a = true)
        // Getters
        .def_property_readonly("symbol", &NumericalRadial::symbol)
        .def_property_readonly("itype", &NumericalRadial::itype)
        .def_property_readonly("l", &NumericalRadial::l)
        .def_property_readonly("nr", &NumericalRadial::nr)
        .def_property_readonly("nk", &NumericalRadial::nk)
        .def_property_readonly("rcut", &NumericalRadial::rcut)
        .def_property_readonly("kcut", &NumericalRadial::kcut)
        .def_property_readonly("rmax", &NumericalRadial::rmax)
        .def_property_readonly("kmax", &NumericalRadial::kmax)
        .def_property_readonly("rgrid", overload_cast_<>()(&NumericalRadial::rgrid, py::const_))
        .def_property_readonly("kgrid", overload_cast_<>()(&NumericalRadial::kgrid, py::const_))
        .def_property_readonly("rvalue", overload_cast_<>()(&NumericalRadial::rvalue, py::const_))
        .def_property_readonly("kvalue", overload_cast_<>()(&NumericalRadial::kvalue, py::const_))
        .def_property_readonly("pr", &NumericalRadial::pr)
        .def_property_readonly("pk", &NumericalRadial::pk)
        .def_property_readonly("is_fft_compliant", &NumericalRadial::is_fft_compliant)
        .def_property_readonly("sbt", &NumericalRadial::sbt)
        .def("rgrid", overload_cast_<const int>()(&NumericalRadial::rgrid, py::const_), "ir"_a)
        .def("kgrid", overload_cast_<const int>()(&NumericalRadial::kgrid, py::const_), "ik"_a)
        .def("rvalue", overload_cast_<const int>()(&NumericalRadial::rvalue, py::const_), "ir"_a)
        .def("kvalue", overload_cast_<const int>()(&NumericalRadial::kvalue, py::const_), "ik"_a);
}