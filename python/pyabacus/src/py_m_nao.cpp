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
        .def(
            "__call__",
            [](RadialCollection& self, const int itype, const int l, const int izeta) -> const NumericalRadial& {
                return self(itype, l, izeta);
            },
            py::return_value_policy::reference_internal,
            "itype"_a,
            "l"_a,
            "izeta"_a)
        // Getters
        .def("symbol", &RadialCollection::symbol, "itype"_a)
        .def_property_readonly("ntype", &RadialCollection::ntype)
        .def("lmax", overload_cast_<const int>()(&RadialCollection::lmax, py::const_), "itype"_a)
        .def_property_readonly("lmax", overload_cast_<>()(&RadialCollection::lmax, py::const_))
        .def("rcut_max", overload_cast_<const int>()(&RadialCollection::rcut_max, py::const_), "itype"_a)
        .def_property_readonly("rcut_max", overload_cast_<>()(&RadialCollection::rcut_max, py::const_))
        .def("nzeta", &RadialCollection::nzeta, "itype"_a, "l"_a)
        .def("nzeta_max", overload_cast_<const int>()(&RadialCollection::nzeta_max, py::const_), "itype"_a)
        .def("nzeta_max", overload_cast_<>()(&RadialCollection::nzeta_max, py::const_))
        .def("nchi", overload_cast_<const int>()(&RadialCollection::nchi, py::const_), "itype"_a)
        .def("nchi", overload_cast_<>()(&RadialCollection::nchi, py::const_));
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
        .def(py::init<>(), R"pbdoc(
    A class that represents a numerical radial function.

    This class is designed to be the container for the radial part of numerical atomic orbitals, Kleinman-Bylander beta functions, and all other similar numerical radial functions in three-dimensional space, each of which is associated with some angular momentum l and whose r and k space values are related by an l-th order spherical Bessel transform.

    A NumericalRadial object can be initialized by "build", which requires the angular momentum, the number of grid points, the grid and the corresponding values. Grid does not have to be uniform. One can initialize the object in either r or k space. After initialization, one can set the
    grid in the other space via set_grid or set_uniform_grid. Values in the other space are automatically computed by a spherical Bessel transform.
        )pbdoc")
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
            R"pbdoc(
    Initializes the object by providing the grid & values in one space.

    Parameters
    ----------
    l : int
        Angular momentum.
    for_r_space : bool
        Specifies whether the input corresponds to r or k space.
    ngrid : int
        Number of grid points.
    grid : array_like
        Grid points, must be positive & strictly increasing.
    value : array_like
        Values on the grid.
    p : float
        Implicit exponent in input values, see pr_ & pk_.
    izeta : int
        Multiplicity index of radial functions of the same itype and l.
    symbol : str
        Chemical symbol.
    itype : int
        Index for the element in calculation.
    init_sbt : bool
        If true, internal SphericalBesselTransformer will be initialized.

    Notes
    -----
    init_sbt is only useful when the internal SphericalBesselTransformer (sbt_) is null-initialized; The function will NOT reset sbt_ if it is already usable.
            )pbdoc",
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
        .def("set_transformer",
             &NumericalRadial::set_transformer,
             R"pbdoc(
    Sets a SphericalBesselTransformer.

    By default, the class uses an internal SphericalBesselTransformer, but one can optionally use a shared one. This could be beneficial when there are a lot of NumericalRadial objects whose grids have the same size.

    Parameters
    ----------
    sbt : SphericalBesselTransformer
        An external transformer.
    update : int
        Specifies whether and how values are recomputed with the new transformer.
        Accepted values are:
        *  0: does not recompute values;
        *  1: calls a forward transform;
        * -1: calls a backward transform.
        )pbdoc",
             "sbt"_a,
             "update"_a = 0)
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
            R"pbdoc(
    Sets up a grid.

    This function can be used to set up the grid which is absent in "build" (in which case values on the new grid are automatically computed by a spherical Bessel transform) or interpolate on an existing grid to a new grid.

    Parameters
    ----------
    for_r_space : bool
        Specifies whether to set grid for the r or k space.
    ngrid : int
        Number of grid points.
    grid : array_like
        Grid points, must be positive & strictly increasing.
    mode : char
        Specifies how values are updated, could be 'i' or 't':
        - 'i': New values are obtained by interpolating and zero-padding
               the existing values from current space. With this option,
               it is an error if the designated space does not have a grid;
        - 't': New values are obtained via transform from the other space.
               With this option, it is an error if the other space does not
               have a grid.
            )pbdoc",
            "for_r_space"_a,
            "ngrid"_a,
            "grid"_a,
            "mode"_a = 'i')
        .def("set_uniform_grid",
             &NumericalRadial::set_uniform_grid,
             R"pbdoc(
    Sets up a uniform grid.

    The functionality of this function is similar to set_grid, except that the new grid is a uniform grid specified by the cutoff and the number of grid points, which are calculated as:

        grid[i] = i * (cutoff / (ngrid - 1))

    Parameters
    ----------
    for_r_space : bool
        Specifies whether to set grid for the r or k space.
    ngrid : int
        Number of grid points.
    cutoff : float
        The maximum value of the grid, which determines the grid spacing.
    enable_fft : bool
        If true, this function will not only set up the grid & values in the designated space, but also sets the grid in the other space such that the r & k grids are FFT-compliant (and updates values via a FFT-based spherical Bessel transform).
     mode : char
        Specifies how values are updated, could be 'i' or 't'.
            )pbdoc",
             "for_r_space"_a,
             "ngrid"_a,
             "cutoff"_a,
             "mode"_a = 'i',
             "enable_fft"_a = false)
        .def(
            "set_value",
            [](NumericalRadial& self, const bool for_r_space, py::array_t<double> value, const int p) {
                py::buffer_info value_info = value.request();
                self.set_value(for_r_space, static_cast<double*>(value_info.ptr), p);
            },
            R"pbdoc(
    Updates values on an existing grid.

    This function does not alter the grid; it merely updates values on the existing grid. The number of values to read from "value" is determined by the current number of points in the r or k space (nr_ or nk_). Values of the other space will be automatically updated if they exist.

    Warning
    -------
    This function does not check the index bound; use with care!
                )pbdoc",
            "for_r_space"_a,
            "value"_a,
            "p"_a)
        .def("wipe", &NumericalRadial::wipe, "r_space"_a = true, "k_space"_a = true)
        .def("normalize",
             &NumericalRadial::normalize,
             R"pbdoc(
    Normalizes the radial function.

    The radial function is normalized such that the integral of the square of the function multiplied by the square of the radial coordinate over the entire space is equal to one:

        ∫ from 0 to +∞ of (x^2 * f(x)^2) dx = 1

    where x is r or k. The integral is evaluated with Simpson's rule. Values in the other space are updated automatically via a spherical Bessel transform.
    )pbdoc",
             "for_r_space"_a = true)
        // Getters
        .def_property_readonly("symbol", &NumericalRadial::symbol)
        .def_property_readonly("itype", &NumericalRadial::itype)
        .def_property_readonly("izeta", &NumericalRadial::izeta)
        .def_property_readonly("l", &NumericalRadial::l)
        .def_property_readonly("nr", &NumericalRadial::nr)
        .def_property_readonly("nk", &NumericalRadial::nk)
        .def_property_readonly("rcut", &NumericalRadial::rcut)
        .def_property_readonly("kcut", &NumericalRadial::kcut)
        .def_property_readonly("rmax", &NumericalRadial::rmax)
        .def_property_readonly("kmax", &NumericalRadial::kmax)
        .def_property_readonly("pr", &NumericalRadial::pr)
        .def_property_readonly("pk", &NumericalRadial::pk)
        .def_property_readonly("sbt", &NumericalRadial::sbt)
        .def_property_readonly("rgrid",
                               [](NumericalRadial& self) {
                                   const double* rgrid = self.rgrid();
                                   return py::array_t<double>({self.nr()}, rgrid);
                               })
        .def_property_readonly("kgrid",
                               [](NumericalRadial& self) {
                                   const double* kgrid = self.kgrid();
                                   return py::array_t<double>({self.nk()}, kgrid);
                               })
        .def_property_readonly("rvalue",
                               [](NumericalRadial& self) {
                                   const double* rvalue = self.rvalue();
                                   return py::array_t<double>({self.nr()}, rvalue);
                               })
        .def_property_readonly("kvalue",
                               [](NumericalRadial& self) {
                                   const double* kvalue = self.kvalue();
                                   return py::array_t<double>({self.nk()}, kvalue);
                               })
        .def_property_readonly("is_fft_compliant", overload_cast_<>()(&NumericalRadial::is_fft_compliant, py::const_));
}