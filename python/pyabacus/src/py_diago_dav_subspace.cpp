#include <complex>
#include <functional>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>

#include "module_hsolver/diago_dav_subspace.h"
#include "module_hsolver/kernels/math_kernel_op.h"
#include "module_base/module_device/types.h"

namespace py = pybind11;
using namespace pybind11::literals;

void bind_diago_dav_subspace(py::module& m)
{
    py::module hsolver = m.def_submodule("hsolver");

    py::class_<hsolver::diag_comm_info>(hsolver, "diag_comm_info")
        .def(py::init<const int, const int>(), "rank"_a, "nproc"_a)
        .def_readonly("rank", &hsolver::diag_comm_info::rank)
        .def_readonly("nproc", &hsolver::diag_comm_info::nproc);

    py::class_<hsolver::Diago_DavSubspace<std::complex<double>, base_device::DEVICE_CPU>>(hsolver, "Diago_DavSubspace")
        .def(py::init<const std::vector<double>, const int, const int,
                    const int, const double, const int, const bool,
                    const hsolver::diag_comm_info>(),
        "precondition"_a, "nband"_a, "nbasis"_a, "david_ndim"_a, "diag_thr"_a, "diag_nmax"_a, "need_subspace"_a, "diag_comm"_a)
        .def("diag",
            [] (hsolver::Diago_DavSubspace<std::complex<double>, base_device::DEVICE_CPU>& self,
                py::array_t<std::complex<double>> h_mat,
                py::array_t<std::complex<double>> psi_in,
                const int psi_in_dmax,
                py::array_t<double> eigenvalue_in,
                const std::vector<bool>& is_occupied,
                const bool& scf_type)
            {
                py::buffer_info h_buf = h_mat.request();
                if (h_buf.ndim != 1)
                {
                    throw std::runtime_error("h_mat must be 1D array representing a matrix");
                }
                py::buffer_info psi_buf = psi_in.request();
                if (psi_buf.ndim != 1)
                {
                    throw std::runtime_error("psi_in must be 1D array representing a matrix");
                }
                py::buffer_info eigenvalue_buf = eigenvalue_in.request();
                if (eigenvalue_buf.ndim != 1)
                {
                    throw std::runtime_error("eigenvalue_in must be 1D array storing eigenvalues");
                }

                std::complex<double>* h_mat_ptr = static_cast<std::complex<double>*>(h_buf.ptr);
                std::complex<double>* psi_in_ptr = static_cast<std::complex<double>*>(psi_buf.ptr);
                double* eigenvalue_in_ptr = static_cast<double*>(eigenvalue_buf.ptr);

                // TODO: Wrap std::function<void(T*, T*, const int, const int, const int, const int)>
                //       to a python callable type
                auto hpsi_func = [h_mat_ptr] (std::complex<double> *hpsi_out,
                            std::complex<double> *psi_in, const int nband_in,
                            const int nbasis_in, const int band_index1,
                            const int band_index2) 
                {
                    const std::complex<double> *one_ = nullptr, *zero_ = nullptr;

                    one_ = new std::complex<double>(1.0, 0.0);
                    zero_ = new std::complex<double>(0.0, 0.0);

                    base_device::DEVICE_CPU *ctx = {};

                    hsolver::gemm_op<std::complex<double>, base_device::DEVICE_CPU>()(
                        ctx, 'N', 'N', 
                        nbasis_in, band_index2 - band_index1 + 1, nbasis_in, 
                        one_, h_mat_ptr, nbasis_in, 
                        psi_in + band_index1 * nbasis_in, nbasis_in,
                        zero_, hpsi_out + band_index1 * nbasis_in, nbasis_in);
                };

                return self.diag(hpsi_func, 
                                 psi_in_ptr, 
                                 psi_in_dmax, 
                                 eigenvalue_in_ptr, 
                                 is_occupied, 
                                 scf_type);
            }
        );
}
