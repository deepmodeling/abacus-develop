#ifndef PYTHON_PYABACUS_SRC_PY_DIAGO_CG_HPP
#define PYTHON_PYABACUS_SRC_PY_DIAGO_CG_HPP

#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <ATen/core/tensor.h>
#include <ATen/core/tensor_types.h>

#include "module_hsolver/diago_cg.h"
#include "module_base/module_device/memory_op.h"

namespace py = pybind11;

namespace py_hsolver
{

class PyDiagoCG
{
    PyDiagoCG() { }
    PyDiagoCG(const PyDiagoCG&) = delete;
    PyDiagoCG& operator=(const PyDiagoCG&) = delete;
    PyDiagoCG(PyDiagoCG&& other)
    {
        psi = other.psi;
        other.psi = nullptr;

        eig = other.eig;
        other.eig = nullptr;
    }

    ~PyDiagoCG() 
    {
        if (psi != nullptr) 
        {
            delete psi;
            psi = nullptr;
        }

        if (eig != nullptr)
        {
            delete eig;
            eig = nullptr;
        }
    }

    void diag(
        std::function<py::array_t<std::complex<double>>(py::array_t<std::complex<double>>)> mm_op,
        bool scf_type
    ) {
        const std::string basis_type = "pw";
        const std::string calculation = scf_type ? "scf" : "nscf";

        auto hpsi_func = [mm_op] (const ct::Tensor& psi_in, ct::Tensor& hpsi_out) {
            const auto ndim = psi_in.shape().ndim();
            REQUIRES_OK(ndim <= 2, "dims of psi_in should be less than or equal to 2");
            const int nvec   = ndim == 1 ? 1 : psi_in.shape().dim_size(0);
            const int ld_psi = ndim == 1 ? psi_in.NumElements() : psi_in.shape().dim_size(1)

            // py::array_t<std::complex<double>> psi({ld_psi, nvec});
            // py::buffer_info psi_buf = psi.request();
            // std::complex<double>* psi_ptr = static_cast<std::complex<double>*>(psi_buf.ptr);
            // std::copy(psi_in, psi_in + nvec * ld_psi, psi_ptr);

            // py::array_t<std::complex<double>> hpsi = mm_op(psi);

            // py::buffer_info hpsi_buf = hpsi.request();
            // std::complex<double>* hpsi_ptr = static_cast<std::complex<double>*>(hpsi_buf.ptr);
            // std::copy(hpsi_ptr, hpsi_ptr + nvec * ld_psi, hpsi_out);
        };

        auto subspace_func = [] (const ct::Tensor& psi_in, ct::Tensor& psi_out) { /*do nothing*/ };

        auto spsi_func = [this] (const ct::Tensor& psi_in, ct::Tensor& spsi_out) {
            const auto ndim = psi_in.shape().ndim();
            REQUIRES_OK(ndim <= 2, "dims of psi_in should be less than or equal to 2");
            const int nrow   = ndim == 1 ? psi_in.NumElements() : psi_in.shape().dim_size(1);
            const int nbands = ndim == 1 ? 1 : psi_in.shape().dim_size(0);
            syncmem_z2z_h2h_op()(
                this->ctx,
                this->ctx,
                spsi_out.data<std::complex<double>>(), 
                psi_in.data<std::complex<double>>(), 
                static_cast<size_t>(nrow * nbands)
            );
        };
    }

private:
    base_device::DEVICE_CPU* ctx = {};

    ct::Tensor* psi = nullptr;
    ct::Tensor* eig = nullptr;

};

} // namespace py_hsolver

#endif