#ifndef PYTHON_PYABACUS_SRC_PY_DIAGO_CG_HPP
#define PYTHON_PYABACUS_SRC_PY_DIAGO_CG_HPP

#include <complex>
#include <functional>

#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <ATen/core/tensor.h>
#include <ATen/core/tensor_map.h>
#include <ATen/core/tensor_types.h>

#include "module_hsolver/diago_cg.h"
#include "module_base/module_device/memory_op.h"

namespace py = pybind11;

namespace py_hsolver
{

class PyDiagoCG
{
public:
    PyDiagoCG(int dim, int num_eigs) : dim{dim}, num_eigs{num_eigs} { }
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

    void init_eig() 
    {
        eig = new ct::Tensor(ct::DataType::DT_DOUBLE, {num_eigs});
        eig->zero();
    }

    py::array_t<double> get_eig() 
    {
        py::array_t<double> eig_out(eig->NumElements());
        py::buffer_info eig_buf = eig_out.request();
        double* eig_out_ptr = static_cast<double*>(eig_buf.ptr);

        if (eig == nullptr) {
            throw std::runtime_error("eig is not initialized");
        }
        double* eig_ptr = eig->data<double>();

        std::copy(eig_ptr, eig_ptr + eig->NumElements(), eig_out_ptr);
        return eig_out;
    }

    void set_psi(py::array_t<std::complex<double>> psi_in)
    {
        py::buffer_info psi_buf = psi_in.request();
        std::complex<double>* psi_ptr = static_cast<std::complex<double>*>(psi_buf.ptr);

        psi = new ct::TensorMap(
            psi_ptr,
            ct::DataType::DT_COMPLEX_DOUBLE,
            ct::DeviceType::CpuDevice,
            ct::TensorShape({num_eigs, dim})
        );
    }

    py::array_t<std::complex<double>> get_psi()
    {
        py::array_t<std::complex<double>> psi_out({num_eigs, dim});
        py::buffer_info psi_buf = psi_out.request();
        std::complex<double>* psi_out_ptr = static_cast<std::complex<double>*>(psi_buf.ptr);

        if (psi == nullptr) {
            throw std::runtime_error("psi is not initialized");
        }
        std::complex<double>* psi_ptr = psi->data<std::complex<double>>();

        std::copy(psi_ptr, psi_ptr + psi->NumElements(), psi_out_ptr);
        return psi_out;
    }

    void set_prec(py::array_t<double> prec_in)
    {
        py::buffer_info prec_buf = prec_in.request();
        double* prec_ptr = static_cast<double*>(prec_buf.ptr);

        prec = new ct::TensorMap(
            prec_ptr,
            ct::DataType::DT_DOUBLE,
            ct::DeviceType::CpuDevice,
            ct::TensorShape({dim})
        );
    }

    void diag(std::function<py::array_t<std::complex<double>>(py::array_t<std::complex<double>>)> mm_op,
              int diag_ndim,
              double tol,
              const std::vector<double>& diag_ethr,
              bool need_subspace,
              bool scf_type,
              int nproc_in_pool = 1
    ) {
        const std::string basis_type = "pw";
        const std::string calculation = scf_type ? "scf" : "nscf";

        auto hpsi_func = [mm_op] (const ct::Tensor& psi_in, ct::Tensor& hpsi_out) {
            const auto ndim = psi_in.shape().ndim();
            REQUIRES_OK(ndim <= 2, "dims of psi_in should be less than or equal to 2");
            const int nvec   = ndim == 1 ? 1 : psi_in.shape().dim_size(0);
            const int ld_psi = ndim == 1 ? psi_in.NumElements() : psi_in.shape().dim_size(1);

            // Note: numpy's py::array_t is row-major, and
            //       our tensor-array is row-major
            py::array_t<std::complex<double>> psi({ld_psi, nvec});
            py::buffer_info psi_buf = psi.request();
            std::complex<double>* psi_ptr = static_cast<std::complex<double>*>(psi_buf.ptr);
            std::copy(psi_in.data<std::complex<double>>(), psi_in.data<std::complex<double>>() + nvec * ld_psi, psi_ptr);

            py::array_t<std::complex<double>> hpsi = mm_op(psi);

            py::buffer_info hpsi_buf = hpsi.request();
            std::complex<double>* hpsi_ptr = static_cast<std::complex<double>*>(hpsi_buf.ptr);
            std::copy(hpsi_ptr, hpsi_ptr + nvec * ld_psi, hpsi_out.data<std::complex<double>>());
        };

        auto subspace_func = [] (const ct::Tensor& psi_in, ct::Tensor& psi_out) { /*do nothing*/ };

        auto spsi_func = [this] (const ct::Tensor& psi_in, ct::Tensor& spsi_out) {
            const auto ndim = psi_in.shape().ndim();
            REQUIRES_OK(ndim <= 2, "dims of psi_in should be less than or equal to 2");
            const int nrow   = ndim == 1 ? psi_in.NumElements() : psi_in.shape().dim_size(1);
            const int nbands = ndim == 1 ? 1 : psi_in.shape().dim_size(0);
            syncmem_z2z_h2h_op()(
                spsi_out.data<std::complex<double>>(), 
                psi_in.data<std::complex<double>>(), 
                static_cast<size_t>(nrow * nbands)
            );
        };

        cg = std::make_unique<hsolver::DiagoCG<std::complex<double>, base_device::DEVICE_CPU>>(
            basis_type,
            calculation,
            need_subspace,
            subspace_func,
            tol,
            diag_ndim,
            nproc_in_pool
        );

        cg->diag(hpsi_func, spsi_func, *psi, *eig, diag_ethr, *prec);
    }

private:
    base_device::DEVICE_CPU* ctx = {};

    int dim;
    int num_eigs;

    ct::Tensor* psi = nullptr;
    ct::Tensor* eig = nullptr;
    ct::Tensor* prec = nullptr;

    std::unique_ptr<hsolver::DiagoCG<std::complex<double>, base_device::DEVICE_CPU>> cg;
};

} // namespace py_hsolver

#endif