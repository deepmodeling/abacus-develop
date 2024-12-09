#ifndef PYTHON_PYABACUS_SRC_PY_DIAGO_DAV_SUBSPACE_HPP
#define PYTHON_PYABACUS_SRC_PY_DIAGO_DAV_SUBSPACE_HPP

#include <complex>
#include <functional>

#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "module_hsolver/diago_dav_subspace.h"

namespace py = pybind11;

namespace py_hsolver
{

class PyDiagoDavSubspace
{
public:
    PyDiagoDavSubspace(int nbasis, int nband) : nbasis(nbasis), nband(nband)
    {
        psi = new std::complex<double>[nbasis * nband];
        eigenvalue = new double[nband];
    }

    PyDiagoDavSubspace(const PyDiagoDavSubspace&) = delete;
    PyDiagoDavSubspace& operator=(const PyDiagoDavSubspace&) = delete;
    PyDiagoDavSubspace(PyDiagoDavSubspace&& other) : nbasis(other.nbasis), nband(other.nband)
    {
        psi = other.psi;
        eigenvalue = other.eigenvalue;

        other.psi = nullptr;
        other.eigenvalue = nullptr;
    }

    ~PyDiagoDavSubspace()
    {
        if (psi != nullptr) 
        { 
            delete[] psi; 
            psi = nullptr; 
        }
        if (eigenvalue != nullptr) 
        { 
            delete[] eigenvalue; 
            eigenvalue = nullptr; 
        }
    }

    void set_psi(py::array_t<std::complex<double>> psi_in)
    {
        assert(psi_in.size() == nbasis * nband);

        for (size_t i = 0; i < nbasis * nband; ++i)
        {
            psi[i] = psi_in.at(i);
        }
    }

    py::array_t<std::complex<double>> get_psi()
    {
        py::array_t<std::complex<double>> psi_out(nband * nbasis);
        py::buffer_info psi_out_buf = psi_out.request();

        std::complex<double>* psi_out_ptr = static_cast<std::complex<double>*>(psi_out_buf.ptr);

        for (size_t i = 0; i < nband * nbasis; ++i)
        {
            psi_out_ptr[i] = psi[i];
        }

        return psi_out;   
    }

    void init_eigenvalue()
    {
        for (size_t i = 0; i < nband; ++i)
        {
            eigenvalue[i] = 0.0;
        }
    }

    py::array_t<double> get_eigenvalue()
    {
        py::array_t<double> eigenvalue_out(nband);
        py::buffer_info eigenvalue_out_buf = eigenvalue_out.request();

        double* eigenvalue_out_ptr = static_cast<double*>(eigenvalue_out_buf.ptr);

        for (size_t i = 0; i < nband; ++i)
        {
            eigenvalue_out_ptr[i] = eigenvalue[i];
        }

        return eigenvalue_out;
    }

    int diag(
        std::function<py::array_t<std::complex<double>>(py::array_t<std::complex<double>>)> mm_op,
        std::vector<double>& precond_vec,
        int dav_ndim,
        double tol,
        int max_iter,
        bool need_subspace,
        std::vector<double>& diag_ethr,
        bool scf_type,
        hsolver::diag_comm_info comm_info
    ) {
        auto hpsi_func = [mm_op] (
            std::complex<double> *psi_in,
            std::complex<double> *hpsi_out, 
            const int ld_psi,
            const int nvec
        ) {
            // Note: numpy's py::array_t is row-major, but
            //       our raw pointer-array is column-major
            py::array_t<std::complex<double>, py::array::f_style> psi({ld_psi, nvec});
            py::buffer_info psi_buf = psi.request();
            std::complex<double>* psi_ptr = static_cast<std::complex<double>*>(psi_buf.ptr);
            std::copy(psi_in, psi_in + nvec * ld_psi, psi_ptr);

            py::array_t<std::complex<double>, py::array::f_style> hpsi = mm_op(psi);

            py::buffer_info hpsi_buf = hpsi.request();
            std::complex<double>* hpsi_ptr = static_cast<std::complex<double>*>(hpsi_buf.ptr);
            std::copy(hpsi_ptr, hpsi_ptr + nvec * ld_psi, hpsi_out);
        };

        obj = std::make_unique<hsolver::Diago_DavSubspace<std::complex<double>, base_device::DEVICE_CPU>>(
            precond_vec, 
            nband, 
            nbasis, 
            dav_ndim, 
            tol, 
            max_iter, 
            need_subspace, 
            comm_info
        );

        return obj->diag(hpsi_func, psi, nbasis, eigenvalue, diag_ethr, scf_type);
    }

private:
    std::complex<double>* psi = nullptr;
    double* eigenvalue = nullptr;

    int nbasis;
    int nband;

    std::unique_ptr<hsolver::Diago_DavSubspace<std::complex<double>, base_device::DEVICE_CPU>> obj;
};

} // namespace py_hsolver

#endif