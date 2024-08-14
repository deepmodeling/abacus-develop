#ifndef PYTHON_PYABACUS_SRC_PY_DIAGO_DAV_SUBSPACE_HPP
#define PYTHON_PYABACUS_SRC_PY_DIAGO_DAV_SUBSPACE_HPP

#include <complex>

#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
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
        py::array_t<std::complex<double>> h_mat,
        std::vector<double> precond_vec,
        int dav_ndim,
        double tol,
        int max_iter,
        bool need_subspace,
        std::vector<bool> is_occupied,
        bool scf_type,
        hsolver::diag_comm_info comm_info
    ) {
        py::buffer_info h_mat_buf = h_mat.request();

        if (h_mat_buf.ndim != 1)
        {
            throw std::runtime_error("h_mat must be a 1D array representing a 2D matrix");
        }

        std::complex<double>* h_mat_ptr = static_cast<std::complex<double>*>(h_mat_buf.ptr);

        // TODO: Wrap std::function<void(T*, T*, const int, const int, const int, const int)>
        //       to a python callable type
        auto hpsi_func = [h_mat_ptr] (std::complex<double> *hpsi_out,
                    std::complex<double> *psi_in, const int nband_in,
                    const int nbasis_in, const int band_index1,
                    const int band_index2) 
        {
            const std::complex<double> one(1.0, 0.0);
            const std::complex<double> zero(0.0, 0.0);

            base_device::DEVICE_CPU *ctx = {};

            hsolver::gemm_op<std::complex<double>, base_device::DEVICE_CPU>()(
                ctx, 'N', 'N', 
                nbasis_in, band_index2 - band_index1 + 1, nbasis_in, 
                &one, h_mat_ptr, nbasis_in, 
                psi_in + band_index1 * nbasis_in, nbasis_in,
                &zero, hpsi_out + band_index1 * nbasis_in, nbasis_in);
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

        return obj->diag(hpsi_func, psi, nbasis, eigenvalue, is_occupied, scf_type);
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