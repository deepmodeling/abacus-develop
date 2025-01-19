#include "propagator.h"

#include "module_base/lapack_connector.h"
#include "module_base/module_container/ATen/kernels/blas.h"
#include "module_base/module_container/ATen/kernels/lapack.h"
#include "module_base/module_container/ATen/kernels/memory.h" // memory operations (Tensor)
#include "module_base/module_device/memory_op.h"              // memory operations
#include "module_base/scalapack_connector.h"
#include "module_parameter/parameter.h"

#include <complex>
#include <iostream>

namespace module_tddft
{
Propagator::~Propagator()
{
}
#ifdef __MPI

inline int globalIndex(int localindex, int nblk, int nprocs, int myproc)
{
    int iblock, gIndex;
    iblock = localindex / nblk;
    gIndex = (iblock * nprocs + myproc) * nblk + localindex % nblk;
    return gIndex;
}

void Propagator::compute_propagator(const int nlocal,
                                    const std::complex<double>* Stmp,
                                    const std::complex<double>* Htmp,
                                    const std::complex<double>* H_laststep,
                                    std::complex<double>* U_operator,
                                    const int print_matrix) const
{
    int tag;
    switch (ptype)
    {
    case 0:
        compute_propagator_cn2(nlocal, Stmp, Htmp, U_operator, print_matrix);
        break;

    case 1:
        tag = 1;
        compute_propagator_taylor(nlocal, Stmp, Htmp, U_operator, print_matrix, tag);
        break;

    case 2:
        compute_propagator_etrs(nlocal, Stmp, Htmp, H_laststep, U_operator, print_matrix);
        break;

    default:
        ModuleBase::WARNING_QUIT("Propagator::compute_propagator", "Method of RT-TDDFT propagator is wrong!");
        break;
    }
}

template <typename Device>
void Propagator::compute_propagator_tensor(const int nlocal,
                                           const ct::Tensor& Stmp,
                                           const ct::Tensor& Htmp,
                                           const ct::Tensor& H_laststep,
                                           ct::Tensor& U_operator,
                                           const int print_matrix,
                                           const bool use_lapack) const
{
    int tag;
    switch (ptype)
    {
    case 0:
        if (!use_lapack)
        {
            compute_propagator_cn2_tensor(nlocal, Stmp, Htmp, U_operator, print_matrix);
        }
        else
        {
            int myid, root_proc = 0;
            MPI_Comm_rank(MPI_COMM_WORLD, &myid);
            if (myid == root_proc)
            {
                compute_propagator_cn2_tensor_lapack<Device>(nlocal, Stmp, Htmp, U_operator, print_matrix);
            }
        }
        break;

    default:
        ModuleBase::WARNING_QUIT("Propagator::compute_propagator_tensor",
                                 "The Tensor-based RT-TDDFT propagator currently supports Crank–Nicolson method only!");
        break;
    }
}

void Propagator::compute_propagator_cn2(const int nlocal,
                                        const std::complex<double>* Stmp,
                                        const std::complex<double>* Htmp,
                                        std::complex<double>* U_operator,
                                        const int print_matrix) const
{
    // (1) copy Htmp to Numerator & Denominator
    std::complex<double>* Numerator = new std::complex<double>[this->ParaV->nloc];
    ModuleBase::GlobalFunc::ZEROS(Numerator, this->ParaV->nloc);
    BlasConnector::copy(this->ParaV->nloc, Htmp, 1, Numerator, 1);

    std::complex<double>* Denominator = new std::complex<double>[this->ParaV->nloc];
    ModuleBase::GlobalFunc::ZEROS(Denominator, this->ParaV->nloc);
    BlasConnector::copy(this->ParaV->nloc, Htmp, 1, Denominator, 1);

    if (print_matrix)
    {
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << " S matrix :" << std::endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                GlobalV::ofs_running << Stmp[i * this->ParaV->ncol + j].real() << "+"
                                     << Stmp[i * this->ParaV->ncol + j].imag() << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << " H matrix :" << std::endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                GlobalV::ofs_running << Numerator[i * this->ParaV->ncol + j].real() << "+"
                                     << Numerator[i * this->ParaV->ncol + j].imag() << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
    }

    // ->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // (2) compute Numerator & Denominator by GEADD
    // Numerator = Stmp - i*para * Htmp;     beta1 = - para = -0.25 * this->dt
    // Denominator = Stmp + i*para * Htmp;   beta2 = para = 0.25 * this->dt
    std::complex<double> alpha = {1.0, 0.0};
    std::complex<double> beta1 = {0.0, -0.25 * this->dt};
    std::complex<double> beta2 = {0.0, 0.25 * this->dt};

    ScalapackConnector::geadd('N',
                              nlocal,
                              nlocal,
                              alpha,
                              Stmp,
                              1,
                              1,
                              this->ParaV->desc,
                              beta1,
                              Numerator,
                              1,
                              1,
                              this->ParaV->desc);
    ScalapackConnector::geadd('N',
                              nlocal,
                              nlocal,
                              alpha,
                              Stmp,
                              1,
                              1,
                              this->ParaV->desc,
                              beta2,
                              Denominator,
                              1,
                              1,
                              this->ParaV->desc);

    if (print_matrix)
    {
        GlobalV::ofs_running << " beta=" << beta1 << std::endl;
        GlobalV::ofs_running << " fenmu:" << std::endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                GlobalV::ofs_running << Denominator[i * this->ParaV->ncol + j].real() << "+"
                                     << Denominator[i * this->ParaV->ncol + j].imag() << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
    }

    //->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // (3) Next, invert Denominator
    // What is the size of ipiv exactly? Need to check ScaLAPACK documentation!
    // But anyway, not this->ParaV->nloc
    int* ipiv = new int[this->ParaV->nrow + this->ParaV->nb];
    ModuleBase::GlobalFunc::ZEROS(ipiv, this->ParaV->nrow + this->ParaV->nb);
    int info = 0;
    // (3.1) compute ipiv
    ScalapackConnector::getrf(nlocal, nlocal, Denominator, 1, 1, this->ParaV->desc, ipiv, &info);

    // Print ipiv
    if (print_matrix)
    {
        GlobalV::ofs_running << " this->ParaV->nloc = " << this->ParaV->nloc << std::endl;
        GlobalV::ofs_running << " this->ParaV->nrow = " << this->ParaV->nrow << std::endl;
        GlobalV::ofs_running << " this->ParaV->ncol = " << this->ParaV->ncol << std::endl;
        GlobalV::ofs_running << " this->ParaV->nb = " << this->ParaV->nb << std::endl;
        GlobalV::ofs_running << " this->ParaV->get_block_size() = " << this->ParaV->get_block_size() << std::endl;
        GlobalV::ofs_running << " nlocal = " << nlocal << std::endl;
        GlobalV::ofs_running << " ipiv:" << std::endl;
        for (int i = 0; i < this->ParaV->nloc; i++)
        {
            GlobalV::ofs_running << ipiv[i] << " ";
        }
        GlobalV::ofs_running << std::endl;
    }

    int lwork = -1;
    int liwotk = -1;
    std::vector<std::complex<double>> work(1, 0);
    std::vector<int> iwork(1, 0);
    // (3.2) compute work
    ScalapackConnector::getri(nlocal,
                              Denominator,
                              1,
                              1,
                              this->ParaV->desc,
                              ipiv,
                              work.data(),
                              &lwork,
                              iwork.data(),
                              &liwotk,
                              &info);
    lwork = work[0].real();
    work.resize(lwork, 0);
    liwotk = iwork[0];
    iwork.resize(liwotk, 0);
    // (3.3) compute inverse matrix of Denominator
    ScalapackConnector::getri(nlocal,
                              Denominator,
                              1,
                              1,
                              this->ParaV->desc,
                              ipiv,
                              work.data(),
                              &lwork,
                              iwork.data(),
                              &liwotk,
                              &info);
    assert(0 == info);

    //->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    // (4) U_operator = Denominator * Numerator;
    ScalapackConnector::gemm('N',
                             'N',
                             nlocal,
                             nlocal,
                             nlocal,
                             1.0,
                             Denominator,
                             1,
                             1,
                             this->ParaV->desc,
                             Numerator,
                             1,
                             1,
                             this->ParaV->desc,
                             0.0,
                             U_operator,
                             1,
                             1,
                             this->ParaV->desc);

    if (print_matrix)
    {
        GlobalV::ofs_running << " fenmu^-1:" << std::endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                GlobalV::ofs_running << Denominator[i * this->ParaV->ncol + j].real() << "+"
                                     << Denominator[i * this->ParaV->ncol + j].imag() << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << " fenzi:" << std::endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                GlobalV::ofs_running << Numerator[i * this->ParaV->ncol + j].real() << "+"
                                     << Numerator[i * this->ParaV->ncol + j].imag() << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << " U operator:" << std::endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                double aa, bb;
                aa = U_operator[i * this->ParaV->ncol + j].real();
                bb = U_operator[i * this->ParaV->ncol + j].imag();
                if (std::abs(aa) < 1e-8)
                {
                    aa = 0.0;
                }
                if (std::abs(bb) < 1e-8)
                {
                    bb = 0.0;
                }
                GlobalV::ofs_running << aa << "+" << bb << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
    }

    delete[] Numerator;
    delete[] Denominator;
    delete[] ipiv;
}

void Propagator::compute_propagator_cn2_tensor(const int nlocal,
                                               const ct::Tensor& Stmp,
                                               const ct::Tensor& Htmp,
                                               ct::Tensor& U_operator,
                                               const int print_matrix) const
{
    // (1) copy Htmp to Numerator & Denominator
    ct::Tensor Numerator(ct::DataType::DT_COMPLEX_DOUBLE,
                         ct::DeviceType::CpuDevice,
                         ct::TensorShape({this->ParaV->nloc}));
    Numerator.zero();
    BlasConnector::copy(this->ParaV->nloc,
                        Htmp.data<std::complex<double>>(),
                        1,
                        Numerator.data<std::complex<double>>(),
                        1);

    ct::Tensor Denominator(ct::DataType::DT_COMPLEX_DOUBLE,
                           ct::DeviceType::CpuDevice,
                           ct::TensorShape({this->ParaV->nloc}));
    Denominator.zero();
    BlasConnector::copy(this->ParaV->nloc,
                        Htmp.data<std::complex<double>>(),
                        1,
                        Denominator.data<std::complex<double>>(),
                        1);

    if (print_matrix)
    {
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << " S matrix :" << std::endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                GlobalV::ofs_running << Stmp.data<std::complex<double>>()[i * this->ParaV->ncol + j].real() << "+"
                                     << Stmp.data<std::complex<double>>()[i * this->ParaV->ncol + j].imag() << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << " H matrix :" << std::endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                GlobalV::ofs_running << Numerator.data<std::complex<double>>()[i * this->ParaV->ncol + j].real() << "+"
                                     << Numerator.data<std::complex<double>>()[i * this->ParaV->ncol + j].imag()
                                     << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
    }

    // ->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // (2) compute Numerator & Denominator by GEADD
    // Numerator = Stmp - i*para * Htmp;     beta1 = - para = -0.25 * this->dt
    // Denominator = Stmp + i*para * Htmp;   beta2 = para = 0.25 * this->dt
    std::complex<double> alpha = {1.0, 0.0};
    std::complex<double> beta1 = {0.0, -0.25 * this->dt};
    std::complex<double> beta2 = {0.0, 0.25 * this->dt};

    ScalapackConnector::geadd('N',
                              nlocal,
                              nlocal,
                              alpha,
                              Stmp.data<std::complex<double>>(),
                              1,
                              1,
                              this->ParaV->desc,
                              beta1,
                              Numerator.data<std::complex<double>>(),
                              1,
                              1,
                              this->ParaV->desc);
    ScalapackConnector::geadd('N',
                              nlocal,
                              nlocal,
                              alpha,
                              Stmp.data<std::complex<double>>(),
                              1,
                              1,
                              this->ParaV->desc,
                              beta2,
                              Denominator.data<std::complex<double>>(),
                              1,
                              1,
                              this->ParaV->desc);

    if (print_matrix)
    {
        GlobalV::ofs_running << " beta=" << beta1 << std::endl;
        GlobalV::ofs_running << " fenmu:" << std::endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                GlobalV::ofs_running << Denominator.data<std::complex<double>>()[i * this->ParaV->ncol + j].real()
                                     << "+"
                                     << Denominator.data<std::complex<double>>()[i * this->ParaV->ncol + j].imag()
                                     << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
    }

    //->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // (3) Next, invert Denominator
    ct::Tensor ipiv(ct::DataType::DT_INT,
                    ct::DeviceType::CpuDevice,
                    ct::TensorShape({this->ParaV->nrow + this->ParaV->nb}));
    ipiv.zero();
    int info = 0;
    // (3.1) compute ipiv
    ScalapackConnector::getrf(nlocal,
                              nlocal,
                              Denominator.data<std::complex<double>>(),
                              1,
                              1,
                              this->ParaV->desc,
                              ipiv.data<int>(),
                              &info);

    // Print ipiv
    if (print_matrix)
    {
        GlobalV::ofs_running << " this->ParaV->nloc = " << this->ParaV->nloc << std::endl;
        GlobalV::ofs_running << " this->ParaV->nrow = " << this->ParaV->nrow << std::endl;
        GlobalV::ofs_running << " this->ParaV->ncol = " << this->ParaV->ncol << std::endl;
        GlobalV::ofs_running << " this->ParaV->nb = " << this->ParaV->nb << std::endl;
        GlobalV::ofs_running << " this->ParaV->get_block_size() = " << this->ParaV->get_block_size() << std::endl;
        GlobalV::ofs_running << " nlocal = " << nlocal << std::endl;
        GlobalV::ofs_running << " ipiv:" << std::endl;
        for (int i = 0; i < this->ParaV->nloc; i++)
        {
            GlobalV::ofs_running << ipiv.data<int>()[i] << " ";
        }
        GlobalV::ofs_running << std::endl;
    }

    int lwork = -1;
    int liwotk = -1;
    ct::Tensor work(ct::DataType::DT_COMPLEX_DOUBLE, ct::DeviceType::CpuDevice, ct::TensorShape({1}));
    ct::Tensor iwork(ct::DataType::DT_INT, ct::DeviceType::CpuDevice, ct::TensorShape({1}));
    // (3.2) compute work
    ScalapackConnector::getri(nlocal,
                              Denominator.data<std::complex<double>>(),
                              1,
                              1,
                              this->ParaV->desc,
                              ipiv.data<int>(),
                              work.data<std::complex<double>>(),
                              &lwork,
                              iwork.data<int>(),
                              &liwotk,
                              &info);
    lwork = work.data<std::complex<double>>()[0].real();
    work.resize(ct::TensorShape({lwork}));
    liwotk = iwork.data<int>()[0];
    iwork.resize(ct::TensorShape({liwotk}));
    // (3.3) compute inverse matrix of Denominator
    ScalapackConnector::getri(nlocal,
                              Denominator.data<std::complex<double>>(),
                              1,
                              1,
                              this->ParaV->desc,
                              ipiv.data<int>(),
                              work.data<std::complex<double>>(),
                              &lwork,
                              iwork.data<int>(),
                              &liwotk,
                              &info);
    assert(0 == info);

    //->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    // (4) U_operator = Denominator * Numerator;
    ScalapackConnector::gemm('N',
                             'N',
                             nlocal,
                             nlocal,
                             nlocal,
                             1.0,
                             Denominator.data<std::complex<double>>(),
                             1,
                             1,
                             this->ParaV->desc,
                             Numerator.data<std::complex<double>>(),
                             1,
                             1,
                             this->ParaV->desc,
                             0.0,
                             U_operator.data<std::complex<double>>(),
                             1,
                             1,
                             this->ParaV->desc);

    if (print_matrix)
    {
        GlobalV::ofs_running << " fenmu^-1:" << std::endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                GlobalV::ofs_running << Denominator.data<std::complex<double>>()[i * this->ParaV->ncol + j].real()
                                     << "+"
                                     << Denominator.data<std::complex<double>>()[i * this->ParaV->ncol + j].imag()
                                     << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << " fenzi:" << std::endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                GlobalV::ofs_running << Numerator.data<std::complex<double>>()[i * this->ParaV->ncol + j].real() << "+"
                                     << Numerator.data<std::complex<double>>()[i * this->ParaV->ncol + j].imag()
                                     << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << " U operator:" << std::endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                double aa, bb;
                aa = U_operator.data<std::complex<double>>()[i * this->ParaV->ncol + j].real();
                bb = U_operator.data<std::complex<double>>()[i * this->ParaV->ncol + j].imag();
                if (std::abs(aa) < 1e-8)
                {
                    aa = 0.0;
                }
                if (std::abs(bb) < 1e-8)
                {
                    bb = 0.0;
                }
                GlobalV::ofs_running << aa << "+" << bb << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
    }
}

//------------------------ Utility function ------------------------//
// Auxiliary function: process non-complex types, return value 1.0
template <typename T>
inline T init_value(typename std::enable_if<!std::is_same<T, std::complex<float>>::value
                                            && !std::is_same<T, std::complex<double>>::value>::type* = nullptr)
{
    return T(1.0);
}

// Auxiliary function: process complex types, return value 1.0 + 0.0i
template <typename T>
inline T init_value(typename std::enable_if<std::is_same<T, std::complex<float>>::value
                                            || std::is_same<T, std::complex<double>>::value>::type* = nullptr)
{
    return T(1.0, 0.0);
}

// Create an identity matrix of size n×n
template <typename T>
ct::Tensor create_identity_matrix(const int n, ct::DeviceType device = ct::DeviceType::CpuDevice)
{
    // Choose the data type of the Tensor
    ct::DataType data_type;
    if (std::is_same<T, float>::value)
    {
        data_type = ct::DataType::DT_FLOAT;
    }
    else if (std::is_same<T, double>::value)
    {
        data_type = ct::DataType::DT_DOUBLE;
    }
    else if (std::is_same<T, std::complex<float>>::value)
    {
        data_type = ct::DataType::DT_COMPLEX;
    }
    else if (std::is_same<T, std::complex<double>>::value)
    {
        data_type = ct::DataType::DT_COMPLEX_DOUBLE;
    }
    else
    {
        static_assert(std::is_same<T, float>::value || std::is_same<T, double>::value
                          || std::is_same<T, std::complex<float>>::value
                          || std::is_same<T, std::complex<double>>::value,
                      "Unsupported data type!");
    }

    ct::Tensor tensor(data_type, device, ct::TensorShape({n, n}));
    tensor.zero();

    // Set the diagonal elements to 1
    if (device == ct::DeviceType::CpuDevice)
    {
        // For CPU, we can directly access the data
        T* data_ptr = tensor.data<T>();
        for (int i = 0; i < n; ++i)
        {
            data_ptr[i * n + i] = init_value<T>();
        }
    }
    else if (device == ct::DeviceType::GpuDevice)
    {
        // For GPU, we need to use a kernel to set the diagonal elements
        T* data_ptr = tensor.data<T>();
        for (int i = 0; i < n; ++i)
        {
            T value = init_value<T>();
            ct::kernels::set_memory<T, ct::DEVICE_GPU>()(data_ptr + i * n + i, value, 1);
        }
    }

    return tensor;
}
//------------------------ Utility function ------------------------//

template <typename Device>
void Propagator::compute_propagator_cn2_tensor_lapack(const int nlocal,
                                                      const ct::Tensor& Stmp,
                                                      const ct::Tensor& Htmp,
                                                      ct::Tensor& U_operator,
                                                      const int print_matrix) const
{
    // ct_device_type = ct::DeviceType::CpuDevice or ct::DeviceType::GpuDevice
    ct::DeviceType ct_device_type = ct::DeviceTypeToEnum<Device>::value;
    // ct_Device = ct::DEVICE_CPU or ct::DEVICE_GPU
    using ct_Device = typename ct::PsiToContainer<Device>::type;

    // (1) copy Htmp to Numerator & Denominator
    ct::Tensor Numerator(ct::DataType::DT_COMPLEX_DOUBLE, ct_device_type, ct::TensorShape({nlocal * nlocal}));
    Numerator.zero();
    base_device::memory::synchronize_memory_op<std::complex<double>, Device, Device>()(
        Numerator.data<std::complex<double>>(),
        Htmp.data<std::complex<double>>(),
        nlocal * nlocal);

    ct::Tensor Denominator(ct::DataType::DT_COMPLEX_DOUBLE, ct_device_type, ct::TensorShape({nlocal * nlocal}));
    Denominator.zero();
    base_device::memory::synchronize_memory_op<std::complex<double>, Device, Device>()(
        Denominator.data<std::complex<double>>(),
        Htmp.data<std::complex<double>>(),
        nlocal * nlocal);

    if (print_matrix)
    {
        ct::Tensor Stmp_cpu = Stmp.to_device<ct::DEVICE_CPU>();
        ct::Tensor Numerator_cpu = Numerator.to_device<ct::DEVICE_CPU>();

        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << " S matrix :" << std::endl;
        for (int i = 0; i < nlocal; i++)
        {
            for (int j = 0; j < nlocal; j++)
            {
                GlobalV::ofs_running << Stmp_cpu.data<std::complex<double>>()[i * nlocal + j].real() << "+"
                                     << Stmp_cpu.data<std::complex<double>>()[i * nlocal + j].imag() << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << " H matrix :" << std::endl;
        for (int i = 0; i < nlocal; i++)
        {
            for (int j = 0; j < nlocal; j++)
            {
                GlobalV::ofs_running << Numerator_cpu.data<std::complex<double>>()[i * nlocal + j].real() << "+"
                                     << Numerator_cpu.data<std::complex<double>>()[i * nlocal + j].imag() << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
    }

    // ->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // (2) compute Numerator & Denominator by GEADD
    // Numerator = Stmp - i*para * Htmp;     beta1 = - para = -0.25 * this->dt
    // Denominator = Stmp + i*para * Htmp;   beta2 = para = 0.25 * this->dt
    std::complex<double> one = {1.0, 0.0};
    std::complex<double> beta1 = {0.0, -0.25 * this->dt};
    std::complex<double> beta2 = {0.0, 0.25 * this->dt};

    // Numerator = -i*para * Htmp
    ct::kernels::blas_scal<std::complex<double>, ct_Device>()(nlocal * nlocal,
                                                              &beta1,
                                                              Numerator.data<std::complex<double>>(),
                                                              1);
    // Numerator = Stmp + (-i*para * Htmp)
    ct::kernels::blas_axpy<std::complex<double>, ct_Device>()(nlocal * nlocal,
                                                              &one,
                                                              Stmp.data<std::complex<double>>(),
                                                              1,
                                                              Numerator.data<std::complex<double>>(),
                                                              1);
    // Denominator = i*para * Htmp
    ct::kernels::blas_scal<std::complex<double>, ct_Device>()(nlocal * nlocal,
                                                              &beta2,
                                                              Denominator.data<std::complex<double>>(),
                                                              1);
    // Denominator = Stmp + (i*para * Htmp)
    ct::kernels::blas_axpy<std::complex<double>, ct_Device>()(nlocal * nlocal,
                                                              &one,
                                                              Stmp.data<std::complex<double>>(),
                                                              1,
                                                              Denominator.data<std::complex<double>>(),
                                                              1);

    if (print_matrix)
    {
        ct::Tensor Denominator_cpu = Denominator.to_device<ct::DEVICE_CPU>();

        GlobalV::ofs_running << " beta=" << beta1 << std::endl;
        GlobalV::ofs_running << " fenmu:" << std::endl;
        for (int i = 0; i < nlocal; i++)
        {
            for (int j = 0; j < nlocal; j++)
            {
                GlobalV::ofs_running << Denominator_cpu.data<std::complex<double>>()[i * nlocal + j].real() << "+"
                                     << Denominator_cpu.data<std::complex<double>>()[i * nlocal + j].imag() << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
    }

    //->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // (3) Next, invert Denominator
    ct::Tensor ipiv(ct::DataType::DT_INT, ct_device_type, ct::TensorShape({nlocal}));
    ipiv.zero();
    // (3.1) compute ipiv
    ct::kernels::lapack_getrf<std::complex<double>, ct_Device>()(nlocal,
                                                                 nlocal,
                                                                 Denominator.data<std::complex<double>>(),
                                                                 nlocal,
                                                                 ipiv.data<int>());

    // Print ipiv
    if (print_matrix)
    {
        ct::Tensor ipiv_cpu = ipiv.to_device<ct::DEVICE_CPU>();

        GlobalV::ofs_running << " this->ParaV->nloc = " << this->ParaV->nloc << std::endl;
        GlobalV::ofs_running << " this->ParaV->nrow = " << this->ParaV->nrow << std::endl;
        GlobalV::ofs_running << " this->ParaV->ncol = " << this->ParaV->ncol << std::endl;
        GlobalV::ofs_running << " this->ParaV->nb = " << this->ParaV->nb << std::endl;
        GlobalV::ofs_running << " this->ParaV->get_block_size() = " << this->ParaV->get_block_size() << std::endl;
        GlobalV::ofs_running << " nlocal = " << nlocal << std::endl;
        GlobalV::ofs_running << " ipiv:" << std::endl;
        for (int i = 0; i < nlocal; i++)
        {
            GlobalV::ofs_running << ipiv_cpu.data<int>()[i] << " ";
        }
        GlobalV::ofs_running << std::endl;
    }

    // (3.2) compute inverse matrix of Denominator
    ct::Tensor Denominator_inv = create_identity_matrix<std::complex<double>>(nlocal, ct_device_type);
    ct::kernels::lapack_getrs<std::complex<double>, ct_Device>()('N',
                                                                 nlocal,
                                                                 nlocal,
                                                                 Denominator.data<std::complex<double>>(),
                                                                 nlocal,
                                                                 ipiv.data<int>(),
                                                                 Denominator_inv.data<std::complex<double>>(),
                                                                 nlocal);

    //->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    // (4) U_operator = Denominator_inv * Numerator;
    std::complex<double> one_gemm = {1.0, 0.0};
    std::complex<double> zero_gemm = {0.0, 0.0};
    ct::kernels::blas_gemm<std::complex<double>, ct_Device>()('N',
                                                              'N',
                                                              nlocal,
                                                              nlocal,
                                                              nlocal,
                                                              &one_gemm,
                                                              Denominator_inv.data<std::complex<double>>(),
                                                              nlocal,
                                                              Numerator.data<std::complex<double>>(),
                                                              nlocal,
                                                              &zero_gemm,
                                                              U_operator.data<std::complex<double>>(),
                                                              nlocal);

    if (print_matrix)
    {
        ct::Tensor Denominator_inv_cpu = Denominator_inv.to_device<ct::DEVICE_CPU>();
        ct::Tensor Numerator_cpu = Numerator.to_device<ct::DEVICE_CPU>();
        ct::Tensor U_operator_cpu = U_operator.to_device<ct::DEVICE_CPU>();

        GlobalV::ofs_running << " fenmu^-1:" << std::endl;
        for (int i = 0; i < nlocal; i++)
        {
            for (int j = 0; j < nlocal; j++)
            {
                GlobalV::ofs_running << Denominator_inv_cpu.data<std::complex<double>>()[i * nlocal + j].real() << "+"
                                     << Denominator_inv_cpu.data<std::complex<double>>()[i * nlocal + j].imag() << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << " fenzi:" << std::endl;
        for (int i = 0; i < nlocal; i++)
        {
            for (int j = 0; j < nlocal; j++)
            {
                GlobalV::ofs_running << Numerator_cpu.data<std::complex<double>>()[i * nlocal + j].real() << "+"
                                     << Numerator_cpu.data<std::complex<double>>()[i * nlocal + j].imag() << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << " U operator:" << std::endl;
        for (int i = 0; i < nlocal; i++)
        {
            for (int j = 0; j < nlocal; j++)
            {
                double aa, bb;
                aa = U_operator_cpu.data<std::complex<double>>()[i * nlocal + j].real();
                bb = U_operator_cpu.data<std::complex<double>>()[i * nlocal + j].imag();
                if (std::abs(aa) < 1e-8)
                {
                    aa = 0.0;
                }
                if (std::abs(bb) < 1e-8)
                {
                    bb = 0.0;
                }
                GlobalV::ofs_running << aa << "+" << bb << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
    }
}

void Propagator::compute_propagator_taylor(const int nlocal,
                                           const std::complex<double>* Stmp,
                                           const std::complex<double>* Htmp,
                                           std::complex<double>* U_operator,
                                           const int print_matrix,
                                           const int tag) const
{
    ModuleBase::GlobalFunc::ZEROS(U_operator, this->ParaV->nloc);
    std::complex<double>* A_matrix = new std::complex<double>[this->ParaV->nloc];
    ModuleBase::GlobalFunc::ZEROS(A_matrix, this->ParaV->nloc);
    std::complex<double>* rank0 = new std::complex<double>[this->ParaV->nloc];
    ModuleBase::GlobalFunc::ZEROS(rank0, this->ParaV->nloc);
    std::complex<double>* rank2 = new std::complex<double>[this->ParaV->nloc];
    ModuleBase::GlobalFunc::ZEROS(rank2, this->ParaV->nloc);
    std::complex<double>* rank3 = new std::complex<double>[this->ParaV->nloc];
    ModuleBase::GlobalFunc::ZEROS(rank3, this->ParaV->nloc);
    std::complex<double>* rank4 = new std::complex<double>[this->ParaV->nloc];
    ModuleBase::GlobalFunc::ZEROS(rank4, this->ParaV->nloc);
    std::complex<double>* tmp1 = new std::complex<double>[this->ParaV->nloc];
    ModuleBase::GlobalFunc::ZEROS(tmp1, this->ParaV->nloc);
    std::complex<double>* tmp2 = new std::complex<double>[this->ParaV->nloc];
    ModuleBase::GlobalFunc::ZEROS(tmp2, this->ParaV->nloc);
    std::complex<double>* Sinv = new std::complex<double>[this->ParaV->nloc];
    ModuleBase::GlobalFunc::ZEROS(Sinv, this->ParaV->nloc);
    BlasConnector::copy(this->ParaV->nloc, Stmp, 1, Sinv, 1);

    if (print_matrix)
    {
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << " S matrix :" << std::endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                GlobalV::ofs_running << Stmp[i * this->ParaV->ncol + j].real() << "+"
                                     << Stmp[i * this->ParaV->ncol + j].imag() << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << " H matrix :" << std::endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                GlobalV::ofs_running << Htmp[i * this->ParaV->ncol + j].real() << "+"
                                     << Htmp[i * this->ParaV->ncol + j].imag() << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
    }

    // set rank0
    int info;
    int naroc[2]; // maximum number of row or column

    for (int iprow = 0; iprow < this->ParaV->dim0; ++iprow)
    {
        for (int ipcol = 0; ipcol < this->ParaV->dim1; ++ipcol)
        {
            if (iprow == ParaV->coord[0] && ipcol == ParaV->coord[1])
            {
                naroc[0] = this->ParaV->nrow;
                naroc[1] = this->ParaV->ncol;
                for (int j = 0; j < naroc[1]; ++j)
                {
                    int igcol = globalIndex(j, this->ParaV->nb, this->ParaV->dim1, ipcol);
                    if (igcol >= nlocal)
                    {
                        continue;
                    }
                    for (int i = 0; i < naroc[0]; ++i)
                    {
                        int igrow = globalIndex(i, this->ParaV->nb, this->ParaV->dim0, iprow);
                        if (igrow >= nlocal)
                        {
                            continue;
                        }
                        if (igcol == igrow)
                        {
                            rank0[j * naroc[0] + i] = {1.0, 0.0};
                        }
                        else
                        {
                            rank0[j * naroc[0] + i] = {0.0, 0.0};
                        }
                    }
                }
            }
        } // loop ipcol
    } // loop iprow

    std::complex<double> beta = {0.0, -0.5 * this->dt / tag}; // for ETRS tag=2 , for taylor tag=1

    //->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // invert Stmp
    int* ipiv = new int[this->ParaV->nloc];
    // (3.1) compute ipiv
    ScalapackConnector::getrf(nlocal, nlocal, Sinv, 1, 1, this->ParaV->desc, ipiv, &info);
    int lwork = -1;
    int liwotk = -1;
    std::vector<std::complex<double>> work(1, 0);
    std::vector<int> iwork(1, 0);
    // (3.2) compute work
    ScalapackConnector::getri(nlocal,
                              Sinv,
                              1,
                              1,
                              this->ParaV->desc,
                              ipiv,
                              work.data(),
                              &lwork,
                              iwork.data(),
                              &liwotk,
                              &info);
    lwork = work[0].real();
    work.resize(lwork, 0);
    liwotk = iwork[0];
    iwork.resize(liwotk, 0);
    ScalapackConnector::getri(nlocal,
                              Sinv,
                              1,
                              1,
                              this->ParaV->desc,
                              ipiv,
                              work.data(),
                              &lwork,
                              iwork.data(),
                              &liwotk,
                              &info);
    assert(0 == info);

    //  A_matrix = - idt S^-1 H ;
    ScalapackConnector::gemm('N',
                             'N',
                             nlocal,
                             nlocal,
                             nlocal,
                             beta,
                             Sinv,
                             1,
                             1,
                             this->ParaV->desc,
                             Htmp,
                             1,
                             1,
                             this->ParaV->desc,
                             0.0,
                             U_operator,
                             1,
                             1,
                             this->ParaV->desc);

    //  rank2 = A^2 ;
    ScalapackConnector::gemm('N',
                             'N',
                             nlocal,
                             nlocal,
                             nlocal,
                             1.0,
                             U_operator,
                             1,
                             1,
                             this->ParaV->desc,
                             U_operator,
                             1,
                             1,
                             this->ParaV->desc,
                             0.0,
                             rank2,
                             1,
                             1,
                             this->ParaV->desc);

    //  rank3 = A^3 ;
    ScalapackConnector::gemm('N',
                             'N',
                             nlocal,
                             nlocal,
                             nlocal,
                             1.0,
                             U_operator,
                             1,
                             1,
                             this->ParaV->desc,
                             rank2,
                             1,
                             1,
                             this->ParaV->desc,
                             0.0,
                             rank3,
                             1,
                             1,
                             this->ParaV->desc);

    //  rank4 = A^4 ;
    ScalapackConnector::gemm('N',
                             'N',
                             nlocal,
                             nlocal,
                             nlocal,
                             1.0,
                             U_operator,
                             1,
                             1,
                             this->ParaV->desc,
                             rank3,
                             1,
                             1,
                             this->ParaV->desc,
                             0.0,
                             rank4,
                             1,
                             1,
                             this->ParaV->desc);

    std::complex<double> p1 = {1.0, 0.0};
    std::complex<double> p2 = {1.0 / 2.0, 0.0};
    std::complex<double> p3 = {1.0 / 6.0, 0.0};
    std::complex<double> p4 = {1.0 / 24.0, 0.0};

    ScalapackConnector::geadd('N',
                              nlocal,
                              nlocal,
                              p1,
                              rank0,
                              1,
                              1,
                              this->ParaV->desc,
                              p1,
                              U_operator,
                              1,
                              1,
                              this->ParaV->desc);

    ScalapackConnector::geadd('N',
                              nlocal,
                              nlocal,
                              p2,
                              rank2,
                              1,
                              1,
                              this->ParaV->desc,
                              p1,
                              U_operator,
                              1,
                              1,
                              this->ParaV->desc);

    ScalapackConnector::geadd('N',
                              nlocal,
                              nlocal,
                              p3,
                              rank3,
                              1,
                              1,
                              this->ParaV->desc,
                              p1,
                              U_operator,
                              1,
                              1,
                              this->ParaV->desc);

    ScalapackConnector::geadd('N',
                              nlocal,
                              nlocal,
                              p4,
                              rank4,
                              1,
                              1,
                              this->ParaV->desc,
                              p1,
                              U_operator,
                              1,
                              1,
                              this->ParaV->desc);

    if (print_matrix)
    {
        GlobalV::ofs_running << " A_matrix:" << std::endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                GlobalV::ofs_running << A_matrix[i * this->ParaV->ncol + j].real() << "+"
                                     << A_matrix[i * this->ParaV->ncol + j].imag() << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << " U operator:" << std::endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                double aa, bb;
                aa = U_operator[i * this->ParaV->ncol + j].real();
                bb = U_operator[i * this->ParaV->ncol + j].imag();
                if (std::abs(aa) < 1e-8)
                {
                    aa = 0.0;
                }
                if (std::abs(bb) < 1e-8)
                {
                    bb = 0.0;
                }
                GlobalV::ofs_running << aa << "+" << bb << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
    }
    delete[] A_matrix;
    delete[] rank0;
    delete[] rank2;
    delete[] rank3;
    delete[] rank4;
    delete[] tmp1;
    delete[] tmp2;
    delete[] Sinv;
    delete[] ipiv;
}

void Propagator::compute_propagator_etrs(const int nlocal,
                                         const std::complex<double>* Stmp,
                                         const std::complex<double>* Htmp,
                                         const std::complex<double>* H_laststep,
                                         std::complex<double>* U_operator,
                                         const int print_matrix) const
{
    std::vector<std::complex<double>> U1(this->ParaV->nloc);
    std::vector<std::complex<double>> U2(this->ParaV->nloc);
    int tag = 2;
    compute_propagator_taylor(nlocal, Stmp, Htmp, U1.data(), print_matrix, tag);
    compute_propagator_taylor(nlocal, Stmp, H_laststep, U2.data(), print_matrix, tag);
    ScalapackConnector::gemm('N',
                             'N',
                             nlocal,
                             nlocal,
                             nlocal,
                             1.0,
                             U1.data(),
                             1,
                             1,
                             this->ParaV->desc,
                             U2.data(),
                             1,
                             1,
                             this->ParaV->desc,
                             0.0,
                             U_operator,
                             1,
                             1,
                             this->ParaV->desc);
}

// Explicit instantiation of template functions
template void Propagator::compute_propagator_tensor<base_device::DEVICE_CPU>(const int nlocal,
                                                                             const ct::Tensor& Stmp,
                                                                             const ct::Tensor& Htmp,
                                                                             const ct::Tensor& H_laststep,
                                                                             ct::Tensor& U_operator,
                                                                             const int print_matrix,
                                                                             const bool use_lapack) const;
template void Propagator::compute_propagator_cn2_tensor_lapack<base_device::DEVICE_CPU>(const int nlocal,
                                                                                        const ct::Tensor& Stmp,
                                                                                        const ct::Tensor& Htmp,
                                                                                        ct::Tensor& U_operator,
                                                                                        const int print_matrix) const;
#if ((defined __CUDA) /* || (defined __ROCM) */)
template void Propagator::compute_propagator_tensor<base_device::DEVICE_GPU>(const int nlocal,
                                                                             const ct::Tensor& Stmp,
                                                                             const ct::Tensor& Htmp,
                                                                             const ct::Tensor& H_laststep,
                                                                             ct::Tensor& U_operator,
                                                                             const int print_matrix,
                                                                             const bool use_lapack) const;
template void Propagator::compute_propagator_cn2_tensor_lapack<base_device::DEVICE_GPU>(const int nlocal,
                                                                                        const ct::Tensor& Stmp,
                                                                                        const ct::Tensor& Htmp,
                                                                                        ct::Tensor& U_operator,
                                                                                        const int print_matrix) const;
#endif // __CUDA
#endif // __MPI
} // namespace module_tddft
