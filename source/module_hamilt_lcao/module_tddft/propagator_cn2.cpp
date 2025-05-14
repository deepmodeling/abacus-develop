#include "module_base/lapack_connector.h"
#include "module_base/module_container/ATen/kernels/blas.h"
#include "module_base/module_container/ATen/kernels/lapack.h"
#include "module_base/module_container/ATen/kernels/memory.h" // memory operations (Tensor)
#include "module_base/module_device/memory_op.h"              // memory operations
#include "module_base/scalapack_connector.h"
#include "module_parameter/parameter.h"
#include "propagator.h"

#include <complex>
#include <iostream>

namespace module_tddft
{
#ifdef __MPI
void Propagator::compute_propagator_cn2(const int nlocal,
                                        const std::complex<double>* Stmp,
                                        const std::complex<double>* Htmp,
                                        std::complex<double>* U_operator,
                                        std::ofstream& ofs_running,
                                        const int print_matrix) const
{
    assert(this->ParaV->nloc > 0);

    // (1) copy Htmp to Numerator & Denominator
    std::complex<double>* Numerator = new std::complex<double>[this->ParaV->nloc];
    ModuleBase::GlobalFunc::ZEROS(Numerator, this->ParaV->nloc);
    BlasConnector::copy(this->ParaV->nloc, Htmp, 1, Numerator, 1);

    std::complex<double>* Denominator = new std::complex<double>[this->ParaV->nloc];
    ModuleBase::GlobalFunc::ZEROS(Denominator, this->ParaV->nloc);
    BlasConnector::copy(this->ParaV->nloc, Htmp, 1, Denominator, 1);

    if (print_matrix)
    {
        ofs_running << std::endl;
        ofs_running << " S matrix :" << std::endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            const int in = i * this->ParaV->ncol;
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                ofs_running << Stmp[in + j].real() << "+" << Stmp[in + j].imag() << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
        ofs_running << std::endl;
        ofs_running << " H matrix :" << std::endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            const int in = i * this->ParaV->ncol;
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                ofs_running << Numerator[in + j].real() << "+" << Numerator[in + j].imag() << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
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
        ofs_running << " beta=" << beta1 << std::endl;
        ofs_running << " fenmu:" << std::endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            const int in = i * this->ParaV->ncol;
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                ofs_running << Denominator[in + j].real() << "+" << Denominator[in + j].imag() << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
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
        ofs_running << " this->ParaV->nloc = " << this->ParaV->nloc << std::endl;
        ofs_running << " this->ParaV->nrow = " << this->ParaV->nrow << std::endl;
        ofs_running << " this->ParaV->ncol = " << this->ParaV->ncol << std::endl;
        ofs_running << " this->ParaV->nb = " << this->ParaV->nb << std::endl;
        ofs_running << " this->ParaV->get_block_size() = " << this->ParaV->get_block_size() << std::endl;
        ofs_running << " nlocal = " << nlocal << std::endl;
        ofs_running << " ipiv:" << std::endl;
        for (int i = 0; i < this->ParaV->nloc; i++)
        {
            ofs_running << ipiv[i] << " ";
        }
        ofs_running << std::endl;
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
        ofs_running << " fenmu^-1:" << std::endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            const int in = i * this->ParaV->ncol;
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                ofs_running << Denominator[in + j].real() << "+" << Denominator[in + j].imag() << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
        ofs_running << " fenzi:" << std::endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            const int in = i * this->ParaV->ncol;
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                ofs_running << Numerator[in + j].real() << "+" << Numerator[in + j].imag() << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
        ofs_running << " U operator:" << std::endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            const int in = i * this->ParaV->ncol;
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                double aa = U_operator[in + j].real();
                double bb = U_operator[in + j].imag();
                if (std::abs(aa) < 1e-8)
                {
                    aa = 0.0;
                }
                if (std::abs(bb) < 1e-8)
                {
                    bb = 0.0;
                }
                ofs_running << aa << "+" << bb << "i ";
            }
            ofs_running << std::endl;
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
                                               std::ofstream& ofs_running,
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
        ofs_running << std::endl;
        ofs_running << " S matrix :" << std::endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            const int in = i * this->ParaV->ncol;
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                ofs_running << Stmp.data<std::complex<double>>()[in + j].real() << "+"
                            << Stmp.data<std::complex<double>>()[in + j].imag() << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
        ofs_running << std::endl;
        ofs_running << " H matrix :" << std::endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            const int in = i * this->ParaV->ncol;
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                ofs_running << Numerator.data<std::complex<double>>()[in + j].real() << "+"
                            << Numerator.data<std::complex<double>>()[in + j].imag() << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
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
        ofs_running << " beta=" << beta1 << std::endl;
        ofs_running << " fenmu:" << std::endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            const int in = i * this->ParaV->ncol;
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                ofs_running << Denominator.data<std::complex<double>>()[in + j].real() << "+"
                            << Denominator.data<std::complex<double>>()[in + j].imag() << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
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
        ofs_running << " this->ParaV->nloc = " << this->ParaV->nloc << std::endl;
        ofs_running << " this->ParaV->nrow = " << this->ParaV->nrow << std::endl;
        ofs_running << " this->ParaV->ncol = " << this->ParaV->ncol << std::endl;
        ofs_running << " this->ParaV->nb = " << this->ParaV->nb << std::endl;
        ofs_running << " this->ParaV->get_block_size() = " << this->ParaV->get_block_size() << std::endl;
        ofs_running << " nlocal = " << nlocal << std::endl;
        ofs_running << " ipiv:" << std::endl;
        for (int i = 0; i < this->ParaV->nloc; i++)
        {
            ofs_running << ipiv.data<int>()[i] << " ";
        }
        ofs_running << std::endl;
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
        ofs_running << " fenmu^-1:" << std::endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            const int in = i * this->ParaV->ncol;
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                ofs_running << Denominator.data<std::complex<double>>()[in + j].real() << "+"
                            << Denominator.data<std::complex<double>>()[in + j].imag() << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
        ofs_running << " fenzi:" << std::endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            const int in = i * this->ParaV->ncol;
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                ofs_running << Numerator.data<std::complex<double>>()[in + j].real() << "+"
                            << Numerator.data<std::complex<double>>()[in + j].imag() << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
        ofs_running << " U operator:" << std::endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            const int in = i * this->ParaV->ncol;
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                double aa = U_operator.data<std::complex<double>>()[in + j].real();
                double bb = U_operator.data<std::complex<double>>()[in + j].imag();
                if (std::abs(aa) < 1e-8)
                {
                    aa = 0.0;
                }
                if (std::abs(bb) < 1e-8)
                {
                    bb = 0.0;
                }
                ofs_running << aa << "+" << bb << "i ";
            }
            ofs_running << std::endl;
        }
    }
}

template <typename Device>
void Propagator::compute_propagator_cn2_tensor_lapack(const int nlocal,
                                                      const ct::Tensor& Stmp,
                                                      const ct::Tensor& Htmp,
                                                      ct::Tensor& U_operator,
                                                      std::ofstream& ofs_running,
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

        ofs_running << std::endl;
        ofs_running << " S matrix :" << std::endl;
        for (int i = 0; i < nlocal; i++)
        {
            const int in = i * nlocal;
            for (int j = 0; j < nlocal; j++)
            {
                ofs_running << Stmp_cpu.data<std::complex<double>>()[in + j].real() << "+"
                            << Stmp_cpu.data<std::complex<double>>()[in + j].imag() << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
        ofs_running << std::endl;
        ofs_running << " H matrix :" << std::endl;
        for (int i = 0; i < nlocal; i++)
        {
            const int in = i * nlocal;
            for (int j = 0; j < nlocal; j++)
            {
                ofs_running << Numerator_cpu.data<std::complex<double>>()[in + j].real() << "+"
                            << Numerator_cpu.data<std::complex<double>>()[in + j].imag() << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
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

        ofs_running << " beta=" << beta1 << std::endl;
        ofs_running << " fenmu:" << std::endl;
        for (int i = 0; i < nlocal; i++)
        {
            const int in = i * nlocal;
            for (int j = 0; j < nlocal; j++)
            {
                ofs_running << Denominator_cpu.data<std::complex<double>>()[in + j].real() << "+"
                            << Denominator_cpu.data<std::complex<double>>()[in + j].imag() << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
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

        ofs_running << " ipiv:" << std::endl;
        for (int i = 0; i < nlocal; i++)
        {
            ofs_running << ipiv_cpu.data<int>()[i] << " ";
        }
        ofs_running << std::endl;
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

        ofs_running << " fenmu^-1:" << std::endl;
        for (int i = 0; i < nlocal; i++)
        {
            const int in = i * nlocal;
            for (int j = 0; j < nlocal; j++)
            {
                ofs_running << Denominator_inv_cpu.data<std::complex<double>>()[in + j].real() << "+"
                            << Denominator_inv_cpu.data<std::complex<double>>()[in + j].imag() << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
        ofs_running << " fenzi:" << std::endl;
        for (int i = 0; i < nlocal; i++)
        {
            const int in = i * nlocal;
            for (int j = 0; j < nlocal; j++)
            {
                ofs_running << Numerator_cpu.data<std::complex<double>>()[in + j].real() << "+"
                            << Numerator_cpu.data<std::complex<double>>()[in + j].imag() << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
        ofs_running << " U operator:" << std::endl;
        for (int i = 0; i < nlocal; i++)
        {
            const int in = i * nlocal;
            for (int j = 0; j < nlocal; j++)
            {
                double aa = U_operator_cpu.data<std::complex<double>>()[in + j].real();
                double bb = U_operator_cpu.data<std::complex<double>>()[in + j].imag();
                if (std::abs(aa) < 1e-8)
                {
                    aa = 0.0;
                }
                if (std::abs(bb) < 1e-8)
                {
                    bb = 0.0;
                }
                ofs_running << aa << "+" << bb << "i ";
            }
            ofs_running << std::endl;
        }
    }
}

// Explicit instantiation of template functions
template void Propagator::compute_propagator_cn2_tensor_lapack<base_device::DEVICE_CPU>(const int nlocal,
                                                                                        const ct::Tensor& Stmp,
                                                                                        const ct::Tensor& Htmp,
                                                                                        ct::Tensor& U_operator,
                                                                                        std::ofstream& ofs_running,
                                                                                        const int print_matrix) const;
#if ((defined __CUDA) /* || (defined __ROCM) */)
template void Propagator::compute_propagator_cn2_tensor_lapack<base_device::DEVICE_GPU>(const int nlocal,
                                                                                        const ct::Tensor& Stmp,
                                                                                        const ct::Tensor& Htmp,
                                                                                        ct::Tensor& U_operator,
                                                                                        std::ofstream& ofs_running,
                                                                                        const int print_matrix) const;
#endif // __CUDA
#endif // __MPI
} // namespace module_tddft
