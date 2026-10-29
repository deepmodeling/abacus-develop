#include "module_basis/module_pw/pw_basis_k.h"
#include "module_hamilt_pw/hamilt_pwdft/operator_pw/meta_pw.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_hamilt_general/module_xc/xc_functional.h"

#include <cuda_runtime.h>
#include <mpi.h>
#include <cmath>
#include <complex>
#include <iostream>
#include <type_traits>
#include <vector>

void check(cudaError_t status)
{
    if (status != cudaSuccess) { std::cerr << cudaGetErrorString(status) << '\n'; std::exit(1); }
}

template <class Real>
size_t unused_allocation()
{
    using Complex = std::complex<Real>;
    using Device = base_device::DEVICE_GPU;
    ModulePW::PW_Basis_K basis("gpu", "double");
    basis.nmaxgr = 1 << 22;
    const int spin = 0;
    size_t before, after, total, released;
    XC_Functional::set_xc_type("PBE");
    check(cudaMemGetInfo(&before, &total));
    {
        hamilt::Meta<hamilt::OperatorPW<Complex, Device>> op(1.0, &spin, nullptr, 0, 0, &basis);
        op.act(1, 1, 1, nullptr, nullptr, 0);
        check(cudaMemGetInfo(&after, &total));
    }
    check(cudaMemGetInfo(&released, &total));
    if (released != before) { std::cerr << "Meta object retained memory\n"; std::exit(2); }
    std::cout << "element_bytes=" << sizeof(Complex) << " unused_allocated_bytes=" << before - after << std::endl;
    return before - after;
}

template <class Device>
void check_meta_action()
{
    using Complex = std::complex<double>;
    const Device* device = nullptr;
    const base_device::DEVICE_CPU* cpu = nullptr;
    const bool gpu = std::is_same<Device, base_device::DEVICE_GPU>::value;
    XC_Functional::set_xc_type("SCAN");
    ModulePW::PW_Basis_K basis(gpu ? "gpu" : "cpu", "double");
    basis.initmpi(1, 0, MPI_COMM_WORLD);
    basis.initgrids(10.0, ModuleBase::Matrix3(1, 0, 0, 0, 1, 0, 0, 0, 1), 12, 14, 16);
    ModuleBase::Vector3<double> point(0.2, -0.15, 0.1);
    basis.initparameters(false, 4.0, 1, &point);
    basis.setuptransform();
    basis.collect_local_pw();
    if (gpu) { basis.get_ig2ixyz_k(); }
    const int stride = basis.npwk_max + 7, bands = 2, spin = 0;
    const double tpiba = 2 * std::acos(-1.0) / 10.0, strength = 3.25;
    std::vector<Complex> input(bands * stride), initial(bands * stride), actual(bands * stride);
    std::vector<double> potential(basis.nrxx, strength);
    for (int i = 0; i < bands * stride; ++i)
    {
        input[i] = Complex(std::sin(0.13 * i), std::cos(0.17 * i));
        initial[i] = Complex(std::cos(0.11 * i), std::sin(0.19 * i));
    }
    Complex *in = nullptr, *out = nullptr;
    double* vk = nullptr;
    base_device::memory::resize_memory_op<Complex, Device>()(device, in, input.size());
    base_device::memory::resize_memory_op<Complex, Device>()(device, out, input.size());
    base_device::memory::resize_memory_op<double, Device>()(device, vk, potential.size());
    base_device::memory::synchronize_memory_op<Complex, Device, base_device::DEVICE_CPU>()(device, cpu, in, input.data(), input.size());
    base_device::memory::synchronize_memory_op<Complex, Device, base_device::DEVICE_CPU>()(device, cpu, out, initial.data(), initial.size());
    base_device::memory::synchronize_memory_op<double, Device, base_device::DEVICE_CPU>()(device, cpu, vk, potential.data(), potential.size());
    double error = 0;
    {
        hamilt::Meta<hamilt::OperatorPW<Complex, Device>> op(tpiba, &spin, vk, 1, basis.nrxx, &basis);
        op.init(0);
        for (int repeat = 1; repeat <= 2; ++repeat)
        {
            op.act(bands, stride, 1, in, out, basis.npwk[0]);
            base_device::memory::synchronize_memory_op<Complex, base_device::DEVICE_CPU, Device>()(cpu, device, actual.data(), out, actual.size());
            for (int band = 0; band < bands; ++band)
            {
                for (int ig = 0; ig < stride; ++ig)
                {
                    const int i = band * stride + ig;
                    Complex expected = initial[i];
                    if (ig < basis.npwk[0])
                    {
                        expected += double(repeat) * strength * tpiba * tpiba
                                  * (basis.getgcar(0, ig) + basis.kvec_c[0]).norm2() * input[i];
                    }
                    if (!std::isfinite(actual[i].real()) || !std::isfinite(actual[i].imag())) { std::exit(3); }
                    error = std::max(error, std::abs(actual[i] - expected) / std::max(1.0, std::abs(expected)));
                }
            }
        }
    }
    base_device::memory::synchronize_memory_op<Complex, base_device::DEVICE_CPU, Device>()(cpu, device, actual.data(), in, actual.size());
    if (actual != input || error > 1e-12) { std::cerr << "Meta action error " << error << '\n'; std::exit(4); }
    base_device::memory::delete_memory_op<Complex, Device>()(device, in);
    base_device::memory::delete_memory_op<Complex, Device>()(device, out);
    base_device::memory::delete_memory_op<double, Device>()(device, vk);
    std::cout << "gpu=" << gpu << " meta_relative_error=" << error << std::endl;
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    GlobalV::NSPIN = GlobalV::NPOL = 1;
    GlobalV::BASIS_TYPE = "pw";
    check(cudaFree(nullptr));
    const size_t unused = unused_allocation<float>() + unused_allocation<double>();
    check_meta_action<base_device::DEVICE_CPU>();
    check_meta_action<base_device::DEVICE_GPU>();
    MPI_Finalize();
    if (unused != 0) { std::cerr << "PBE allocated unused Meta-GGA storage\n"; return 5; }
}
