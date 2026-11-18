#include "../libxc_abacus.h"
#include "source_base/module_device/device.h"
#include "gtest/gtest.h"

#include <algorithm>
#include <cmath>
#include <tuple>
#include <vector>

#if defined(__CUDA) && defined(ABACUS_LIBXC_CUDA)
namespace
{
class LibxcCudaPotential : public testing::TestWithParam<std::tuple<int, int, int>>
{
};

TEST_P(LibxcCudaPotential, MatchesHostAcrossBatchBoundaries)
{
    if (!base_device::information::probe_gpu_availability())
    {
        GTEST_SKIP() << "CUDA-enabled ABACUS and Libxc are required.";
    }
    const std::tuple<int, int, int> parameters = GetParam();
    const int id = std::get<0>(parameters);
    const int spin = std::get<1>(parameters);
    const int points = std::get<2>(parameters);
    const int nsigma = spin == 1 ? 1 : 3;
    const std::vector<int> ids = {id};
    std::vector<xc_func_type> host = XC_Functional_Libxc::init_func(ids, spin, 0.31, 0.09);
    std::vector<xc_func_type> device = XC_Functional_Libxc::init_func(ids, spin, 0.31, 0.09, true);
    const bool is_lda = host[0].info->family == XC_FAMILY_LDA;
    const std::size_t rho_size = static_cast<std::size_t>(points) * spin;
    const std::size_t sigma_size = static_cast<std::size_t>(points) * nsigma;
    std::vector<double> rho(rho_size);
    std::vector<double> sigma(sigma_size);
    for (int i = 0; i < points; ++i)
    {
        // Include vacuum, threshold, ordinary and strongly polarized densities.
        const double density = i % 17 == 0 ? 0.0 : (i % 19 == 0 ? 1.0e-12 : 0.05 + 0.001 * (i % 997));
        rho[i * spin] = density;
        if (spin == 2)
        {
            rho[i * spin + 1] = density * (i % 23 == 0 ? 1.0e-8 : 0.4);
        }
        for (int j = 0; j < nsigma; ++j)
        {
            sigma[i * nsigma + j] = density * density * (j == 1 ? 0.1 : 0.3);
        }
    }
    std::vector<double> exc_host(points, 23.0);
    std::vector<double> vrho_host(rho_size, 23.0);
    std::vector<double> vsigma_host(sigma_size, 23.0);
    std::vector<double> exc_device(points, 23.0);
    std::vector<double> vrho_device(rho_size, 23.0);
    std::vector<double> vsigma_device(sigma_size, 23.0);
    const double* rho_data = rho.data();
    const double* sigma_data = sigma.data();
    double* exc_host_data = exc_host.data();
    double* vrho_host_data = vrho_host.data();
    double* vsigma_host_data = vsigma_host.data();
    double* exc_device_data = exc_device.data();
    double* vrho_device_data = vrho_device.data();
    double* vsigma_device_data = vsigma_device.data();
    xc_func_set_dens_threshold(&host[0], 1.0e-10);
    xc_func_set_dens_threshold(&device[0], 1.0e-10);
    if (is_lda)
    {
        xc_lda_exc_vxc(&host[0], points, rho_data, exc_host_data, vrho_host_data);
    }
    else
    {
        xc_gga_exc_vxc(&host[0], points, rho_data, sigma_data, exc_host_data, vrho_host_data, vsigma_host_data);
    }
    XC_Functional_Libxc::xc_cuda_exc_vxc(device[0], points, rho_data, sigma_data,
                                        exc_device_data, vrho_device_data, vsigma_device_data);
    for (std::size_t i = 0; i < exc_host.size(); ++i)
    {
        const double tolerance = 2.0e-11 * std::max(1.0, std::abs(exc_host[i]));
        ASSERT_NEAR(exc_device[i], exc_host[i], tolerance) << i;
    }
    for (std::size_t i = 0; i < vrho_host.size(); ++i)
    {
        const double tolerance = 2.0e-10 * std::max(1.0, std::abs(vrho_host[i]));
        ASSERT_NEAR(vrho_device[i], vrho_host[i], tolerance) << i;
    }
    if (!is_lda)
    {
        for (std::size_t i = 0; i < vsigma_host.size(); ++i)
        {
            const double tolerance = 2.0e-10 * std::max(1.0, std::abs(vsigma_host[i]));
            ASSERT_NEAR(vsigma_device[i], vsigma_host[i], tolerance) << i;
        }
    }
    XC_Functional_Libxc::finish_func(host);
    XC_Functional_Libxc::finish_func(device);
}

INSTANTIATE_TEST_SUITE_P(LdaGgaHse,
                        LibxcCudaPotential,
                        testing::Combine(testing::Values(XC_LDA_X, XC_GGA_C_PBE, XC_HYB_GGA_XC_HSE06),
                                         testing::Values(XC_UNPOLARIZED, XC_POLARIZED),
                                         testing::Values(1, 65536, 65537, 131075)));
} // namespace
#endif
