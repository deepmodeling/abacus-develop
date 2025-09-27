#include "cuda_runtime.h"
#include "fftw3.h"
#include "source_base/module_device/device.h"
#include "source_base/vector3.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_basis/module_pw/pw_basis_k.h"

#include <complex>
#include <gtest/gtest.h>
#include <typeinfo>

using namespace std;
template <typename DataType, typename DeviceType>
struct TypePair
{
    using T = DataType;
    using Device = DeviceType;
};

template <typename TypePair>
class PW_BASIS_K_GPU_TEST : public ::testing::Test
{
  public:
    using T = typename TypePair::T;
    using Device = typename TypePair::Device;
    ModulePW::PW_Basis pwtest;
    complex<T>* d_rhog = nullptr;
    complex<T>* d_rhogr = nullptr;
    complex<T>* d_rhogout = nullptr;
    complex<T>* d_rhor = nullptr;
    complex<T>* tmp = nullptr;
    complex<T>* h_rhog = nullptr;
    complex<T>* h_rhogout = nullptr;
    complex<T>* h_rhor = nullptr;
    void init(ModulePW::PW_Basis& pwtest)
    {
        ModuleBase::Matrix3 latvec(1, 1, 0, 0, 1, 1, 0, 0, 2);
        T wfcecut;
        T lat0 = 2.2;

        bool gamma_only = false;
        wfcecut = 18;
        gamma_only = false;
        int distribution_type = 1;
        bool xprime = false;
        const int nks = 1;
        // init
        const int mypool = 0;
        const int key = 1;
        const int nproc_in_pool = 1;
        const int rank_in_pool = 0;
        MPI_Comm POOL_WORLD;
        MPI_Comm_split(MPI_COMM_WORLD, mypool, key, &POOL_WORLD);
        pwtest.initmpi(nproc_in_pool, rank_in_pool, POOL_WORLD);
        pwtest.initgrids(lat0, latvec, wfcecut);
        pwtest.initparameters(gamma_only, wfcecut, distribution_type, xprime);
        pwtest.setuptransform();
        pwtest.collect_local_pw();

        const int npw = pwtest.npw;
        const int nrxx = pwtest.nrxx;
        const int nmaxgr = pwtest.nmaxgr;
        const int nx = pwtest.nx;
        const int ny = pwtest.ny;
        const int nz = pwtest.nz;
        const int nplane = pwtest.nplane;

        const T tpiba2 = ModuleBase::TWO_PI * ModuleBase::TWO_PI / lat0 / lat0;
        const T ggecut = wfcecut / tpiba2;
        ModuleBase::Matrix3 GT, G, GGT;
        GT = latvec.Inverse();
        G = GT.Transpose();
        GGT = G * GT;
        tmp = new complex<T>[nx * ny * nz];
            if (rank_in_pool == 0)
            {
                for (int ix = 0; ix < nx; ++ix)
                {
                    const T vx = ix - int(nx / 2);
                    for (int iy = 0; iy < ny; ++iy)
                    {
                        const int offset = (ix * ny + iy) * nz;
                        const T vy = iy - int(ny / 2);
                        for (int iz = 0; iz < nz; ++iz)
                        {
                            tmp[offset + iz] = 0.0;
                            T vz = iz - int(nz / 2);
                            ModuleBase::Vector3<double> v(vx, vy, vz);
                            T modulus = v * (GGT * v);
                            if (modulus <= ggecut)
                            {
                                tmp[offset + iz] = 1.0 / (modulus + 1);
                                if (vy > 0)
                                {
                                    tmp[offset + iz]
                                        += std::complex<T>(0, 1.0) / (std::abs(static_cast<T>(v.x) + 1) + 1);
                                }
                                else if (vy < 0)
                                {
                                    tmp[offset + iz]
                                        -= std::complex<T>(0, 1.0) / (std::abs(-static_cast<T>(v.x) + 1) + 1);
                                }
                            }
                        }
                    }
                }
                if (typeid(T) == typeid(double))
                {
                    fftw_plan pp = fftw_plan_dft_3d(nx,
                                                    ny,
                                                    nz,
                                                    (fftw_complex*)tmp,
                                                    (fftw_complex*)tmp,
                                                    FFTW_BACKWARD,
                                                    FFTW_ESTIMATE);
                    fftw_execute(pp);
                    fftw_destroy_plan(pp);
                }
                else if (typeid(T) == typeid(float))
                {
                    fftwf_plan pp = fftwf_plan_dft_3d(nx,
                                                      ny,
                                                      nz,
                                                      (fftwf_complex*)tmp,
                                                      (fftwf_complex*)tmp,
                                                      FFTW_BACKWARD,
                                                      FFTW_ESTIMATE);
                    fftwf_execute(pp);
                    fftwf_destroy_plan(pp);
                }
                ModuleBase::Vector3<T> delta_g(T(int(nx / 2)) / nx, T(int(ny / 2)) / ny, T(int(nz / 2)) / nz);
                for (int ixy = 0; ixy < nx * ny; ++ixy)
                {
                    const int ix = ixy / ny;
                    const int iy = ixy % ny;
                    for (int iz = 0; iz < nz; ++iz)
                    {
                        ModuleBase::Vector3<T> real_r(ix, iy, iz);
                        T phase_im = -delta_g * real_r;
                        complex<T> phase(0, ModuleBase::TWO_PI * phase_im);
                        tmp[ixy * nz + iz] *= exp(phase);
                    }
                }

                h_rhog = new complex<T>[npw];
                h_rhogout = new complex<T>[npw];
                for (int ig = 0; ig < npw; ++ig)
                {
                    h_rhog[ig] = 1.0 / (pwtest.gg[ig] + 1);
                   
                    if (pwtest.gdirect[ig].y > 0)
                    {
                        h_rhog[ig] += std::complex<float>(0, 1.0) / (std::abs(float(pwtest.gdirect[ig].x) + 1) + 1);
                    }
                    else if (pwtest.gdirect[ig].y < 0)
                    {
                        h_rhog[ig] -= std::complex<float>(0, 1.0) / (std::abs(float(-pwtest.gdirect[ig].x) + 1) + 1);
                    }
                }

                cudaMalloc((void**)&d_rhog, npw * sizeof(complex<T>));
                cudaMalloc((void**)&d_rhor, nrxx * sizeof(complex<T>));
                cudaMemcpy(d_rhog, h_rhog, npw * sizeof(complex<T>), cudaMemcpyHostToDevice);

                h_rhor = new complex<T>[nrxx];

                pwtest.recip_to_real<std::complex<T>, std::complex<T>,base_device::DEVICE_GPU>(d_rhog, d_rhor);
                cudaMemcpy(h_rhor, d_rhor, nrxx * sizeof(complex<T>), cudaMemcpyDeviceToHost);

                pwtest.real_to_recip<std::complex<T>, std::complex<T>,base_device::DEVICE_GPU>(d_rhor, d_rhog);
                cudaMemcpy(h_rhogout, d_rhog, npw * sizeof(complex<T>), cudaMemcpyDeviceToHost);
            }
    }
    ModulePW::PW_Basis* access_pw()
    {
        return &pwtest;
    }
    void TearDown() override
    {
        delete[] h_rhog;
        delete[] h_rhogout;
        delete[] h_rhor;
        delete[] tmp;
        cudaFree(d_rhog);
        cudaFree(d_rhogr);
        cudaFree(d_rhogout);
        cudaFree(d_rhor);
    }
};

using MixedTypes = ::testing::Types<TypePair<float, base_device::DEVICE_GPU>, 
                                    TypePair<double, base_device::DEVICE_GPU> >;

TYPED_TEST_CASE(PW_BASIS_K_GPU_TEST, MixedTypes);

TYPED_TEST(PW_BASIS_K_GPU_TEST, Mixing)
{
    using T = typename TestFixture::T;
    using Device = typename TestFixture::Device;
    ModulePW::PW_Basis pwtest;
    pwtest.set_device("gpu");
    pwtest.set_precision("mixing");
    pwtest.fft_bundle.setfft("gpu", "mixing");
    this->init(pwtest);
    int startiz = pwtest.startz_current;
    const int nx = pwtest.nx;
    const int ny = pwtest.ny;
    const int nz = pwtest.nz;
    const int nplane = pwtest.nplane;
    const int npw = pwtest.npw;
    for (int ixy = 0; ixy < nx * ny; ++ixy)
    {
        const int offset = ixy * nz + startiz;
        const int startz = ixy * nplane;
        for (int iz = 0; iz < nplane; ++iz)
        {
            EXPECT_NEAR(this->tmp[offset + iz].real(), this->h_rhor[startz + iz].real(), 1e-4);
        }
    }
    for (int ig = 0; ig < npw; ++ig)
    {
        EXPECT_NEAR(this->h_rhog[ig].real(), this->h_rhogout[ig].real(), 1e-4);
        EXPECT_NEAR(this->h_rhog[ig].imag(), this->h_rhogout[ig].imag(), 1e-4);
    }
}

TYPED_TEST(PW_BASIS_K_GPU_TEST, FloatDouble)
{
    using T = typename TestFixture::T;
    using Device = typename TestFixture::Device;
    ModulePW::PW_Basis pwtest;
    pwtest.set_device("gpu");
    pwtest.set_precision("mixing");
    if (typeid(T) == typeid(float))
    {
        pwtest.fft_bundle.setfft("gpu", "single");
    }
    else if (typeid(T) == typeid(double))
    {
        pwtest.fft_bundle.setfft("gpu", "double");
    }
    else
    {
        cout << "Error: Unsupported type" << endl;
        return;
    }
    this->init(pwtest);
    int startiz = pwtest.startz_current;
    const int nx = pwtest.nx;
    const int ny = pwtest.ny;
    const int nz = pwtest.nz;
    const int nplane = pwtest.nplane;
    const int npw = pwtest.npw;
    for (int ixy = 0; ixy < nx * ny; ++ixy)
    {
        const int offset = ixy * nz + startiz;
        const int startz = ixy * nplane;
        for (int iz = 0; iz < nplane; ++iz)
        {
            EXPECT_NEAR(this->tmp[offset + iz].real(), this->h_rhor[startz + iz].real(), 1e-4);
        }
    }

    for (int ig = 0; ig < npw; ++ig)
    {
        EXPECT_NEAR(this->h_rhog[ig].real(), this->h_rhogout[ig].real(), 1e-4);
        EXPECT_NEAR(this->h_rhog[ig].imag(), this->h_rhogout[ig].imag(), 1e-4);
    }
}
