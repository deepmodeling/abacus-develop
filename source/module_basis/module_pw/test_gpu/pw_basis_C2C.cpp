//---------------------------------------------
// TEST for CUFFT
//---------------------------------------------
#include "../pw_basis.h"
#include "cuda_runtime.h"
#include "fftw3.h"
#include "module_base/constants.h"
#include "module_base/global_function.h"
#include "pw_test.h"

using namespace std;
TEST_F(PWTEST, recip_to_real_C2C_double)
{
    cout << "dividemthd 1, gamma_only: off, check fft between double and complex" << endl;
    ModulePW::PW_Basis pwtest("gpu", precision_flag);
    pwtest.fft_bundle.setfft("gpu", "double");
    ModuleBase::Matrix3 latvec(1, 1, 0, 0, 1, 1, 0, 0, 2);
    double wfcecut;
    double lat0 = 2.2;
    bool gamma_only = false;
    wfcecut = 18;
    gamma_only = false;
    int distribution_type = 1;
    bool xprime = false;

    // init
    pwtest.initmpi(1,0,MPI_COMM_NULL);
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

    const double tpiba2 = ModuleBase::TWO_PI * ModuleBase::TWO_PI / lat0 / lat0;
    const double ggecut = wfcecut / tpiba2;
    ModuleBase::Matrix3 GT, G, GGT;
    GT = latvec.Inverse();
    G = GT.Transpose();
    GGT = G * GT;
    complex<double>* tmp = new complex<double>[nx * ny * nz];
    if (rank_in_pool == 0)
    {
        for (int ix = 0; ix < nx; ++ix)
        {
            const double vx = ix - int(nx / 2);
            for (int iy = 0; iy < ny; ++iy)
            {
                const int offset = (ix * ny + iy) * nz;
                const double vy = iy - int(ny / 2);
                for (int iz = 0; iz < nz; ++iz)
                {
                    tmp[offset + iz] = 0.0;
                    const double vz = iz - int(nz / 2);
                    ModuleBase::Vector3<double> v(vx, vy, vz);
                    double modulus = v * (GGT * v);
                    if (modulus <= ggecut)
                    {
                        tmp[offset + iz] = 1.0 / (modulus + 1);
                        if (vy > 0)
                        {
                            tmp[offset + iz] += ModuleBase::IMAG_UNIT / (std::abs(v.x + 1) + 1);
                        }
                        else if (vy < 0)
                        {
                            tmp[offset + iz] -= ModuleBase::IMAG_UNIT / (std::abs(-v.x + 1) + 1);
                        }
                    }
                }
            }
        }
        fftw_plan pp
            = fftw_plan_dft_3d(nx, ny, nz, (fftw_complex*)tmp, (fftw_complex*)tmp, FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(pp);
        fftw_destroy_plan(pp);

        ModuleBase::Vector3<double> delta_g(double(int(nx / 2)) / nx,
                                            double(int(ny / 2)) / ny,
                                            double(int(nz / 2)) / nz);
        for (int ixy = 0; ixy < nx * ny; ++ixy)
        {
            const int ix = ixy / ny;
            const int iy = ixy % ny;
            for (int iz = 0; iz < nz; ++iz)
            {
                ModuleBase::Vector3<double> real_r(ix, iy, iz);
                double phase_im = -delta_g * real_r;
                complex<double> phase(0, ModuleBase::TWO_PI * phase_im);
                tmp[ixy * nz + iz] *= exp(phase);
            }
        }
    }
    // const int size = nx * ny * nz;
    complex<double>* h_rhog = new complex<double>[npw];
    complex<double>* h_rhogout = new complex<double>[npw];
    complex<double>* d_rhog = nullptr;
    complex<double>* d_rhogr = nullptr;
    complex<double>* d_rhogout = nullptr;
    cudaMalloc((void**)&d_rhog, npw * sizeof(complex<double>));
    cudaMalloc((void**)&d_rhogr, npw * sizeof(complex<double>));
    cudaMalloc((void**)&d_rhogout, npw * sizeof(complex<double>));

    for (int ig = 0; ig < npw; ++ig)
    {
        h_rhog[ig] = 1.0 / (pwtest.gg[ig] + 1);
        if (pwtest.gdirect[ig].y > 0)
        {
            h_rhog[ig] += ModuleBase::IMAG_UNIT / (std::abs(pwtest.gdirect[ig].x + 1) + 1);
        }
        else if (pwtest.gdirect[ig].y < 0)
        {
            h_rhog[ig] -= ModuleBase::IMAG_UNIT / (std::abs(-pwtest.gdirect[ig].x + 1) + 1);
        }
    }
    cudaMemcpy(d_rhog, h_rhog, npw * sizeof(complex<double>), cudaMemcpyHostToDevice);

    std::complex<double>* h_rhor = new std::complex<double>[nrxx];
    std::complex<double>* d_rhor = nullptr;
    cudaMalloc((void**)&d_rhor, nrxx * sizeof(std::complex<double>));
    pwtest.recip_to_real<std::complex<double>, std::complex<double>, base_device::DEVICE_GPU>(d_rhog, d_rhor);
    cudaMemcpy(h_rhor, d_rhor, nrxx * sizeof(std::complex<double>), cudaMemcpyDeviceToHost);

    int startiz = pwtest.startz_current;
    for (int ixy = 0; ixy < nx * ny; ++ixy)
    {
        const int offset = ixy * nz + startiz;
        const int startz = ixy * nplane ;
        for (int iz = 0; iz < nplane; ++iz)
        {
            EXPECT_NEAR(tmp[offset + iz].real(), h_rhor[startz + iz].real(), 1e-6);
            EXPECT_NEAR(tmp[offset + iz].imag(), h_rhor[startz + iz].imag(), 1e-6);
        }
    }

    pwtest.real_to_recip<std::complex<double>,std::complex<double>,base_device::DEVICE_GPU>(d_rhor,d_rhog);
    cudaMemcpy(h_rhogout,d_rhog,npw * sizeof(complex<double>),cudaMemcpyDeviceToHost);
    for (int ig = 0; ig < npw; ++ig)
    {
        EXPECT_NEAR(h_rhog[ig].real(), h_rhogout[ig].real(), 1e-6);
        EXPECT_NEAR(h_rhog[ig].imag(), h_rhogout[ig].imag(), 1e-6);
    }
    delete[] h_rhog;
    delete[] h_rhogout;
    delete[] h_rhor;
    delete[] tmp;
    cudaFree(d_rhog);
    cudaFree(d_rhogr);
    cudaFree(d_rhogout);
    cudaFree(d_rhor);
}

TEST_F(PWTEST, recip_to_real_C2C_float)
{
    cout << "dividemthd 1, gamma_only: off, check fft between double and complex" << endl;
    ModulePW::PW_Basis pwtest("gpu", precision_flag);
    pwtest.fft_bundle.setfft("gpu", "single");
    ModuleBase::Matrix3 latvec(1, 1, 0, 0, 1, 1, 0, 0, 2);
    double wfcecut = 18;
    double lat0 = 2.2;
    bool gamma_only = false;
    gamma_only = false;
    int distribution_type = 1;
    bool xprime = false;

    pwtest.initmpi(1,0,MPI_COMM_NULL);
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
    const double tpiba2 = ModuleBase::TWO_PI * ModuleBase::TWO_PI / lat0 / lat0;
    const double ggecut = wfcecut / tpiba2;
    ModuleBase::Matrix3 GT = latvec.Inverse();
    ModuleBase::Matrix3 G = GT.Transpose();
    ModuleBase::Matrix3 GGT = G * GT;
    complex<float>* tmp = new complex<float>[nx * ny * nz];
    if (rank_in_pool == 0)
    {
        for (int ix = 0; ix < nx; ++ix)
        {
            const float vx = ix - int(nx / 2);
            for (int iy = 0; iy < ny; ++iy)
            {
                const float vy = iy - int(ny / 2);
                const int offset = (ix * ny + iy) * nz;
                for (int iz = 0; iz < nz; ++iz)
                {
                    tmp[offset+ iz] = 0.0;
                    float vz = iz - int(nz / 2);
                    ModuleBase::Vector3<double> v(vx, vy, vz);
                    float modulus = v * (GGT * v);
                    if (modulus <= ggecut)
                    {
                        tmp[offset+ iz] = 1.0 / (modulus + 1);
                        if (vy > 0)
                            tmp[offset+ iz] += std::complex<float>(0, 1.0) / (std::abs(vx + 1) + 1);
                        else if (vy < 0)
                            tmp[offset+ iz] -= std::complex<float>(0, 1.0) / (std::abs(-vx + 1) + 1);
                    }
                }
            }
        }
        fftwf_plan pp
            = fftwf_plan_dft_3d(nx, ny, nz, (fftwf_complex*)tmp, (fftwf_complex*)tmp, FFTW_BACKWARD, FFTW_ESTIMATE);
        fftwf_execute(pp);
        fftwf_destroy_plan(pp);

        ModuleBase::Vector3<float> delta_g(float(int(nx / 2)) / nx, float(int(ny / 2)) / ny, float(int(nz / 2)) / nz);
        for (int ixy = 0; ixy < nx * ny; ++ixy)
        {
            const int ix = ixy / ny;
            const int iy = ixy % ny;
            for (int iz = 0; iz < nz; ++iz)
            {
                ModuleBase::Vector3<float> real_r(ix, iy, iz);
                float phase_im = -delta_g * real_r;
                complex<float> phase(0, ModuleBase::TWO_PI * phase_im);
                tmp[ixy * nz + iz] *= exp(phase);
            }
        }
    }
    // const int size = nx * ny * nz;
    complex<float>* h_rhog = new complex<float>[npw];
    complex<float>* h_rhogout = new complex<float>[npw];
    complex<float>* d_rhog = nullptr;
    complex<float>* d_rhogr = nullptr;
    complex<float>* d_rhogout = nullptr;
    cudaMalloc((void**)&d_rhog, npw * sizeof(complex<float>));
    cudaMalloc((void**)&d_rhogr, npw * sizeof(complex<float>));
    cudaMalloc((void**)&d_rhogout, npw * sizeof(complex<float>));

    for (int ig = 0; ig < npw; ++ig)
    {
        h_rhog[ig] = 1.0 / (pwtest.gg[ig] + 1);
        if (pwtest.gdirect[ig].y > 0)
        {
            h_rhog[ig] += ModuleBase::IMAG_UNIT / (std::abs(pwtest.gdirect[ig].x + 1) + 1);
        }
        else if (pwtest.gdirect[ig].y < 0)
        {
            h_rhog[ig] -= ModuleBase::IMAG_UNIT / (std::abs(-pwtest.gdirect[ig].x + 1) + 1);
        }
    }
    cudaMemcpy(d_rhog, h_rhog, npw * sizeof(complex<float>), cudaMemcpyHostToDevice);
    cudaMemcpy(d_rhogout, h_rhogout, npw * sizeof(complex<float>), cudaMemcpyHostToDevice);

    std::complex<float>* h_rhor = new std::complex<float>[nrxx];
    std::complex<float>* d_rhor = nullptr;
    cudaMalloc((void**)&d_rhor, nrxx * sizeof(std::complex<float>));
    pwtest.recip_to_real<std::complex<float>, std::complex<float>, base_device::DEVICE_GPU>(d_rhog, d_rhor);
    cudaMemcpy(h_rhor, d_rhor, nrxx * sizeof(std::complex<float>), cudaMemcpyDeviceToHost);

    int startiz = pwtest.startz_current;
    for (int ixy = 0; ixy < nx * ny; ++ixy)
    {
        const int offset = ixy * nz + startiz;
        const int startz = ixy * nplane;
        for (int iz = 0; iz < nplane; ++iz)
        {
            EXPECT_NEAR(tmp[offset + iz].real(), h_rhor[startz + iz].real(), 1e-4);
            EXPECT_NEAR(tmp[offset + iz].imag(), h_rhor[startz + iz].imag(), 1e-4);
        }
    }

    pwtest.real_to_recip<std::complex<float>,std::complex<float>,base_device::DEVICE_GPU>(d_rhor,d_rhog);
    cudaMemcpy(h_rhogout,d_rhog,npw * sizeof(complex<float>),cudaMemcpyDeviceToHost);
    for (int ig = 0; ig < npw; ++ig)
    {
        EXPECT_NEAR(h_rhog[ig].real(), h_rhogout[ig].real(), 1e-6);
        EXPECT_NEAR(h_rhog[ig].imag(), h_rhogout[ig].imag(), 1e-6);
    }

    delete[] h_rhog;
    delete[] h_rhogout;
    delete[] h_rhor;
    delete[] tmp;
    cudaFree(d_rhog);
    cudaFree(d_rhogr);
    cudaFree(d_rhogout);
    cudaFree(d_rhor);
}
