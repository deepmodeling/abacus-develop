//---------------------------------------------
// TEST for FFT
//---------------------------------------------
#include "../pw_basis_k.h"
#include "cuda_runtime.h"
#include "module_base/constants.h"
#include "module_base/global_function.h"
#include "pw_test.h"
using namespace std;
TEST_F(PWTEST, pw_basis_k_C2C_double)
{
    cout << "dividemthd 1, gamma_only: on, xprime: false, gamma kpoint, check fft" << endl;
    ModulePW::PW_Basis_K pwtest("gpu", "double");
    ModuleBase::Matrix3 latvec(1, 0.3, 0, 0, 2, 0, 0, 0, 2);
    double wfcecut =10;
    double lat0 = 2.7;
    bool gamma_only;
    ModuleBase::Vector3<double>* kvec_d;
    int nks = 1;
    kvec_d = new ModuleBase::Vector3<double>[nks];
    kvec_d[0].set(0, 0, 0);
    wfcecut = 10;
    gamma_only = false;
    int distribution_type = 1;
    bool xprime = false;
    //--------------------------------------------------
    // init //real parameter
    const int mypool = 0;
    const int key = 1;
    const int nproc_in_pool = 1;
    const int rank_in_pool = 0;
    MPI_Comm POOL_WORLD;
    MPI_Comm_split(MPI_COMM_WORLD,mypool,key,&POOL_WORLD);
    pwtest.initmpi(nproc_in_pool, rank_in_pool, POOL_WORLD);
    pwtest.initgrids(lat0, latvec, 4 * wfcecut);
    pwtest.initparameters(gamma_only, wfcecut, nks, kvec_d, distribution_type, xprime);
    pwtest.setuptransform();
    pwtest.collect_local_pw();

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
    double* rhor = new double[nrxx];
    for (int ik = 0; ik < nks; ++ik)
    {
        int npwk = pwtest.npwk[ik];
        if (rank_in_pool == 0)
        {
            ModuleBase::Vector3<double> kk = kvec_d[ik];
            for (int ix = 0; ix < nx; ++ix)
            {
                const double vx = ix - int(nx / 2);
                for (int iy = 0; iy < ny; ++iy)
                {
                    const double vy = iy - int(ny / 2);
                    const int offset = (ix * ny + iy) * nz;
                    for (int iz = 0; iz < nz; ++iz)
                    {
                        tmp[offset+ iz] = 0.0;
                        double vz = iz - int(nz / 2);
                        ModuleBase::Vector3<double> v(vx, vy, vz);
                        double modulusgk = (v + kk) * (GGT * (v + kk));
                        if (modulusgk <= ggecut)
                        {
                            tmp[offset+ iz] = 1.0 / (modulusgk + 1);
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
        complex<double>* h_rhog = new complex<double>[npwk];
        complex<double>* h_rhogout = new complex<double>[npwk];
        complex<double>* h_rhor = new complex<double>[nrxx];
        for (int ig = 0; ig < npwk; ++ig)
        {
            h_rhog[ig] = 1.0 / (pwtest.getgk2(ik, ig) + 1);
            ModuleBase::Vector3<double> f = pwtest.getgdirect(ik, ig);
            if (f.y > 0)
            {
                h_rhog[ig] += ModuleBase::IMAG_UNIT / (std::abs(f.x + 1) + 1);
            }
            else if (f.y < 0)
            {
                h_rhog[ig] -= ModuleBase::IMAG_UNIT / (std::abs(-f.x + 1) + 1);
            }
        }
        complex<double>* d_rhog = nullptr;
        complex<double>* d_rhor = nullptr;
        cudaMalloc((void**)&d_rhog, npwk * sizeof(complex<double>));
        cudaMalloc((void**)&d_rhor, nrxx * sizeof(complex<double>));
        cudaMemcpy(d_rhog, h_rhog, npwk * sizeof(complex<double>), cudaMemcpyHostToDevice);
        pwtest.recip_to_real<std::complex<double>, base_device::DEVICE_GPU>(d_rhog,
                                                                            d_rhor,
                                                                            ik); 
        cudaMemcpy(h_rhor, d_rhor, nrxx * sizeof(complex<double>), cudaMemcpyDeviceToHost);
        int startiz = pwtest.startz_current;
        for (int ixy = 0; ixy < nx * ny; ++ixy)
        {
            const int offset = ixy * nz + startiz;
            const int startz = ixy * nplane;
            for (int iz = 0; iz < nplane; ++iz)
            {
                EXPECT_NEAR(tmp[offset + iz].real(), h_rhor[startz + iz].real(), 1e-6);
                EXPECT_NEAR(tmp[offset + iz].imag(), h_rhor[startz + iz].imag(), 1e-6);
            }
        }

        pwtest.real_to_recip<std::complex<double>,base_device::DEVICE_GPU>(d_rhor,d_rhog,ik);
        cudaMemcpy(h_rhogout,d_rhog,npwk * sizeof(complex<double>),cudaMemcpyDeviceToHost);
        for (int ig = 0; ig < npwk; ++ig)
        {
            EXPECT_NEAR(h_rhog[ig].real(), h_rhogout[ig].real(), 1e-6);
            EXPECT_NEAR(h_rhog[ig].imag(), h_rhogout[ig].imag(), 1e-6);
        }
        delete[] h_rhog;
        delete[] h_rhor;
        delete[] h_rhogout;
        cudaFree(d_rhor);
        cudaFree(d_rhog);
    }
    delete[] tmp;
    delete[] rhor;
    delete[] kvec_d;
    fftw_cleanup();
}

TEST_F(PWTEST, pw_basis_k_C2C_float)
{
    cout << "dividemthd 1, gamma_only: on, xprime: false, gamma kpoint, check fft" << endl;
    ModulePW::PW_Basis_K pwtest("gpu", "single");
    ModuleBase::Matrix3 latvec(1, 0.3, 0, 0, 2, 0, 0, 0, 2);
    double wfcecut =10;
    double lat0 = 2.7;
    bool gamma_only;
    ModuleBase::Vector3<double>* kvec_d;
    int nks = 1;
    kvec_d = new ModuleBase::Vector3<double>[nks];
    kvec_d[0].set(0, 0, 0);
    wfcecut = 10;
    gamma_only = false;
    int distribution_type = 1;
    bool xprime = false;
    //--------------------------------------------------
    // init //real parameter
    const int mypool = 0;
    const int key = 1;
    const int nproc_in_pool = 1;
    const int rank_in_pool = 0;
    MPI_Comm POOL_WORLD;
    MPI_Comm_split(MPI_COMM_WORLD,mypool,key,&POOL_WORLD);
    pwtest.initmpi(nproc_in_pool, rank_in_pool, POOL_WORLD);

    pwtest.initgrids(lat0, latvec, 4 * wfcecut);
    pwtest.initparameters(gamma_only, wfcecut, nks, kvec_d, distribution_type, xprime);
    pwtest.setuptransform();
    pwtest.collect_local_pw();

    const int nrxx = pwtest.nrxx;
    const int nmaxgr = pwtest.nmaxgr;
    const int nx = pwtest.nx;
    const int ny = pwtest.ny;
    const int nz = pwtest.nz;
    const int nplane = pwtest.nplane;
    const float tpiba2 = ModuleBase::TWO_PI * ModuleBase::TWO_PI / lat0 / lat0;
    const float ggecut = wfcecut / tpiba2;
    ModuleBase::Matrix3 GT, G, GGT;
    GT = latvec.Inverse();
    G = GT.Transpose();
    GGT = G * GT;
    complex<float>* tmp = new complex<float>[nx * ny * nz];
    for (int ik = 0; ik < nks; ++ik)
    {
        int npwk = pwtest.npwk[ik];
        if (rank_in_pool == 0)
        {
            ModuleBase::Vector3<double> kk = kvec_d[ik];
            for (int ix = 0; ix < nx; ++ix)
            {
                const float vx = ix - int(nx / 2);
                for (int iy = 0; iy < ny; ++iy)
                {
                    const int offset = (ix * ny + iy) * nz;
                    const float vy = iy - int(ny / 2);
                    for (int iz = 0; iz < nz; ++iz)
                    {
                        tmp[offset+ iz] = 0.0;
                        double vz = iz - int(nz / 2);
                        ModuleBase::Vector3<double> v(vx, vy, vz);
                        float modulusgk = float((v + kk) * (GGT * (v + kk)));
                        if (modulusgk <= ggecut)
                        {
                            tmp[offset+ iz] = float(1.0 / (modulusgk + 1));
                            if (vy > 0)
                            {
                                tmp[offset+ iz] += std::complex<float>(0,1.0)  / (std::abs(float(v.x) + 1) + 1);
                            }
                            else if (vy < 0)
                            {
                                tmp[offset+ iz] -= std::complex<float>(0,1.0)  / (std::abs(float(-v.x) + 1) + 1);
                            }
                        }
                    }
                }
            }
            fftwf_plan pp
                = fftwf_plan_dft_3d(nx, ny, nz, (fftwf_complex*)tmp, (fftwf_complex*)tmp, FFTW_BACKWARD, FFTW_ESTIMATE);
            fftwf_execute(pp);
            fftwf_destroy_plan(pp);

            ModuleBase::Vector3<float> delta_g(float(int(nx / 2)) / nx,
                                                float(int(ny / 2)) / ny,
                                                float(int(nz / 2)) / nz);
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
        complex<float>* h_rhog = new complex<float>[npwk];
        complex<float>* h_rhogout = new complex<float>[npwk];
        complex<float>* h_rhor = new complex<float>[nrxx];
        for (int ig = 0; ig < npwk; ++ig)
        {
            h_rhog[ig] = float(1.0 / (pwtest.getgk2(ik, ig) + 1));
            ModuleBase::Vector3<double> f = pwtest.getgdirect(ik, ig);
            if (f.y > 0)
            {
                h_rhog[ig] += std::complex<float>(0,1.0) / (std::abs(float(f.x) + 1) + 1);
            }
            else if (f.y < 0)
            {
                h_rhog[ig] -= std::complex<float>(0,1.0)  / (std::abs(float(-f.x) + 1) + 1);
            }
        }
        complex<float>* d_rhog = nullptr;
        complex<float>* d_rhor = nullptr;
        cudaMalloc((void**)&d_rhog, npwk * sizeof(complex<float>));
        cudaMalloc((void**)&d_rhor, nrxx * sizeof(complex<float>));
        cudaMemcpy(d_rhog, h_rhog, npwk * sizeof(complex<float>), cudaMemcpyHostToDevice);
        pwtest.recip_to_real<std::complex<float>, base_device::DEVICE_GPU>(d_rhog,
                                                                           d_rhor,
                                                                           ik); 
        cudaMemcpy(h_rhor, d_rhor, nrxx * sizeof(complex<float>), cudaMemcpyDeviceToHost);
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

        pwtest.real_to_recip<std::complex<float>,base_device::DEVICE_GPU>(d_rhor,d_rhog,ik);
        cudaMemcpy(h_rhogout,d_rhog,npwk * sizeof(complex<float>),cudaMemcpyDeviceToHost);
        for (int ig = 0; ig < npwk; ++ig)
        {
            EXPECT_NEAR(h_rhog[ig].real(), h_rhogout[ig].real(), 1e-4);
            EXPECT_NEAR(h_rhog[ig].imag(), h_rhogout[ig].imag(), 1e-4);
        }
        delete[] h_rhog;
        delete[] h_rhor;
        delete[] h_rhogout;
        cudaFree(d_rhor);
        cudaFree(d_rhog);
    }
    delete[] tmp;
    delete[] kvec_d;
    fftw_cleanup();
}
