#include "../pw_basis.h"
#ifdef __MPI
#include "source_base/parallel_global.h"
#include "mpi.h"
#include "test_tool.h"
#endif
#include "source_base/global_function.h"
#include "pw_test.h"

extern int nproc_in_pool, rank_in_pool;

namespace
{
class PW_Basis_Comm_Accessor : public ModulePW::PW_Basis
{
  public:
    PW_Basis_Comm_Accessor(const std::string& device_, const std::string& precision_)
        : ModulePW::PW_Basis(device_, precision_)
    {
    }

    using ModulePW::PW_Basis::gatherp_scatters;
    using ModulePW::PW_Basis::gathers_scatterp;
};

template <typename BasisType>
void zero_complex_buffer(std::complex<double>* data, const int size)
{
    for (int i = 0; i < size; ++i)
    {
        data[i] = std::complex<double>(0.0, 0.0);
    }
}

template <typename BasisType>
void fill_plane_major_sticks(const BasisType& pw, std::complex<double>* plane)
{
    zero_complex_buffer<BasisType>(plane, pw.nrxx);
    for (int istot = 0; istot < pw.nstot; ++istot)
    {
        const int ixy = pw.istot2ixy[istot];
        for (int iz = 0; iz < pw.nplane; ++iz)
        {
            const int gz = pw.startz_current + iz;
            const double real = (rank_in_pool + 1) * 1000.0 + istot * 10.0 + gz;
            const double imag = (ixy + 1) * 0.25 + gz * 0.5;
            plane[ixy * pw.nplane + iz] = std::complex<double>(real, imag);
        }
    }
}

template <typename BasisType>
void expect_plane_major_equal(const BasisType& pw,
                              const std::complex<double>* expected,
                              const std::complex<double>* actual)
{
    for (int ir = 0; ir < pw.nrxx; ++ir)
    {
        EXPECT_DOUBLE_EQ(expected[ir].real(), actual[ir].real());
        EXPECT_DOUBLE_EQ(expected[ir].imag(), actual[ir].imag());
    }
}

template <typename BasisType>
int comm_roundtrip_work_size(const BasisType& pw)
{
    const int gather_size = pw.nst * pw.nz;
    const int scatter_recv_size = pw.startr[pw.poolnproc - 1] + pw.numr[pw.poolnproc - 1];
    return std::max(gather_size, scatter_recv_size);
}

template <typename BasisType>
void expect_stick_major_equal(const BasisType& pw, const std::complex<double>* sticks)
{
    int istot0 = 0;
    for (int ip = 0; ip < pw.poolrank; ++ip)
    {
        istot0 += pw.nst_per[ip];
    }

    for (int is = 0; is < pw.nst; ++is)
    {
        const int global_istot = istot0 + is;
        const int ixy = pw.is2fftixy[is];
        EXPECT_EQ(ixy, pw.istot2ixy[global_istot]);
        for (int iz = 0; iz < pw.nz; ++iz)
        {
            int owner = -1;
            for (int ip = 0; ip < pw.poolnproc; ++ip)
            {
                if (iz >= pw.startz[ip] && iz < pw.startz[ip] + pw.numz[ip])
                {
                    owner = ip;
                    break;
                }
            }
            EXPECT_GE(owner, 0);
            const double real = (owner + 1) * 1000.0 + global_istot * 10.0 + iz;
            const double imag = (ixy + 1) * 0.25 + iz * 0.5;
            EXPECT_DOUBLE_EQ(real, sticks[is * pw.nz + iz].real());
            EXPECT_DOUBLE_EQ(imag, sticks[is * pw.nz + iz].imag());
        }
    }
}

bool case_has_zero_plane_stress(const int nx, const int ny, const int nz)
{
    PW_Basis_Comm_Accessor pwtest(device_flag, precision_flag);
    ModuleBase::Matrix3 latvec(1, 0, 0, 0, 1, 0, 0, 0, 1);
    const double lat0 = 4.0;
    const double wfcecut = 20.0;

#ifdef __MPI
    pwtest.initmpi(nproc_in_pool, rank_in_pool, POOL_WORLD);
#endif
    pwtest.initgrids(lat0, latvec, nx, ny, nz);
    pwtest.initparameters(false, wfcecut, 1, true);
    pwtest.setuptransform();

#ifdef __MPI
    const int local_stress = (pwtest.nplane == 0 && pwtest.nst > 0) ? 1 : 0;
    int any_stress = 0;
    MPI_Allreduce(&local_stress, &any_stress, 1, MPI_INT, MPI_MAX, POOL_WORLD);
    return any_stress == 1;
#else
    return false;
#endif
}

void run_comm_roundtrip_case(const int nx, const int ny, const int nz)
{
    PW_Basis_Comm_Accessor pwtest(device_flag, precision_flag);
    ModuleBase::Matrix3 latvec(1, 0, 0, 0, 1, 0, 0, 0, 1);
    const double lat0 = 4.0;
    const double wfcecut = 20.0;

#ifdef __MPI
    pwtest.initmpi(nproc_in_pool, rank_in_pool, POOL_WORLD);
#endif
    pwtest.initgrids(lat0, latvec, nx, ny, nz);
    pwtest.initparameters(false, wfcecut, 1, true);
    pwtest.setuptransform();

    std::complex<double>* plane_in = new std::complex<double>[pwtest.nrxx];
    std::complex<double>* plane_ref = new std::complex<double>[pwtest.nrxx];
    std::complex<double>* plane_out = new std::complex<double>[pwtest.nrxx];
    const int sticks_work_size = comm_roundtrip_work_size(pwtest);
    std::complex<double>* sticks = new std::complex<double>[sticks_work_size];

    fill_plane_major_sticks(pwtest, plane_in);
    for (int ir = 0; ir < pwtest.nrxx; ++ir)
    {
        plane_ref[ir] = plane_in[ir];
    }
    zero_complex_buffer<PW_Basis_Comm_Accessor>(plane_out, pwtest.nrxx);
    zero_complex_buffer<PW_Basis_Comm_Accessor>(sticks, sticks_work_size);

    pwtest.gatherp_scatters(plane_in, sticks);
    expect_stick_major_equal(pwtest, sticks);
    pwtest.gathers_scatterp(sticks, plane_out);

    expect_plane_major_equal(pwtest, plane_ref, plane_out);

    delete[] plane_in;
    delete[] plane_ref;
    delete[] plane_out;
    delete[] sticks;
}
} // namespace

TEST_F(PWTEST, test_comm_roundtrip_pw_basis)
{
    run_comm_roundtrip_case(10, 10, 10);
}

TEST_F(PWTEST, test_comm_roundtrip_pw_basis_zero_plane_pressure)
{
    const int candidate_cases[][3] = {
        {10, 10, 2},
        {16, 16, 2},
        {20, 20, 2},
        {24, 24, 2},
        {20, 20, 3},
        {24, 24, 3},
        {32, 16, 2},
        {32, 32, 2},
    };

    for (const auto& candidate_case : candidate_cases)
    {
        const int nx = candidate_case[0];
        const int ny = candidate_case[1];
        const int nz = candidate_case[2];
        if (case_has_zero_plane_stress(nx, ny, nz))
        {
            run_comm_roundtrip_case(nx, ny, nz);
            return;
        }
    }

    GTEST_SKIP() << "No zero-plane/stick stress layout found for the current MPI decomposition.";
}
