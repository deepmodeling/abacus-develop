#ifdef __MPI
#include "source_base/parallel_device.h"

#include "source_base/parallel_cell.h"
#include "source_base/parallel_comm.h"

#include "gtest/gtest.h"
#include <complex>
#include <vector>

namespace
{

template <typename T>
void exercise_cpu_point()
{
    const int count = 2;
    T values[count] = {static_cast<T>(1), static_cast<T>(2)};
    T temporary[count] = {};
    Parallel_Common::object_cpu_point<T, base_device::DEVICE_CPU> point;

    EXPECT_EQ(point.get_buffer(values, count), values);
    EXPECT_EQ(point.get(values, count), values);
    EXPECT_EQ(point.get_buffer(values, count, temporary), values);
    EXPECT_EQ(point.get(values, count, temporary), values);
    point.sync_h2d(values, temporary, count);
    point.sync_d2h(temporary, values, count);
    point.del(values);
}

template <typename T>
void exercise_gpu_staging_stubs()
{
    const int count = 2;
    T values[count] = {static_cast<T>(1), static_cast<T>(2)};
    T temporary[count] = {};
    Parallel_Common::object_cpu_point<T, base_device::DEVICE_GPU> point;

    EXPECT_EQ(point.get_buffer(values, count, temporary), temporary);
    EXPECT_EQ(point.get(values, count, temporary), temporary);
    point.sync_h2d(values, temporary, count);
    point.sync_d2h(temporary, values, count);
    point.del(temporary);

    T* allocated = point.get_buffer(values, count);
    EXPECT_NE(allocated, nullptr);
    point.del(allocated);
}

template <typename T>
void exercise_mpi_wrappers(const ModuleBase::CommunicationDomain& domain)
{
    MPI_Comm communicator = domain.communicator();
    const int rank = domain.rank();
    MPICommGroup world_group(communicator);
    const int size = world_group.gsize;
    const int count = 2;

    T sent[count] = {static_cast<T>(rank + 1), static_cast<T>(rank + 2)};
    T received[count] = {};
    MPI_Status status;
    MPI_Request request;
    Parallel_Common::send_dev<T, base_device::DEVICE_CPU>(sent, count, MPI_PROC_NULL, 0, communicator);
    Parallel_Common::isend_dev<T, base_device::DEVICE_CPU>(sent,
                                                           count,
                                                           MPI_PROC_NULL,
                                                           0,
                                                           communicator,
                                                           &request,
                                                           nullptr);
    Parallel_Common::recv_dev<T, base_device::DEVICE_CPU>(received, count, MPI_PROC_NULL, 0, communicator, &status);

    T broadcast[count] = {};
    if (rank == 0)
    {
        broadcast[0] = static_cast<T>(3);
        broadcast[1] = static_cast<T>(5);
    }
    Parallel_Common::bcast_dev<T, base_device::DEVICE_CPU>(broadcast, count, communicator, 0);
    EXPECT_EQ(broadcast[0], static_cast<T>(3));
    EXPECT_EQ(broadcast[1], static_cast<T>(5));

    T reduced[count] = {static_cast<T>(rank + 1), static_cast<T>(2 * (rank + 1))};
    Parallel_Common::reduce_dev<T, base_device::DEVICE_CPU>(reduced, count, communicator);
    const T sum = static_cast<T>(size * (size + 1) / 2);
    EXPECT_EQ(reduced[0], sum);
    EXPECT_EQ(reduced[1], static_cast<T>(2) * sum);

    const T gathered_value = static_cast<T>(rank + 1);
    std::vector<T> gathered(size);
    std::vector<int> receive_counts(size, 1);
    std::vector<int> displacements(size);
    for (int index = 0; index < size; ++index)
    {
        displacements[index] = index;
    }
    Parallel_Common::gatherv_dev<T, base_device::DEVICE_CPU>(&gathered_value,
                                                             1,
                                                             gathered.data(),
                                                             receive_counts.data(),
                                                             displacements.data(),
                                                             communicator);
    for (int index = 0; index < size; ++index)
    {
        EXPECT_EQ(gathered[index], static_cast<T>(index + 1));
    }
}

TEST(ParallelDevice, CoversCpuPointSpecializations)
{
    exercise_cpu_point<float>();
    exercise_cpu_point<double>();
    exercise_cpu_point<std::complex<float>>();
    exercise_cpu_point<std::complex<double>>();
}

TEST(ParallelDevice, CoversGpuStagingWithoutAccelerator)
{
    exercise_gpu_staging_stubs<float>();
    exercise_gpu_staging_stubs<double>();
    exercise_gpu_staging_stubs<std::complex<float>>();
    exercise_gpu_staging_stubs<std::complex<double>>();
}

TEST(ParallelDevice, CoversMpiTypeOverloads)
{
    const ModuleBase::CommunicationDomain domain = ModuleBase::world_communication_domain();
    exercise_mpi_wrappers<float>(domain);
    exercise_mpi_wrappers<double>(domain);
    exercise_mpi_wrappers<std::complex<float>>(domain);
    exercise_mpi_wrappers<std::complex<double>>(domain);
}

} // namespace
#endif
