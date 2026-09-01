#include <gtest/gtest.h>

#include <cstdint>

#include <ATen/core/tensor.h>
#include <ATen/core/tensor_map.h>
#include <ATen/kernels/memory.h>
#include <base/utils/gtest.h>

namespace container {
namespace kernels {

template <typename T>
class MemoryTest : public testing::Test {
public:
    MemoryTest() = default;
    ~MemoryTest() override = default;
};

TYPED_TEST_SUITE(MemoryTest, base::utils::Types);

TYPED_TEST(MemoryTest, ResizeAndSynchronizeMemory) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    kernels::resize_memory<Type, Device> resizeMemory;
    kernels::synchronize_memory<Type, DEVICE_CPU, Device> syncMemoryDeviceToHost;
    kernels::synchronize_memory<Type, Device, DEVICE_CPU> syncMemoryHostToDevice;

    Tensor A = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0)}).to_device<Device>());

    Type* d_B = nullptr;
    resizeMemory(d_B, 3, "B");
    Tensor B = std::move(TensorMap(d_B, A.data_type(), A.device_type(), {3}));
    B.zero();

    syncMemoryDeviceToHost(B.data<Type>(), A.data<Type>(), 3);
    EXPECT_EQ(A, B);
    
    A.zero();
    syncMemoryHostToDevice(A.data<Type>(), B.data<Type>(), 3);
    EXPECT_EQ(A, B);
}

TYPED_TEST(MemoryTest, SetMemory) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    kernels::set_memory<Type, Device> setMemory;

    Tensor A = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0)}).to_device<Device>());
    Tensor B = A;
    
    A.zero();
    setMemory(B.data<Type>(), 0, 3);
    EXPECT_EQ(A, B);
}

TYPED_TEST(MemoryTest, CastAndDeleteMemory) {
    using Type = std::complex<typename GetTypeReal<typename std::tuple_element<0, decltype(TypeParam())>::type>::type>;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    kernels::delete_memory<std::complex<float>, DEVICE_CPU> deleteMemory;
    kernels::cast_memory<std::complex<float>, Type, DEVICE_CPU, Device> castMemory_H2D_D2S;

    Tensor A = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0)}).to_device<Device>());
    Tensor B = A.to_device<DEVICE_CPU>().cast<std::complex<float>>();
    
    auto * d_A = (std::complex<float>*)malloc(sizeof(std::complex<float>) * 3);
    castMemory_H2D_D2S(d_A, A.data<Type>(), 3);
    Tensor C = std::move(TensorMap(d_A, B.data_type(), B.device_type(), {3}));

    EXPECT_EQ(B, C);
    deleteMemory(d_A);
}

template <typename T>
void test_cpu_memory_operations(const T& value)
{
    T* data = nullptr;
    kernels::resize_memory<T, DEVICE_CPU>()(data, 6, "typed cpu buffer");
    ASSERT_NE(data, nullptr);

    kernels::set_memory<T, DEVICE_CPU>()(data, value, 6);
    for (int i = 0; i < 6; ++i)
    {
        EXPECT_EQ(data[i], value);
    }

    const T source[4] = {T(1), T(2), T(3), T(4)};
    kernels::synchronize_memory<T, DEVICE_CPU, DEVICE_CPU>()(data, source, 4);
    for (int i = 0; i < 4; ++i)
    {
        EXPECT_EQ(data[i], source[i]);
    }

    const std::vector<int64_t> input_shape = {2, 2};
    const std::vector<int64_t> output_shape = {2, 3};
    kernels::set_memory<T, DEVICE_CPU>()(data, T(), 6);
    kernels::synchronize_memory_stride<T, DEVICE_CPU, DEVICE_CPU>()(data, source, output_shape, input_shape);
    EXPECT_EQ(data[0], source[0]);
    EXPECT_EQ(data[1], source[1]);
    EXPECT_EQ(data[2], T());
    EXPECT_EQ(data[3], source[2]);
    EXPECT_EQ(data[4], source[3]);
    EXPECT_EQ(data[5], T());

    kernels::delete_memory<T, DEVICE_CPU>()(data);
}

TEST(MemoryTestCPU, IntegerAndComplexInstantiations)
{
    test_cpu_memory_operations<int>(7);
    test_cpu_memory_operations<int64_t>(9);
    test_cpu_memory_operations<float>(1.5F);
    test_cpu_memory_operations<double>(2.5);
    test_cpu_memory_operations<std::complex<float>>(std::complex<float>(1.0F, -2.0F));
    test_cpu_memory_operations<std::complex<double>>(std::complex<double>(2.0, -3.0));
}

template <typename Output, typename Input>
void test_cpu_cast(const Input& input, const Output& expected)
{
    const Input source[1] = {input};
    Output result[1] = {Output()};
    kernels::cast_memory<Output, Input, DEVICE_CPU, DEVICE_CPU>()(result, source, 1);
    EXPECT_EQ(result[0], expected);
}

TEST(MemoryTestCPU, CoversAllCastInstantiations)
{
    test_cpu_cast<float, float>(1.25F, 1.25F);
    test_cpu_cast<double, double>(2.5, 2.5);
    test_cpu_cast<float, double>(3.5, 3.5F);
    test_cpu_cast<double, float>(4.5F, 4.5);
    test_cpu_cast<std::complex<float>, std::complex<float>>({1.0F, 2.0F}, {1.0F, 2.0F});
    test_cpu_cast<std::complex<double>, std::complex<double>>({2.0, 3.0}, {2.0, 3.0});
    test_cpu_cast<std::complex<float>, std::complex<double>>({3.0, 4.0}, {3.0F, 4.0F});
    test_cpu_cast<std::complex<double>, std::complex<float>>({4.0F, 5.0F}, {4.0, 5.0});
}

#if !(defined(__CUDA) || defined(__ROCM))
template <typename T>
void test_gpu_placeholder_operations()
{
    T* data = nullptr;
    const T* source = nullptr;
    kernels::resize_memory<T, DEVICE_GPU>()(data, 0, "gpu placeholder");
    EXPECT_EQ(data, nullptr);
    kernels::set_memory<T, DEVICE_GPU>()(data, T(), 0);
    kernels::synchronize_memory<T, DEVICE_GPU, DEVICE_GPU>()(data, source, 0);
    kernels::synchronize_memory<T, DEVICE_GPU, DEVICE_CPU>()(data, source, 0);
    kernels::synchronize_memory<T, DEVICE_CPU, DEVICE_GPU>()(data, source, 0);
    kernels::delete_memory<T, DEVICE_GPU>()(data);
}

template <typename Output, typename Input>
void test_gpu_placeholder_casts()
{
    Output* output = nullptr;
    const Input* input = nullptr;
    kernels::cast_memory<Output, Input, DEVICE_GPU, DEVICE_GPU>()(output, input, 0);
    kernels::cast_memory<Output, Input, DEVICE_GPU, DEVICE_CPU>()(output, input, 0);
    kernels::cast_memory<Output, Input, DEVICE_CPU, DEVICE_GPU>()(output, input, 0);
}

TEST(MemoryTestGPUPlaceholder, CoversNoAcceleratorInstantiations)
{
    test_gpu_placeholder_operations<int>();
    test_gpu_placeholder_operations<int64_t>();
    test_gpu_placeholder_operations<float>();
    test_gpu_placeholder_operations<double>();
    test_gpu_placeholder_operations<std::complex<float>>();
    test_gpu_placeholder_operations<std::complex<double>>();

    test_gpu_placeholder_casts<float, float>();
    test_gpu_placeholder_casts<double, double>();
    test_gpu_placeholder_casts<float, double>();
    test_gpu_placeholder_casts<double, float>();
    test_gpu_placeholder_casts<std::complex<float>, std::complex<float>>();
    test_gpu_placeholder_casts<std::complex<double>, std::complex<double>>();
    test_gpu_placeholder_casts<std::complex<float>, std::complex<double>>();
    test_gpu_placeholder_casts<std::complex<double>, std::complex<float>>();
}
#endif

} // namespace op
} // namespace container
