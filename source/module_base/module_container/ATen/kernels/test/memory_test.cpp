#include <gtest/gtest.h>

#include <ATen/core/tensor.h>
#include <ATen/core/tensor_map.h>
#include <ATen/kernels/memory.h>
#include <test/test_utils.h>

namespace container {
namespace kernels {

template <typename T>
class MemoryTest : public testing::Test {
public:
    MemoryTest() = default;
    ~MemoryTest() override = default;
};

TYPED_TEST_SUITE(MemoryTest, test_utils::Types);

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

} // namespace op
} // namespace container
