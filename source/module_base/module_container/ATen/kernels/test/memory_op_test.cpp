#include <gtest/gtest.h>

#include <ATen/core/tensor.h>
#include <ATen/core/tensor_map.h>
#include <ATen/kernels/memory_op.h>
#include <ATen/kernels/test/op_test_utils.h>

namespace container {
namespace op {

template <typename T>
class MemoryOpTest : public testing::Test {
public:
    MemoryOpTest() {
        test_utils::init_blas_handle();
        test_utils::init_cusolver_handle();
    }
    ~MemoryOpTest() override {
        test_utils::delete_blas_handle();
        test_utils::delete_cusolver_handle();
    }
};

TYPED_TEST_SUITE(MemoryOpTest, test_utils::Types);

TYPED_TEST(MemoryOpTest, ResizeAndSynchronizeMemory) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    op::resize_memory_op<Type, Device> resizeMemory;
    op::synchronize_memory_op<Type, DEVICE_CPU, Device> syncMemoryDeviceToHost;
    op::synchronize_memory_op<Type, Device, DEVICE_CPU> syncMemoryHostToDevice;

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

TYPED_TEST(MemoryOpTest, SetMemory) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    op::set_memory_op<Type, Device> setMemory;

    Tensor A = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0)}).to_device<Device>());
    Tensor B = A;
    
    A.zero();
    setMemory(B.data<Type>(), 0, 3);
    EXPECT_EQ(A, B);
}

TYPED_TEST(MemoryOpTest, CastAndDeleteMemory) {
    using Type = std::complex<typename GetTypeReal<typename std::tuple_element<0, decltype(TypeParam())>::type>::type>;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    op::delete_memory_op<std::complex<float>, DEVICE_CPU> deleteMemory;
    op::cast_memory_op<std::complex<float>, Type, DEVICE_CPU, Device> castMemory_H2D_D2S;

    Tensor A = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0)}).to_device<Device>());
    Tensor B = A.to_device<DEVICE_CPU>().cast<std::complex<float>>();
    
    std::complex<float> * d_A = (std::complex<float>*)malloc(sizeof(std::complex<float>) * 3);
    castMemory_H2D_D2S(d_A, A.data<Type>(), 3);
    Tensor C = std::move(TensorMap(d_A, B.data_type(), B.device_type(), {3}));

    EXPECT_EQ(B, C);
    deleteMemory(d_A);
}

} // namespace op
} // namespace container
