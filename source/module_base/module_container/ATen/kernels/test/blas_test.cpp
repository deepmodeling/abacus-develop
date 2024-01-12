#include <gtest/gtest.h>

#include <ATen/core/tensor.h>
#include <ATen/kernels/blas.h>
#include <base/utils/gtest.h>

namespace container {
namespace kernels {

template <typename T>
class BlasTest : public testing::Test {
public:
    BlasTest() {
        base::utils::init_blas_handle();
    }
    ~BlasTest() override {
        base::utils::delete_blas_handle();
    }
};

TYPED_TEST_SUITE(BlasTest, base::utils::Types);

TYPED_TEST(BlasTest, Dot) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    blas_dot<Type, Device> dotCalculator;

    const int n = 3;
    const Tensor x = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0)}).to_device<Device>());
    const Tensor y = std::move(Tensor({static_cast<Type>(4.0), static_cast<Type>(5.0), static_cast<Type>(6.0)}).to_device<Device>());
    
    Type result = {};
    dotCalculator(n, x.data<Type>(), 1, y.data<Type>(), 1, &result);
    const Type expected = static_cast<Type>(32.0);

    EXPECT_EQ(result, expected);
}

TYPED_TEST(BlasTest, Scal) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    blas_scal<Type, Device> scalCalculator;

    const int n = 3;
    const Type alpha = static_cast<Type>(2.0);
    Tensor x = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0)}).to_device<Device>());
    
    scalCalculator(n, &alpha, x.data<Type>(), 1);
    const Tensor expected = std::move(Tensor({static_cast<Type>(2.0), static_cast<Type>(4.0), static_cast<Type>(6.0)}).to_device<Device>());

    EXPECT_EQ(x, expected);
}


TYPED_TEST(BlasTest, Axpy) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    blas_axpy<Type, Device> axpyCalculator;

    const int n = 3;
    const Type alpha = static_cast<Type>(2.0);
    const Tensor x = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0)}).to_device<Device>());
    Tensor       y = std::move(Tensor({static_cast<Type>(4.0), static_cast<Type>(5.0), static_cast<Type>(6.0)}).to_device<Device>());
    
    axpyCalculator(n, &alpha, x.data<Type>(), 1, y.data<Type>(), 1);
    const Tensor expected = std::move(Tensor({static_cast<Type>(6.0), static_cast<Type>(9.0), static_cast<Type>(12.0)}).to_device<Device>());

    EXPECT_EQ(y, expected);
}


TYPED_TEST(BlasTest, Gemv) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    blas_gemv<Type, Device> gemvCalculator;

    const char trans = 'N';
    const int m = 3;
    const int n = 2;
    const Type alpha = static_cast<Type>(2.0);
    const Type beta  = static_cast<Type>(3.0);
    const Tensor A = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0), 
                                       static_cast<Type>(4.0), static_cast<Type>(5.0), static_cast<Type>(6.0)}).to_device<Device>());
    const Tensor x = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0)}).to_device<Device>());
    Tensor       y = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0)}).to_device<Device>());
    
    gemvCalculator(trans, m, n, &alpha, A.data<Type>(), m, x.data<Type>(), 1, &beta, y.data<Type>(), 1);
    const Tensor expected = std::move(Tensor({static_cast<Type>(21.0), static_cast<Type>(30.0), static_cast<Type>(39.0)}).to_device<Device>());

    EXPECT_EQ(y, expected);
}


TYPED_TEST(BlasTest, GemvBatched) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    blas_gemv<Type, Device> gemvCalculator;
    blas_gemv_batched<Type, Device> gemvBatchedCalculator;

    const char trans = 'N';
    const int m = 3;
    const int n = 2;
    const int batch_size = 2;
    const Type alpha = static_cast<Type>(2.0);
    const Type beta  = static_cast<Type>(3.0);

    std::vector<Type*> A = {};
    std::vector<Type*> x = {};
    std::vector<Type*> y = {};

    const Tensor _A = std::move(Tensor({
        static_cast<Type>(1.0), static_cast<Type>(2.0), 
        static_cast<Type>(3.0), static_cast<Type>(4.0), 
        static_cast<Type>(5.0), static_cast<Type>(6.0),
        
        static_cast<Type>(7.0), static_cast<Type>(8.0),
        static_cast<Type>(9.0), static_cast<Type>(10.0),
        static_cast<Type>(11.0),static_cast<Type>(12.0)}).to_device<Device>());
    
    A.push_back(_A.data<Type>());
    A.push_back(_A.data<Type>() + m * n);

    const Tensor _x = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0)}).to_device<Device>());
    x.push_back(_x.data<Type>());
    x.push_back(_x.data<Type>());

    Tensor _y1 = std::move(Tensor({static_cast<Type>(4.0), static_cast<Type>(5.0), static_cast<Type>(6.0),
                                   static_cast<Type>(7.0), static_cast<Type>(8.0), static_cast<Type>(9.0)}).to_device<Device>());
    Tensor _y2 = _y1;
    y.push_back(_y1.data<Type>());
    y.push_back(_y1.data<Type>() + m);

    gemvBatchedCalculator(trans, m, n, &alpha, A.data(), m, x.data(), 1, &beta, y.data(), 1, batch_size);

    for (int ii = 0; ii < batch_size; ++ii) {
        gemvCalculator(trans, m, n, &alpha, A[ii], m, x[ii], 1, &beta, _y2.data<Type>() + ii * m, 1);
    }

    EXPECT_EQ(_y1, _y2);
}


TYPED_TEST(BlasTest, GemvBatchedStrided) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    blas_gemv<Type, Device> gemvCalculator;
    blas_gemv_batched_strided<Type, Device> gemvBatchedStridedCalculator;

    const char trans = 'N';
    const int m = 3;
    const int n = 2;
    const int batch_size = 2;
    const Type alpha = static_cast<Type>(2.0);
    const Type beta  = static_cast<Type>(3.0);

    std::vector<Type*> A = {};
    std::vector<Type*> x = {};
    std::vector<Type*> y = {};

    const Tensor _A = std::move(Tensor({
        static_cast<Type>(1.0), static_cast<Type>(2.0), 
        static_cast<Type>(3.0), static_cast<Type>(4.0), 
        static_cast<Type>(5.0), static_cast<Type>(6.0),
        
        static_cast<Type>(7.0), static_cast<Type>(8.0),
        static_cast<Type>(9.0), static_cast<Type>(10.0),
        static_cast<Type>(11.0),static_cast<Type>(12.0)}).to_device<Device>());
    
    A.push_back(_A.data<Type>());
    A.push_back(_A.data<Type>() + m * n);

    const Tensor _x = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0)}).to_device<Device>());
    x.push_back(_x.data<Type>());
    x.push_back(_x.data<Type>());

    Tensor _y1 = std::move(Tensor({static_cast<Type>(4.0), static_cast<Type>(5.0), static_cast<Type>(6.0),
                                   static_cast<Type>(7.0), static_cast<Type>(8.0), static_cast<Type>(9.0)}).to_device<Device>());
    Tensor _y2 = _y1;
    y.push_back(_y1.data<Type>());
    y.push_back(_y1.data<Type>() + m);

    gemvBatchedStridedCalculator(trans, m, n, &alpha, A[0], m, m * n, x[0], 1, 0, &beta, y[0], 1, m, batch_size);

    for (int ii = 0; ii < batch_size; ++ii) {
        gemvCalculator(trans, m, n, &alpha, A[ii], m, x[ii], 1, &beta, _y2.data<Type>() + ii * m, 1);
    }
    EXPECT_EQ(_y1, _y2);
}


TYPED_TEST(BlasTest, Gemm) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    blas_gemm<Type, Device> gemmCalculator;

    const char trans = 'N';
    const int m = 3;
    const int n = 2;
    const Type alpha = static_cast<Type>(2.0);
    const Type beta  = static_cast<Type>(3.0);
    const Tensor A = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0), 
                                       static_cast<Type>(4.0), static_cast<Type>(5.0), static_cast<Type>(6.0)}).to_device<Device>());
    const Tensor x = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0)}).to_device<Device>());
    Tensor       y = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0), static_cast<Type>(3.0)}).to_device<Device>());
    
    gemmCalculator(trans, trans, m, 1, n, &alpha, A.data<Type>(), m, x.data<Type>(), n, &beta, y.data<Type>(), m);
    const Tensor expected = std::move(Tensor({static_cast<Type>(21.0), static_cast<Type>(30.0), static_cast<Type>(39.0)}).to_device<Device>());

    EXPECT_EQ(y, expected);
}


TYPED_TEST(BlasTest, GemmBatched) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    blas_gemv_batched<Type, Device> gemvBatchedCalculator;
    blas_gemm_batched<Type, Device> gemmBatchedCalculator;

    const char trans = 'N';
    const int m = 3;
    const int n = 2;
    const int batch_size = 2;
    const Type alpha = static_cast<Type>(2.0);
    const Type beta  = static_cast<Type>(3.0);

    std::vector<Type*> A = {};
    std::vector<Type*> x = {};
    std::vector<Type*> y1 = {};
    std::vector<Type*> y2 = {};

    const Tensor _A = std::move(Tensor({
        static_cast<Type>(1.0), static_cast<Type>(2.0), 
        static_cast<Type>(3.0), static_cast<Type>(4.0), 
        static_cast<Type>(5.0), static_cast<Type>(6.0),
        
        static_cast<Type>(7.0), static_cast<Type>(8.0),
        static_cast<Type>(9.0), static_cast<Type>(10.0),
        static_cast<Type>(11.0),static_cast<Type>(12.0)}).to_device<Device>());
    
    A.push_back(_A.data<Type>());
    A.push_back(_A.data<Type>() + m * n);

    const Tensor _x = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0)}).to_device<Device>());
    x.push_back(_x.data<Type>());
    x.push_back(_x.data<Type>());

    Tensor _y1 = std::move(Tensor({static_cast<Type>(4.0), static_cast<Type>(5.0), static_cast<Type>(6.0),
                                   static_cast<Type>(7.0), static_cast<Type>(8.0), static_cast<Type>(9.0)}).to_device<Device>());
    Tensor _y2 = _y1;
    y1.push_back(_y1.data<Type>());
    y1.push_back(_y1.data<Type>() + m);
    y2.push_back(_y2.data<Type>());
    y2.push_back(_y2.data<Type>() + m);

    gemvBatchedCalculator(trans, m, n, &alpha, A.data(), m, x.data(), 1, &beta, y1.data(), 1, batch_size);
    gemmBatchedCalculator(trans, trans, m, 1, n, &alpha, A.data(), m, x.data(), n, &beta, y2.data(), m, batch_size);

    EXPECT_EQ(_y1, _y2);
}


TYPED_TEST(BlasTest, GemmBatchedStrided) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;

    blas_gemv_batched_strided<Type, Device> gemvBatchedStridedCalculator;
    blas_gemm_batched_strided<Type, Device> gemmBatchedStridedCalculator;

    const char trans = 'N';
    const int m = 3;
    const int n = 2;
    const int batch_size = 2;
    const Type alpha = static_cast<Type>(2.0);
    const Type beta  = static_cast<Type>(3.0);

    std::vector<Type*> A = {};
    std::vector<Type*> x = {};
    std::vector<Type*> y1 = {};
    std::vector<Type*> y2 = {};

    const Tensor _A = std::move(Tensor({
        static_cast<Type>(1.0), static_cast<Type>(2.0), 
        static_cast<Type>(3.0), static_cast<Type>(4.0), 
        static_cast<Type>(5.0), static_cast<Type>(6.0),
        
        static_cast<Type>(7.0), static_cast<Type>(8.0),
        static_cast<Type>(9.0), static_cast<Type>(10.0),
        static_cast<Type>(11.0),static_cast<Type>(12.0)}).to_device<Device>());
    
    A.push_back(_A.data<Type>());
    A.push_back(_A.data<Type>() + m * n);

    const Tensor _x = std::move(Tensor({static_cast<Type>(1.0), static_cast<Type>(2.0)}).to_device<Device>());
    x.push_back(_x.data<Type>());
    x.push_back(_x.data<Type>());

    Tensor _y1 = std::move(Tensor({static_cast<Type>(4.0), static_cast<Type>(5.0), static_cast<Type>(6.0),
                                   static_cast<Type>(7.0), static_cast<Type>(8.0), static_cast<Type>(9.0)}).to_device<Device>());
    Tensor _y2 = _y1;
    y1.push_back(_y1.data<Type>());
    y1.push_back(_y1.data<Type>() + m);
    y2.push_back(_y2.data<Type>());
    y2.push_back(_y2.data<Type>() + m);

    gemvBatchedStridedCalculator(trans, m, n, &alpha, A[0], m, m * n, x[0], 1, 0, &beta, y1[0], 1, m, batch_size);
    gemmBatchedStridedCalculator(trans, trans, m, 1, n, &alpha, A[0], m, m * n, x[0], n, 0, &beta, y2[0], m, m, batch_size);

    EXPECT_EQ(_y1, _y2);
}

} // namespace op
} // namespace container

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
