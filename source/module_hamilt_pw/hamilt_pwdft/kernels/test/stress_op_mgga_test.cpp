#include <module_hamilt_pw/hamilt_pwdft/kernels/stress_op.h>

#include <base/utils/gtest.h>
#include <ATen/core/tensor.h>

namespace hamilt {

template<typename T>
class StressMggaTest : public testing::Test {
public:
    StressMggaTest() = default;

    ~StressMggaTest() override = default;
};

TYPED_TEST_SUITE(StressMggaTest, base::utils::ComplexTypes);

TYPED_TEST(StressMggaTest, cal_stress_mgga_op) {
    using Type = typename std::tuple_element<0, decltype(TypeParam())>::type;
    using Device = typename std::tuple_element<1, decltype(TypeParam())>::type;
    using Real = typename GetTypeReal<Type>::type;
    
    const int nrxx = 4;
    const int nspin = 0;
    const Real w1 = 0.00031685874122631746;

    ct::Tensor gradwfc = std::move(ct::Tensor(
               {static_cast<Type>( 0.0000001909712), static_cast<Type>(-0.0000001687233), static_cast<Type>(-0.0000002392604),
                static_cast<Type>( 0.0000000154653), static_cast<Type>( 0.0000004211339), static_cast<Type>( 0.0000007242030),
                static_cast<Type>( 0.0000007559914), static_cast<Type>( 0.0000005232479), static_cast<Type>( 0.0000001903915),
                static_cast<Type>(-0.0000000207916), static_cast<Type>( 0.0000000366767), static_cast<Type>( 0.0000003353377)}).to_device<Device>());

    ct::Tensor crosstaus = std::move(ct::Tensor(
    ct::DataTypeToEnum<Real>::value, ct::DeviceTypeToEnum<Device>::value, {nrxx * 6}));
    crosstaus.zero();
    ct::Tensor expected_crosstaus = crosstaus;
    expected_crosstaus.zero();

    auto cal_stress_mgga_solver = cal_stress_mgga_op<Type, typename ct::ContainerToPsi<Device>::type>();
    cal_stress_mgga_solver(nspin, nrxx, w1, gradwfc.data<Type>(), crosstaus.data<Real>());

    EXPECT_EQ(crosstaus, expected_crosstaus);
}

} // namespace hamilt