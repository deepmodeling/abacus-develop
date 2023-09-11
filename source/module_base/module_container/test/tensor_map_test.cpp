#include <vector>
#include <gtest/gtest.h>

#include <ATen/core/tensor_map.h>

TEST(TensorMap, Constructor) {
    // Test reference constructor
    std::vector<float> vec{1.0, 2.0, 3.0};
    container::TensorMap t4(&vec[0], container::DataType::DT_FLOAT, container::DeviceType::CpuDevice, container::TensorShape({1, 3}));
    EXPECT_EQ(t4.data_type(), container::DataType::DT_FLOAT);
    EXPECT_EQ(t4.device_type(), container::DeviceType::CpuDevice);
    EXPECT_EQ(t4.shape().dims(), std::vector<int64_t>({1, 3}));
    EXPECT_EQ(t4.NumElements(), 3);
    EXPECT_EQ(t4.data(), vec.data());
}

TEST(TensorMap, Resize) {
    // Test reference constructor
    std::vector<float> vec{1.0, 2.0, 3.0};
    container::TensorMap t4(&vec[0], container::DataType::DT_FLOAT, container::DeviceType::CpuDevice, container::TensorShape({1, 3}));
    EXPECT_EQ(t4.data_type(), container::DataType::DT_FLOAT);
    EXPECT_EQ(t4.device_type(), container::DeviceType::CpuDevice);
    EXPECT_EQ(t4.shape().dims(), std::vector<int64_t>({1, 3}));
    EXPECT_EQ(t4.NumElements(), 3);
    EXPECT_EQ(t4.data(), vec.data());

    EXPECT_THROW(t4.resize({2, 2});, std::logic_error);
}