#include <vector>
#include <gtest/gtest.h>

#include <ATen/core/tensor_map.h>

TEST(TensorMap, Constructor) {
    // Test reference constructor
    std::vector<float> vec{1.0, 2.0, 3.0};
    container::TensorMap t1(&vec[0], container::DataType::DT_FLOAT, container::DeviceType::CpuDevice, container::TensorShape({1, 3}));
    EXPECT_EQ(t1.data_type(), container::DataType::DT_FLOAT);
    EXPECT_EQ(t1.device_type(), container::DeviceType::CpuDevice);
    EXPECT_EQ(t1.shape().dims(), std::vector<int64_t>({1, 3}));
    EXPECT_EQ(t1.NumElements(), 3);
    EXPECT_EQ(t1.data(), vec.data());

    container::TensorMap t2(&vec[0], t1, container::TensorShape({1, 3}));
    EXPECT_EQ(t2.data_type(), container::DataType::DT_FLOAT);
    EXPECT_EQ(t2.device_type(), container::DeviceType::CpuDevice);
    EXPECT_EQ(t2.shape().dims(), std::vector<int64_t>({1, 3}));
    EXPECT_EQ(t2.NumElements(), 3);
    EXPECT_EQ(t2.data(), vec.data());

    container::TensorMap t3(&vec[0], t2);
    EXPECT_EQ(t3.data_type(), container::DataType::DT_FLOAT);
    EXPECT_EQ(t3.device_type(), container::DeviceType::CpuDevice);
    EXPECT_EQ(t3.shape().dims(), std::vector<int64_t>({1, 3}));
    EXPECT_EQ(t3.NumElements(), 3);
    EXPECT_EQ(t3.data(), vec.data());
}