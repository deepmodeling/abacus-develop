#include <vector>
#include <gtest/gtest.h>

#include <ATen/core/tensor.h>
#include <ATen/core/tensor_map.h>


TEST(Tensor, Constructor) {
    // Test constructor with default allocator
    container::Tensor t1(container::DataType::DT_FLOAT, container::TensorShape({2, 3}));
    EXPECT_EQ(t1.data_type(), container::DataType::DT_FLOAT);
    EXPECT_EQ(t1.device_type(), container::DeviceType::CpuDevice);
    // EXPECT_EQ(t1.shape().dims(), std::vector<int64_t>({2, 3}));
    EXPECT_EQ(t1.NumElements(), 6);

#if __CUDA || __ROCM
    // Test constructor with specified device type
    container::Tensor t2(container::DataType::DT_DOUBLE, container::DeviceType::GpuDevice, container::TensorShape({3, 4}));
    EXPECT_EQ(t2.data_type(), container::DataType::DT_DOUBLE);
    EXPECT_EQ(t2.device_type(), container::DeviceType::GpuDevice);
    // EXPECT_EQ(t2.shape().dims(), std::vector<int64_t>({3, 4}));
    EXPECT_EQ(t2.NumElements(), 12);
#endif

    // Test copy constructor
    container::Tensor t3(t1);
    EXPECT_EQ(t3.data_type(), container::DataType::DT_FLOAT);
    EXPECT_EQ(t3.device_type(), container::DeviceType::CpuDevice);
    // EXPECT_EQ(t3.shape().dims(), std::vector<int64_t>({2, 3}));
    EXPECT_EQ(t3.NumElements(), 6);
    EXPECT_NE(t3.data(), t1.data());

    // Test reference constructor
    std::vector<float> vec{1.0, 2.0, 3.0};
    container::TensorMap t4(&vec[0], container::DataType::DT_FLOAT, container::DeviceType::CpuDevice, container::TensorShape({1, 3}));
    EXPECT_EQ(t4.data_type(), container::DataType::DT_FLOAT);
    EXPECT_EQ(t4.device_type(), container::DeviceType::CpuDevice);
    // EXPECT_EQ(t4.shape().dims(), std::vector<int64_t>({1, 3}));
    EXPECT_EQ(t4.NumElements(), 3);
    EXPECT_EQ(t4.data(), vec.data());
}


TEST(Tensor, GetDataPointer) {
    // Create a 1x1 float tensor with data [1.0, 2.0, 3.0, 4.0].
    container::Tensor t1(container::DataType::DT_INT, container::TensorShape({1, 1}));
    container::Tensor t2(container::DataType::DT_INT64, container::TensorShape({1, 1}));
    container::Tensor t3(container::DataType::DT_FLOAT, container::TensorShape({1, 1}));
    container::Tensor t4(container::DataType::DT_DOUBLE, container::TensorShape({1, 1}));
    container::Tensor t5(container::DataType::DT_COMPLEX, container::TensorShape({1, 1}));
    container::Tensor t6(container::DataType::DT_COMPLEX_DOUBLE, container::TensorShape({1, 1}));
    t1.data<int>()[0] = 1;
    t2.data<int64_t>()[0] = 1;
    t3.data<float>()[0] = 1.0f;
    t4.data<double>()[0] = 1.0f;
    t5.data<std::complex<float>>()[0] = {1.0f, 0.0f};
    t6.data<std::complex<double>>()[0] = {1.0f, 0.0f};
    // Get a pointer to the data buffer.
    void* ptr1 = t1.data();
    void* ptr2 = t2.data();
    void* ptr3 = t3.data();
    void* ptr4 = t4.data();
    void* ptr5 = t5.data();
    void* ptr6 = t6.data();
    // Ensure that the returned pointer is not null and points to the expected data.
    EXPECT_NE(ptr1, nullptr);
    EXPECT_NE(ptr2, nullptr);
    EXPECT_NE(ptr3, nullptr);
    EXPECT_NE(ptr4, nullptr);
    EXPECT_NE(ptr5, nullptr);
    EXPECT_NE(ptr6, nullptr);
    EXPECT_EQ(static_cast<int*>(ptr1)[0], 1);
    EXPECT_EQ(static_cast<int64_t*>(ptr2)[0], 1);
    EXPECT_EQ(static_cast<float*>(ptr3)[0], 1.0f);
    EXPECT_EQ(static_cast<double*>(ptr4)[0], 1.0f);

    EXPECT_EQ(static_cast<std::complex<float>*>(ptr5)[0].real(), 1.0);
    EXPECT_EQ(static_cast<std::complex<float>*>(ptr5)[0].imag(), 0.0);
    EXPECT_EQ(static_cast<std::complex<double>*>(ptr6)[0].real(), 1.0);
    EXPECT_EQ(static_cast<std::complex<double>*>(ptr6)[0].imag(), 0.0);
}


TEST(Tensor, GetDataPointerDeathTest) {
    // Try to get a typed pointer with a type that does not match the tensor's data type.
    // This should cause an error message to be printed and the program to exit with failure.
    container::Tensor tensor(container::DataType::DT_FLOAT, container::TensorShape({1, 1}));
    // Verify that requesting data with an unsupported data type causes the program to exit.
    ASSERT_EXIT(
      tensor.data<int>(), // Unsupported data type
      ::testing::ExitedWithCode(EXIT_FAILURE),
      "Tensor data type does not match requested type."
    );
}

TEST(Tensor, SizeOfType) {
    // Test DT_FLOAT
    EXPECT_EQ(container::Tensor::SizeOfType(container::DataType::DT_FLOAT), sizeof(float));

    // Test DT_INT
    EXPECT_EQ(container::Tensor::SizeOfType(container::DataType::DT_INT), sizeof(int32_t));

    // Test DT_INT64
    EXPECT_EQ(container::Tensor::SizeOfType(container::DataType::DT_INT64), sizeof(int64_t));

    // Test DT_DOUBLE
    EXPECT_EQ(container::Tensor::SizeOfType(container::DataType::DT_DOUBLE), sizeof(double));

    // Test DT_COMPLEX
    EXPECT_EQ(container::Tensor::SizeOfType(container::DataType::DT_COMPLEX), sizeof(std::complex<float>));

    // Test DT_COMPLEX_DOUBLE
    EXPECT_EQ(container::Tensor::SizeOfType(container::DataType::DT_COMPLEX_DOUBLE), sizeof(std::complex<double>));

}

TEST(Tensor, SizeOfTypeDeathTest) {
    // Verify that requesting data with an unsupported data type causes the program to exit.
    ASSERT_EXIT(
      container::Tensor::SizeOfType(container::DataType::DT_INVALID),
      ::testing::ExitedWithCode(EXIT_FAILURE),
      "Unsupported data type!"
    );
}

TEST(Tensor, ToDeviceAndSetZero) {
    // Create tensor on CPU
    container::Tensor tensor(container::DataType::DT_FLOAT, container::DeviceType::CpuDevice, {2, 3});

    // Set zero test
    tensor.zero();

    // Create tensor on GPU
    container::Tensor cpu_tensor = tensor.to_device<container::DEVICE_CPU>();

    // Check device type
    EXPECT_EQ(cpu_tensor.device_type(), container::DeviceType::CpuDevice);

    // Check data type
    EXPECT_EQ(cpu_tensor.data_type(), container::DataType::DT_FLOAT);

    // Check shape
    EXPECT_EQ(cpu_tensor.shape(), tensor.shape());

    // Check data
    for (int ii = 0; ii < cpu_tensor.NumElements(); ii++) {
        EXPECT_EQ(cpu_tensor.data<float>()[ii], 0.0);
    }
}

TEST(Tensor, Cast) {
    // Create a tensor object with float data type and device CPU
    container::Tensor t(container::DataType::DT_COMPLEX_DOUBLE, container::DeviceType::CpuDevice, {2, 3});
    t.data<std::complex<double>>()[0] = {1.0, 0.0};
    t.data<std::complex<double>>()[1] = {2.0, 0.0};
    t.data<std::complex<double>>()[2] = {3.0, 0.0};
    t.data<std::complex<double>>()[3] = {4.0, 0.0};
    t.data<std::complex<double>>()[4] = {5.0, 0.0};
    t.data<std::complex<double>>()[5] = {6.0, 0.0};

    // Cast the tensor to integer data type
    container::Tensor t_float = t.cast<std::complex<float>>();

    // Check that the data type and device of the output tensor are correct
    EXPECT_EQ(t_float.data_type(), container::DataType::DT_COMPLEX);
    EXPECT_EQ(t_float.device_type(), container::DeviceType::CpuDevice);

    // Check that the shape of the output tensor is correct
    EXPECT_EQ(t_float.shape().dims(), std::vector<int64_t>({2, 3}));

    // Check that the data of the output tensor is correct
    EXPECT_EQ(t_float.data<std::complex<float>>()[0].real(), 1.0);
    EXPECT_EQ(t_float.data<std::complex<float>>()[0].imag(), 0.0);
    EXPECT_EQ(t_float.data<std::complex<float>>()[1].real(), 2.0);
    EXPECT_EQ(t_float.data<std::complex<float>>()[1].imag(), 0.0);
    EXPECT_EQ(t_float.data<std::complex<float>>()[2].real(), 3.0);
    EXPECT_EQ(t_float.data<std::complex<float>>()[2].imag(), 0.0);
    EXPECT_EQ(t_float.data<std::complex<float>>()[3].real(), 4.0);
    EXPECT_EQ(t_float.data<std::complex<float>>()[3].imag(), 0.0);
    EXPECT_EQ(t_float.data<std::complex<float>>()[4].real(), 5.0);
    EXPECT_EQ(t_float.data<std::complex<float>>()[4].imag(), 0.0);
    EXPECT_EQ(t_float.data<std::complex<float>>()[5].real(), 6.0);
    EXPECT_EQ(t_float.data<std::complex<float>>()[5].imag(), 0.0);
}

// Tests the reshape() function of the Tensor class.
TEST(Tensor, Reshape) {
    container::Tensor t(container::DataType::DT_FLOAT, container::DeviceType::CpuDevice, {2, 3, 4});
    container::TensorShape new_shape({-1, 8});
    ASSERT_NO_THROW(t.reshape(new_shape));
    EXPECT_EQ(t.shape().ndim(), 2);
    EXPECT_EQ(t.shape().dim_size(0), 3);
    EXPECT_EQ(t.shape().dim_size(1), 8);

    container::Tensor t1(container::DataType::DT_FLOAT, container::DeviceType::CpuDevice, {2, 3, 4});
    container::TensorShape new_shape1({2, 3, 4});
    ASSERT_NO_THROW(t1.reshape(new_shape1));
    EXPECT_EQ(t1.shape(), new_shape1);
}

TEST(Tensor, GetValueAndInnerMostPtr) {
    container::Tensor t(container::DataType::DT_INT, container::DeviceType::CpuDevice, {2, 2, 4});
    std::vector<int> vec = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
    EXPECT_EQ(t.shape().NumElements(), vec.size());
    memcpy(t.data<int>(), vec.data(), sizeof(int) * vec.size());
    EXPECT_EQ(t.get_value<int>(0, 0, 1), 2);
    EXPECT_EQ(t.get_value<int>(1, 1, 2), 15);
    t.reshape({4, 4});

    // check the inner_most_ptr meshod
    auto row_ptr = t.inner_most_ptr<int>(2);
    EXPECT_EQ(row_ptr[0], 9);
    EXPECT_EQ(row_ptr[3], 12);
}

TEST(Tensor, ReshapeDeathTest) {
    container::Tensor t(container::DataType::DT_FLOAT, container::DeviceType::CpuDevice, {2, 3, 4});
    container::TensorShape new_shape({-1, 8});
    new_shape.set_dim_size(1, -2);
    EXPECT_THROW(t.reshape(new_shape), std::invalid_argument);

    container::Tensor t1(container::DataType::DT_FLOAT, container::DeviceType::CpuDevice, {2, 3, 4});
    container::TensorShape new_shape1({-1, 5});
    EXPECT_THROW(t1.reshape(new_shape1), std::invalid_argument);

    container::Tensor t2(container::DataType::DT_FLOAT, container::DeviceType::CpuDevice, {2, 3, 4});
    container::TensorShape new_shape2({-1, -1});
    EXPECT_THROW(t2.reshape(new_shape2), std::invalid_argument);

    container::Tensor t3(container::DataType::DT_FLOAT, container::DeviceType::CpuDevice, {2, 3, 4});
    container::TensorShape new_shape3({2, 7});
    EXPECT_THROW(t3.reshape(new_shape3), std::invalid_argument);
}

// Tests the slice() function of the Tensor class.
TEST(Tensor, Slice) {
    container::Tensor t(container::DataType::DT_FLOAT, container::DeviceType::CpuDevice, {2, 3});
    // fill with test data
    for (int i = 0; i < t.NumElements(); ++i) {
        t.data<float>()[i] = i;
    }
    // test the slice() function
    container::Tensor output = t.slice({0, 0}, {2, 2});
    EXPECT_EQ(output.shape().ndim(), 2);
    EXPECT_EQ(output.shape().dim_size(0), 2);
    EXPECT_EQ(output.shape().dim_size(1), 2);
    EXPECT_EQ(output.data<float>()[0], 0.0f);
    EXPECT_EQ(output.data<float>()[1], 1.0f);
    EXPECT_EQ(output.data<float>()[2], 3.0f);
    EXPECT_EQ(output.data<float>()[3], 4.0f);

    // test error handling
    container::Tensor t2(container::DataType::DT_FLOAT, container::DeviceType::CpuDevice, {2, 3});
    EXPECT_THROW(t2.slice({-1, 0}, {2, 2}), std::invalid_argument);
    EXPECT_THROW(t2.slice({0, 0}, {2, 4}), std::invalid_argument);
    EXPECT_THROW(t2.slice({0, 0, 0}, {2, 4, 3}), std::invalid_argument);
    EXPECT_THROW(t2.slice({0, 0, 0, 0}, {2, 4, 3, 6}), std::invalid_argument);

    container::Tensor t3(container::DataType::DT_FLOAT, container::DeviceType::CpuDevice, {3});
    // fill with test data
    for (int i = 0; i < t3.NumElements(); ++i) {
        t3.data<float>()[i] = i;
    }
    // test the slice() function
    container::Tensor output3 = t3.slice({0}, {1});
    EXPECT_EQ(output3.shape().ndim(), 1);
    EXPECT_EQ(output3.shape().dim_size(0), 1);
    EXPECT_EQ(output3.data<float>()[0], 0.0f);

    container::Tensor t4(container::DataType::DT_FLOAT, container::DeviceType::CpuDevice, {3});
    // fill with test data
    for (int i = 0; i < t4.NumElements(); ++i) {
        t4.data<float>()[i] = i;
    }
    t4.reshape({1, 1, 3});
    // test the slice() function
    container::Tensor output4 = t4.slice({0, 0, 0}, {1, 1, 2});
    EXPECT_EQ(output4.shape().ndim(), 3);
    EXPECT_EQ(output4.shape().dim_size(2), 2);
    EXPECT_EQ(output4.data<float>()[0], 0.0f);
}

TEST(Tensor, Buffer) {
    // create a tensor of shape (2, 3)
    container::TensorShape shape({2, 3});
    container::Tensor tensor(container::DataType::DT_FLOAT, container::DeviceType::CpuDevice, shape);

    // fill the tensor with some values
    auto* data_ptr = tensor.data<float>();
    for (int i = 0; i < tensor.NumElements(); i++) {
        data_ptr[i] = static_cast<float>(i);
    }

    // get the tensor buffer
    const container::TensorBuffer& buffer = tensor.buffer();

    // check if the data pointer is the same as the tensor data pointer
    assert(buffer.data() == static_cast<void*>(data_ptr));
}

TEST(Tensor, Resize) {
    container::Tensor t1(container::DataType::DT_FLOAT, container::TensorShape({2, 2}));
    const float* data_ptr1 = t1.data<float>();

    container::TensorShape new_shape({3, 3});
    t1.resize(new_shape);

    // Check if the data type remains the same after resize
    EXPECT_EQ(t1.data_type(), container::DataType::DT_FLOAT);

    // Check if the shape of the tensor object is updated
    EXPECT_EQ(t1.shape(), new_shape);

    // Check if the data buffer of the tensor object is reallocated
    EXPECT_NE(t1.data<float>(), data_ptr1);

    // Check if the data buffer is correctly zeroed
    const float* data_ptr2 = t1.data<float>();
    for (int i = 0; i < new_shape.NumElements(); ++i) {
        EXPECT_FLOAT_EQ(data_ptr2[i], 0.0);
    }

    container::Tensor t2(container::DataType::DT_FLOAT, container::DeviceType::CpuDevice, {2, 3, 4});
    container::TensorShape new_shape2({2, 3, 4});
    ASSERT_NO_THROW(t2.resize(new_shape2));
}

TEST(Tensor, GetAllocatorDeathTest) {
    container::Tensor t1(container::DataType::DT_FLOAT, container::TensorShape({2, 2}));
    ASSERT_EXIT(
      container::Allocator* alloc = container::Tensor::GetAllocator(container::DeviceType::UnKnown),
      ::testing::ExitedWithCode(EXIT_FAILURE),
      "Tensor device type unknown does not match requested type."
    );
}

TEST(Tensor, OutputOperator) {
    // Create a tensor of shape [2, 2] with random values
    const int64_t num_elements = 4;
    int* data1 = new int[num_elements];
    int64_t* data2 = new int64_t[num_elements];
    float* data3 = new float[num_elements];
    double* data4 = new double[num_elements];
    std::complex<float>* data5 = new std::complex<float>[num_elements];
    std::complex<double>* data6 = new std::complex<double>[num_elements];
    for (int ii = 0; ii < num_elements; ++ii) {
        data1[ii] = static_cast<int>(ii);
        data2[ii] = static_cast<int64_t>(ii);
        data3[ii] = static_cast<float>(ii);
        data4[ii] = static_cast<double>(ii);
        data5[ii] = std::complex<float>{static_cast<float>(ii), static_cast<float>(ii)};
        data6[ii] = std::complex<double>{static_cast<double>(ii), static_cast<double>(ii)};
    }
    const container::TensorShape shape({2, 2});
    const container::TensorMap t1(data1, container::DataType::DT_INT, container::DeviceType::CpuDevice, shape);
    const container::TensorMap t2(data2, container::DataType::DT_INT64, container::DeviceType::CpuDevice, shape);
    const container::TensorMap t3(data3, container::DataType::DT_FLOAT, container::DeviceType::CpuDevice, shape);
    const container::TensorMap t4(data4, container::DataType::DT_DOUBLE, container::DeviceType::CpuDevice, shape);
    const container::TensorMap t5(data5, container::DataType::DT_COMPLEX, container::DeviceType::CpuDevice, shape);
    const container::TensorMap t6(data6, container::DataType::DT_COMPLEX_DOUBLE, container::DeviceType::CpuDevice, shape);
    // Test if the output operator produces the expected output
    std::ostringstream oss;
    oss << t1 << t2 << t3 << t4 << t5 << t6;
    const std::string expected_output = "Tensor(shape=[2,2], data_type=int32, device_type=cpu, owns_memory=0, buffer=";
    EXPECT_TRUE(oss.str().find(expected_output) == 0);

    delete[] data1;
}
