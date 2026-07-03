#ifndef DEEPKS_TEST_RUNNER_H_
#define DEEPKS_TEST_RUNNER_H_

#include "deepks_test.h"

#include <gtest/gtest.h>
#include <torch/torch.h>
#include <vector>

namespace DeepksTestRunner
{
template <typename T>
void build_phialpha(test_deepks<T>& test)
{
    test.check_phialpha();
}

template <typename T>
void build_pdm(test_deepks<T>& test)
{
    build_phialpha(test);
    if (testing::Test::HasFatalFailure())
    {
        return;
    }
    test.check_pdm();
}

template <typename T>
void build_descriptor(test_deepks<T>& test, std::vector<torch::Tensor>& descriptor)
{
    build_pdm(test);
    if (testing::Test::HasFatalFailure())
    {
        return;
    }
    test.check_descriptor(descriptor);
}

template <typename T>
void build_edelta(test_deepks<T>& test, std::vector<torch::Tensor>& descriptor)
{
    build_descriptor(test, descriptor);
    if (testing::Test::HasFatalFailure())
    {
        return;
    }
    test.check_edelta(descriptor);
}
} // namespace DeepksTestRunner

#endif
