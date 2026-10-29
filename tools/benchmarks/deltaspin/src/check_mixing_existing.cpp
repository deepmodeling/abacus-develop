#include "module_base/module_mixing/test/mixing_test.cpp"

template <class T>
void check_scratch(int length)
{
    Base_Mixing::Broyden_Mixing mixer(2, 0.2);
    Base_Mixing::Mixing_Data data;
    mixer.init_mixing_data(data, length, sizeof(T));
    std::vector<T> input(length), output(length), actual(length);
    for (int i = 0; i < length; ++i)
    {
        input[i] = T(0.01 * (i % 7));
        output[i] = T(0.03 * (i % 11));
    }
    auto screen = [length](T* residual)
    {
        for (int i = 0; i < length; ++i) { residual[i] *= 0.5; }
    };
    auto partial_mix = [length](T* mixed, const T* in, const T* residual)
    {
        for (int i = 0; i < length; i += 3) { mixed[i] = in[i] + 0.2 * residual[i]; }
    };
    mixer.push_data(data, input.data(), output.data(), screen, partial_mix, true);
    mixer.cal_coef(data, [length](T* a, T* b)
    {
        double value = 0;
        for (int i = 0; i < length; ++i) { value += std::real(std::conj(a[i]) * b[i]); }
        return value;
    });
    mixer.mix_data(data, actual.data());
    double error = 0;
    for (int i = 0; i < length; ++i)
    {
        const T expected = i % 3 == 0 ? input[i] + 0.2 * (0.5 * (output[i] - input[i])) : T{};
        error = std::max(error, std::abs(actual[i] - expected));
    }
    EXPECT_EQ(error, 0.0);
}

TEST(MixingScratch, PreservesScreenedResidualAndZeroInitializedOutput)
{
    omp_set_num_threads(12);
    for (int length : {7, 1048593})
    {
        check_scratch<double>(length);
        check_scratch<std::complex<double>>(length);
    }
}

int main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
