#include "deepks_test.h"

#include <cerrno>
#include <complex>
#include <cstdlib>
#include <gtest/gtest.h>
#include <string>

namespace
{
bool parse_double_token(const std::string& token, double* value)
{
    char* end = nullptr;
    errno = 0;
    const double parsed = std::strtod(token.c_str(), &end);
    if (end == token.c_str() || *end != '\0' || errno == ERANGE)
    {
        return false;
    }
    *value = parsed;
    return true;
}

bool parse_complex_token(const std::string& token, double* real, double* imag)
{
    if (token.size() < 5 || token[0] != '(' || token[token.size() - 1] != ')')
    {
        return false;
    }

    const std::string value = token.substr(1, token.size() - 2);
    const std::string::size_type comma = value.find(',');
    if (comma == std::string::npos || comma == 0 || comma == value.size() - 1)
    {
        return false;
    }

    return parse_double_token(value.substr(0, comma), real) && parse_double_token(value.substr(comma + 1), imag);
}
} // namespace

namespace Test_Deepks
{
Grid_Driver GridD(false, false);
}

template <typename T>
test_deepks<T>::test_deepks()
{
}

template <typename T>
test_deepks<T>::~test_deepks()
{
    delete this->p_elec_DM;
}

template <typename T>
void test_deepks<T>::check_dstable()
{
    // OGT.talpha.print_Table_DSR(ORB);
    // this->assert_file_matches_reference("S_I_mu_alpha.dat", "S_I_mu_alpha_ref.dat");
}

template <typename T>
void test_deepks<T>::assert_file_matches_reference(const std::string& actual_file, const std::string& reference_file)
{
    SCOPED_TRACE("Comparing " + actual_file + " with " + reference_file);
    std::ifstream actual(actual_file.c_str());
    std::ifstream reference(reference_file.c_str());
    const double test_thr = 1e-8;

    ASSERT_TRUE(actual.is_open()) << "Cannot open actual file " << actual_file;
    ASSERT_TRUE(reference.is_open()) << "Cannot open reference file " << reference_file;

    std::string actual_word;
    std::string reference_word;
    int entry = 0;
    while (actual >> actual_word)
    {
        ASSERT_TRUE(reference >> reference_word)
            << reference_file << " has fewer entries than " << actual_file << " at entry " << entry;

        double actual_num = 0.0;
        double reference_num = 0.0;
        const bool actual_is_number = parse_double_token(actual_word, &actual_num);
        const bool reference_is_number = parse_double_token(reference_word, &reference_num);
        double actual_real = 0.0;
        double actual_imag = 0.0;
        double reference_real = 0.0;
        double reference_imag = 0.0;
        const bool actual_is_complex = parse_complex_token(actual_word, &actual_real, &actual_imag);
        const bool reference_is_complex = parse_complex_token(reference_word, &reference_real, &reference_imag);

        if (actual_is_number || reference_is_number)
        {
            ASSERT_TRUE(actual_is_number) << "Cannot parse actual numeric entry " << entry << ": " << actual_word;
            ASSERT_TRUE(reference_is_number)
                << "Cannot parse reference numeric entry " << entry << ": " << reference_word;
            EXPECT_NEAR(actual_num, reference_num, test_thr) << "numeric mismatch at entry " << entry;
        }
        else if (actual_is_complex || reference_is_complex)
        {
            ASSERT_TRUE(actual_is_complex) << "Cannot parse actual complex entry " << entry << ": " << actual_word;
            ASSERT_TRUE(reference_is_complex)
                << "Cannot parse reference complex entry " << entry << ": " << reference_word;
            EXPECT_NEAR(actual_real, reference_real, test_thr) << "complex real mismatch at entry " << entry;
            EXPECT_NEAR(actual_imag, reference_imag, test_thr) << "complex imag mismatch at entry " << entry;
        }
        else
        {
            EXPECT_EQ(actual_word, reference_word) << "text mismatch at entry " << entry;
        }
        ++entry;
    }
    EXPECT_FALSE(reference >> reference_word)
        << reference_file << " has more entries than " << actual_file << " starting with " << reference_word;
}

template class test_deepks<double>;
template class test_deepks<std::complex<double>>;
