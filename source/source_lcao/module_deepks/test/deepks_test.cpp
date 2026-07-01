#include "deepks_test.h"

#include <complex>
#include <gtest/gtest.h>

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
        if ((actual_word[0] - '0' >= 0 && actual_word[0] - '0' < 10) || actual_word[0] == '-')
        {
            const double num1 = std::stod(actual_word);
            const double num2 = std::stod(reference_word);
            EXPECT_NEAR(num1, num2, test_thr) << "numeric mismatch at entry " << entry;
        }
        else if (actual_word[0] == '(' && actual_word[actual_word.size() - 1] == ')' && reference_word[0] == '('
                 && reference_word[reference_word.size() - 1] == ')')
        {
            const std::string actual_str = actual_word.substr(1, actual_word.size() - 2);
            const std::string reference_str = reference_word.substr(1, reference_word.size() - 2);
            const double actual_real = std::stod(actual_str.substr(0, actual_str.find(',')));
            const double actual_imag = std::stod(actual_str.substr(actual_str.find(',') + 1));
            const double reference_real = std::stod(reference_str.substr(0, reference_str.find(',')));
            const double reference_imag = std::stod(reference_str.substr(reference_str.find(',') + 1));
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
