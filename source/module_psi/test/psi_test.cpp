#include "module_psi/psi.h"

#include <gtest/gtest.h>

class TestPsi : public ::testing::Test
{
  public:
    const int ink = 1;
    const int inbands = 4;
    const int inbasis = 10;
    int ngk[4] = {10, 10, 10, 10};

    const psi::Psi<std::complex<double>>* psi_object31 = new psi::Psi<std::complex<double>>(ink, inbands, inbasis, &ngk[0]);
    const psi::Psi<double>* psi_object32 = new psi::Psi<double>(ink, inbands, inbasis, &ngk[0]);
    const psi::Psi<std::complex<float>>* psi_object33 = new psi::Psi<std::complex<float>>(ink, inbands, inbasis, &ngk[0]);
    const psi::Psi<float>* psi_object34 = new psi::Psi<float>(ink, inbands, inbasis, &ngk[0]);

    psi::Psi<std::complex<double>>* psi_object4 = new psi::Psi<std::complex<double>>(*psi_object31, ink, 0);
    psi::Psi<std::complex<double>>* psi_object5 = new psi::Psi<std::complex<double>>(psi_object31->get_pointer(), *psi_object31, ink, 0);
};

TEST_F(TestPsi, get_val)
{
    psi::Psi<std::complex<double>>* psi_object11 = new psi::Psi<std::complex<double>>();
    psi::Psi<double>* psi_object12 = new psi::Psi<double>();
    psi::Psi<std::complex<float>>* psi_object13 = new psi::Psi<std::complex<float>>();
    psi::Psi<float>* psi_object14 = new psi::Psi<float>();

    EXPECT_EQ(psi_object11->get_nk(), 1);
    EXPECT_EQ(psi_object11->get_nbands(), 1);
    EXPECT_EQ(psi_object11->get_nbasis(), 1);
    EXPECT_EQ(psi_object11->get_current_k(), 0);
    EXPECT_EQ(psi_object11->get_current_b(), 0);
    EXPECT_EQ(psi_object11->get_current_nbas(), 1);
    EXPECT_EQ(psi_object11->get_k_first(), true);
    EXPECT_EQ(psi_object11->get_psi_bias(), 0);

    EXPECT_EQ(psi_object12->get_nk(), 1);
    EXPECT_EQ(psi_object12->get_nbands(), 1);
    EXPECT_EQ(psi_object12->get_nbasis(), 1);
    EXPECT_EQ(psi_object12->get_current_k(), 0);
    EXPECT_EQ(psi_object12->get_current_b(), 0);
    EXPECT_EQ(psi_object12->get_current_nbas(), 1);
    EXPECT_EQ(psi_object12->get_k_first(), true);
    EXPECT_EQ(psi_object12->get_psi_bias(), 0);

    EXPECT_EQ(psi_object13->get_nk(), 1);
    EXPECT_EQ(psi_object13->get_nbands(), 1);
    EXPECT_EQ(psi_object13->get_nbasis(), 1);
    EXPECT_EQ(psi_object13->get_current_k(), 0);
    EXPECT_EQ(psi_object13->get_current_b(), 0);
    EXPECT_EQ(psi_object13->get_current_nbas(), 1);
    EXPECT_EQ(psi_object13->get_k_first(), true);
    EXPECT_EQ(psi_object13->get_psi_bias(), 0);

    EXPECT_EQ(psi_object14->get_nk(), 1);
    EXPECT_EQ(psi_object14->get_nbands(), 1);
    EXPECT_EQ(psi_object14->get_nbasis(), 1);
    EXPECT_EQ(psi_object14->get_current_k(), 0);
    EXPECT_EQ(psi_object14->get_current_b(), 0);
    EXPECT_EQ(psi_object14->get_current_nbas(), 1);
    EXPECT_EQ(psi_object14->get_k_first(), true);
    EXPECT_EQ(psi_object14->get_psi_bias(), 0);
}

TEST_F(TestPsi, get_ngk)
{
    psi::Psi<std::complex<double>>* psi_object21 = new psi::Psi<std::complex<double>>(&ngk[0]);
    psi::Psi<double>* psi_object22 = new psi::Psi<double>(&ngk[0]);
    psi::Psi<std::complex<float>>* psi_object23 = new psi::Psi<std::complex<float>>(&ngk[0]);
    psi::Psi<float>* psi_object24 = new psi::Psi<float>(&ngk[0]);

    EXPECT_EQ(psi_object21->get_ngk(2), ngk[2]);
    EXPECT_EQ(psi_object21->get_ngk_pointer()[0], ngk[0]);

    EXPECT_EQ(psi_object22->get_ngk(2), ngk[2]);
    EXPECT_EQ(psi_object22->get_ngk_pointer()[0], ngk[0]);

    EXPECT_EQ(psi_object23->get_ngk(2), ngk[2]);
    EXPECT_EQ(psi_object23->get_ngk_pointer()[0], ngk[0]);

    EXPECT_EQ(psi_object24->get_ngk(2), ngk[2]);
    EXPECT_EQ(psi_object24->get_ngk_pointer()[0], ngk[0]);
}

TEST_F(TestPsi, get_pointer_op_zero_complex_double)
{
    for (int i = 0; i < ink; i++)
    {
        psi_object31->fix_k(i);
        for (int j = 0; j < inbands; j++)
        {
            for (int k = 0; k< inbasis; k++)
            {
                psi_object31->get_pointer()[j * inbasis + k].real(j * inbasis + k);
                psi_object31->get_pointer()[j * inbasis + k].imag(j * inbasis + k);
            }
        }
    }
    // fix_k
    psi_object31->fix_k(0);
    // get_pointer
    EXPECT_EQ(psi_object31->get_pointer(2)[0].real(), 20);
    EXPECT_EQ(psi_object31->get_pointer(2)[0].imag(), 20);
    // op
    EXPECT_EQ((*psi_object31)(0, 1, 2).real(), 12);
    EXPECT_EQ((*psi_object31)(0, 1, 2).real(), 12);
    EXPECT_EQ((*psi_object31)(2, 2).real(), 22);
    EXPECT_EQ((*psi_object31)(2, 2).imag(), 22);
    EXPECT_EQ((*psi_object31)(2).real(), 2);
    EXPECT_EQ((*psi_object31)(2).imag(), 2);
    // zero
    psi::Psi<std::complex<double>>* psi_object6 = new psi::Psi<std::complex<double>>(*psi_object31);
    psi_object6->zero_out();
    EXPECT_EQ((*psi_object6)(0, 1, 2).real(), 0);
    EXPECT_EQ((*psi_object6)(0, 1, 2).imag(), 0);
    delete psi_object6;

    // cover all lines in fix_k func
    psi_object31->fix_k(1);
    EXPECT_EQ(psi_object31->get_psi_bias(), 0);
    psi::Psi<std::complex<double>>* psi_temp = new psi::Psi<std::complex<double>>(ink, inbands, inbasis);
    psi_temp->fix_k(0);
    EXPECT_EQ(psi_object31->get_current_nbas(), inbasis);
    delete psi_temp;
}

TEST_F(TestPsi, get_pointer_op_zero_complex_float)
{
    for (int i = 0; i < ink; i++)
    {
        psi_object33->fix_k(i);
        for (int j = 0; j < inbands; j++)
        {
            for (int k = 0; k< inbasis; k++)
            {
                psi_object33->get_pointer()[j * inbasis + k].real(j * inbasis + k);
                psi_object33->get_pointer()[j * inbasis + k].imag(j * inbasis + k);
            }
        }
    }
    // fix_k
    psi_object33->fix_k(0);
    // get_pointer
    EXPECT_EQ(psi_object33->get_pointer(2)[0].real(), 20);
    EXPECT_EQ(psi_object33->get_pointer(2)[0].imag(), 20);
    // op
    EXPECT_EQ((*psi_object33)(0, 1, 2).real(), 12);
    EXPECT_EQ((*psi_object33)(0, 1, 2).real(), 12);
    EXPECT_EQ((*psi_object33)(2, 2).real(), 22);
    EXPECT_EQ((*psi_object33)(2, 2).imag(), 22);
    EXPECT_EQ((*psi_object33)(2).real(), 2);
    EXPECT_EQ((*psi_object33)(2).imag(), 2);
    // zero
    psi::Psi<std::complex<float>>* psi_object6 = new psi::Psi<std::complex<float>>(*psi_object33);
    psi_object6->zero_out();
    EXPECT_EQ((*psi_object6)(0, 1, 2).real(), 0);
    EXPECT_EQ((*psi_object6)(0, 1, 2).imag(), 0);
    delete psi_object6;
}

TEST_F(TestPsi, get_pointer_op_zero_double)
{
    for (int i = 0; i < ink; i++)
    {
        psi_object32->fix_k(i);
        for (int j = 0; j < inbands; j++)
        {
            for (int k = 0; k< inbasis; k++)
            {
                psi_object32->get_pointer()[j * inbasis + k] = j * inbasis + k;
            }
        }
    }
    // fix_k
    psi_object32->fix_k(0);
    // get_pointer
    EXPECT_EQ(psi_object32->get_pointer(2)[0], 20);
    // op
    EXPECT_EQ((*psi_object32)(0, 1, 2), 12);
    EXPECT_EQ((*psi_object32)(2, 2), 22);
    EXPECT_EQ((*psi_object32)(2), 2);
    // zero
    psi::Psi<double>* psi_object6 = new psi::Psi<double>(*psi_object32);
    psi_object6->zero_out();
    EXPECT_EQ((*psi_object6)(0, 1, 2), 0);
    delete psi_object6;
}

TEST_F(TestPsi, get_pointer_op_zero_float)
{
    for (int i = 0; i < ink; i++)
    {
        psi_object34->fix_k(i);
        for (int j = 0; j < inbands; j++)
        {
            for (int k = 0; k< inbasis; k++)
            {
                psi_object34->get_pointer()[j * inbasis + k] = j * inbasis + k;
            }
        }
    }
    // fix_k
    psi_object34->fix_k(0);
    // get_pointer
    EXPECT_EQ(psi_object34->get_pointer(2)[0], 20);
    // op
    EXPECT_EQ((*psi_object34)(0, 1, 2), 12);
    EXPECT_EQ((*psi_object34)(2, 2), 22);
    EXPECT_EQ((*psi_object34)(2), 2);
    // zero
    psi::Psi<float>* psi_object6 = new psi::Psi<float>(*psi_object34);
    psi_object6->zero_out();
    EXPECT_EQ((*psi_object6)(0, 1, 2), 0);
    delete psi_object6;
}


TEST_F(TestPsi, size)
{
    EXPECT_EQ(psi_object31->size(), ink * inbands * inbasis);
    EXPECT_EQ(psi_object32->size(), ink * inbands * inbasis);
    EXPECT_EQ(psi_object33->size(), ink * inbands * inbasis);
    EXPECT_EQ(psi_object34->size(), ink * inbands * inbasis);

    // cover all lines in size func
    psi::Psi<std::complex<double>>* psi_object1 = new psi::Psi<std::complex<double>>();
    EXPECT_EQ(psi_object1->size(), 0);
}

TEST_F(TestPsi, range)
{
    psi::Range range1(1);// k_first = 1;index_1 = 0;range_1 = 1;range_2 = 1;
    psi::Range range2(0,1,0,0);
    EXPECT_EQ(range1.range_1, 1);
    EXPECT_EQ(range1.range_2, 1);
    EXPECT_EQ(range2.k_first, 0);
    EXPECT_EQ(range2.index_1, 1);
    EXPECT_EQ(range2.range_1, 0);
    EXPECT_EQ(range2.range_2, 0);

    int num1 = std::get<1>(psi_object31->to_range(range1));
    EXPECT_EQ(num1, 1);
    int num2 = std::get<1>(psi_object31->to_range(range2));
    EXPECT_EQ(num2, 0);

    num1 = std::get<1>(psi_object32->to_range(range1));
    EXPECT_EQ(num1, 1);
    num2 = std::get<1>(psi_object32->to_range(range2));
    EXPECT_EQ(num2, 0);

    num1 = std::get<1>(psi_object33->to_range(range1));
    EXPECT_EQ(num1, 1);
    num2 = std::get<1>(psi_object33->to_range(range2));
    EXPECT_EQ(num2, 0);

    num1 = std::get<1>(psi_object34->to_range(range1));
    EXPECT_EQ(num1, 1);
    num2 = std::get<1>(psi_object34->to_range(range2));
    EXPECT_EQ(num2, 0);
}



#if __UT_USE_CUDA || __UT_USE_ROCM
TEST_F(TestPsi, Range)
{
    psi::Psi<std::complex<double>, psi::DEVICE_GPU>* psi_object3_gpu1 = new psi::Psi<std::complex<double>, psi::DEVICE_GPU>(*psi_object31);
    delete psi_object3_gpu1;

    psi::Psi<double, psi::DEVICE_GPU>* psi_object3_gpu2 = new psi::Psi<double, psi::DEVICE_GPU>(*psi_object32);
    delete psi_object3_gpu2;

    psi::Psi<std::complex<float>, psi::DEVICE_GPU>* psi_object3_gpu3 = new psi::Psi<std::complex<float>, psi::DEVICE_GPU>(*psi_object33);
    delete psi_object3_gpu3;

    psi::Psi<float, psi::DEVICE_GPU>* psi_object3_gpu4 = new psi::Psi<float, psi::DEVICE_GPU>(*psi_object34);
    delete psi_object3_gpu4;
}
#endif // __UT_USE_CUDA || __UT_USE_ROCM
