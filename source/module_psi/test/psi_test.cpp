#include "module_psi/psi.h"

#include <gtest/gtest.h>

class TestPsi : public ::testing::Test
{
  public:
      const int ink = 2;
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
    psi_object31->fix_k(2);
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

TEST_F(TestPsi, band_first)
{
    const psi::Psi<std::complex<double>>* psi_band_c64 = new psi::Psi<std::complex<double>>(ink, inbands, inbasis, &ngk[0], false);
    const psi::Psi<double>* psi_band_64 = new psi::Psi<double>(ink, inbands, inbasis, &ngk[0], false);
    const psi::Psi<std::complex<float>>* psi_band_c32 = new psi::Psi<std::complex<float>>(ink, inbands, inbasis, &ngk[0], false);
    const psi::Psi<float>* psi_band_32 = new psi::Psi<float>(ink, inbands, inbasis, &ngk[0], false);

    // set values: cover 4 different cases
    for (int ib = 0;ib < inbands;++ib)
    {
        psi_band_c64->fix_b(ib); // 1. fix_b, fix_k, (ibasis)
        psi_band_64->fix_b(ib);// 2. fix_kb, (ibasis)
        psi_band_c32->fix_b(ib);// 3. fix_b, (ik, ibasis)
        EXPECT_EQ(psi_band_c64->get_current_b(), ib);
        EXPECT_EQ(psi_band_64->get_current_b(), ib);
        EXPECT_EQ(psi_band_c32->get_current_b(), ib);
        EXPECT_EQ(psi_band_c64->get_psi_bias(), ib * ink * inbasis);
        EXPECT_EQ(psi_band_64->get_psi_bias(), ib * ink * inbasis);
        EXPECT_EQ(psi_band_c32->get_psi_bias(), ib * ink * inbasis);
        for (int ik = 0;ik < ink;++ik)
        {
            psi_band_c64->fix_k(ik);
            psi_band_64->fix_kb(ik, ib);
            EXPECT_EQ(psi_band_c64->get_current_k(), ik);
            EXPECT_EQ(psi_band_64->get_current_k(), ik);
            EXPECT_EQ(psi_band_c64->get_psi_bias(), (ib * ink + ik) * inbasis);
            EXPECT_EQ(psi_band_64->get_psi_bias(), (ib * ink + ik) * inbasis);
            for (int ibas = 0;ibas < inbasis;++ibas)
            {
                int index = ((ib * ink) + ik) * inbasis + ibas;
                (*psi_band_c64)(ibas).real(index);
                (*psi_band_64)(ibas) = index;
                (*psi_band_c32)(ik, ibas).real(index);
                (*psi_band_32)(ib, ik, ibas) = index; //4. no fix, (ib, ik, ibasis)
            }
        }
    }

    // get_pointer and operator() (using different fix from setter)
    psi_band_c64->fix_k(1);
    EXPECT_EQ(psi_band_c64->get_pointer()[1].real(), 71);
    EXPECT_EQ((*psi_band_c64)(2).real(), 72);
    psi_band_64->fix_b(1);
    EXPECT_EQ(psi_band_64->get_pointer()[1], 21);
    EXPECT_EQ(psi_band_64->get_pointer(1)[1], 31);
    EXPECT_EQ((*psi_band_64)(0, 2), 22);
    EXPECT_EQ((*psi_band_64)(1, 2), 32);
    psi_band_c32->fix_b(0);
    EXPECT_EQ(psi_band_c32->get_pointer()[1].real(), 1);
    EXPECT_EQ((*psi_band_c32)(1, 1, 2).real(), 32);
    psi_band_32->fix_kb(0, 2);
    EXPECT_EQ(psi_band_32->get_pointer()[1], 41);
    EXPECT_EQ(psi_band_32->get_pointer(1)[1], 51);
    EXPECT_EQ((*psi_band_32)(2), 42);

    // range
    psi::Range b2_k11(0, 2, 1, 1);
    psi::Range b13(0, -1, 1, 3);
    psi::Range illegal_kfirst(1, 2, 1, 1);
    psi::Range illegal_index1(0, 4, 2, 1);
    psi::Range illegal_range1(0, -1, 3, 1);
    psi::Range illegal_range2(0, 2, 1, 3);
    EXPECT_EQ(std::get<0>(psi_band_c64->to_range(b2_k11))[1].real(), 51);
    EXPECT_EQ(std::get<1>(psi_band_c64->to_range(b2_k11)), 1);
    EXPECT_EQ(std::get<0>(psi_band_64->to_range(b13))[50], 70);
    EXPECT_EQ(std::get<1>(psi_band_64->to_range(b13)), 3);
    EXPECT_EQ(std::get<0>(psi_band_c32->to_range(illegal_kfirst)), nullptr);
    EXPECT_EQ(std::get<1>(psi_band_c32->to_range(illegal_index1)), 0);
    EXPECT_EQ(std::get<0>(psi_band_32->to_range(illegal_range1)), nullptr);
    EXPECT_EQ(std::get<1>(psi_band_32->to_range(illegal_range2)), 0);

    // pointer constructor
    // band-first to k-first
    psi::Psi<float> psi_band_32_k(psi_band_32->get_pointer(), psi_band_32->get_nk(), psi_band_32->get_nbands(), psi_band_32->get_nbasis(), psi_band_32->get_ngk_pointer(), true);
    // k-first to band-first
    psi::Psi<float> psi_band_32_b(psi_band_32_k.get_pointer(), psi_band_32_k.get_nk(), psi_band_32_k.get_nbands(), psi_band_32_k.get_nbasis(), psi_band_32_k.get_ngk_pointer(), false);
    EXPECT_EQ(psi_band_32_k.get_nk(), ink);
    EXPECT_EQ(psi_band_32_k.get_nbands(), inbands);
    EXPECT_EQ(psi_band_32_k.get_nbasis(), inbasis);
    EXPECT_EQ(psi_band_32_b.get_nk(), ink);
    EXPECT_EQ(psi_band_32_b.get_nbands(), inbands);
    EXPECT_EQ(psi_band_32_b.get_nbasis(), inbasis);
    for (int ik = 0;ik < ink;++ik)
        for (int ib = 0;ib < inbands;++ib)
        {
            psi_band_32->fix_kb(ik, ib);
            psi_band_32_k.fix_kb(ik, ib);
            psi_band_32_b.fix_kb(ik, ib);
            EXPECT_EQ(psi_band_32->get_psi_bias(), (ib * ink + ik) * inbasis);
            EXPECT_EQ(psi_band_32_k.get_psi_bias(), (ik * inbands + ib) * inbasis);
            EXPECT_EQ(psi_band_32_b.get_psi_bias(), (ib * ink + ik) * inbasis);
        }

    delete psi_band_c64;
    delete psi_band_64;
    delete psi_band_c32;
    delete psi_band_32;
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
