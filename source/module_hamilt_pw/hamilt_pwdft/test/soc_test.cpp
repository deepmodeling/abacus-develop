#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <string>
#include <cmath>
#include <complex>

/************************************************
 *  unit test of class Soc and Fcoef
 ***********************************************/

/**
 * - Tested Functions:
 *   - Fcoef::create to create a 5 dimensional array of complex numbers
 *   - Soc::set_fcoef to set the fcoef array
 *   - Soc::spinor to calculate the spinor
 *   - Soc::rot_ylm to calculate the rotation matrix
 *   - Soc::sph_ind to calculate the m index of the spherical harmonics
*/

//compare two complex by using EXPECT_DOUBLE_EQ()
void EXPECT_COMPLEX_DOUBLE_EQ(const std::complex<double>& a,const std::complex<double>& b)
{
    EXPECT_DOUBLE_EQ(a.real(),b.real());
    EXPECT_DOUBLE_EQ(a.imag(),b.imag());
}

#define private public
#include "module_hamilt_pw/hamilt_pwdft/soc.h"
class FcoefTest : public testing::Test
{
protected:
                  Fcoef fcoef;
                  std::string output;
};

TEST_F(FcoefTest, Create)
{
                  fcoef.create(2, 3, 4);
                  EXPECT_EQ(fcoef.ind1, 2);
                  EXPECT_EQ(fcoef.ind4, 3);
                  EXPECT_EQ(fcoef.ind5, 4);
                  EXPECT_NE(fcoef.p, nullptr);
                  EXPECT_DOUBLE_EQ(fcoef(0,0,0,0,0).real(),0.0);
                  EXPECT_DOUBLE_EQ(fcoef(0,0,0,0,0).imag(),0.0);
                  const Fcoef &fcoef_const = fcoef;
                  EXPECT_DOUBLE_EQ(fcoef_const(0,0,0,0,0).real(),0.0);
                  EXPECT_DOUBLE_EQ(fcoef_const(0,0,0,0,0).imag(),0.0);
}

TEST_F(FcoefTest, CreateWarning)
{
                  fcoef.create(0, 2, 2);
                  EXPECT_EQ(fcoef.ind1, 1);
                  EXPECT_EQ(fcoef.ind4, 1);
                  EXPECT_EQ(fcoef.ind5, 1);
                  EXPECT_EQ(fcoef.p, nullptr);
}

class SocTest : public testing::Test
{
protected:
                  Soc soc;
                  Fcoef fcoef;
                  std::string output;
};

TEST_F(SocTest, Spinor)
{
                  EXPECT_DOUBLE_EQ(soc.spinor(0, 0.5, 0, 0), 1.0);
                  EXPECT_DOUBLE_EQ(soc.spinor(0, 0.5, 0, 1), 0.0);
                  EXPECT_NEAR(soc.spinor(1, 0.5, 0, 0), std::sqrt(2./3.),1e-15);
                  EXPECT_NEAR(soc.spinor(1, 0.5, 0, 1), -std::sqrt(1./3.), 1e-15);
                  EXPECT_NEAR(soc.spinor(1, 0.5, 1, 0), std::sqrt(1./3.), 1e-15);
                  EXPECT_NEAR(soc.spinor(1, 0.5, 1, 1), -std::sqrt(2./3.),1e-15);
                  EXPECT_DOUBLE_EQ(soc.spinor(1, 0.5, -1, 1), 0.0);
                  //warning message
                  testing::internal::CaptureStdout();
                  EXPECT_EXIT(soc.spinor(2, 0.5, 0, 0), ::testing::ExitedWithCode(0),"");
                  output = testing::internal::GetCapturedStdout();
                  EXPECT_THAT(output, testing::HasSubstr("j and l not compatible"));
                  //warning message
                  testing::internal::CaptureStdout();
                  EXPECT_EXIT(soc.spinor(2, 0.5, 0, 2), ::testing::ExitedWithCode(0),"");
                  output = testing::internal::GetCapturedStdout();
                  EXPECT_THAT(output, testing::HasSubstr("spin direction unknown"));
                  //warning message
                  testing::internal::CaptureStdout();
                  EXPECT_EXIT(soc.spinor(2, 0.5, 3, 0), ::testing::ExitedWithCode(0),"");
                  output = testing::internal::GetCapturedStdout();
                  EXPECT_THAT(output, testing::HasSubstr("m not allowed"));
}

TEST_F(SocTest, SphInd)
{
                  EXPECT_EQ(soc.sph_ind(0, 0.5, 0, 0), 0);
                  EXPECT_EQ(soc.sph_ind(0, 0.5, 0, 1), 0);
                  EXPECT_EQ(soc.sph_ind(1, 0.5, 0, 0), -1);
                  EXPECT_EQ(soc.sph_ind(1, 0.5, 0, 1), 0);
                  EXPECT_EQ(soc.sph_ind(1, 0.5, -1, 0), 0);
                  //warning message
                  testing::internal::CaptureStdout();
                  EXPECT_EXIT(soc.sph_ind(2, 0.5, 0, 0), ::testing::ExitedWithCode(0),"");
                  output = testing::internal::GetCapturedStdout();
                  EXPECT_THAT(output, testing::HasSubstr("l and j not suitable"));
                  //warning message
                  testing::internal::CaptureStdout();
                  EXPECT_EXIT(soc.sph_ind(2, 0.5, 0, 2), ::testing::ExitedWithCode(0),"");
                  output = testing::internal::GetCapturedStdout();
                  EXPECT_THAT(output, testing::HasSubstr("spin must be 0 1"));
                  //warning message
                  testing::internal::CaptureStdout();
                  EXPECT_EXIT(soc.sph_ind(2, 0.5, 3, 0), ::testing::ExitedWithCode(0),"");
                  output = testing::internal::GetCapturedStdout();
                  EXPECT_THAT(output, testing::HasSubstr("m not allowed"));
}

TEST_F(SocTest, RotYlm)
{
                  soc.rot_ylm(0);
                  EXPECT_NE(soc.p_rot, nullptr);
                  EXPECT_COMPLEX_DOUBLE_EQ(soc.p_rot[0],std::complex<double>(1.0, 0.0));
                  soc.rot_ylm(2);
                  int l = 2;
                  int l2p1 = 2*l + 1;
                  for (int i = 1; i < 2*l + 1 ; i += 2)
                  {
                                    int m = (i+1)/2;
                                    int n = l-m;
                                    EXPECT_COMPLEX_DOUBLE_EQ(soc.p_rot[l2p1*i + n], std::complex<double>(pow(-1.0,m)/sqrt(2), 0.0));
                                    EXPECT_COMPLEX_DOUBLE_EQ(soc.p_rot[l2p1*(i+1) + n], std::complex<double>(0.0,-pow(-1.0,m)/sqrt(2)));
                                    n = l+m;
                                    EXPECT_COMPLEX_DOUBLE_EQ(soc.p_rot[l2p1*i + n], std::complex<double>(1.0/sqrt(2),0.0));
                                    EXPECT_COMPLEX_DOUBLE_EQ(soc.p_rot[l2p1*(i+1) + n], std::complex<double>(0.0,1.0/sqrt(2)));
                  }
                  const std::complex<double> &rot = soc.rotylm(2, 1);
                  EXPECT_DOUBLE_EQ(rot.real(),0.0);
                  EXPECT_DOUBLE_EQ(rot.imag(),1.0/sqrt(2.0));
}

TEST_F(SocTest, SetFcoef)
{
                  soc.fcoef.create(1,1,1);
                  soc.rot_ylm(0);
                  soc.set_fcoef(0, 0, 0, 0, 0, 0, 0.5, 0.5, 0, 0, 0);
                  EXPECT_COMPLEX_DOUBLE_EQ(soc.fcoef(0, 0, 0, 0, 0), std::complex<double>(1.0, 0.0));
}
#undef private
