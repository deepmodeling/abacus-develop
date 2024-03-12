#include"../complexmatrix.h"
#include"../matrix.h"
#include"gtest/gtest.h"
#include "gmock/gmock.h"

/************************************************
*  unit test of class ComplexMatrix and related functions
***********************************************/

/**
 * Tested functions of class ComplexMatrix
 *  - constructor
 *      - ComplexMatrix()
 *      - ComplexMatrix(const int nrows,const int ncols,const bool flag_zero=true)
 *      - ComplexMatrix(const ComplexMatrix &m1)
 *      - ComplexMatrix(ComplexMatrix && m1)
 *      - explicit ComplexMatrix(const matrix &m)
 *  - operator "=" : a ComplexMatrix and the rvalue of a ComplexMatrix
 *  - operator "()": access the element. the normal function and the the "const" function
 *  - operator "*=" "+=" "-="
 *  - real()
 *  - zero_out();
 *  - set_as_identity_matrix()
 *  - print():Output the elements of this complex matrix greater than threshold.
 *  - checkreal()
 *
 * Tested relative functions
 *  - operator "+" "-" "*" between two ComplexMatrix
 *  - operator "*" between a ComplexMatrix and double or complex, and reverse.
 *  - trace()
 *  - abs2_row()
 *  - abs2_column()
 *  - abs2():
 *  - transpose()
 *  - conj()
 *  - scale_accumulate():
 *  - scaled_sum():
 *
 */

//a mock function of WARNING_QUIT, to avoid the uncorrected call by matrix.cpp at line 37.
namespace ModuleBase
{
    void WARNING_QUIT(const std::string &file,const std::string &description) {exit(1);}
}

inline void EXPECT_COMPLEX_EQ(const std::complex<double>& a,const std::complex<double>& b)
{
    EXPECT_DOUBLE_EQ(a.real(),b.real());
    EXPECT_DOUBLE_EQ(a.imag(),b.imag());
}

class ComplexMatrixTest : public testing::Test
{
    protected:
    ModuleBase::ComplexMatrix cm22, cm23, cm32, cm33;
    std::complex<double> comzero {0.0,0.0};

    void SetUp()
    {
        cm22.create(2,2);
        cm22(0,0) = std::complex<double>{1.0,2.0};
        cm22(0,1) = std::complex<double>{2.0,3.0};
        cm22(1,0) = std::complex<double>{3.0,4.0};
        cm22(1,1) = std::complex<double>{4.0,5.0};

        cm23.create(2,3);
        cm23(0,0) = std::complex<double>{1.0,2.0};
        cm23(0,1) = std::complex<double>{2.0,3.0};
        cm23(0,2) = std::complex<double>{3.0,4.0};
        cm23(1,0) = std::complex<double>{4.0,5.0};
        cm23(1,1) = std::complex<double>{5.0,6.0};
        cm23(1,2) = std::complex<double>{6.0,7.0};
    }
};

TEST(ComplexMatrix,Constructor)
{
    ModuleBase::ComplexMatrix cm;
    EXPECT_EQ(cm.nr,0);
    EXPECT_EQ(cm.nc,0);
    EXPECT_EQ(cm.size,0);
}

TEST(ComplexMatrix,ConstructorNrNc)
{
    ModuleBase::ComplexMatrix cm(3,4);
    EXPECT_EQ(cm.nr,3);
    EXPECT_EQ(cm.nc,4);
    EXPECT_EQ(cm.size,12);
}

TEST_F(ComplexMatrixTest,ConstructorCM)
{
    ModuleBase::ComplexMatrix cm1(cm22);
    EXPECT_EQ(cm1.nr,cm22.nr);
    EXPECT_EQ(cm1.nc,cm22.nc);
    EXPECT_EQ(cm1.size,cm22.size);
    for(int i=0; i<cm1.size; ++i)
    {
        EXPECT_COMPLEX_EQ(cm1.c[i],cm22.c[i]);
    }
}

TEST_F(ComplexMatrixTest,ConstructorCMrvalue)
{
    ModuleBase::ComplexMatrix cm2(cm22);
    ModuleBase::ComplexMatrix cm1(std::move(cm22));
    EXPECT_EQ(cm1.nr,cm2.nr);
    EXPECT_EQ(cm1.nc,cm2.nc);
    EXPECT_EQ(cm1.size,cm2.size);
    for(int i=0; i<cm1.size; ++i)
    {
        EXPECT_COMPLEX_EQ(cm1.c[i],cm2.c[i]);
    }

    EXPECT_EQ(cm22.nr,0);
    EXPECT_EQ(cm22.nc,0);
    EXPECT_EQ(cm22.size,0);
    EXPECT_EQ(cm22.c,nullptr);
}

TEST(ComplexMatrix,ConstrucotrMatrix)
{
    ModuleBase::matrix m(2,2);
    m(0,0) = 1.0; m(0,1) = 2.0;
    m(1,0) = 3.0; m(1,1) = 4.0;

    ModuleBase::ComplexMatrix cm(m);
    EXPECT_EQ(cm.nr,m.nr);
    EXPECT_EQ(cm.nc,m.nc);
    EXPECT_EQ(cm.size,m.nr * m.nc);

    EXPECT_COMPLEX_EQ(cm(0,0),std::complex<double>{1.0,0.0});
    EXPECT_COMPLEX_EQ(cm(0,1),std::complex<double>{2.0,0.0});
    EXPECT_COMPLEX_EQ(cm(1,0),std::complex<double>{3.0,0.0});
    EXPECT_COMPLEX_EQ(cm(1,1),std::complex<double>{4.0,0.0});
}

TEST(ComplexMatrix,Create)
{
    ModuleBase::ComplexMatrix cm;
    cm.create(111,222);
    EXPECT_EQ(cm.nr,111);
    EXPECT_EQ(cm.nc,222);
    EXPECT_EQ(cm.size,111*222);
}

TEST_F(ComplexMatrixTest,OperatorEqual)
{
    ModuleBase::ComplexMatrix cm;
    cm = cm22;
    EXPECT_EQ(cm.nr,cm22.nr);
    EXPECT_EQ(cm.nc,cm22.nc);
    EXPECT_EQ(cm.size,cm22.size);

    for(int i=0; i<cm.size; ++i)
    {
        EXPECT_COMPLEX_EQ(cm.c[i],cm22.c[i]);
    }
}

TEST_F(ComplexMatrixTest,OperatorEqualrvalue)
{
    ModuleBase::ComplexMatrix cm;
    cm = 2.0*cm22;
    EXPECT_EQ(cm.nr,cm22.nr);
    EXPECT_EQ(cm.nc,cm22.nc);
    EXPECT_EQ(cm.size,cm22.size);

    for(int i=0; i<cm.size; ++i)
    {
        EXPECT_COMPLEX_EQ(cm.c[i],cm22.c[i]*2.0);
    }
}

TEST_F(ComplexMatrixTest,OperatorParentheses)
{
    EXPECT_EQ(cm22(0,0),cm22.c[0]);
}

TEST_F(ComplexMatrixTest,OperatorParenthesesConst)
{
    const ModuleBase::ComplexMatrix cm(cm22);
    EXPECT_COMPLEX_EQ((cm)(0,0),cm22(0,0));
}


TEST_F(ComplexMatrixTest,OperatorMultiEqual)
{
    ModuleBase::ComplexMatrix cm(cm22);
    std::complex<double> com{1.0,2.0};
    cm *= com;
    EXPECT_EQ(cm.nr,cm22.nr);
    EXPECT_EQ(cm.nc,cm22.nc);
    EXPECT_EQ(cm.size,cm22.size);

    for(int i=0; i<cm.size; ++i)
    {
        EXPECT_COMPLEX_EQ(cm.c[i],cm22.c[i]*com);
    }
}

TEST_F(ComplexMatrixTest,OperatorPlusEqual)
{
    ModuleBase::ComplexMatrix cm(cm22);
    cm += cm22;
    EXPECT_EQ(cm.nr,cm22.nr);
    EXPECT_EQ(cm.nc,cm22.nc);
    EXPECT_EQ(cm.size,cm22.size);

    for(int i=0; i<cm.size; ++i)
    {
        EXPECT_COMPLEX_EQ(cm.c[i],cm22.c[i] + cm22.c[i]);
    }

    //EXPECT_DEATH(cm22 += cm23,"");
}

TEST_F(ComplexMatrixTest,OperatorMinusEqual)
{
    ModuleBase::ComplexMatrix cm(3.0*cm22);
    cm -= cm22;
    EXPECT_EQ(cm.nr,cm22.nr);
    EXPECT_EQ(cm.nc,cm22.nc);
    EXPECT_EQ(cm.size,cm22.size);

    for(int i=0; i<cm.size; ++i)
    {
        EXPECT_COMPLEX_EQ(cm.c[i],3.0*cm22.c[i] - cm22.c[i]);
    }

    //EXPECT_DEATH(cm22 -= cm23,"");
}

TEST_F(ComplexMatrixTest,real)
{
    ModuleBase::matrix m = cm22.real();
    EXPECT_EQ(m.nr,cm22.nr);
    EXPECT_EQ(m.nc,cm22.nc);
    EXPECT_COMPLEX_EQ(m(0,0),cm22(0,0).real());
    EXPECT_COMPLEX_EQ(m(0,1),cm22(0,1).real());
    EXPECT_COMPLEX_EQ(m(1,0),cm22(1,0).real());
    EXPECT_COMPLEX_EQ(m(1,1),cm22(1,1).real());
}

TEST_F(ComplexMatrixTest,ZeroOut)
{
    cm22.zero_out();
    for(int i=0; i<cm22.size; ++i)
    {
        EXPECT_COMPLEX_EQ(cm22.c[i],std::complex<double>{0.0,0.0});
    }
}


TEST_F(ComplexMatrixTest,SetAsIdentityMatrix)
{
    cm22.set_as_identity_matrix();
    EXPECT_COMPLEX_EQ(cm22(0,0),std::complex<double> {1.0,0.0});
    EXPECT_COMPLEX_EQ(cm22(0,1),std::complex<double> {0.0,0.0});
    EXPECT_COMPLEX_EQ(cm22(1,0),std::complex<double> {0.0,0.0});
    EXPECT_COMPLEX_EQ(cm22(1,1),std::complex<double> {1.0,0.0});
}

TEST(ComplexMatrix,CheckReal)
{
    ModuleBase::ComplexMatrix cm22(2,2);
    cm22(0,0) = std::complex<double> {0.0,0.0};
    cm22(0,1) = std::complex<double> {1.0,0.0};
    cm22(1,0) = std::complex<double> {2.0,0.0};
    cm22(1,1) = std::complex<double> {3.0,0.0};
    EXPECT_TRUE(cm22.checkreal());

    cm22(0,0) = std::complex<double> {0.0,0.01};
    EXPECT_FALSE(cm22.checkreal());
}

TEST_F(ComplexMatrixTest,OperatorPlus)
{
    ModuleBase::ComplexMatrix cm1(cm22),cm2;
    cm2 = cm1 + cm22;
    EXPECT_EQ(cm2.nr,cm22.nr);
    EXPECT_EQ(cm2.nc,cm22.nc);
    EXPECT_EQ(cm2.size,cm22.size);

    for(int i=0; i<cm2.size; ++i)
    {
        EXPECT_COMPLEX_EQ(cm2.c[i],cm22.c[i] + cm22.c[i]);
    }

    //EXPECT_DEATH(cm22 + cm23,"");
}

TEST_F(ComplexMatrixTest,OperatorMinus)
{
    ModuleBase::ComplexMatrix cm1(2.0*cm22),cm2;
    cm2 = cm1 - cm22;
    EXPECT_EQ(cm2.nr,cm22.nr);
    EXPECT_EQ(cm2.nc,cm22.nc);
    EXPECT_EQ(cm2.size,cm22.size);

    for(int i=0; i<cm2.size; ++i)
    {
        EXPECT_COMPLEX_EQ(cm2.c[i],cm22.c[i]);
    }

    //EXPECT_DEATH(cm22 + cm23,"");
}

TEST_F(ComplexMatrixTest,OperatorMultMatrix)
{
    ModuleBase::ComplexMatrix cm23(2,3),cm32(3,2),cm22,cm33;
    cm23(0,0)=std::complex<double>{1.0,1.0};
    cm23(0,1)=std::complex<double>{2.0,0.0};
    cm23(0,2)=std::complex<double>{3.0,-1.0};
    cm23(1,0)=std::complex<double>{4.0,-2.0};
    cm23(1,1)=std::complex<double>{5.0,-3.0};
    cm23(1,2)=std::complex<double>{6.0,-4.0};

    cm32(0,0)=std::complex<double>{-11.0,11.0};
    cm32(0,1)=std::complex<double>{-12.0,12.0};
    cm32(1,0)=std::complex<double>{-13.0,13.0};
    cm32(1,1)=std::complex<double>{-14.0,14.0};
    cm32(2,0)=std::complex<double>{-15.0,15.0};
    cm32(2,1)=std::complex<double>{-16.0,16.0};

    cm22 = cm23 * cm32;
    EXPECT_EQ(cm22.nr,2);
    EXPECT_EQ(cm22.nc,2);
    EXPECT_EQ(cm22.size,4);
    EXPECT_COMPLEX_EQ(cm22(0,0),std::complex<double>{-78.0,86.0});
    EXPECT_COMPLEX_EQ(cm22(0,1),std::complex<double>{-84.0,92.0});
    EXPECT_COMPLEX_EQ(cm22(1,0),std::complex<double>{-78.0,320.0});
    EXPECT_COMPLEX_EQ(cm22(1,1),std::complex<double>{-84.0,344.0});

    cm33 = cm32 * cm23;
    EXPECT_EQ(cm33.nr,3);
    EXPECT_EQ(cm33.nc,3);
    EXPECT_EQ(cm33.size,9);
    EXPECT_COMPLEX_EQ(cm33(0,0),std::complex<double>{-46.0,72.0  });
    EXPECT_COMPLEX_EQ(cm33(0,1),std::complex<double>{-46.0,118.0 });
    EXPECT_COMPLEX_EQ(cm33(0,2),std::complex<double>{-46.0,164.0 });
    EXPECT_COMPLEX_EQ(cm33(1,0),std::complex<double>{-54.0,84.0  });
    EXPECT_COMPLEX_EQ(cm33(1,1),std::complex<double>{-54.0,138.0 });
    EXPECT_COMPLEX_EQ(cm33(1,2),std::complex<double>{-54.0,192.0 });
    EXPECT_COMPLEX_EQ(cm33(2,0),std::complex<double>{-62.0,96.0  });
    EXPECT_COMPLEX_EQ(cm33(2,1),std::complex<double>{-62.0,158.0 });
    EXPECT_COMPLEX_EQ(cm33(2,2),std::complex<double>{-62.0,220.0 });

    EXPECT_DEATH(cm22 * cm32,"");
}

TEST_F(ComplexMatrixTest,OperatorMultDouble)
{
    ModuleBase::ComplexMatrix cm2,cm3;

    cm2 = cm22 * 2.0;
    EXPECT_EQ(cm2.nr,cm22.nr);
    EXPECT_EQ(cm2.nc,cm22.nc);
    EXPECT_EQ(cm2.size,cm22.size);
    for(int i=0; i<cm2.size; ++i)
    {
        EXPECT_COMPLEX_EQ(cm2.c[i],2.0 * cm22.c[i]);
    }

    cm3 = 3.0 * cm22;
    EXPECT_EQ(cm3.nr,cm22.nr);
    EXPECT_EQ(cm3.nc,cm22.nc);
    EXPECT_EQ(cm3.size,cm22.size);
    for(int i=0; i<cm3.size; ++i)
    {
        EXPECT_COMPLEX_EQ(cm3.c[i],3.0 * cm22.c[i]);
    }
}

TEST_F(ComplexMatrixTest,OperatorMultComplex)
{
    ModuleBase::ComplexMatrix cm2,cm3;
    std::complex<double> com {2.0,3.0};

    cm2 = cm22 * com;
    EXPECT_EQ(cm2.nr,cm22.nr);
    EXPECT_EQ(cm2.nc,cm22.nc);
    EXPECT_EQ(cm2.size,cm22.size);
    for(int i=0; i<cm2.size; ++i)
    {
        EXPECT_COMPLEX_EQ(cm2.c[i],com * cm22.c[i]);
    }

    cm3 = com * cm22;
    EXPECT_EQ(cm3.nr,cm22.nr);
    EXPECT_EQ(cm3.nc,cm22.nc);
    EXPECT_EQ(cm3.size,cm22.size);
    for(int i=0; i<cm3.size; ++i)
    {
        EXPECT_COMPLEX_EQ(cm3.c[i],com * cm22.c[i]);
    }
}

TEST_F(ComplexMatrixTest,Trace)
{
    EXPECT_COMPLEX_EQ(trace(cm22),std::complex<double>{5.0,7.0});
}

TEST_F(ComplexMatrixTest,abs2Row)
{
    EXPECT_EQ(abs2_row(cm22,0),18.0);
    EXPECT_EQ(abs2_row(cm22,1),66.0);
}

TEST_F(ComplexMatrixTest,abs2Column)
{
    EXPECT_EQ(abs2_column(cm22,0),30.0);
    EXPECT_EQ(abs2_column(cm22,1),54.0);
}

TEST_F(ComplexMatrixTest,abs2)
{
    EXPECT_EQ(abs2(cm22),84.0);
}

TEST_F(ComplexMatrixTest,abs2arraymatrix)
{
    ModuleBase::ComplexMatrix **m;
    m = new ModuleBase::ComplexMatrix*[2];
    m[0] = &cm22;
    m[1] = &cm23;
    EXPECT_EQ(abs2(1,m),84.0);
    EXPECT_EQ(abs2(2,m),314.0);
    delete [] m;
}

TEST_F(ComplexMatrixTest,transpose)
{
    ModuleBase::ComplexMatrix m(transpose(cm22,false));
    EXPECT_COMPLEX_EQ(m(0,0),cm22(0,0));
    EXPECT_COMPLEX_EQ(m(0,1),cm22(1,0));
    EXPECT_COMPLEX_EQ(m(1,0),cm22(0,1));
    EXPECT_COMPLEX_EQ(m(1,1),cm22(1,1));

    ModuleBase::ComplexMatrix m1(transpose(cm22,true));
    EXPECT_COMPLEX_EQ(m1(0,0),conj(cm22(0,0)));
    EXPECT_COMPLEX_EQ(m1(0,1),conj(cm22(1,0)));
    EXPECT_COMPLEX_EQ(m1(1,0),conj(cm22(0,1)));
    EXPECT_COMPLEX_EQ(m1(1,1),conj(cm22(1,1)));
}

TEST_F(ComplexMatrixTest,conj)
{
    ModuleBase::ComplexMatrix m = conj(cm22);
    EXPECT_COMPLEX_EQ(m(0,0),std::conj(cm22(0,0)));
    EXPECT_COMPLEX_EQ(m(0,1),std::conj(cm22(0,1)));
    EXPECT_COMPLEX_EQ(m(1,0),std::conj(cm22(1,0)));
    EXPECT_COMPLEX_EQ(m(1,1),std::conj(cm22(1,1)));
}

TEST_F(ComplexMatrixTest,ScaleAccumulate)
{
    std::complex<double> com1{2.0,2.0};
    ModuleBase::ComplexMatrix m(2.0*cm22),m1;
    m1 = m;
    ModuleBase::scale_accumulate(com1,cm22,m);
    for(int i=0; i<m.size; ++i)
    {
        EXPECT_COMPLEX_EQ(m.c[i],m1.c[i] + com1 * cm22.c[i]);
    }

    EXPECT_DEATH(ModuleBase::scale_accumulate(com1,cm22,cm23),"");
}

TEST_F(ComplexMatrixTest,ScaleAccumulateArray)
{
    std::complex<double> com1{2.0,2.0};
    ModuleBase::ComplexMatrix cm1(cm22),cm2(2.0*cm22),cm3(3.0*cm22),cm4(4.0*cm22);
    ModuleBase::ComplexMatrix **cmout;
    ModuleBase::ComplexMatrix **cmin;
    cmout = new ModuleBase::ComplexMatrix*[2];
    cmin  = new ModuleBase::ComplexMatrix*[2];

    cmout[0] = &cm1;
    cmout[1] = &cm2;
    cmin[0] = &cm3;
    cmin[1] = &cm4;
    int size = cm22.size;

    ModuleBase::scale_accumulate(2,com1,cmin,cmout);

    for(int i=0; i<size; ++i)
    {
        EXPECT_COMPLEX_EQ((*cmout[0]).c[i], cm22.c[i] + com1 * (*cmin[0]).c[i] );
        EXPECT_COMPLEX_EQ((*cmout[1]).c[i], 2.0 * cm22.c[i] + com1 * (*cmin[1]).c[i] );
    }

    delete [] cmout;
    delete [] cmin;
}


TEST_F(ComplexMatrixTest,ScaleSum)
{
    std::complex<double> com1{2.0,2.0},com2{-1.0,-2.0};
    ModuleBase::ComplexMatrix cm1(1.1*cm22),cm2(2.5*cm22),cmout(2,2);

    ModuleBase::scaled_sum(com1,cm1,com2,cm2,cmout);
    for(int i=0; i<cmout.size; ++i)
    {
        EXPECT_COMPLEX_EQ(cmout.c[i], com1 * cm1.c[i] + com2 * cm2.c[i]);

    }

    EXPECT_DEATH(ModuleBase::scaled_sum(com1,cm22,com2,cm23,cmout),"");
}



TEST_F(ComplexMatrixTest,ScaleSumArray)
{
    std::complex<double> com1{2.0,2.0},com2{-1.0,-2.0};
    ModuleBase::ComplexMatrix cm1(cm22),cm2(2.0*cm22),cm3(3.0*cm22),cm4(4.0*cm22);
    ModuleBase::ComplexMatrix cm5(2,2),cm6(2,2);
    ModuleBase::ComplexMatrix **cmout;
    ModuleBase::ComplexMatrix **cmin1;
    ModuleBase::ComplexMatrix **cmin2;
    cmout = new ModuleBase::ComplexMatrix*[2];
    cmin1 = new ModuleBase::ComplexMatrix*[2];
    cmin2 = new ModuleBase::ComplexMatrix*[2];

    cmin1[0] = &cm1;
    cmin1[1] = &cm2;
    cmin2[0] = &cm3;
    cmin2[1] = &cm4;
    cmout[0] = &cm5;
    cmout[1] = &cm6;

    int size = cm22.size;

    ModuleBase::scaled_sum(2,com1,cmin1,com2,cmin2,cmout);

    for(int i=0; i<size; ++i)
    {
        EXPECT_COMPLEX_EQ((*cmout[0]).c[i], com1 * (*cmin1[0]).c[i] + com2 * (*cmin2[0]).c[i] );
        EXPECT_COMPLEX_EQ((*cmout[1]).c[i], com1 * (*cmin1[1]).c[i] + com2 * (*cmin2[1]).c[i] );
    }

    delete [] cmout;
    delete [] cmin1;
    delete [] cmin2;
}

TEST_F(ComplexMatrixTest,print)
{
   std::ifstream ifs;
   std::ofstream ofs;
   ofs.open("printtest1.log");
   cm22.print(ofs,1e-10,1e-10);
   ofs.close();
   ifs.open("printtest1.log");
   std::string output;
   getline(ifs,output);
   EXPECT_THAT(output,testing::HasSubstr("(1,2)\t(2,3)\t"));
   getline(ifs,output);
   EXPECT_THAT(output,testing::HasSubstr("(3,4)\t(4,5)\t"));
   ifs.close();
   remove("printtest1.log");
// The condition of  std::abs(data)>threshold_abs && std::imag(data)) <= threshold_imag
   ofs.open("printtest2.log");
   cm22.print(ofs,1e-10,2);
   ofs.close();
   ifs.open("printtest2.log");
   getline(ifs,output);
   EXPECT_THAT(output,testing::HasSubstr("1\t(2,3)\t"));
   getline(ifs,output);
   EXPECT_THAT(output,testing::HasSubstr("(3,4)\t(4,5)\t"));
   ifs.close();
   remove("printtest2.log");
// The condition of  std::abs(data)<threshold_abs && std::imag(data)) > threshold_imag
   ofs.open("printtest3.log");
   cm22.print(ofs,3,1e-10);
   ofs.close();
   ifs.open("printtest3.log");
   getline(ifs,output);
   EXPECT_THAT(output,testing::HasSubstr("0\t(2,3)\t"));
   getline(ifs,output);
   EXPECT_THAT(output,testing::HasSubstr("(3,4)\t(4,5)\t"));
   ifs.close();
   remove("printtest3.log");
}
