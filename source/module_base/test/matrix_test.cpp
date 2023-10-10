#include"../matrix.h"
#include"gtest/gtest.h"
#include "gmock/gmock.h"

/************************************************
*  unit test of class matrix and related functions
***********************************************/

/**
 *  - Tested functions of class matrix:
 *      - constructor:
 *          - constructed by nrow and ncloumn
 *          - constructed by a matrix
 *          - constructed by the rvalue of a matrix
 *      - function create
 *      - operator "=": assigned by a matrix or the rvalue of a matrix
 *      - operator "()": access the element
 *      - operator "*=", "+=", "-="
 *      - function trace_on
 *      - function zero_out
 *      - function fill_out(fill a matrix with a const double)
 *      - function max/min/absmax
 *      - function norm
 *      - function print(print the element which is larger than threshold)
 *      - function reshape(change the index of the array)
 * 
 * - Tested functions related to class matrix
 *      - operator "+", "-", "*" between two matrixs
 *      - operator "*" between a double and a matrix, and reverse.
 *      - function transpose
 *      - function trace_on
 *      - function mdot
 *      - function matrixAlloc
 */

//a mock function of WARNING_QUIT, to avoid the uncorrected call by matrix.cpp at line 37.
namespace ModuleBase
{
    void matrixAlloc();
}

class matrixTest : public testing::Test
{
    protected:
    ModuleBase::matrix m23a,m33a,m33b,m33c,m34a,m34b;

    void SetUp()
    {
        m23a.create(2,3);
        for (int i=1;i<=6;++i) {m23a.c[i-1] = i*1.0;}

        m33a.create(3,3);
        for (int i=1;i<=9;++i) {m33a.c[i-1] = i*1.0;}

        m33b.create(3,3);
        for (int i=1;i<=9;++i) {m33b.c[i-1] = i*11.1;}

        m33c.create(3,3,true);
        m34a.create(3,4,true);
        m34b.create(3,4,true);
    }

};

TEST(matrix,ConstructorNrNc)
{
    ModuleBase::matrix m(3,4,true);
    EXPECT_EQ(m.nr,3);
    EXPECT_EQ(m.nc,4);
    EXPECT_DOUBLE_EQ(m(0,0),0.0);   
}

TEST_F(matrixTest,ConstructorMatrix)
{
    ModuleBase::matrix m(m33a);
    int mnr = m.nr;
    EXPECT_EQ(mnr,m33a.nr);
    EXPECT_EQ(m.nc,m33a.nc);
    for (int i=0;i<9;++i)
    {
        EXPECT_DOUBLE_EQ(m.c[i],m33a.c[i]);
    }    
}

TEST_F(matrixTest,ConstructorMtrixRValue)
{

    ModuleBase::matrix m(3.0*m33a);
    EXPECT_EQ(m.nr,m33a.nr);
    EXPECT_EQ(m.nc,m33a.nc);
    for (int i=0;i<9;++i)
    {
        EXPECT_DOUBLE_EQ(m.c[i],m33a.c[i] * 3.0);
    }    
}

TEST_F(matrixTest,Create)
{
    m33a.create(13,14,true);
    EXPECT_EQ(m33a.nr,13);
    EXPECT_EQ(m33a.nc,14);
    for(int i=0;i<13*14;++i)
    {
        EXPECT_DOUBLE_EQ(m33a.c[i],0.0);
    }
}

TEST_F(matrixTest,OperatorEqualMatrix)
{
    ModuleBase::matrix m;
    m = m33a;
    EXPECT_EQ(m.nr,m33a.nr);
    EXPECT_EQ(m.nc,m33a.nc);
    for (int i=0;i<9;++i)
    {
        EXPECT_DOUBLE_EQ(m.c[i],m33a.c[i]);
    }   

    m23a = m33a;
    EXPECT_EQ(m23a.nr,m33a.nr);
    EXPECT_EQ(m23a.nc,m33a.nc); 
}

TEST_F(matrixTest,OperatorEqualMatrixRvalue)
{
    ModuleBase::matrix m;
    m = 3.0 * m33a;
    EXPECT_EQ(m.nr,m33a.nr);
    EXPECT_EQ(m.nc,m33a.nc);
    for (int i=0;i<9;++i)
    {
        EXPECT_DOUBLE_EQ(m.c[i],m33a.c[i] * 3.0);
    }    
}

TEST_F(matrixTest,OperatorParentheses)
{
    //EXPECT_DEATH(m33a(3,3),"");
    //EXPECT_DEATH(m33a(-1,0),"");
    m33a(0,0) = 1.1;
    EXPECT_DOUBLE_EQ(m33a(0,0),1.1);
}

TEST_F(matrixTest,OperatorMultiplyEqual)
{
    m33b = m33a;
    m33a *= 11.1;
    for (int i=0;i<m33a.nr;++i)
    {
        for (int j=0;j<m33a.nc;++j)
        {
            EXPECT_DOUBLE_EQ(m33a(i,j),m33b(i,j)*11.1);
        }
    }
}

TEST_F(matrixTest,OperatorPlusEqual)
{
    EXPECT_DEATH(m33a += m34a,"");

    m33c = m33a;
    m33a += m33b;
    for (int i=0;i<m33a.nr;++i)
    {
        for (int j=0;j<m33a.nc;++j)
        {
            EXPECT_DOUBLE_EQ(m33a(i,j),m33b(i,j)+m33c(i,j));
        }
    }
}

TEST_F(matrixTest,OperatorMinusEqual)
{
    EXPECT_DEATH(m33a -= m34a,"");
    
    m33c = m33a;
    m33a -= m33b;
    for (int i=0;i<m33a.nr;++i)
    {
        for (int j=0;j<m33a.nc;++j)
        {
            EXPECT_DOUBLE_EQ(m33a(i,j),m33c(i,j)-m33b(i,j));
        }
    }
}

TEST_F(matrixTest,ClassMatrixTraceOn)
{
    m33a(0,0) = 1.1;
    m33a(1,1) = 2.2;
    m33a(2,2) = 5.5;
    EXPECT_DOUBLE_EQ(m33a.trace_on(),8.8);
}

TEST_F(matrixTest,ClassMatrixZeroOut)
{
    m33a.zero_out();
    for (int i=0;i<m33a.nr;++i)
    {
        for (int j=0;j<m33a.nc;++j)
        {
            EXPECT_DOUBLE_EQ(m33a(i,j),0.0);
        }
    }
}

TEST_F(matrixTest,ClassMatrixMaxMinAbsmax)
{
    m33a(1,0) = 9999.9;
    m33a(0,0) = -999999.9;
    EXPECT_DOUBLE_EQ(m33a.max(),9999.9);
    EXPECT_DOUBLE_EQ(m33a.min(),-999999.9);
    EXPECT_DOUBLE_EQ(m33a.absmax(),999999.9);
}

TEST_F(matrixTest,ClassMatrixNorm)
{
    EXPECT_NEAR(m33a.norm(),16.881943016134133728,1E-12);
    EXPECT_NEAR(m33b.norm(),187.38956747908889611,1E-12);
}

TEST_F(matrixTest,OperatorPlus)
{
    m33c = m33a + m33b;
    for (int i=0;i<m33a.nr;++i)
    {
        for (int j=0;j<m33a.nc;++j)
        {
            EXPECT_DOUBLE_EQ(m33c(i,j),m33a(i,j) + m33b(i,j));
        }
    }
    //m33c = m33a + m34a;
    EXPECT_DEATH(m33a+m34a,"");
}


TEST_F(matrixTest,OperatorMinus)
{
    m33c = m33a - m33b;
    for (int i=0;i<m33a.nr;++i)
    {
        for (int j=0;j<m33a.nc;++j)
        {
            EXPECT_DOUBLE_EQ(m33c(i,j),m33a(i,j) - m33b(i,j));
        }
    }
    EXPECT_DEATH(m33a-m34a,"");
}

TEST_F(matrixTest,OperatorMultiplyTwoMatrix)
{
    EXPECT_DEATH(m34a*m33a,"");

    m33a(0,0) = 1.0; m33a(0,1) = 2.0; m33a(0,2) = 3.0;
    m33a(1,0) = 4.0; m33a(1,1) = 5.0; m33a(1,2) = 6.0;
    m33a(2,0) = 7.0; m33a(2,1) = 8.0; m33a(2,2) = 9.0;

    m33b(0,0) = -3.0; m33b(0,1) = -2.0; m33b(0,2) = -1.0;
    m33b(1,0) = -6.0; m33b(1,1) = -5.0; m33b(1,2) = -4.0;
    m33b(2,0) = -9.0; m33b(2,1) = -8.0; m33b(2,2) = -7.0;

    m33c = m33a * m33b;
    EXPECT_DOUBLE_EQ(m33c(0,0),-42.0);
    EXPECT_DOUBLE_EQ(m33c(0,1),-36.0);
    EXPECT_DOUBLE_EQ(m33c(0,2),-30.0);
    EXPECT_DOUBLE_EQ(m33c(1,0),-96.0);
    EXPECT_DOUBLE_EQ(m33c(1,1),-81.0);
    EXPECT_DOUBLE_EQ(m33c(1,2),-66.0);
    EXPECT_DOUBLE_EQ(m33c(2,0),-150.0);
    EXPECT_DOUBLE_EQ(m33c(2,1),-126.0);
    EXPECT_DOUBLE_EQ(m33c(2,2),-102.0);

    m33c = m33b * m33a;
    EXPECT_DOUBLE_EQ(m33c(0,0),-18.0);
    EXPECT_DOUBLE_EQ(m33c(0,1),-24.0);
    EXPECT_DOUBLE_EQ(m33c(0,2),-30.0);
    EXPECT_DOUBLE_EQ(m33c(1,0),-54.0);
    EXPECT_DOUBLE_EQ(m33c(1,1),-69.0);
    EXPECT_DOUBLE_EQ(m33c(1,2),-84.0);
    EXPECT_DOUBLE_EQ(m33c(2,0),-90.0);
    EXPECT_DOUBLE_EQ(m33c(2,1),-114.0);
    EXPECT_DOUBLE_EQ(m33c(2,2),-138.0);
}

TEST_F(matrixTest,OperatorMultiplyDouble)
{
    m33b = 3.0 * m33a;
    m33c = m33a * 3.0;
    for(int i=1;i<=9;++i)
    {
        EXPECT_DOUBLE_EQ(m33b.c[i-1],m33a.c[i-1] * 3.0);
        EXPECT_DOUBLE_EQ(m33c.c[i-1],m33a.c[i-1] * 3.0);
    }
}

TEST_F(matrixTest,Transpose)
{
    m23a(0,0) = 1.0; m23a(0,1) = 2.0; m23a(0,2) = 3.0;
    m23a(1,0) = 4.0; m23a(1,1) = 5.0; m23a(1,2) = 6.0;
    ModuleBase::matrix m32a = transpose(m23a);

    EXPECT_DOUBLE_EQ(m32a(0,0),1.0);
    EXPECT_DOUBLE_EQ(m32a(0,1),4.0);
    EXPECT_DOUBLE_EQ(m32a(1,0),2.0);
    EXPECT_DOUBLE_EQ(m32a(1,1),5.0);
    EXPECT_DOUBLE_EQ(m32a(2,0),3.0);
    EXPECT_DOUBLE_EQ(m32a(2,1),6.0);
}

TEST_F(matrixTest,TraceOn)
{
    EXPECT_DEATH(trace_on(m23a,m33a),"");
    EXPECT_DEATH(trace_on(m34a,m33a),"");

    ModuleBase::matrix m32(3,2), m23(2,3);
    m23(0,0) = 1.0; m23(0,1) = 2.0; m23(0,2) = 3.0;
    m23(1,0) = 4.0; m23(1,1) = 5.0; m23(1,2) = 6.0;

    m32(0,0) = 12.0; m32(0,1) = 23.0;
    m32(1,0) = 34.0; m32(1,1) = 45.0;
    m32(2,0) = 56.0; m32(2,1) = 67.0;

    EXPECT_DOUBLE_EQ(trace_on(m23,m32),967.0);
    EXPECT_DOUBLE_EQ(trace_on(m32,m23),967.0);
}

TEST_F(matrixTest,MDot)
{
    EXPECT_DEATH(mdot(m23a,m33a),"");
    EXPECT_DEATH(mdot(m34a,m33a),"");
    EXPECT_DOUBLE_EQ(mdot(m33a,m33b),3163.5);
}

TEST_F(matrixTest,Alloc)
{
    std::string output;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(ModuleBase::matrixAlloc(), ::testing::ExitedWithCode(0),"");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("Allocation error for Matrix"));
}

TEST_F(matrixTest,Fillout)
{
    double k=2.4;
    m33a.fill_out(k);
    for (int i=0;i<m33a.nr;++i)
    {
        for (int j=0;j<m33a.nc;++j)
        {
            EXPECT_DOUBLE_EQ(m33a(i,j),k);
        }
    }
}


TEST_F(matrixTest,Print)
{
    std::string output;
    const double threshold=4.0;
    testing::internal::CaptureStdout();
    //EXPECT_EXIT(m33a.print(std::cout,threshold),::testing::ExitedWithCode(0),"");
    m33a.print(std::cout,threshold);
    output  = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output,testing::HasSubstr("0\t0\t0\t\n0\t5\t6\t\n7\t8\t9\t\n"));
}

TEST_F(matrixTest,Reshape)
{
    const int nr_new=3;
    const int nc_new=2;
    const bool flag_zero=true;
    m23a.reshape(nr_new,nc_new,flag_zero);
    ASSERT_DEATH(m33a.reshape(nr_new,nc_new),"");
    EXPECT_DOUBLE_EQ(m23a.nr,nr_new);
    EXPECT_DOUBLE_EQ(m23a.nc,nc_new);
    for (int i=0;i<m23a.nr;++i)
    {
        for (int j=0;j<m23a.nc;++j)
        {
            EXPECT_DOUBLE_EQ(m23a(i,j),0.0);
        }
    }
}