#include"../complexarray.h"
#include"gtest/gtest.h"
#include"gmock/gmock.h"

/************************************************
*  unit test of class ComplexArray and related functions
***********************************************/

/**
* - Tested functions of class ComplexArray:
*   - constructor:
*       - ComplexArray(const int bnd1=0, const int bnd2=1, const int bnd3=1, const int bnd4=1)
*       - ComplexArray(const ComplexArray &cd)
*       - ComplexArray(ComplexArray &&cd)
*   - operator "=":
*       - assign a complex to all elements of a ComplexArray
*       - assign a ComplexArray to another ComplexArray
*       - rvalue reference to a ComplexArray
*
*   - operator "+":
*       - one ComplexArray plus another ComplexArray that has
*         the same dimension.
*       - throw error when one ComplexArray plus another ComplexArray 
*         that has different dimension.
*
*   - operator "+=":
*       - one ComplexArray plus another ComplexArray, and assign the 
*         value to first ComplexArray
*       - throw error when one ComplexArray plus another ComplexArray
*         that has different dimension.
*
*   - operator "-":
*       - one ComplexArray minus another ComplexArray that has
*         the same dimension.
*       - throw error when one ComplexArray minus another ComplexArray
*         that has different dimension.
*
*   - operator "-=":
*       - one ComplexArray minus another ComplexArray, and assign the
*         value to first ComplexArray
*       - throw error when one ComplexArray minus another ComplexArray
*         that has different dimension.
*
*   - operator "*":
*       - one ComplexArray is multiplied by a double
*       - one ComplexArray is multiplied by a complex
*       - one ComplexArray is multiplied by another ComplexArray that has same dimension
*       - throw error when one ComplexArray is miltiplied by another ComplexArray
*         that has different dimension.
*
*   - operator "*=":
*       similar as "*"
*
*   - operator "==":
*       judge if two ComplexArray is equal, all element is equal.
*
*   - operator "!=":
*       judge if two ComplexArray is not equal
*
*   - oprator "()":
*       - access the element
*
*   - function zero_out, negate, randomize, getBound/Size, create.
*
*   Test relative functions:
*   - overloading of operator "*". a double/complex multiply a Class ComplexArray.
*   - functon abs2(), return the sum of squares of all elements 
*   - function dot()
*   - function scale_accumulate()
*   - function scaled_sum()
*   - function point_mult()
*   - function complexArrayxAlloc
*   - overloading of the function scaled_sum(). Does cd3 = c1*cd1 + c2*cd2. c1 and c2 are complex numbers.
*   - overloading of operator "()". The operator is effectively overload by "const". 
*/

//compare two complex by using EXPECT_DOUBLE_EQ()
void EXPECT_COMPLEX_EQUAL(const std::complex<double>& a,const std::complex<double>& b)
{
    EXPECT_DOUBLE_EQ(a.real(),b.real());
    EXPECT_DOUBLE_EQ(a.imag(),b.imag());
}

namespace ModuleBase
{
	void complexArrayxAlloc();
}

class ComplexArray_test : public testing::Test
{
    protected:
        ModuleBase::ComplexArray a2,a4,b2,b4,c2,c4,d2;
        std::complex<double> com1 {1.0,2.0};
        std::complex<double> com2 {3.0,4.0};
        std::complex<double> com3 {-2.0,-3.0};
        std::complex<double> comzero {0.0,0.0};

        void SetUp()
        {
            a2 = ModuleBase::ComplexArray(2,1,1,1);  // if define this class as a matrix
            b2 = ModuleBase::ComplexArray(2,1,1,1);  // of 4 dimesions,
            c2 = ModuleBase::ComplexArray(2,1,1,1);  // is a2 +/-/* d2 allowed or not???
            d2 = ModuleBase::ComplexArray(1,2,1,1);  // it does not matter, this situation will not appear in ABACUS
            a4 = ModuleBase::ComplexArray(2,2,1,1);
            b4 = ModuleBase::ComplexArray(2,2,1,1);
            c4 = ModuleBase::ComplexArray(2,2,1,1);
        }

};


TEST(ComplexArray, constructor_bnd1234_test)
{
    ModuleBase::ComplexArray a(1,2,3,4);
    ASSERT_EQ(a.getSize(),24);
    ASSERT_EQ(a.getBound1(),1);
    ASSERT_EQ(a.getBound2(),2);
    ASSERT_EQ(a.getBound3(),3);
    ASSERT_EQ(a.getBound4(),4);
    /*
    ASSERT_DEATH(ModuleBase::ComplexArray a(0,2,3,4),"");
    ASSERT_DEATH(ModuleBase::ComplexArray a(1,0,3,4),"");
    ASSERT_DEATH(ModuleBase::ComplexArray a(1,2,0,4),"");
    ASSERT_DEATH(ModuleBase::ComplexArray a(1,2,3,0),"");
    */
}

TEST(ComplexArray, constructor_copy_test)
{
    ModuleBase::ComplexArray a(1,2,3,4);
    ModuleBase::ComplexArray b(a);
    ASSERT_EQ(b.getSize(),24);
    ASSERT_EQ(b.getBound1(),1);
    ASSERT_EQ(b.getBound2(),2);
    ASSERT_EQ(b.getBound3(),3);
    ASSERT_EQ(b.getBound4(),4);
       
    ASSERT_EQ(a.getSize(),24);
    ASSERT_EQ(a.getBound1(),1);
    ASSERT_EQ(a.getBound2(),2);
    ASSERT_EQ(a.getBound3(),3);
    ASSERT_EQ(a.getBound4(),4);

}

TEST(ComplexArray, constructor_rvalue_test)
{
    ModuleBase::ComplexArray b(std::move(ModuleBase::ComplexArray(1,2,3,4)));
    ASSERT_EQ(b.getSize(),24);
    ASSERT_EQ(b.getBound1(),1);
    ASSERT_EQ(b.getBound2(),2);
    ASSERT_EQ(b.getBound3(),3);
    ASSERT_EQ(b.getBound4(),4);    
}


TEST_F(ComplexArray_test,operator_equal_complex)
{   
    a2 = com1;
    for(int i = 0;i<a2.getSize();++i)
    {
        EXPECT_COMPLEX_EQUAL(a2.ptr[i],com1);
    }
}


TEST_F(ComplexArray_test,operator_equal_ComplexArray)
{   
    a2 = com1;
    ModuleBase::ComplexArray b = a2;
    for (int i = 0;i<a2.getSize();++i)
    {
        EXPECT_COMPLEX_EQUAL(b.ptr[i],a2.ptr[i]);
    }
}

TEST_F(ComplexArray_test,operator_equal_ComplexArray_rvalue)
{   
    a2 = com1;
    ModuleBase::ComplexArray b = a2 * 1.0;

     for (int i = 0;i<a2.getSize();++i)
    {
        EXPECT_COMPLEX_EQUAL(b.ptr[i],com1);
    }   

}

TEST_F(ComplexArray_test,operator_plus)
{
    a2 = com1;
    b2 = com2;
    c2 = a2 + b2;
    for(int i = 0;i<a2.getSize();++i) 
    {
        EXPECT_COMPLEX_EQUAL(c2.ptr[i],(com1+com2));
    }
    EXPECT_DEATH(a2+a4,"");
    //EXPECT_DEATH(a2+d2,"");
}

TEST_F(ComplexArray_test,operator_plus_equal)
{
    a2 = com1;
    b2 = com2;
    a2 += b2;
    for(int i = 0;i<a2.getSize();++i) 
    {
        EXPECT_COMPLEX_EQUAL(a2.ptr[i],(com1+com2));
    }
    EXPECT_DEATH(a2+=a4,"");
    //EXPECT_DEATH(a2+=d2,"");
}

TEST_F(ComplexArray_test,operator_minus)
{
    a2 = com1;
    b2 = com2;
    c2 = a2 - b2;
    for(int i = 0;i<a2.getSize();++i) 
    {
        EXPECT_COMPLEX_EQUAL(c2.ptr[i],(com1-com2));
    }
    EXPECT_DEATH(a2-a4,"");
    //EXPECT_DEATH(a2-d2,"");
}

TEST_F(ComplexArray_test,operator_minus_equal)
{      
    a2 = com1;
    b2 = com2;
    a2 -= b2;
    for(int i = 0;i<a2.getSize();++i) 
    {
        EXPECT_COMPLEX_EQUAL(a2.ptr[i],(com1-com2));
    }
    EXPECT_DEATH(a2-=a4,"");
    //EXPECT_DEATH(a2-=d2,"");
}


TEST_F(ComplexArray_test,operator_multiply_double)
{
    a2 = com1;
    c2 = a2 * 2.0;
    for (int i = 0;i<a2.getSize();++i)
    {
        EXPECT_COMPLEX_EQUAL(c2.ptr[i],com1*2.0);
    }  
}

TEST_F(ComplexArray_test,operator_multiply_complex)
{
    a2 = com1;
    c2 = a2 * com2;
    for (int i = 0;i<a2.getSize();++i)
    {
        EXPECT_COMPLEX_EQUAL(c2.ptr[i],com1*com2);
    }
}

TEST_F(ComplexArray_test,operator_multiply_equal_double)
{
    a2 = com1;
    a2 *= 3.0;
    for (int i = 0;i<a2.getSize();++i)
    {
        EXPECT_COMPLEX_EQUAL(a2.ptr[i],com1*3.0);
    }
}

TEST_F(ComplexArray_test,operator_multiply_equal_complex)
{
    a2 = com1;
    a2 *= com2;
    for (int i = 0;i<a2.getSize();++i)
    {
        EXPECT_COMPLEX_EQUAL(a2.ptr[i],com1*com2);
    }
}

TEST_F(ComplexArray_test,operator_multiply_equal_ComplexArray)
{
    a2 = com1;
    b2 = com2;
    
    a2 *= b2;
    for (int i = 0;i<a2.getSize();++i)
    {
        EXPECT_COMPLEX_EQUAL(a2.ptr[i],com1*com2);
    }

    EXPECT_DEATH(a2*=a4,"");
}

TEST_F(ComplexArray_test,operator_equalequal)
{
    a2 = com1;
    b2 = com1;
    a4 = com1;
    EXPECT_TRUE(a2 == b2);
    EXPECT_FALSE(a2 == a4);
}

TEST_F(ComplexArray_test,operator_notequal)
{
    a2 = com1;
    b2 = com1;
    a4 = com1;
    EXPECT_TRUE(a2 != a4);
    EXPECT_FALSE(a2 != b2);
}

TEST_F(ComplexArray_test,operator_parentheses)
{
    a2 = com1;
    c2 = a2;
    EXPECT_COMPLEX_EQUAL(c2(0,0,0,0), com1);

    c2(1,0,0,0) = com2;
    EXPECT_NE(c2,a2);

    EXPECT_DEATH(a2(1,1,1,0),"");
}

TEST_F(ComplexArray_test,zero_out)
{
    a2 = com1;
    a2.zero_out();
    for (int i=0;i<a2.getSize();++i)
    {
        EXPECT_COMPLEX_EQUAL(a2.ptr[i],comzero);
    }
}

TEST_F(ComplexArray_test,negate)
{
    a2 = com1;
    a2.negate();
    for (int i=0;i<a2.getSize();++i)
    {
        EXPECT_COMPLEX_EQUAL(a2.ptr[i],com1*-1.0);
    }
}

TEST(ComplexArray,randomize)
{
    ModuleBase::ComplexArray a(10,10,10,10);
    a.randomize();

    bool allequal = true;
    for (int i=1;i<5;++i)
    {
        allequal = allequal && (a.ptr[i] == a.ptr[0]);
    }
    EXPECT_FALSE(allequal);

    for (int i=0;i<a.getSize();++i)
    {
        EXPECT_TRUE((a.ptr[i].real() >= -0.5) && (a.ptr[i].imag() >= -0.5));
        EXPECT_TRUE((a.ptr[i].real() < 0.5) && (a.ptr[i].imag() < 0.5));
    }

}

TEST(ComplexArray,getBoundSize)
{
    ModuleBase::ComplexArray a(1,2,3,4);
    EXPECT_EQ(a.getSize(),24);
    EXPECT_EQ(a.getBound1(),1);
    EXPECT_EQ(a.getBound2(),2);
    EXPECT_EQ(a.getBound3(),3);
    EXPECT_EQ(a.getBound4(),4);
}


TEST_F(ComplexArray_test,create)
{
    a2.create(2,3,4,5);
    EXPECT_EQ(a2.getSize(),120);
}

TEST_F(ComplexArray_test,operator_double_multiply)
{
    a2 = com1;
    c2 = 2.0 * a2 ;
    for (int i = 0;i<a2.getSize();++i)
    {
        EXPECT_COMPLEX_EQUAL(c2.ptr[i],2.0 * com1);
    }  
}

TEST_F(ComplexArray_test,operator_complex_multiply)
{
    a2 = com1;
    c2 = com2 * a2 ;
    for (int i = 0;i<a2.getSize();++i)
    {
        EXPECT_COMPLEX_EQUAL(c2.ptr[i],com2 * com1);
    }  
}

TEST(ComplexArray,abs2)
{
    ModuleBase::ComplexArray a(2,2,1,1);
    a(0,0,0,0) = std::complex<double>{0.0,0.0};
    a(1,0,0,0) = std::complex<double>{1.0,1.0};
    a(0,1,0,0) = std::complex<double>{2.0,2.0};
    a(1,1,0,0) = std::complex<double>{3.0,3.0};
    EXPECT_DOUBLE_EQ(abs2(a),28.0);
}

TEST(ComplexArray,dot)
{
    ModuleBase::ComplexArray a(2,2,1,1);
    a(0,0,0,0) = std::complex<double>{0.0,1.0};
    a(1,0,0,0) = std::complex<double>{1.0,0.0};
    a(0,1,0,0) = std::complex<double>{2.0,3.0};
    a(1,1,0,0) = std::complex<double>{3.0,2.0};

    ModuleBase::ComplexArray b(2,2,1,1);
    b(0,0,0,0) = std::complex<double>{1.0,0.0};
    b(1,0,0,0) = std::complex<double>{2.0,-1.0};
    b(0,1,0,0) = std::complex<double>{3.0,-2.0};
    b(1,1,0,0) = std::complex<double>{4.0,-3.0};
    std::complex<double> expectab {8.0,-32.0};
    std::complex<double> expectba {8.0,32.0};
    EXPECT_COMPLEX_EQUAL(dot(a,b),expectab);
    EXPECT_COMPLEX_EQUAL(dot(b,a),expectba);
}

TEST_F(ComplexArray_test,scale_accumulate_double)
{
    a2 = com1;
    b2 = com2;
    scale_accumulate(0.3,a2,b2);
    for (int i = 0;i<b2.getSize();++i)
    {
        EXPECT_COMPLEX_EQUAL(b2.ptr[i],0.3*com1+com2);
    }  
}

TEST_F(ComplexArray_test,scale_accumulate_complex)
{
    a2 = com1;
    b2 = com2;
    scale_accumulate(com3,a2,b2);
    for (int i = 0;i<b2.getSize();++i)
    {
        EXPECT_COMPLEX_EQUAL(b2.ptr[i],com3*com1+com2);
    }  
}

TEST_F(ComplexArray_test,scaled_sum)
{
    a2 = com1;
    b2 = com2;
    scaled_sum(0.3,a2,0.4,b2,c2);
    for (int i = 0;i<c2.getSize();++i)
    {
        EXPECT_COMPLEX_EQUAL(c2.ptr[i],0.3*com1 + 0.4*com2);
    }  
}

TEST_F(ComplexArray_test,point_mult)
{
    a2 = com1;
    b2 = com2;
    point_mult(a2,b2,c2);
    for (int i = 0;i<c2.getSize();++i)
    {
        EXPECT_COMPLEX_EQUAL(c2.ptr[i],com1*com2);
    }  
}

TEST_F(ComplexArray_test,xAlloc)
{
	std::string output;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(ModuleBase::complexArrayxAlloc(),::testing::ExitedWithCode(0),"");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("Allocation error for complexArray"));
}

TEST_F(ComplexArray_test,scaled_sum_2)
{
    a2 = com1;
    b2 = com2;
    scaled_sum(com1,a2,com2,b2,c2);
    for (int i = 0;i<c2.getSize();++i)
    {
        EXPECT_COMPLEX_EQUAL(c2.ptr[i],com1*com1 + com2*com2);
    }  
    scaled_sum(comzero,a2,comzero,b2,c2);
    for (int i = 0;i<c2.getSize();++i)
    {
        EXPECT_COMPLEX_EQUAL(c2.ptr[i],0);
    } 
    scaled_sum(com3,a2,com3,b2,c2);
    for (int i = 0;i<c2.getSize();++i)
    {
        EXPECT_COMPLEX_EQUAL(c2.ptr[i],com3*com1 + com3*com2);
    } 
}

TEST_F(ComplexArray_test,operator_parentheses_const)
{
    a2=com1;
    const ModuleBase::ComplexArray ca2=a2;
    EXPECT_COMPLEX_EQUAL(ca2(0,0,0,0), com1);
}
