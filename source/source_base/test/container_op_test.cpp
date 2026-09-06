#include "../container_operator.h"
#include "gtest/gtest.h"
#include <iostream>

/************************************************
 *  unit test of container operator
 ***********************************************/

/**
 * - Tested Functions:
 *   - VectorPlus
 *     - "+" operator for vectors
 *   - VectorMinus
 *     - "-" operator for vectors
 *   - VectorMultiply
 *     - "*" operator for scalar*vector
 *   - MapPlus
 *     - "+" operator for maps
 *   - MapMinus
 *     - "-" operator for maps
 *   - MapMultiply
 *     - "*" operator for scalar*map
 */

TEST(ContainerOperator,VectorPlus)
{
	std::vector<double> a(20,10.0);
	std::vector<double> b(20,1.0);
	std::vector<double> c(20);
	c = a+b;
	for (int i=0;i<c.size();i++)
		EXPECT_EQ(c[i],11.0);
}

TEST(ContainerOperator,VectorMinus)
{
	std::vector<double> a(20,10.0);
	std::vector<double> b(20,1.0);
	std::vector<double> c(20);
	c = a-b;
	for (int i=0;i<c.size();i++)
		EXPECT_EQ(c[i],9.0);
}

TEST(ContainerOperator,VectorMultiply)
{
	std::vector<double> a(20,10.0);
	double b = 2.0;
	std::vector<double> c(20);
	c = b*a;
	for (int i=0;i<c.size();i++)
		EXPECT_EQ(c[i],20.0);
}

TEST(ContainerOperator,VectorLengthCheck)
{
	std::vector<double> a(20,10.0);
	std::vector<double> b(19,1.0);
	std::vector<double> c(20);
	EXPECT_DEATH(c=a+b,"");
	EXPECT_DEATH(c=a-b,"");
}

TEST(ContainerOperator,MapLengthCheck)
{
	std::map<int,double> a;
	std::map<int,double> b;
	std::map<int,double> c;
	for (int i=0;i<10;i++)
	{
		a.insert(std::pair<int,double> (i, i*2.0));
		if (i<9) b.insert(std::pair<int,double> (i, i*3.0));
	}
	EXPECT_DEATH(c=a+b,"");
	EXPECT_DEATH(c=a-b,"");
}

TEST(ContainerOperator,MapPlus)
{
	std::map<int,double> a;
	std::map<int,double> b;
	std::map<int,double> c;
	for (int i=0;i<10;i++)
	{
		a.insert(std::pair<int,double> (i, i*2.0));
		b.insert(std::pair<int,double> (i, i*3.0));
	}
	c = a+b;
	for (int i=0;i<10;i++)
	{
		//std::cout << c[i] << std::endl;
		EXPECT_EQ(c[i],i*5.0);
	}
}

TEST(ContainerOperator,MapMinus)
{
	std::map<int,double> a;
	std::map<int,double> b;
	std::map<int,double> c;
	for (int i=0;i<10;i++)
	{
		a.insert(std::pair<int,double> (i, i*4.0));
		b.insert(std::pair<int,double> (i, i*2.0));
	}
	c = a-b;
	for (int i=0;i<10;i++)
	{
		//std::cout << c[i] << std::endl;
		EXPECT_EQ(c[i],i*2.0);
	}
}

TEST(ContainerOperator,MapMultiply)
{
	std::map<int,double> a;
	double b = 3.0;
	std::map<int,double> c;
	for (int i=0;i<10;i++)
	{
		a.insert(std::pair<int,double> (i, i*4.0));
	}
	c = b*a;
	for (int i=0;i<10;i++)
	{
		//std::cout << c[i] << std::endl;
		EXPECT_EQ(c[i],i*4.0*3.0);
	}
}
