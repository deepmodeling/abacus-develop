#include "../matrix3.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <random>
#include <vector>

/************************************************
 *  unit test of class Matrix3
 ***********************************************/

/**
 * - Tested Functions:
 *   - Construct
 *     - two ways of constructing a 3x3 matrix
 *   - Identity
 *     - set a 3x3 matrix to identity matrix
 *   - Zero
 *     - set all elements of a 3x3 matrix to zero
 *   - Det
 *     - calculate the determinant of 3x3 matrix
 *   - Transpose
 *     - do the transpose of 3x3 matrix
 *   - Inverse
 *     - do the inverse of 3x3 matrix
 *   - Assignment
 *     - overload assignment operator "=" for 3x3 matrix
 *   - AddEqual
 *     - overload operator "+=" for 3x3 matrix
 *   - MinusEqual
 *     - overload operator "-=" for 3x3 matrix
 *   - MultiplyEqual
 *     - overload operator "*=" for (3x3 matrix) * scalar
 *   - OverEqual
 *     - overload operator "/=" for (3x3 matrix) / scaler
 *   - Print
 *     - print a 3x3 matrix
 *   - MaddM
 *     - overload operator "+"  for two 3x3 matrices
 *   - MminusM
 *     - overload operator "-"  for two 3x3 matrices
 *   - MoverNum
 *     - overload operator "/"  for a (3x3 matrix)/(scalar)
 *   - MmultiplyM
 *     - overload operator "*"  for (3x3 matrix)*(3x3 matrix)
 *   - MmultiplyNum
 *     - overload operator "*"  for (3x3 matrix)*(scalar)
 *     - and (scalar)*(3x3 matrix)
 *   - MmultiplyV
 *     - overload operator "*"  for (3x3 matrix)*(Vector3)
 *   - VmultiplyM
 *     - overload operator "*"  for (Vector3)*(3x3 matrix)
 *   - MeqM
 *     - overload operator "==" to assert
 *     - the equality between two 3x3 matrices
 *   - MneM
 *     - overload operator "!=" to assert
 *     - the inequality between two 3x3 matrices
 *   - ToMatrix
 *     - change the form of a 3x3 matrix from that of class Matrix3
 *     - to that of class matrix
 *	 -Template V Multiply M
 *	   -tamplate classes for (Vector3)*(3x3 matrix)
 */

class Matrix3Test : public testing::Test
{
protected:
	ModuleBase::Matrix3 matrix_a, matrix_a1, matrix_b;
	ModuleBase::Matrix3 get_random_matrix3()
	{
		std::vector<double> v(9);
		for (auto &i : v)
		{
			i = std::rand();
		}
		auto matrix_a = ModuleBase::Matrix3(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8]);
		return matrix_a;
	}
	// for capturing stdout
	std::string output;
};

TEST_F(Matrix3Test, Construct)
{
	// method 1
	ModuleBase::Matrix3 ma;
	EXPECT_EQ(ma.e11,1); EXPECT_EQ(ma.e12,0); EXPECT_EQ(ma.e13,0);
	EXPECT_EQ(ma.e21,0); EXPECT_EQ(ma.e22,1); EXPECT_EQ(ma.e23,0);
	EXPECT_EQ(ma.e31,0); EXPECT_EQ(ma.e32,0); EXPECT_EQ(ma.e33,1);
	// method 2
	ModuleBase::Matrix3 mb(1,2,3,4,5,6,7,8,9);
	EXPECT_EQ(mb.e11,1); EXPECT_EQ(mb.e12,2); EXPECT_EQ(mb.e13,3);
	EXPECT_EQ(mb.e21,4); EXPECT_EQ(mb.e22,5); EXPECT_EQ(mb.e23,6);
	EXPECT_EQ(mb.e31,7); EXPECT_EQ(mb.e32,8); EXPECT_EQ(mb.e33,9);
}

TEST_F(Matrix3Test, Idenity)
{
	ModuleBase::Matrix3 mb(1,2,3,4,5,6,7,8,9);
	mb.Identity();
	EXPECT_EQ(mb.e11,1); EXPECT_EQ(mb.e12,0); EXPECT_EQ(mb.e13,0);
	EXPECT_EQ(mb.e21,0); EXPECT_EQ(mb.e22,1); EXPECT_EQ(mb.e23,0);
	EXPECT_EQ(mb.e31,0); EXPECT_EQ(mb.e32,0); EXPECT_EQ(mb.e33,1);
}

TEST_F(Matrix3Test, Zero)
{
	ModuleBase::Matrix3 ma;
	ma.Zero();
	EXPECT_EQ(ma.e11,0); EXPECT_EQ(ma.e12,0); EXPECT_EQ(ma.e13,0);
	EXPECT_EQ(ma.e21,0); EXPECT_EQ(ma.e22,0); EXPECT_EQ(ma.e23,0);
	EXPECT_EQ(ma.e31,0); EXPECT_EQ(ma.e32,0); EXPECT_EQ(ma.e33,0);
}

TEST_F(Matrix3Test, Det)
{
	ModuleBase::Matrix3 ma;
	ma = get_random_matrix3();
	double determinant = ma.e11 * ma.e22 * ma.e33
				- ma.e11 * ma.e32 * ma.e23
				- ma.e12 * ma.e21 * ma.e33
				+ ma.e12 * ma.e31 * ma.e23
				+ ma.e13 * ma.e21 * ma.e32
				- ma.e13 * ma.e22 * ma.e31;
	EXPECT_DOUBLE_EQ(ma.Det(),determinant);
}

TEST_F(Matrix3Test, Transpose)
{
	ModuleBase::Matrix3 ma, mb;
	ma = get_random_matrix3();
	mb = ma.Transpose();
	EXPECT_EQ(ma.e11,mb.e11); EXPECT_EQ(ma.e12,mb.e21); EXPECT_EQ(ma.e13,mb.e31);
	EXPECT_EQ(ma.e21,mb.e12); EXPECT_EQ(ma.e22,mb.e22); EXPECT_EQ(ma.e23,mb.e32);
	EXPECT_EQ(ma.e31,mb.e13); EXPECT_EQ(ma.e32,mb.e23); EXPECT_EQ(ma.e33,mb.e33);
}

TEST_F(Matrix3Test, Inverse)
{
	ModuleBase::Matrix3 ma, mb;
	ma = get_random_matrix3();
	mb = ma.Inverse();
	EXPECT_NEAR( (ma.e11*mb.e11 + ma.e12*mb.e21 + ma.e13*mb.e31), 1.0, 1e-15);
	EXPECT_NEAR( (ma.e11*mb.e12 + ma.e12*mb.e22 + ma.e13*mb.e32), 0.0, 1e-15);
	EXPECT_NEAR( (ma.e11*mb.e13 + ma.e12*mb.e23 + ma.e13*mb.e33), 0.0, 1e-15);
	EXPECT_NEAR( (ma.e21*mb.e11 + ma.e22*mb.e21 + ma.e23*mb.e31), 0.0, 1e-15);
	EXPECT_NEAR( (ma.e21*mb.e12 + ma.e22*mb.e22 + ma.e23*mb.e32), 1.0, 1e-15);
	EXPECT_NEAR( (ma.e21*mb.e13 + ma.e22*mb.e23 + ma.e23*mb.e33), 0.0, 1e-15);
	EXPECT_NEAR( (ma.e31*mb.e11 + ma.e32*mb.e21 + ma.e33*mb.e31), 0.0, 1e-15);
	EXPECT_NEAR( (ma.e31*mb.e12 + ma.e32*mb.e22 + ma.e33*mb.e32), 0.0, 1e-15);
	EXPECT_NEAR( (ma.e31*mb.e13 + ma.e32*mb.e23 + ma.e33*mb.e33), 1.0, 1e-15);
}

TEST_F(Matrix3Test, Assignment)
{
	ModuleBase::Matrix3 ma, mb;
	ma = get_random_matrix3();
	mb = ma;
	EXPECT_EQ(ma.e11,mb.e11); EXPECT_EQ(ma.e12,mb.e12); EXPECT_EQ(ma.e13,mb.e13);
	EXPECT_EQ(ma.e21,mb.e21); EXPECT_EQ(ma.e22,mb.e22); EXPECT_EQ(ma.e23,mb.e23);
	EXPECT_EQ(ma.e31,mb.e31); EXPECT_EQ(ma.e32,mb.e32); EXPECT_EQ(ma.e33,mb.e33);
}

TEST_F(Matrix3Test, AddEqual)
{
	ModuleBase::Matrix3 ma, mb;
	ma = get_random_matrix3();
	mb += ma;
	EXPECT_EQ(ma.e11+1.0,mb.e11); EXPECT_EQ(ma.e12,mb.e12); EXPECT_EQ(ma.e13,mb.e13);
	EXPECT_EQ(ma.e21,mb.e21); EXPECT_EQ(ma.e22+1.0,mb.e22); EXPECT_EQ(ma.e23,mb.e23);
	EXPECT_EQ(ma.e31,mb.e31); EXPECT_EQ(ma.e32,mb.e32); EXPECT_EQ(ma.e33+1.0,mb.e33);
}

TEST_F(Matrix3Test, MinusEqual)
{
	ModuleBase::Matrix3 ma, mb;
	ma = get_random_matrix3();
	mb -= ma;
	EXPECT_EQ(1.0-ma.e11,mb.e11); EXPECT_EQ(-ma.e12,mb.e12); EXPECT_EQ(-ma.e13,mb.e13);
	EXPECT_EQ(-ma.e21,mb.e21); EXPECT_EQ(1.0-ma.e22,mb.e22); EXPECT_EQ(-ma.e23,mb.e23);
	EXPECT_EQ(-ma.e31,mb.e31); EXPECT_EQ(-ma.e32,mb.e32); EXPECT_EQ(1.0-ma.e33,mb.e33);
}

TEST_F(Matrix3Test, MultiplyEqual)
{
	ModuleBase::Matrix3 ma, mb;
	ma = get_random_matrix3();
	mb = ma;
	mb *= 3.0;
	EXPECT_EQ(ma.e11*3.0,mb.e11); EXPECT_EQ(ma.e12*3.0,mb.e12); EXPECT_EQ(ma.e13*3.0,mb.e13);
	EXPECT_EQ(ma.e21*3.0,mb.e21); EXPECT_EQ(ma.e22*3.0,mb.e22); EXPECT_EQ(ma.e23*3.0,mb.e23);
	EXPECT_EQ(ma.e31*3.0,mb.e31); EXPECT_EQ(ma.e32*3.0,mb.e32); EXPECT_EQ(ma.e33*3.0,mb.e33);
}

TEST_F(Matrix3Test, OverEqual)
{
	ModuleBase::Matrix3 ma, mb;
	ma = get_random_matrix3();
	mb = ma;
	mb /= 3.0;
	EXPECT_EQ(ma.e11/3.0,mb.e11); EXPECT_EQ(ma.e12/3.0,mb.e12); EXPECT_EQ(ma.e13/3.0,mb.e13);
	EXPECT_EQ(ma.e21/3.0,mb.e21); EXPECT_EQ(ma.e22/3.0,mb.e22); EXPECT_EQ(ma.e23/3.0,mb.e23);
	EXPECT_EQ(ma.e31/3.0,mb.e31); EXPECT_EQ(ma.e32/3.0,mb.e32); EXPECT_EQ(ma.e33/3.0,mb.e33);
}

TEST_F(Matrix3Test, Print)
{
	ModuleBase::Matrix3 ma;
	ma = get_random_matrix3();
	testing::internal::CaptureStdout();
	ma.print();
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("e"));
}

TEST_F(Matrix3Test, MaddM)
{
	ModuleBase::Matrix3 ma, mb, mc;
	ma = get_random_matrix3();
	mb = get_random_matrix3();
	mc = ma + mb;
	EXPECT_EQ(ma.e11+mb.e11, mc.e11);
	EXPECT_EQ(ma.e12+mb.e12, mc.e12);
	EXPECT_EQ(ma.e13+mb.e13, mc.e13);
	EXPECT_EQ(ma.e21+mb.e21, mc.e21);
	EXPECT_EQ(ma.e22+mb.e22, mc.e22);
	EXPECT_EQ(ma.e23+mb.e23, mc.e23);
	EXPECT_EQ(ma.e31+mb.e31, mc.e31);
	EXPECT_EQ(ma.e32+mb.e32, mc.e32);
	EXPECT_EQ(ma.e33+mb.e33, mc.e33);
}

TEST_F(Matrix3Test, MminusM)
{
	ModuleBase::Matrix3 ma, mb, mc;
	ma = get_random_matrix3();
	mb = get_random_matrix3();
	mc = ma - mb;
	EXPECT_EQ(ma.e11-mb.e11, mc.e11);
	EXPECT_EQ(ma.e12-mb.e12, mc.e12);
	EXPECT_EQ(ma.e13-mb.e13, mc.e13);
	EXPECT_EQ(ma.e21-mb.e21, mc.e21);
	EXPECT_EQ(ma.e22-mb.e22, mc.e22);
	EXPECT_EQ(ma.e23-mb.e23, mc.e23);
	EXPECT_EQ(ma.e31-mb.e31, mc.e31);
	EXPECT_EQ(ma.e32-mb.e32, mc.e32);
	EXPECT_EQ(ma.e33-mb.e33, mc.e33);
}

TEST_F(Matrix3Test, MoverNum)
{
	ModuleBase::Matrix3 ma, mb;
	ma = get_random_matrix3();
	mb = ma/3.0;
	EXPECT_EQ(ma.e11/3.0,mb.e11); EXPECT_EQ(ma.e12/3.0,mb.e12); EXPECT_EQ(ma.e13/3.0,mb.e13);
	EXPECT_EQ(ma.e21/3.0,mb.e21); EXPECT_EQ(ma.e22/3.0,mb.e22); EXPECT_EQ(ma.e23/3.0,mb.e23);
	EXPECT_EQ(ma.e31/3.0,mb.e31); EXPECT_EQ(ma.e32/3.0,mb.e32); EXPECT_EQ(ma.e33/3.0,mb.e33);
}

TEST_F(Matrix3Test, MmultiplyM)
{
	ModuleBase::Matrix3 ma, mb, mc;
	ma = get_random_matrix3();
	mb = get_random_matrix3();
	mc = ma * mb;
	EXPECT_EQ( (ma.e11*mb.e11 + ma.e12*mb.e21 + ma.e13*mb.e31), mc.e11);
	EXPECT_EQ( (ma.e11*mb.e12 + ma.e12*mb.e22 + ma.e13*mb.e32), mc.e12);
	EXPECT_EQ( (ma.e11*mb.e13 + ma.e12*mb.e23 + ma.e13*mb.e33), mc.e13);
	EXPECT_EQ( (ma.e21*mb.e11 + ma.e22*mb.e21 + ma.e23*mb.e31), mc.e21);
	EXPECT_EQ( (ma.e21*mb.e12 + ma.e22*mb.e22 + ma.e23*mb.e32), mc.e22);
	EXPECT_EQ( (ma.e21*mb.e13 + ma.e22*mb.e23 + ma.e23*mb.e33), mc.e23);
	EXPECT_EQ( (ma.e31*mb.e11 + ma.e32*mb.e21 + ma.e33*mb.e31), mc.e31);
	EXPECT_EQ( (ma.e31*mb.e12 + ma.e32*mb.e22 + ma.e33*mb.e32), mc.e32);
	EXPECT_EQ( (ma.e31*mb.e13 + ma.e32*mb.e23 + ma.e33*mb.e33), mc.e33);
}

TEST_F(Matrix3Test, MmultiplyNum)
{
	ModuleBase::Matrix3 ma, mb, mc;
	ma = get_random_matrix3();
	mb = ma;
	mc = ma;
	mb = ma*3.0;
	mc = 3.0*ma;
	EXPECT_EQ(ma.e11*3.0,mb.e11); EXPECT_EQ(ma.e12*3.0,mb.e12); EXPECT_EQ(ma.e13*3.0,mb.e13);
	EXPECT_EQ(ma.e21*3.0,mb.e21); EXPECT_EQ(ma.e22*3.0,mb.e22); EXPECT_EQ(ma.e23*3.0,mb.e23);
	EXPECT_EQ(ma.e31*3.0,mb.e31); EXPECT_EQ(ma.e32*3.0,mb.e32); EXPECT_EQ(ma.e33*3.0,mb.e33);
	EXPECT_EQ(ma.e11*3.0,mc.e11); EXPECT_EQ(ma.e12*3.0,mc.e12); EXPECT_EQ(ma.e13*3.0,mc.e13);
	EXPECT_EQ(ma.e21*3.0,mc.e21); EXPECT_EQ(ma.e22*3.0,mc.e22); EXPECT_EQ(ma.e23*3.0,mc.e23);
	EXPECT_EQ(ma.e31*3.0,mc.e31); EXPECT_EQ(ma.e32*3.0,mc.e32); EXPECT_EQ(ma.e33*3.0,mc.e33);
}

TEST_F(Matrix3Test, MmultiplyV)
{
	ModuleBase::Matrix3 ma;
	ModuleBase::Vector3<double> u(3.0,4.0,5.0);
	ModuleBase::Vector3<double> v;
	ma = get_random_matrix3();
	v = ma * u;
	EXPECT_EQ(v.x, u.x*ma.e11+u.y*ma.e12+u.z*ma.e13);
	EXPECT_EQ(v.y, u.x*ma.e21+u.y*ma.e22+u.z*ma.e23);
	EXPECT_EQ(v.z, u.x*ma.e31+u.y*ma.e32+u.z*ma.e33);
}

TEST_F(Matrix3Test, VmultiplyM)
{
	ModuleBase::Matrix3 ma;
	ModuleBase::Vector3<double> u(3.0,4.0,5.0);
	ModuleBase::Vector3<double> v;
	ma = get_random_matrix3();
	v = u*ma;
	EXPECT_EQ(v.x, u.x*ma.e11+u.y*ma.e21+u.z*ma.e31);
	EXPECT_EQ(v.y, u.x*ma.e12+u.y*ma.e22+u.z*ma.e32);
	EXPECT_EQ(v.z, u.x*ma.e13+u.y*ma.e23+u.z*ma.e33);
}

TEST_F(Matrix3Test, MeqM)
{
	matrix_a = get_random_matrix3();
	matrix_a1 = matrix_a;
	matrix_b = get_random_matrix3();
	EXPECT_TRUE(matrix_a == matrix_a1);
	EXPECT_FALSE(matrix_a == matrix_b);
}

TEST_F(Matrix3Test, MneM)
{
	matrix_a = get_random_matrix3();
	matrix_a1 = matrix_a;
	matrix_b = get_random_matrix3();
	EXPECT_FALSE(matrix_a != matrix_a1);
	EXPECT_TRUE(matrix_a != matrix_b);
}

TEST_F(Matrix3Test, ToMatrix)
{
	ModuleBase::Matrix3 ma;
	ModuleBase::matrix mb;
	ma = get_random_matrix3();
	mb = ma.to_matrix();
	EXPECT_EQ(ma.e11,mb(0,0)); EXPECT_EQ(ma.e12,mb(0,1)); EXPECT_EQ(ma.e13,mb(0,2));
	EXPECT_EQ(ma.e21,mb(1,0)); EXPECT_EQ(ma.e22,mb(1,1)); EXPECT_EQ(ma.e23,mb(1,2));
	EXPECT_EQ(ma.e31,mb(2,0)); EXPECT_EQ(ma.e32,mb(2,1)); EXPECT_EQ(ma.e33,mb(2,2));
}

TEST_F(Matrix3Test,TemplateVectorMultiplyMatrix)
{
	ModuleBase::Vector3<int> ui(1,2,3);
	ModuleBase::Vector3<double> ud(1.0,2.0,3.0);
	ModuleBase::Vector3<float> uf(1.0,2.0,3.0);
	ModuleBase::Vector3<double> vi,vd,vf;
	ModuleBase::Matrix3 ma;
	ma = get_random_matrix3();
	vi=ui*ma;
	vd=ud*ma;
	vf=uf*ma;
	EXPECT_DOUBLE_EQ(vi.x,ui.x * ma.e11 + ui.y * ma.e21 + ui.z * ma.e31);
	EXPECT_DOUBLE_EQ(vi.y,ui.x * ma.e12 + ui.y * ma.e22 + ui.z * ma.e32);
	EXPECT_DOUBLE_EQ(vi.z,ui.x * ma.e13 + ui.y * ma.e23 + ui.z * ma.e33);
	EXPECT_DOUBLE_EQ(vd.x,ud.x * ma.e11 + ud.y * ma.e21 + ud.z * ma.e31);
	EXPECT_DOUBLE_EQ(vd.y,ud.x * ma.e12 + ud.y * ma.e22 + ud.z * ma.e32);
	EXPECT_DOUBLE_EQ(vd.z,ud.x * ma.e13 + ud.y * ma.e23 + ud.z * ma.e33);
	EXPECT_DOUBLE_EQ(vf.x,uf.x * ma.e11 + uf.y * ma.e21 + uf.z * ma.e31);
	EXPECT_DOUBLE_EQ(vf.y,uf.x * ma.e12 + uf.y * ma.e22 + uf.z * ma.e32);
	EXPECT_DOUBLE_EQ(vf.z,uf.x * ma.e13 + uf.y * ma.e23 + uf.z * ma.e33);
}