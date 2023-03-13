#include "../vector3.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"

/************************************************
 *  unit test of class Vector3
 ***********************************************/

/**
 * - Tested Functions:
 *   - Construct
 *     - two ways of constructing a 3d vector
 *   - Set
 *     - set a 3d vector
 *   - Equal
 *     - overload operator "=" for 3d vector
 *   - equal
 *     - overload operator "=" for scalar
 *   - PlusEqual
 *     - overload operator "+=" for 3d vector
 *   - MinusEqual
 *     - overload operator "-=" for 3d vector
 *   - MultiplyEqual
 *     - overload operator "*=" for (3d vector) * scalar
 *   - OverEqual
 *     - overload operator "/=" for (3d vector) / scalar
 *   - Negative
 *     - overload operator "-" to get - Vector3
 *   - Reverse
 *     - same as negative
 *   - Access
 *     - access elements by using "[]"
 *   - ConstAccess
 *     - access elements by using "[]" through pinters
 *     - withough chaning element values
 *   - VectorPlus
 *     - overload operator "+"  for two 3d vectors
 *   - VectorMinus
 *     - overload operator "-"  for two 3d vectors
 *   - Norm2
 *     - get the square of norm of a 3d vector
 *   - Norm
 *     - get the norm of a 3d vector
 *   - Normalize
 *     - normalize a 3d vector
 *   - VmultiplyV
 *     - overload operator "*" to calculate
 *     - the dot product of two 3d vectors
 *   - VdotV
 *     - dot product of two 3d vectors
 *   - VmultiplyNum
 *     - overload operator "*" to calculate
 *     - the product of a 3d vector with a scalar
 *     - of the product of a scalar with a 3d vector
 *   - VoverNum
 *     - overload operator "/" to calculate
 *     -  a 3d vector over a scalar
 *   - OperatorCaret
 *     - overload operator "^" to calculate
 *     - the cross product of two 3d vectors
 *   - VeqV
 *     - overload operator "==" to assert
 *     - the equality between two 3d vectors
 *   - VneV
 *     - overload operator "!=" to assert
 *     - the inequality between two 3d vectors
 *   - StdOutV
 *     - overload operator "<<" to print out
 *     - a 3d vectors on standard output
 *   - PrintV
 *     - print a 3d vectors on standard output
 *     - with formats
 */

class Vector3Test : public testing::Test
{
protected:
	double da = 3.0;
	double db = 4.0;
	double dc = 5.0;
	int    ia = 3;
	int    ib = 4;
	int    ic = 5;
	float  fa = 3.0;
	float  fb = 4.0;
	float  fc = 5.0;
	// for capturing stdout
	std::string output;
};

TEST_F(Vector3Test,Construct)
{
	// double Vector3
	ModuleBase::Vector3<double> u (da,db,dc);
	ModuleBase::Vector3<double> up (u);
	EXPECT_EQ(u.x,3.0);
	EXPECT_EQ(u.y,4.0);
	EXPECT_EQ(u.z,5.0);
	EXPECT_EQ(up.x,3.0);
	EXPECT_EQ(up.y,4.0);
	EXPECT_EQ(up.z,5.0);
	// float Vector3
	ModuleBase::Vector3<float> v (fa,fb,fc);
	ModuleBase::Vector3<float> vp (v);
	EXPECT_EQ(v.x,3.0);
	EXPECT_EQ(v.y,4.0);
	EXPECT_EQ(v.z,5.0);
	EXPECT_EQ(vp.x,3.0);
	EXPECT_EQ(vp.y,4.0);
	EXPECT_EQ(vp.z,5.0);
	// int Vector3
	ModuleBase::Vector3<int> w (ia,ib,ic);
	ModuleBase::Vector3<int> wp (w);
	EXPECT_EQ(w.x,3);
	EXPECT_EQ(w.y,4);
	EXPECT_EQ(w.z,5);
	EXPECT_EQ(wp.x,3);
	EXPECT_EQ(wp.y,4);
	EXPECT_EQ(wp.z,5);
}

TEST_F(Vector3Test,Set)
{
	// double Vector3
	ModuleBase::Vector3<double> u;
	u.set(da,db,dc);
	EXPECT_EQ(u.x,3.0);
	EXPECT_EQ(u.y,4.0);
	EXPECT_EQ(u.z,5.0);
	// float Vector3
	ModuleBase::Vector3<float> v;
	v.set(fa,fb,fc);
	EXPECT_EQ(v.x,3.0);
	EXPECT_EQ(v.y,4.0);
	EXPECT_EQ(v.z,5.0);
	// int Vector3
	ModuleBase::Vector3<int> w;
	w.set(ia,ib,ic);
	EXPECT_EQ(w.x,3);
	EXPECT_EQ(w.y,4);
	EXPECT_EQ(w.z,5);
}

TEST_F(Vector3Test,Equal)
{
	// double Vector3
	ModuleBase::Vector3<double> u, up;
	u.set(da,db,dc);
	up = u;
	EXPECT_EQ(up.x,3.0);
	EXPECT_EQ(up.y,4.0);
	EXPECT_EQ(up.z,5.0);
	// float Vector3
	ModuleBase::Vector3<float> v, vp;
	v.set(fa,fb,fc);
	vp = v;
	EXPECT_EQ(vp.x,3.0);
	EXPECT_EQ(vp.y,4.0);
	EXPECT_EQ(vp.z,5.0);
	// int Vector3
	ModuleBase::Vector3<int> w, wp;
	w.set(ia,ib,ic);
	wp = w;
	EXPECT_EQ(wp.x,3);
	EXPECT_EQ(wp.y,4);
	EXPECT_EQ(wp.z,5);
}

TEST_F(Vector3Test,equal)
{
	// double Vector3
	ModuleBase::Vector3<double> u;
	u.set(da,db,dc);
	u = 2;
	EXPECT_EQ(u.x,2.0);
	EXPECT_EQ(u.y,2.0);
	EXPECT_EQ(u.z,2.0);
	// float Vector3
	ModuleBase::Vector3<float> v;
	v.set(fa,fb,fc);
	v = 2;
	EXPECT_EQ(v.x,2.0);
	EXPECT_EQ(v.y,2.0);
	EXPECT_EQ(v.z,2.0);
	// int Vector3
	ModuleBase::Vector3<int> w;
	w.set(ia,ib,ic);
	w = 2;
	EXPECT_EQ(w.x,2);
	EXPECT_EQ(w.y,2);
	EXPECT_EQ(w.z,2);
}

TEST_F(Vector3Test,PlusEqual)
{
	// double Vector3
	ModuleBase::Vector3<double> u, up;
	u.set(da,db,dc);
	up.set(da,db,dc);
	up += u;
	EXPECT_EQ(up.x,6.0);
	EXPECT_EQ(up.y,8.0);
	EXPECT_EQ(up.z,10.0);
	// float Vector3
	ModuleBase::Vector3<float> v, vp;
	v.set(fa,fb,fc);
	vp.set(fa,fb,fc);
	vp += v;
	EXPECT_EQ(vp.x,6.0);
	EXPECT_EQ(vp.y,8.0);
	EXPECT_EQ(vp.z,10.0);
	// int Vector3
	ModuleBase::Vector3<int> w, wp;
	w.set(ia,ib,ic);
	wp.set(ia,ib,ic);
	wp += w;
	EXPECT_EQ(wp.x,6);
	EXPECT_EQ(wp.y,8);
	EXPECT_EQ(wp.z,10);
}

TEST_F(Vector3Test,MinusEqual)
{
	// double Vector3
	ModuleBase::Vector3<double> u, up;
	u.set(da,db,dc);
	up.set(3*da,3*db,3*dc);
	up -= u;
	EXPECT_EQ(up.x,6.0);
	EXPECT_EQ(up.y,8.0);
	EXPECT_EQ(up.z,10.0);
	// float Vector3
	ModuleBase::Vector3<float> v, vp;
	v.set(fa,fb,fc);
	vp.set(3*fa,3*fb,3*fc);
	vp -= v;
	EXPECT_EQ(vp.x,6.0);
	EXPECT_EQ(vp.y,8.0);
	EXPECT_EQ(vp.z,10.0);
	// int Vector3
	ModuleBase::Vector3<int> w, wp;
	w.set(ia,ib,ic);
	wp.set(3*ia,3*ib,3*ic);
	wp -= w;
	EXPECT_EQ(wp.x,6);
	EXPECT_EQ(wp.y,8);
	EXPECT_EQ(wp.z,10);
}

TEST_F(Vector3Test,MultiplyEqual)
{
	// double Vector3
	ModuleBase::Vector3<double> u;
	u.set(da,db,dc);
	u *= 2;
	EXPECT_EQ(u.x,6.0);
	EXPECT_EQ(u.y,8.0);
	EXPECT_EQ(u.z,10.0);
	// float Vector3
	ModuleBase::Vector3<float> v;
	v.set(fa,fb,fc);
	v *= 2;
	EXPECT_EQ(v.x,6.0);
	EXPECT_EQ(v.y,8.0);
	EXPECT_EQ(v.z,10.0);
	// int Vector3
	ModuleBase::Vector3<int> w;
	w.set(ia,ib,ic);
	w *= 2;
	EXPECT_EQ(w.x,6);
	EXPECT_EQ(w.y,8);
	EXPECT_EQ(w.z,10);
}

TEST_F(Vector3Test,OverEqual)
{
	// double Vector3
	ModuleBase::Vector3<double> u;
	u.set(4*da,4*db,4*dc);
	u /= 2;
	EXPECT_EQ(u.x,6.0);
	EXPECT_EQ(u.y,8.0);
	EXPECT_EQ(u.z,10.0);
	// float Vector3
	ModuleBase::Vector3<float> v;
	v.set(4*fa,4*fb,4*fc);
	v /= 2;
	EXPECT_EQ(v.x,6.0);
	EXPECT_EQ(v.y,8.0);
	EXPECT_EQ(v.z,10.0);
	// int Vector3
	ModuleBase::Vector3<int> w;
	w.set(4*ia,4*ib,4*ic);
	w /= 2;
	EXPECT_EQ(w.x,6);
	EXPECT_EQ(w.y,8);
	EXPECT_EQ(w.z,10);
}

TEST_F(Vector3Test,Negative)
{
	// double Vector3
	ModuleBase::Vector3<double> u, up;
	u.set(da,db,dc);
	up = -u;
	EXPECT_EQ(up.x,-3.0);
	EXPECT_EQ(up.y,-4.0);
	EXPECT_EQ(up.z,-5.0);
	// float Vector3
	ModuleBase::Vector3<float> v, vp;
	v.set(fa,fb,fc);
	vp = -v;
	EXPECT_EQ(vp.x,-3.0);
	EXPECT_EQ(vp.y,-4.0);
	EXPECT_EQ(vp.z,-5.0);
	// int Vector3
	ModuleBase::Vector3<int> w, wp;
	w.set(ia,ib,ic);
	wp = -w;
	EXPECT_EQ(wp.x,-3);
	EXPECT_EQ(wp.y,-4);
	EXPECT_EQ(wp.z,-5);
}

TEST_F(Vector3Test,Access)
{
	// double Vector3
	ModuleBase::Vector3<double> u;
	u.set(da,db,dc);
	EXPECT_EQ(u[0],3.0);
	EXPECT_EQ(u[1],4.0);
	EXPECT_EQ(u[2],5.0);
	// float Vector3
	ModuleBase::Vector3<float> v;
	v.set(fa,fb,fc);
	EXPECT_EQ(v.x,3.0);
	EXPECT_EQ(v.y,4.0);
	EXPECT_EQ(v.z,5.0);
	// int Vector3
	ModuleBase::Vector3<int> w;
	w.set(ia,ib,ic);
	EXPECT_EQ(w.x,3);
	EXPECT_EQ(w.y,4);
	EXPECT_EQ(w.z,5);
}


TEST_F(Vector3Test,ConstAccess)
{
	// double Vector3
	ModuleBase::Vector3<double> u;
	u.set(da,db,dc);
	const ModuleBase::Vector3<double> *up(&u);
	EXPECT_EQ((*up)[0],3.0);
	EXPECT_EQ((*up)[1],4.0);
	EXPECT_EQ((*up)[2],5.0);
	// float Vector3
	ModuleBase::Vector3<float> v;
	const ModuleBase::Vector3<float> *vp(&v);
	v.set(fa,fb,fc);
	EXPECT_EQ((*vp).x,3.0);
	EXPECT_EQ((*vp).y,4.0);
	EXPECT_EQ((*vp).z,5.0);
	// int Vector3
	//ModuleBase::Vector3<int> w;
	//w.set(ia,ib,ic);
	//EXPECT_EQ(w.x,3);
	//EXPECT_EQ(w.y,4);
	//EXPECT_EQ(w.z,5);
}


TEST_F(Vector3Test,Reverse)
{
	// double Vector3
	ModuleBase::Vector3<double> u;
	u.set(da,db,dc);
	u.reverse();
	EXPECT_EQ(u.x,-3.0);
	EXPECT_EQ(u.y,-4.0);
	EXPECT_EQ(u.z,-5.0);
	// float Vector3
	ModuleBase::Vector3<float> v;
	v.set(fa,fb,fc);
	v.reverse();
	EXPECT_EQ(v.x,-3.0);
	EXPECT_EQ(v.y,-4.0);
	EXPECT_EQ(v.z,-5.0);
	// int Vector3
	ModuleBase::Vector3<int> w;
	w.set(ia,ib,ic);
	w.reverse();
	EXPECT_EQ(w.x,-3);
	EXPECT_EQ(w.y,-4);
	EXPECT_EQ(w.z,-5);
}

TEST_F(Vector3Test,VectorPlus)
{
	// double Vector3
	ModuleBase::Vector3<double> u,up,upp;
	u.set(da,db,dc);
	up.set(da,db,dc);
	upp = u + up;
	EXPECT_EQ(upp[0],6.0);
	EXPECT_EQ(upp[1],8.0);
	EXPECT_EQ(upp[2],10.0);
	// float Vector3
	ModuleBase::Vector3<float> v,vp,vpp;
	v.set(fa,fb,fc);
	vp.set(fa,fb,fc);
	vpp = v + vp;
	EXPECT_EQ(vpp.x,6.0);
	EXPECT_EQ(vpp.y,8.0);
	EXPECT_EQ(vpp.z,10.0);
	// int Vector3
	ModuleBase::Vector3<int> w,wp,wpp;
	w.set(ia,ib,ic);
	wp.set(ia,ib,ic);
	wpp = w + wp;
	EXPECT_EQ(wpp.x,6);
	EXPECT_EQ(wpp.y,8);
	EXPECT_EQ(wpp.z,10);
}

TEST_F(Vector3Test,VectorMinus)
{
	// double Vector3
	ModuleBase::Vector3<double> u,up,upp;
	u.set(da,db,dc);
	up.set(2*da,2*db,2*dc);
	upp = u - up;
	EXPECT_EQ(upp[0],-3.0);
	EXPECT_EQ(upp[1],-4.0);
	EXPECT_EQ(upp[2],-5.0);
	// float Vector3
	ModuleBase::Vector3<float> v,vp,vpp;
	v.set(fa,fb,fc);
	vp.set(3*fa,3*fb,3*fc);
	vpp = v - vp;
	EXPECT_EQ(vpp.x,-6.0);
	EXPECT_EQ(vpp.y,-8.0);
	EXPECT_EQ(vpp.z,-10.0);
	// int Vector3
	ModuleBase::Vector3<int> w,wp,wpp;
	w.set(3*ia,3*ib,3*ic);
	wp.set(ia,ib,ic);
	wpp = w - wp;
	EXPECT_EQ(wpp.x,6);
	EXPECT_EQ(wpp.y,8);
	EXPECT_EQ(wpp.z,10);
}

TEST_F(Vector3Test,Norm2)
{
	// double Vector3
	ModuleBase::Vector3<double> u;
	u.set(da,db,dc);
	EXPECT_EQ(u.norm2(),50.0);
	// float Vector3
	ModuleBase::Vector3<float> v;
	v.set(fa,fb,fc);
	EXPECT_EQ(v.norm2(),50.0);
	// int Vector3
	ModuleBase::Vector3<int> w;
	w.set(ia,ib,ic);
	EXPECT_EQ(w.norm2(),50);
}


TEST_F(Vector3Test,Norm)
{
	// double Vector3
	ModuleBase::Vector3<double> u;
	u.set(da,db,dc);
	double nm = u.norm();
	double nm2= sqrt(50.0);
	EXPECT_DOUBLE_EQ(nm,nm2);
	EXPECT_FLOAT_EQ(nm,sqrt(50.0));
	// float Vector3
	ModuleBase::Vector3<float> v;
	v.set(fa,fb,fc);
	float nmp = v.norm();
	float nmp2= sqrt(50.0);
	EXPECT_FLOAT_EQ(nmp,sqrt(50.0));
}


TEST_F(Vector3Test,Normalize)
{
	// double Vector3
	ModuleBase::Vector3<double> u;
	u.set(da,db,dc);
	u.normalize();
	EXPECT_DOUBLE_EQ(u.norm(),1.0);
	// float Vector3
	ModuleBase::Vector3<float> v;
	v.set(fa,fb,fc);
	v.normalize();
	EXPECT_FLOAT_EQ(v.norm(),1.0);
}

TEST_F(Vector3Test,VmultiplyV)
{
	// double Vector3
	ModuleBase::Vector3<double> u,up;
	u.set(da,db,dc);
	up.set(da,db,dc);
	double mpd = u * up;
	EXPECT_EQ(mpd,50.0);
	// float Vector3
	ModuleBase::Vector3<float> v,vp;
	v.set(fa,fb,fc);
	vp.set(fa,fb,fc);
	float mpf = v*vp;
	EXPECT_EQ(mpf,50.0);
	// int Vector3
	ModuleBase::Vector3<int> w,wp;
	w.set(ia,ib,ic);
	wp.set(ia,ib,ic);
	int mpi = w*wp;
	EXPECT_EQ(mpf,50);
}

TEST_F(Vector3Test,VdotV)
{
	// double Vector3
	ModuleBase::Vector3<double> u,up;
	u.set(da,db,dc);
	up.set(da,db,dc);
	double mpd = dot(u,up);
	EXPECT_EQ(mpd,50.0);
	// float Vector3
	ModuleBase::Vector3<float> v,vp;
	v.set(fa,fb,fc);
	vp.set(fa,fb,fc);
	float mpf = dot(v,vp);
	EXPECT_EQ(mpf,50.0);
	// int Vector3
	ModuleBase::Vector3<int> w,wp;
	w.set(ia,ib,ic);
	wp.set(ia,ib,ic);
	int mpi = dot(w,wp);
	EXPECT_EQ(mpf,50);
}

TEST_F(Vector3Test,VmultiplyNum)
{
	// double Vector3
	ModuleBase::Vector3<double> u,up,upp;
	u.set(da,db,dc);
	double s = 3.0;
	up = s*u; upp = u*s;
	EXPECT_EQ(upp[0],up[0]);
	EXPECT_EQ(upp[1],up[1]);
	EXPECT_EQ(upp[2],up[2]);
	EXPECT_EQ(upp[0],9.0);
	EXPECT_EQ(upp[1],12.0);
	EXPECT_EQ(upp[2],15.0);
	// float Vector3
	ModuleBase::Vector3<float> v,vp,vpp;
	v.set(fa,fb,fc);
	float t = 3.0;
	vp = t*v; vpp = v*t;
	EXPECT_EQ(vpp[0],vp[0]);
	EXPECT_EQ(vpp[1],vp[1]);
	EXPECT_EQ(vpp[2],vp[2]);
	EXPECT_EQ(vpp[0],9.0);
	EXPECT_EQ(vpp[1],12.0);
	EXPECT_EQ(vpp[2],15.0);
	// int Vector3
	ModuleBase::Vector3<int> w,wp,wpp;
	w.set(ia,ib,ic);
	int q = 3;
	wp = q*w; wpp = w*q;
	EXPECT_EQ(wpp[0],wp[0]);
	EXPECT_EQ(wpp[1],wp[1]);
	EXPECT_EQ(wpp[2],wp[2]);
	EXPECT_EQ(wpp[0],9.0);
	EXPECT_EQ(wpp[1],12.0);
	EXPECT_EQ(wpp[2],15.0);
}

TEST_F(Vector3Test,VoverNum)
{
	// double Vector3
	ModuleBase::Vector3<double> u,up;
	u.set(2*da,2*db,2*dc);
	double s = 2.0;
	up = u/s;
	EXPECT_EQ(up.x,3.0);
	EXPECT_EQ(up.y,4.0);
	EXPECT_EQ(up.z,5.0);
	// float Vector3
	ModuleBase::Vector3<float> v,vp;
	v.set(2*fa,2*fb,2*fc);
	float t = 2.0;
	vp = v/t;
	EXPECT_EQ(vp.x,3.0);
	EXPECT_EQ(vp.y,4.0);
	EXPECT_EQ(vp.z,5.0);
	// int Vector3
	ModuleBase::Vector3<int> w,wp;
	w.set(2*ia,2*ib,2*ic);
	int q = 2;
	wp = w/q;
	EXPECT_EQ(wp.x,3);
	EXPECT_EQ(wp.y,4);
	EXPECT_EQ(wp.z,5);
}

TEST_F(Vector3Test,OperatorCaret)
{
	// double Vector3
	ModuleBase::Vector3<double> u,up,upp;
	u.set(da,db,dc);
	up.set(da,db,dc);
	upp = u^up;
	EXPECT_EQ(upp.x,u.y*up.z - u.z*up.y);
	EXPECT_EQ(upp.y,u.z*up.x - u.x*up.z);
	EXPECT_EQ(upp.z,u.x*up.y - u.y*up.x);
	// float Vector3
	ModuleBase::Vector3<float> v,vp,vpp;
	v.set(2*fa,2*fb,2*fc);
	vp.set(fa,fb,fc);
	vpp = v^vp;
	EXPECT_EQ(vpp.x,v.y*vp.z - v.z*vp.y);
	EXPECT_EQ(vpp.y,v.z*vp.x - v.x*vp.z);
	EXPECT_EQ(vpp.z,v.x*vp.y - v.y*vp.x);
	// int Vector3
	ModuleBase::Vector3<int> w,wp,wpp;
	w.set(2*ia,2*ib,2*ic);
	wp.set(ia,ib,ic);
	wpp = w^wp;
	EXPECT_EQ(wpp.x,w.y*wp.z - w.z*wp.y);
	EXPECT_EQ(wpp.y,w.z*wp.x - w.x*wp.z);
	EXPECT_EQ(wpp.z,w.x*wp.y - w.y*wp.x);
}

TEST_F(Vector3Test,Cross)
{
	// double Vector3
	ModuleBase::Vector3<double> u,up,upp;
	u.set(da,db,dc);
	up.set(da,db,dc);
	upp = cross(u,up);
	EXPECT_EQ(upp.x,u.y*up.z - u.z*up.y);
	EXPECT_EQ(upp.y,u.z*up.x - u.x*up.z);
	EXPECT_EQ(upp.z,u.x*up.y - u.y*up.x);
	// float Vector3
	ModuleBase::Vector3<float> v,vp,vpp;
	v.set(2*fa,2*fb,2*fc);
	vp.set(fa,fb,fc);
	vpp = cross(v,vp);
	EXPECT_EQ(vpp.x,v.y*vp.z - v.z*vp.y);
	EXPECT_EQ(vpp.y,v.z*vp.x - v.x*vp.z);
	EXPECT_EQ(vpp.z,v.x*vp.y - v.y*vp.x);
	// int Vector3
	ModuleBase::Vector3<int> w,wp,wpp;
	w.set(2*ia,2*ib,2*ic);
	wp.set(ia,ib,ic);
	wpp = cross(w,wp);
	EXPECT_EQ(wpp.x,w.y*wp.z - w.z*wp.y);
	EXPECT_EQ(wpp.y,w.z*wp.x - w.x*wp.z);
	EXPECT_EQ(wpp.z,w.x*wp.y - w.y*wp.x);
}

TEST_F(Vector3Test,VeqV)
{
	// double Vector3
	ModuleBase::Vector3<double> u,up;
	u.set(da,db,dc);
	up.set(da,db,dc);
	EXPECT_TRUE(up == u);
	// float Vector3
	ModuleBase::Vector3<float> v,vp;
	v.set(fa,fb,fc);
	vp.set(fa,fb,fc);
	EXPECT_TRUE(vp == v);
	// int Vector3
	ModuleBase::Vector3<int> w,wp;
	w.set(ia,ib,ic);
	wp.set(ia,ib,ic);
	EXPECT_TRUE(wp == w);
}

TEST_F(Vector3Test,VneV)
{
	// double Vector3
	ModuleBase::Vector3<double> u,up;
	u.set(da,db,dc);
	up.set(da,db,2*dc);
	EXPECT_TRUE(up != u);
	// float Vector3
	ModuleBase::Vector3<float> v,vp;
	v.set(fa,fb,2*fc);
	vp.set(fa,fb,fc);
	EXPECT_TRUE(vp != v);
	// int Vector3
	ModuleBase::Vector3<int> w,wp;
	w.set(ia,ib,2*ic);
	wp.set(ia,ib,ic);
	EXPECT_TRUE(wp != w);
}

TEST_F(Vector3Test,StdOutV)
{
	// double Vector3
	ModuleBase::Vector3<double> u(da,db,dc);
	testing::internal::CaptureStdout();
	std::cout << u << std::endl;
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("("));
	// float Vector3
	ModuleBase::Vector3<float> v(fa,fb,fc);
	testing::internal::CaptureStdout();
	std::cout << v << std::endl;
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr(","));
	// int Vector3
	ModuleBase::Vector3<int> w(ia,ib,ic);
	testing::internal::CaptureStdout();
	std::cout << w << std::endl;
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr(")"));
}

TEST_F(Vector3Test,PrintV)
{
	// double Vector3
	ModuleBase::Vector3<double> u(3.1415926,db,dc);
	testing::internal::CaptureStdout();
	u.print();
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("3.1416"));
	// float Vector3
	ModuleBase::Vector3<float> v(fa,fb,3.14);
	testing::internal::CaptureStdout();
	v.print();
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("3.14"));
	// int Vector3
	ModuleBase::Vector3<int> w(ia,101,ic);
	testing::internal::CaptureStdout();
	w.print();
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("101"));
}

