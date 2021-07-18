//==========================================================
// AUTHOR : Zhengpan , mohan , spshu
// DATE : 2007-9
//==========================================================
#ifndef SYMMETRY_BASIC_H
#define SYMMETRY_BASIC_H
//#include "../src_pw/tools.h"
//#include "symmetry_basic.h"
#include "symm_other.h"
#include "../module_base/mymath.h"
#include "../module_base/ylm.h"
#include "../module_base/matrix3.h"
class Symmetry_Basic
{
	public:

	 Symmetry_Basic();
	~Symmetry_Basic();

	double epsilon;

	void maxmin(const double, const double, const double, double &, double &, double &);
	void maxmin(const double, const double, double &, double &);
	
	// control accuray
	bool equal(const double &m, const double &n)const;
	void check_boundary(double &x)const;
	double get_translation_vector(const double &x1, const double &x2);
	void check_translation(double &x, const double &t);
	double check_diff(const double &x1, const double &x2);
	void change_value(double &x1, double &x2);
	void change_value(int &x1, int &x2);
	
	void shortest_vector(Vector3<double> &t1,Vector3<double> &t2,Vector3<double> &t3);
	void shorter_vector(Vector3<double> &t1,Vector3<double> &t2,double &abs0,double &abs1,bool &flag1);
	void recip(
			const double a,
			const Vector3<double> &a1,
			const Vector3<double> &a2,
			const Vector3<double> &a3,
			Vector3<double> &b1,
			Vector3<double> &b2,
			Vector3<double> &b3
			);
	void veccon(
			double *va,
			double *vb,
			const int num,
			const Vector3<double> &aa1,
			const Vector3<double> &aa2,
			const Vector3<double> &aa3,
			const Vector3<double> &bb1,
			const Vector3<double> &bb2,
			const Vector3<double> &bb3
			);
	void matrigen(Matrix3 *symgen, const int ngen, Matrix3* symop, int &nop);
	void setgroup(Matrix3 *symop, int &nop, const int &ibrav);
	void rotate( Matrix3 &gmatrix, Vector3<double> &gtrans, int i, int j, int k, const int, const int, const int, int&, int&, int&);

	protected:

	string get_brav_name(const int ibrav);
	void pointgroup(const int &nrot,int &pgnumber,string &pgname, const Matrix3* gmatrix);
	void atom_ordering(double *posi, const int natom, int *subindex);

	private:

	void order_atoms(double* pos, const int &nat, const int *index);
	void order_y(double *pos, const int &oldpos, const int &newpos);
	void order_z(double *pos, const int &oldpos, const int &newpos);
};

#endif
