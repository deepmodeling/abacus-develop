//==========================================================
// AUTHOR : Zhengpan , mohan , spshu
// DATE : 2007-9
//==========================================================
#ifndef SYMMETRY_BASIC_H
#define SYMMETRY_BASIC_H
#include "symm_other.h"
#include "module_base/mymath.h"
#include "module_base/ylm.h"
#include "module_base/matrix3.h"
namespace ModuleSymmetry
{
class Symmetry_Basic
{
	public:

	 Symmetry_Basic();
	~Symmetry_Basic();

    double epsilon;
    double epsilon_input;   ///< the input value of symmetry_prec, should not  be changed

	// control accuray
	bool equal(const double &m, const double &n)const;
	void check_boundary(double &x)const;
	double get_translation_vector(const double &x1, const double &x2);
	void check_translation(double &x, const double &t) const;
	double check_diff(const double &x1, const double &x2);
	
	void veccon(
			double *va,
			double *vb,
			const int num,
			const ModuleBase::Vector3<double> &aa1,
			const ModuleBase::Vector3<double> &aa2,
			const ModuleBase::Vector3<double> &aa3,
			const ModuleBase::Vector3<double> &bb1,
			const ModuleBase::Vector3<double> &bb2,
			const ModuleBase::Vector3<double> &bb3
			);
	void matrigen(ModuleBase::Matrix3 *symgen, const int ngen, ModuleBase::Matrix3* symop, int &nop) const;
	void setgroup(ModuleBase::Matrix3 *symop, int &nop, const int &ibrav) const;
	void rotate(
			ModuleBase::Matrix3 &gmatrix, ModuleBase::Vector3<double> &gtrans, 
			int i, int j, int k, const int, const int, const int, int&, int&, int&);
	void test_atom_ordering(double *posi, const int natom, int *subindex) const;

	protected:

	std::string get_brav_name(const int ibrav) const;
	void pointgroup(const int &nrot,int &pgnumber,std::string &pgname, const ModuleBase::Matrix3* gmatrix, std::ofstream &ofs_running);
	void atom_ordering(double *posi, const int natom, int *subindex);
	void atom_ordering_new(double *posi, const int natom, int *subindex) const;

	private:

	void order_atoms(double* pos, const int &nat, const int *index) const;
	void order_y(double *pos, const int &oldpos, const int &newpos);
	void order_z(double *pos, const int &oldpos, const int &newpos);
};

//for test only
extern bool test_brav;

}//end of define namespace

#endif
