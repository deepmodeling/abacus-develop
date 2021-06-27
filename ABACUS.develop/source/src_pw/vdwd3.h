//==========================================================
// AUTHOR : Yuyang Ji
// DATE : 2019-04-22
// UPDATE : 2021-4-19
//==========================================================
#ifndef VDWD3_H
#define VDWD3_H
#include "vdwd3_parameters.h"
#include "module_cell/unitcell_pseudo.h"
#include "module_base/vector3.h"
#include <vector>
#include <cmath>

class Vdwd3
{
public:

	Vdwd3(const UnitCell_pseudo &unit_in, Vdwd3_Parameters &para_in);
	
	void cal_energy();	
	void cal_force();
	void cal_stress();

	const double                       &get_energy()const{ return energy; }
	const vector<Vector3<double>> 	   &get_force() const{ return force;  }
	const Matrix3                      &get_stress()const{ return stress; }

private:

    const UnitCell_pseudo &ucell;
	Vdwd3_Parameters &para;

	double energy = 0;
	vector<Vector3<double>> force;
	Matrix3 stress;

	vector<Vector3<double>> lat;
    vector<int> iz;
	vector<Vector3<double>> xyz;
	vector<int> rep_vdw;
    vector<int> rep_cn;

	void init(const UnitCell_pseudo &ucell);	
	void set_criteria(double &rthr, vector<Vector3<double>> &lat, vector<double> &tau_max);
    vector<double> atomkind(const UnitCell_pseudo &ucell);	
	void getc6(int &iat, int &jat, double &nci, double &ncj, double &c6);
	void pbcncoord(vector<double> &cn);
	void pbcthreebody(vector<int> &iz,  vector<Vector3<double>> &lat, vector<Vector3<double>> &xyz, vector<int> &rep_cn, vector<double> &cc6ab, double &eabc);
	void pbcgdisp(vector<Vector3<double>> &g, matrix &sigma);
	void get_dc6_dcnij(int &mxci, int &mxcj, double &cni, double &cnj, int &izi, int &izj, int &iat, int &jat, double &c6check, double &dc6i, double &dc6j);
	int lin(int &i1, int &i2) {int idum1 = max(i1+1, i2+1); int idum2 = min(i1+1, i2+1); int res = idum2+idum1*(idum1-1)/2-1; return res;}
};

#endif
