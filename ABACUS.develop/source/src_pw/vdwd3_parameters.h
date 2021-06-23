//==========================================================
// AUTHOR:	Yuyang Ji
// DATE:	2021-04-19
//==========================================================
#ifndef VDWD3_PARAMETERS_H
#define VDWD3_PARAMETERS_H

#include <vector>
#include <string>
#include "src_global/vector3.h"
#include "module_base/matrix3.h"
#include "module_cell/unitcell_pseudo.h"


class Vdwd3_Parameters
{
public:
	Vdwd3_Parameters();

    void init_C6();
    void init_r2r4();
    void init_rcov();
    void init_r0ab();

    bool flag_vdwd3;    
    bool abc; //third-order term?
	double rthr2; //R^2 distance neglect threshold (important for speed in case of large systems) (a.u.)
	double cn_thr2; //R^2 distance to cutoff for CN_calculation (a.u.)
    double s6;
	double rs6;
	double s18;
	double rs18;
    string version;
    string model;
	Vector3<int> period;

    // internal used
    int max_elem = 94;
    vector<int> mxc;
    vector<vector<vector<vector<vector<double>>>>> c6ab;
    vector<double> r2r4;
    vector<double> rcov;
    vector<vector<double>> r0ab;
    double k1 = 16.0, k2 = 4.0/3.0, k3 = -4.0;
    double alp6 = 14.0, alp8 = alp6+2, alp10 = alp8+2;

private:
    int limit(int &i);
};

#endif // define VDWD3_PARAMETERS_H 
