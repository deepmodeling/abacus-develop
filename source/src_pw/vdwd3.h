//==========================================================
// AUTHOR : Yuyang Ji
// DATE : 2019-04-22
// UPDATE : 2021-4-19
//==========================================================
#ifndef VDWD3_H
#define VDWD3_H
#include "vdwd3_parameters.h"
#include "../module_cell/unitcell_pseudo.h"
#include "../module_base/vector3.h"
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
	const std::vector<ModuleBase::Vector3<double>> 	   &get_force() const{ return force;  }
	const ModuleBase::Matrix3                      &get_stress()const{ return stress; }

private:

    const UnitCell_pseudo &ucell;
	Vdwd3_Parameters &para;

	double energy = 0;
	std::vector<ModuleBase::Vector3<double>> force;
	ModuleBase::Matrix3 stress;

	std::vector<ModuleBase::Vector3<double>> lat;
    std::vector<int> iz;
	std::vector<ModuleBase::Vector3<double>> xyz;
	std::vector<int> rep_vdw;
    std::vector<int> rep_cn;

	void init(const UnitCell_pseudo &ucell);	

	void set_criteria(
		double &rthr, 
		std::vector<ModuleBase::Vector3<double>> &lat, 
		std::vector<double> &tau_max);

    std::vector<double> atomkind(const UnitCell_pseudo &ucell);	

	void getc6(int &iat, int &jat, double &nci, double &ncj, double &c6);

	void pbcncoord(std::vector<double> &cn);

	void pbcthreebody(
		std::vector<int> &iz,  
		std::vector<ModuleBase::Vector3<double>> &lat, 
		std::vector<ModuleBase::Vector3<double>> &xyz, 
		std::vector<int> &rep_cn, 
		std::vector<double> &cc6ab, 
		double &eabc);

	void pbcgdisp(
		std::vector<ModuleBase::Vector3<double>> &g, 
		ModuleBase::matrix &smearing_sigma);

	void get_dc6_dcnij(
		int &mxci, 
		int &mxcj, 
		double &cni, 
		double &cnj, 
		int &izi, 
		int &izj, 
		int &iat, 
		int &jat, 
		double &c6check, 
		double &dc6i, 
		double &dc6j);

	int lin(int &i1, int &i2) {
		int idum1 = std::max(i1+1, i2+1); 
		int idum2 = std::min(i1+1, i2+1); 
		int res = idum2+idum1*(idum1-1)/2-1; 
		return res;}
};

#endif
