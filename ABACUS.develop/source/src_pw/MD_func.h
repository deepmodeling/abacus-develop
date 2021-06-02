#ifndef MD_FUNC_H
#define MD_FUNC_H

#include "src_pw/unitcell_pseudo.h"


class MD_func
{
    public:
    MD_func(){};
    ~MD_func(){};
	bool RestartMD(const int& numIon, Vector3<double>* vel, int& step_rst);
    void mdRestartOut(const int& step, const int& recordFreq, const int& numIon, Vector3<double>* vel);
	double GetAtomKE(const int& numIon, const Vector3<double>* vel, const double* allmass);
	void InitVelocity(
		const int& numIon, 
		const double& temperature, 
		const double& fundamentalTime, 
		const double* allmass,
		Vector3<double>* vel);
//	void ReadNewTemp(int step);
	string intTurnTostring(long int iter,string path);
	int getMassMbl(const UnitCell_pseudo &unit_in, double* allmass, Vector3<int>* ionmbl);
	void callInteraction_LCAO(const int& numIon, Vector3<double>* force, matrix& stress_lcao);
	void callInteraction_PW(const int& numIon, Vector3<double>* force, matrix& stress_pw);
    void printpos(const string& file, const int& iter, const int& recordFreq, const UnitCell_pseudo& unit_in);
    void scalevel(
		const int& numIon,
		const int& nfrozen,
		const double& temperature,
		Vector3<double>* vel,
		const double* allmass);
	double MAXVALF(const int numIon, const Vector3<double>* force);
	double Conserved(const double KE, const double PE, const int number);
};
#endif
