#ifndef MD_RUN_H
#define MD_RUN_H

#include "MD_func.h"
#include "MD_thermo.h"
#include "MD_parameters.h"

//using namespace std;
class MD_run
{
	public:	
	MD_run(){};
	~MD_run(){};
	void runnvt(int step1);
	void runNVE(int step1); //NVE ensemble MD
	bool runFIRE(int step1); //relax method FIRE

	private:
    MD_func mdf;
    MD_thermo mdt;
	MD_parameters &mdp;
};

#endif