#ifndef MD_PARAMETERS_H
#define MD_PARAMETERS_H

class MD_parameters
{
public:
    MD_parameters()
	{
		mdtype=1;
		NVT_tau=0;
		NVT_control=1;
		dt=1;
		Qmass=1;
		tfirst=0;         //kelvin
		tlast=-1;
		recordFreq=1;
		mdoutputpath="mdoutput";
		rstMD=0;
		fixTemperature=1;
		ediff=1e-4;
		ediffg=1e-3;
		MNHC=4;
	};
    ~MD_parameters(){};

    int rstMD;          // 1 : restart MD, vel.restart and ion.restart files will be read
                        // 0 : not restart from ion and vel files
	int mdtype;         // 1: NVT:Nose Hoover, 2:NVT:velocity scaling 3: NPT, 0: NVE 
    double Qmass;       // Inertia of extended system variable
    double dt ;         // Time increment (hbar/E_hartree)
    double tfirst; //Temperature (in Hartree, 1 Hartree ~ 3E5 K)
    double tlast;
    int recordFreq;     // The period to dump MD information for monitoring and restarting MD
    int NVT_control;
    double NVT_tau;
    int MNHC;
    string mdoutputpath;// output directory of md files: .ion .vel
	double ediff;       //parameter for constrain 
	double ediffg;      //parameter for constrain
    int fixTemperature;

};


#endif