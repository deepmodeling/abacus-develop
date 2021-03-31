#ifndef MD_PARAMETERS_H
#define MD_PARAMETERS_H

class MD_parameters
{
    MD_parameters(){};
    ~MD_parameters(){};

    private:
    int rstMD;          // 1 : restart MD, vel.restart and ion.restart files will be read
                        // 0 : not restart from ion and vel files
	int mdtype;         // 1: NVT:Nose Hoover, 2:NVT:velocity scaling 3: NPT, 0: NVE 
    double Qmass;       // Inertia of extended system variable
    double dt ;         // Time increment (hbar/E_hartree)
    double temperature; //Temperature (in Hartree, 1 Hartree ~ 3E5 K)
    int ntype;          //
    int dumpmdfreq;     // The period to dump MD information for monitoring and restarting MD
    int NVT_control;
    double NVT_tau;
    int MNHC;
    string mdoutputpath;// output directory of md files: .ion .vel
	double ediff;       //parameter for constrain 
	double ediffg;      //parameter for constrain

};

#endif