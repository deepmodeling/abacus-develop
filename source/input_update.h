//==========================================================
// Author: mohan
// DATE : 2011-06-06
//==========================================================
#ifndef UPDATE_INPUT_H
#define UPDATE_INPUT_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <stdio.h>
#include <string.h>

using namespace std;

class Update_input
{
public:
    
	Update_input();
    ~Update_input();
    void init(const std::string &fn);

private:
	
    int cal_force;
    double force_thr;		// threshold of force in unit (Ry/Bohr)
	double force_thr_ev2;	// invalid force threshold, mohan add 2011-04-17
    double scf_thr_rho;				// \sum |rhog_out - rhog_in |^2
    int scf_nmax;				// number of max elec iter
    int relax_nmax;				// number of max ionic iter
    int md_nstep;               // number of md steps
    double mixing_beta;		// 0 : no_mixing
	int printe;
	//int chg_extrap;	// mohan add 2011-03-13, xiaohui modify 2015-02-01
	std::string chg_extrap;	//xiaohui add 2015-02-01
    int out_chg;		// output charge density.
	int out_dm; // output density matrix.
	int out_dos;			// dos calculation. mohan add 20090909
	int	out_wfc_lcao;			// output the wave functions in local basis.
    
	bool Read(const std::string &fn);
#ifdef __MPI
    void Bcast(void);
#endif

	template <class T>
	static void change(std::ofstream &ofs, const std::string &name, T &v1, T &v2)
	{
		ofs << " " << name << " change from " << v1 << " to " << v2 << std::endl;
	}

    template <class T>
    static void read_value(std::ifstream &ifs, T &var)
    {
        ifs >> var;
        ifs.ignore(150, '\n');
        return;
    }
};

#endif //UPDATE_INPUT
