#ifndef EPSILON0_PWSCF_H
#define EPSILON0_PWSCF_H

#include "../src_pw/wavefunc.h"
#include "../src_parallel/parallel_global.h"

class Epsilon0_pwscf
{

public:
        
		 Epsilon0_pwscf();
		 ~Epsilon0_pwscf();
		 
		 bool epsilon;
		 double intersmear;
		 double intrasmear;
		 double domega;
		 int nomega;
		 double shift;
		 bool metalcalc;
		 double smearing_sigma;
		 
		 void Init();
		 void Cal_epsilon0();
		 void Delete();
		 double wgrid(int iw);
		 double focc(int ib, int ik);
		 void Cal_dipole(int ik);
		 
		 double **epsr;
		 double **epsi;
		 
private:
		 
		 bool init_finish;
		 std::complex<double> ***dipole_aux;
		 std::complex<double> ***dipole;
		 
};

namespace GlobalC
{
extern Epsilon0_pwscf epsilon0_pwscf;
}

#endif		 
