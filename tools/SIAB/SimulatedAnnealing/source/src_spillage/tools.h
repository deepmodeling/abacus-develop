#ifndef TOOL_H
#define TOOL_H

#include "common.h"
#include "../src_pw/pw_basis.h"
#include "MultiZeta.h"
#include "read_INPUT.h"
#include "../src_parallel/parallel_kpoints.h" // add 2010-06-12
//plane wave , mohan add 2010-06-14
extern PW_Basis PW;
extern MultiZeta mz;
extern Read_INPUT input;
extern Parallel_Kpoints Pkpoints;
extern ofstream ofs_running;

// parallel information
extern int NPROC; // mohan add 2010-06-12
extern int KPAR;
extern int MY_RANK;
extern int MY_POOL;
extern int NPROC_IN_POOL;
extern int RANK_IN_POOL;

// global + input information.
extern int NKSTOT;
extern int NKS;
extern int NBANDS;
extern int NWFCALL;
extern int NE;

// global information
extern int STRNUM;
extern int NTYPE;
extern string *LABEL;
extern int *NA;
extern double **CARPOSX;
extern double **CARPOSY;
extern double **CARPOSZ;
extern Matrix3 LATVEC;

// plane wave information
extern bool USEPW;
extern double LAT0;

// k points
extern double *CARKX;
extern double *CARKY;
extern double *CARKZ;

//orbital + energy cutoff
extern double ECUT;
extern double ECUT_JLQ;
extern double RCUT;
extern bool SMOOTH; // use to make the second derivative of psi is continues.
extern double SIGMA;
extern double TOLERENCE;
extern int LMAXALL;
extern int LMAXUSED;
extern int NMAXUSED; // mohan 2010-06-17
extern int NCHIUSED;

// band information
extern bool BANDS_CONTROL;
extern int BANDS_START;//mohan add 2010-05-02
extern int BANDS_END;

// other global information
extern int CALSPI; // mohan 2010-06-17
extern int BLOCK_NE; // if ie < BLOCK_NE, we use this Jlq(ie), mohan add 2009-08-26
extern int BLOCK_NE_MIN; // mohan add 2009-08-27
extern bool RESTART;
extern int TEST1;
extern int SCHEME_VALUE;

extern int ACCEPTANCE_STEPS; // mohan add 2009-10-31
extern double ACCEPTANCE_HIGH;
extern double ACCEPTANCE_LOW;

extern double KINETIC_MAX; // mohan add 2009-10-31
extern double KINETIC_DR; // mohan add 2010-04-12
extern int OPTIMIZE_METHOD;// mohan add 2010-04-14

// constant
extern const double PI; // mohan add 2010-04-16
extern const double TWO_PI; // mohan add 2010-06-14
extern const double FOUR_PI;
extern const complex<double> IMAG_UNIT;
extern const complex<double> NEG_IMAG_UNIT;//mohan add 2010-06-14
extern const double SQRT_INVERSE_FOUR_PI;
extern const double PI_HALF;
extern const double SQRT2;

// read wavefunction
extern string WFC_FILE;

void DESTROY();
#endif
