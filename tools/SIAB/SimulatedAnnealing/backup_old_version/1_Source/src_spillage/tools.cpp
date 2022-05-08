#include "tools.h"

//plane wave , mohan add 2010-06-14
PW_Basis PW;
MultiZeta mz;
Read_INPUT input;
Parallel_Kpoints Pkpoints;
ofstream ofs_running;

// parallel information
int NPROC = 1; // mohan add 2010-06-12
int KPAR = 1;
int MY_RANK = 0;
int MY_POOL = 0;
int NPROC_IN_POOL = 1;
int RANK_IN_POOL = 0;

// global + input information
int NKSTOT = -1;
int NKS = -1;
int NBANDS = -1;
int NWFCALL = -1;
int NE = -1;

// global information
int STRNUM = -1;
int NTYPE = -1;
string *LABEL;
int *NA;
double **CARPOSX;
double **CARPOSY;
double **CARPOSZ;
Matrix3 LATVEC;

// plane wave information
bool USEPW = true;
double LAT0 = -1.0;

// k points
double *CARKX;
double *CARKY;
double *CARKZ;

// orbital + energy cutoff
double ECUT = -1.0;
double ECUT_JLQ = -1.0;
double RCUT = -1.0;
bool SMOOTH = false;  
double SIGMA = -0.1;
double TOLERENCE = 1.0e-12;
int LMAXALL = -1;
int LMAXUSED = -1;
int NMAXUSED = -1;
int NCHIUSED = -1;

// band information
bool BANDS_CONTROL=false; // mohan add 2010-05-02
int BANDS_START=0;
int BANDS_END=0; // needed to be checked carefully if not read in. 

// other global information
int CALSPI = 0; // mohan 2010-06-17
int BLOCK_NE = -10000;
int BLOCK_NE_MIN = -10000;
bool RESTART = false;
int TEST1 = false;
int SCHEME_VALUE = 1;

int ACCEPTANCE_STEPS = 100; // mohan add 2009-10-31
double ACCEPTANCE_HIGH = 0.6;
double ACCEPTANCE_LOW = 0.3;

double KINETIC_MAX = 30.0; // (unit?????)
double KINETIC_DR = 0.01; //
int OPTIMIZE_METHOD = 1; // 1: Kin 2: Ecut

// constant
const double PI = 3.1415926535897;
const double TWO_PI = 2.0 * PI;
const double FOUR_PI = 4.0 * PI;
const double SQRT_INVERSE_FOUR_PI = sqrt(1.0/FOUR_PI);
const complex<double> IMAG_UNIT = complex<double>(0, 1.0);
const complex<double> NEG_IMAG_UNIT = complex<double>(0,-1.0);
const double PI_HALF = PI / 2.0;
const double SQRT2 = 1.41421356237309504880;

// read wavefunction
string WFC_FILE = "WAVEFUNC.dat";

void DESTROY(void)
{
	if(NTYPE>0)
	{
		delete[] LABEL;
		delete[] NA;
		for(int i=0; i<NTYPE; i++)
		{
			delete[] CARPOSX[i];
			delete[] CARPOSY[i];
			delete[] CARPOSZ[i];
		}
		delete[] CARPOSX;
		delete[] CARPOSY;
		delete[] CARPOSZ;
	}
	if(NKS>0)
	{
		delete[] CARKX;
		delete[] CARKY;
		delete[] CARKZ;
	}
}



