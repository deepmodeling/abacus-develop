#ifndef CONSTANT_H
#define CONSTANT_H
#include <complex>

//==========================================================
// GLOBAL CONSTANTS 
//==========================================================

//==========================================================
// EXPLAIN : constants
// NAME : pi
// NAME : tpi (two pi)
// NAME : fpi (four pi)
// NAME : sqrtpi ( sqrt pi )
// NAME : sqrtpm1 ( 1/sqrtpi )
//==========================================================
namespace ModuleBase
{
const double PI		  				= 3.14159265358979323846;
const double PI_HALF				= PI / 2.0;
const double TWO_PI					= 2 * PI;
const double FOUR_PI   				= 4.0 * 3.14159265358979323846;
//const double SQRT_PI				= 1.77245385090551602729;
//const double INVERSE_SQRT_PI		= 1.0 / SQRT_PI;
const double INVERSE_FOUR_PI		= 1.0/FOUR_PI;
const double SQRT_INVERSE_FOUR_PI 	= sqrt(INVERSE_FOUR_PI);
const double SQRT2 = 1.41421356237309504880;
//const double SQRT3 = 1.73205080756887729352;

//==========================================================
// EXPLAIN : std::complex constants
//==========================================================
const std::complex<double> ZERO(0.0, 0.0);
const std::complex<double> ONE(1.0, 0.0);
const std::complex<double> NEG_ONE(-1.0, 0.0);
const std::complex<double> IMAG_UNIT(0.0,1.0);
const std::complex<double> NEG_IMAG_UNIT(0.0,-1.0);

//==========================================================
// EXPLAIN : physical constants
//==========================================================
const double K_BOLTZMAN_SI    = 1.3806504e-23;// J K^-1
const double K_BOLTZMAN_AU    = 3.1667e-6;// Hartree K^-1
const double Hartree_to_K     = 3.1577464e5; // Hartree to K
//const double K_BOLTZMAN_RY    = 6.3335e-6;// Rydberg K^-1; mohan add 2010-09-03
//const double K_BOLTZMAN_EV    = 8.6173e-5; // eV; mohan add 2010-09-03
//const double K_BOLTZMAN_M1_AU = 315795.260;// Hartree^-1 K
//const double FACTEM           = 315795.260;// 27.212d0*11605.d0 Hartree^-1 K

//==========================================================
// EXPLAIN : physical constants define the Atomic Units 
//==========================================================
const double BOHR_RADIUS_SI   = 0.529177e-10;	// m
//const double BOHR_RADIUS_CM   = 0.529177e-8;	// cm
const double BOHR_TO_A = 0.5291770;		// angstrom
//const double ELECTRONMASS_SI  = 9.10953e-31;	// kg
//const double ELECTRONMASS_UMA = 5.4858e-4;		// uma

//==========================================================
// EXPLAIN : units conversion factors
//==========================================================
const double ELECTRONVOLT_SI  = 1.6021892e-19;	// J
//const double UMA_SI           = 1.66057e-27;	// Kg
const double ANGSTROM_AU      = 1.8897270;		// au
//const double AU_TO_OHMCMM1    = 46000.00;		// (ohm cm)^-1
//const double AU_KB            = 294210.00;		// Kbar
//const double KB_AU            = 1.00 / 294210.00;// au
//const double AU_GPA           = 29421.00;		// GPa
//const double GPA_AU           = 1.00 / 29421.00;// au
//const double SCMASS           = 1822.890;		// uma to au ( mass of a proton )
//const double UMA_AU           = 1822.890;		// au
//const double AU_TERAHERTZ     = 2.418e-5;		// THz
//const double TERAHERTZ        = 2.418e-5;		// from au to THz
//const double AU_SEC           = 2.4189e-17;		// sec
const double AU_to_FS           = 2.418884326505e-2;		// from a.u. to fs
//const double rhothr = 1.0e-5;					// tolerance
//const double gsmall = 1.0e-12;
const double e2 = 2.0;							// the square of the electron charge
const double DEGSPIN = 2.0;						// the number of spins per level
const double Hartree_to_eV = 27.211396;// slcbb  // 27.21138344; // eV
const double Ry_to_eV = 13.605698;				// 13.60569172; // conversion from Ry to eV
//const double eV_to_kelvin = 1.16044;  			// from ev to Kelvin
const double NA = 6.02214129e23;				// mol
const double EMASS_SI = 9.1093826e-31;         // mass of electron (kg)
const double AU_to_MASS = NA*EMASS_SI*1e3;     // mass a.u. to g/mol

const double HARTREE_SI = 4.35974394e-18; //J
const double RYDBERG_SI =  HARTREE_SI/2.0; //J
//const double RY_TO_KELVIN = RYDBERG_SI / K_BOLTZMAN_SI;

//const double AMCONV = 1.660538782e-27 / 9.10938215e-31 * 0.50; // mass conversion: a.m.u to a.u. (Ry)

//const double uakbar = 147105.0;					// pressure conversion from Ry/(a.u)^3 to K

// zero up to a given accuracy
//const double epsr  = 1.0e-6;
const double threshold_wg  = 1.0e-10;
}

#endif 
