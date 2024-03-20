//==========================================================
// AUTHOR : Lixin He, Fang Wei, Mohan Chen
// DATA : 2006-11 ~ 2008-11
//==========================================================
#ifndef PLANEWAVE_H
#define PLANEWAVE_H

#include "../src_spillage/common.h"
#include "numerical_basis.h"

struct Way2iw // dimension : WFCALL
{
    int type;
	int i; // atom
    int L;
	int m;
};

class PW_Basis
{
public:

	// about lattice.
	Matrix3 G;
	Matrix3 GT;
	Matrix3 GGT;
	double tpiba; // 2pi/lat0
	double tpiba2; // tpiba^2
	
	// about k point
	int* ngk;
	int** igk;
	int** itia2iat;
	int npwx;

	// about fft grid.
	int nx;
	int ny;
	int nz;
	int nxyz;

    double ggpsi; 					// planewave cut off for the wavefunctions, unit (2*PI/lat0)^2
    double ggwfc2;					// ggwav=wfact*ggpsi, default value: wfact=4.0
	double ggchg;					
    
    complex<double>** strucFac;			// StrucFac [ntype,ngmax]

    int ngmw;						//(= ngmax) / num. of G vectors within ggfft
	int ngmw_g;
	int ngmw_start;

	int ngmc;
	int ngmc_g;
	int ngmc_start;
    
	Vector3<double> *gdirect;		//(= *G1d) ; // ig = new Vector igc[ngmc],
    Vector3<double> *gdirect_global;	//(= *G1d) ; // ig = new Vector igc[ngmc],
    Vector3<double> *gcar;   			//G vectors in cartesian corrdinate
    Vector3<double> *gcar_global;   	//G vectors in cartesian corrdinate
    double *gg;         			// modulus (G^2) of G vectors
    double *gg_global;    			// modulus (G^2) of G vectors

    ComplexMatrix eigts1;   		//
    ComplexMatrix eigts2;   		//the phases e^{-iG*tau_s}
    ComplexMatrix eigts3;   		//
    int *ig1;        				//
    int *ig2;        				// the indices of G components
    int *ig3;        				//

    PW_Basis();
    ~PW_Basis();

    void init(void);
    void setup_structure_factor(void); 		// Calculate structur factors
	void get_sk(void);
	complex<double>* get_sk(const int ik, const int it, const int ia)const;
	Vector3<double> get_1qvec_cartesian(const int ik,const int ig)const;

	void table(void);
	Numerical_Basis NBasis;

	void allocate_psi1d(const int &il);
	void allocate_psi3d(const int &level);
	realArray psi1d;
	ComplexMatrix *psi3d;
	matrix *ylm;

	void calculate_psi1d(void);
	// FUNCTION:
	// use psi1d to generate 3D wave functions.
	// psi3d(it, ia, l, n, m) for each k point.
	void calculate_psi3d(const int &ilevel, const int &ik);

	void update_psi1d(const int &il, const int &ic, const int &ie, 
		const double &c4_now, const double &c4_old);

	int nwfc2;

	void update_psi3d( const int &il, const int &ic, const int &ik);

	int kmesh;
	double Dk;
	
	// FNCTION:
	// calculate the overlap between psi(it, ia, l, n, m)
	// be called in SpillageStep::init_QS_matrix
	complex<double> calculateS(const int &iw, const int &iw2, const int &ik);

	void calculate_Jlq(const int &ik, const int &iw, const int &ie);
	complex<double> *Jlq;
	complex<double> calculate_Jlq_Phi(const int &ik, const int &mu);

	complex<double> *save;
	int* posmu;
	int* posnu;

private:
	void set_igk(void);
    void setup_gg(void);
    void setup_FFT_dimension(void); 		// set up FFT dimensions
#ifdef __MPI
    void get_MPI_GVectors(void);
#else
    void get_GVectors(void);
#endif

	// FUNCTION:
	// use interpolation scheme to get the one dimension wave function for
	// each k point
	double Polynomial_Interpolation(const int &it, const int &l, const int &n, const double &gnorm)const;

	Way2iw* iwindex;
	
};
#endif //PlaneWave class
