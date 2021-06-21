#ifndef PLANEWAVE_H
#define PLANEWAVE_H

#include "unitcell.h"
#include "klist.h"
#include "../src_global/complexmatrix.h"
#include "../src_global/vector3.h"
#include "../src_parallel/ft.h"

using namespace std;

class PW_Basis
{

public:
    PW_Basis();
    ~PW_Basis();

    void gen_pw(ofstream &log, const UnitCell &Ucell_in, const K_Vectors &Klist_in);

    void set
    (
        const bool &gamma_only_in,
        const double &ecutwfc_in,
        const double &ecutrho_in,
        const int &nx_in,
        const int &ny_in,
        const int &nz_in,
        const int &ncx_in,
        const int &ncy_in,
        const int &ncz_in,
		const int &bx_in,
		const int &by_in,
		const int &bz_in
    );

	const UnitCell *Ucell; // pointer for UnitCell
    const K_Vectors *Klist; // pointer for K-point list

    FFT FFT_chg; // FFT operations for charge density
	FFT FFT_wfc; // FFT operations for wave functions



//===============================================
// Part 1: Kinetic energy cutoff for plane waves
//===============================================
public:
    bool gamma_only;	// only half G are used.
    double ecutwfc;		// Energy cutoff for wavefunctions.
    double ecutrho; 	// Energy cutoff for charge/potential.

    // ggpsi=2*Ecut*(lat0*lat0)/(4.0*PI*PI);
    // in NCPP,  ggchg=ggfft; in USPP (double grids), ggchg=4*ggfft
    double ggpsi; 		// planewave cut off for the wavefunctions, unit (2*PI/lat0)^2
    double ggwfc;		// cutoff for wave function, ggwfc=4*ggpsi
    double ggwfc2;		// ggwav=wfact*ggpsi, default value: wfact=4.0
    double wfac;		// weighting factor from ggpsi to ggwfc2
    double ggchg;   	// G^2 cutoff for supporting charge density,

private:
	// setup for ggpsi, ggwfc, ggwfc2, ggchg
    void setup_gg(void);



//===============================================
// Part 2: FFT dimensions in real space
//===============================================
public:
	// FFT dimensions for wave functions.
	int nx, ny, nz, nxyz;
	// FFT dimensions for charge/potential.
	int ncx, ncy, ncz, ncxyz;
	// nczp: number of x-y planes in each processor
	// nczp_start: starting nczp index in each processor
	int nczp, nczp_start;
	// nrxx: for parallel,local FFT grid dimension
    int nrxx;
	// nrxx_start: starting nrxx in each processor 
    int nrxx_start;

private:
    void setup_FFT_dimension(void);	// set up FFT dimensions



//===============================================
// Part 3: FFT dimensions in real space
// used in LCAO algorithm: grid integrals
//===============================================
public:
	// combine [bx,by,bz] FFT grids into a big one
	// typical values are bx=2, by=2, bz=2
	// bxyz = bx*by*bz, typical value is 8
    int bx, by, bz, bxyz;
	// nbx=ncx/bx, nby=ncy/by, nbz=ncz/bz, 
	// nbxyz=nbx*nby*nbz
	int nbx, nby, nbz, nbxyz;
	// nbzp: number of x-y planes in each processor 
	// nbzp_start: starting nbzp index in each processor
	int nbzp, nbzp_start;
	// nbxx is local dimension of nbxyz
	// nbxx=nbx*nby*nbzp
	int nbxx;



//===============================================
// Part 4: G vectors in reciprocal FFT box
//===============================================
public:
	// mohan add 2008-4-2
    int ngmc_g; // global ngmc (total number of G vectors)
    int ngmc; // num. of G vectors within ggchg in each proc.
	// map: 1D G vector -> charge in FFT grid 
    int *ig2fftc; // dimension: [ngmc]

	// PW_complement::get_ngmw(ngmc, ggwfc2, gg_global, ngmw);
    int ngmw; // num. of G vectors within ggwfc2 in each proc.
	// map: 1D G vector -> wave function in FFT grid
    int *ig2fftw; // dimension: [ngmw]

	// structure factor (ntype, ngmc)
    ComplexMatrix strucFac;
    void setup_structure_factor(void); 		// Calculate structur factors

private:
#ifdef __MPI
    void divide_fft_grid(void);
    void get_MPI_GVectors(void);
    void columns_and_pw_distribution_2(void);
#else
    void get_GVectors(void);
#endif



//===============================================
// Part 5: G vectors, |G|^2, G index [ngmc] 
//===============================================
public:
    Vector3<double> *gdirect;		//(= *G1d) ; // ig = new Vector igc[ngmc],
    Vector3<double> *gdirect_global;	//(= *G1d) ; // ig = new Vector igc[ngmc],
    // store the 3D G vector coordinates for charge density grid;
    // ig.x, ig.y, ig.z should be integers !!!
    // G vectors are in order of increasing G^2
    // can be shared by charge density/potential and wave functions.

    Vector3<double> *gcar;   			//G vectors in cartesian corrdinate
    Vector3<double> *gcar_global;   	//G vectors in cartesian corrdinate
    //g=ig*G ?? HLX (05-26-06): need to check if this is ok!
    //tau is also defined in Cartesian coordintes unit lat0
    Vector3<double> get_GPlusK_cartesian(const int ik, const int ig) const {
        assert(ig>=0 && ig<this->ngmc && ik>=0 && ik<Klist->nks);
        Vector3<double> g_temp_ = Klist->kvec_c[ik] + this->gcar[ig];
        return g_temp_;
    };
    double get_GPlusK_cartesian_projection(const int ik, const int ig, const int axis) const
    {
        assert(ig >= 0 && ig < this->ngmc && ik >= 0 && ik < Klist->nks && axis >= 0 && axis <= 2);
        Vector3<double> g_temp_ = Klist->kvec_c[ik] + this->gcar[ig];
        if (axis == 0)
        {
            return g_temp_.x;
        }
        else if (axis == 1)
        {
            return g_temp_.y;
        }
        else if (axis == 2)
        {
            return g_temp_.z;
        }
        return 0.0;
    }
    double get_SquareGPlusK_cartesian(const int ik, const int ig) const 
    {
        assert(ig >= 0 && ig < this->ngmc && ik >= 0 && ik < Klist->nks);
        Vector3<double> g_temp_ = Klist->kvec_c[ik] + this->gcar[ig];
        return (g_temp_ * g_temp_);
    };
    Vector3<double> get_G_cartesian(const int ig) const 
    {
        assert(ig>=0 && ig<this->ngmc);
        return this->gcar[ig];
    };
    double get_G_cartesian_projection(const int ig, const int axis) const 
    {
        assert(ig>=0 && ig<this->ngmc && axis>=0 && axis<=2);
        if(axis == 0) 
        {
            return this->gcar[ig].x;
        }
        else if(axis == 1)
        {
            return this->gcar[ig].y;
        }
        else if(axis == 2)
        {
            return this->gcar[ig].z;
        }
        return 0.0;
    }
    double get_NormG_cartesian(const int ig) const
    {
        assert(ig >= 0 && ig < this->ngmc);
        return (this->gcar[ig].x * this->gcar[ig].x + this->gcar[ig].y * this->gcar[ig].y + this->gcar[ig].z * this->gcar[ig].z);
    }


    double *gg;       	// modulus (G^2) of G vectors [ngmc]
    double *gg_global;  // modulus (G^2) of G vectors [ngmc_g]
    //gg[ng]=ig[ng]*GGT*ig[ng]/(lat0*lat0)=g[ng]*g[ng] (/lat0*lat0)
	// gg_global dimension: [cutgg_num_now] (save memory skill is used)

	// index of cartesian coordinates [ngmc]
    int *ig1;	//
    int *ig2;	// the indices of G components
    int *ig3;	//



//===============================================
// Part 6: |G|^2 [nggm] 
//===============================================
public:
    // number of |G| shells (i.e., how many different |G|s)
    int nggm;
    double *ggs;	// store |G|^2 for each shell
    int *ig2ngg;	// dimension [ngmc], index from ngmc to nggm
    int gstart;		// first nonzero g vector

private:
    void get_nggm(const int ngmc_local);



//===============================================
// Part 7: phase of e^{-iG*tau_s}
//===============================================
public:
	// phase of e^{-iG*tau_s}
    ComplexMatrix eigts1; // dimension: [Ucell->nat, 2*this->ncx + 1] 
    ComplexMatrix eigts2; // dimension: [Ucell->nat, 2*this->ncy + 1] 
    ComplexMatrix eigts3; // dimension: [Ucell->nat, 2*this->ncz + 1]



//===============================================
// Part 8: update_gvectors 
// LiuXh add 20180515
//===============================================
public:
    double *gg_global0;
    int *cutgg_num_table;
    int ggchg_time_global;

    void update_gvectors(ofstream &log, const UnitCell &Ucell_in);
};
#endif //PlaneWave class
