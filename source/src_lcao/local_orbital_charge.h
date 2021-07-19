#ifndef LOCAL_ORBITAL_CHARGE
#define LOCAL_ORBITAL_CHARGE

#include "../src_pw/tools.h"
#include "grid_technique.h"
#include "wfc_dm_2d.h"

class Local_Orbital_Charge
{

	public:

	Local_Orbital_Charge();
	~Local_Orbital_Charge();

	// mohan added 2021-02-08
	void allocate_dm_wfc(const Grid_Technique &gt);
	// sum bands to compute the electron charge density
	void sum_bands(void);

	//-----------------
	// in DM_gamma.h
	//-----------------
	void allocate_gamma(const Grid_Technique &gt);

	void gamma_file(const Grid_Technique &gt);
	void cal_dk_gamma_from_2D_pub(void);


	//-----------------
	// in DM_k.h
	//-----------------
	void allocate_DM_k(void);
	

	// liaochen modify on 2010-3-23 
	// change its state from private to public
	double*** DM;	
	double** DM_R;

	// whether to printout density matrix
	int out_dm; // output density matrix or not.

	void write_dm(const int &is, const int &iter, const string &fn, const int &precision);

	void read_dm(const int &is, const string &fn);
	
	Wfc_Dm_2d wfc_dm_2d;		// Peize Lin test 2019-01-16

private:

	// whether the DM array has been allocated
	bool init_DM;
	// whether the DM(R) array has been allocated
	bool init_DM_R;

	void cal_dk_gamma(void);

	void cal_dk_k(const Grid_Technique &gt);

	// mohan add 2010-09-06
	int lgd_last;// sub-FFT-mesh orbitals number in previous step.
	int lgd_now;// sub-FFT-mesh orbitals number in this step.

	int nnrg_last;// sub-FFT-mesh orbtials number in previous step, with k.
	int nnrg_now; // sub-FFT-mesh orbitals number in this step, with k.

	// add by yshen on 9/22/2014
	// these variables are memory pool for DM series matrixes, 
	// so that these matrixes will be storaged continuously in the memory.
	double **DM_pool;
	
	// Buffer parameters for tranforming 2D block-cyclic distributed DM matrix 
	// to grid distributed DM matrix
    int *sender_2D_index;
    int sender_size;
    int *sender_size_process;
    int *sender_displacement_process;
    double* sender_buffer;

    int *receiver_local_index;
    int receiver_size;
    int *receiver_size_process;
    int *receiver_displacement_process;
    double* receiver_buffer;

    int setAlltoallvParameter(MPI_Comm comm_2D, int blacs_ctxt, int nblk);

    void cal_dk_gamma_from_2D(void);
};

#endif
