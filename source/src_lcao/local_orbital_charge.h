#ifndef LOCAL_ORBITAL_CHARGE
#define LOCAL_ORBITAL_CHARGE

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/matrix.h"
#include "../module_base/complexmatrix.h"
#include "../src_parallel/parallel_common.h"
#include "../module_gint/grid_technique.h"
#include "src_lcao/record_adj.h"
#include "src_lcao/local_orbital_wfc.h"
#include "src_lcao/LCAO_hamilt.h"
#include "module_psi/psi.h"

class Local_Orbital_Charge
{
	public:

	Local_Orbital_Charge();
	~Local_Orbital_Charge();

	// mohan added 2021-02-08
    void allocate_dm_wfc(const int& lgd,
        Local_Orbital_wfc &lowf,
        psi::Psi<double>* psid,
        psi::Psi<std::complex<double>>* psi);
    // sum bands to compute the electron charge density
	void sum_bands(LCAO_Hamilt &UHM);

	//-----------------
	// in DM_gamma.cpp
	//-----------------
	void allocate_gamma(const int &lgd, psi::Psi<double>* psid);

    void gamma_file(psi::Psi<double>* psid, Local_Orbital_wfc &lowf);
    void cal_dk_gamma_from_2D_pub(void);
    //transformation from 2d block to grid, only gamma_only used it now
    //template<typename T>
    void dm2dToGrid(const psi::Psi<double>& dm2d, double** dm_grid);

	//-----------------
	// in DM_k.cpp
	//-----------------
	void allocate_DM_k(void);

	// liaochen modify on 2010-3-23 
	// change its state from private to public
	double*** DM;	
	double** DM_R;

	// whether to printout density matrix
    static int out_dm; // output density matrix or not.

	void write_dm(const int &is, const int &iter, const std::string &fn, const int &precision);

	void read_dm(const int &is, const std::string &fn);

    
    //-----------------
	// in dm_2d.cpp
    //-----------------
    // dm stands for density matrix
    std::vector<ModuleBase::matrix> dm_gamma;			// dm_gamma[is](iw1,iw2);
    std::vector<ModuleBase::ComplexMatrix> dm_k;		// dm_k[ik](iw1,iw2);

    // use the original formula (Hamiltonian matrix) to calculate energy density matrix	
    std::vector<ModuleBase::ComplexMatrix> edm_k_tddft;

    void init_dm_2d(void);
    
    // dm = wfc.T * wg * wfc.conj(); used in gamma_only
    void cal_dm(const ModuleBase::matrix& wg,   // wg(ik,ib), cal all dm 
        std::vector<ModuleBase::matrix>& wfc_gamma,
        std::vector<ModuleBase::matrix>& dm_gamma);

    // in multi-k,  it is dm(k)
    void cal_dm(const ModuleBase::matrix& wg,    // wg(ik,ib), cal all dm 
        std::vector<ModuleBase::ComplexMatrix>& wfc_k,
        std::vector<ModuleBase::ComplexMatrix>& dm_k);

    // dm(R) = wfc.T * wg * wfc.conj()*kphase, only used in multi-k 
    void cal_dm_R(
        std::vector<ModuleBase::ComplexMatrix>& dm_k,
        Record_adj& ra,
        double** dm2d);     //output, dm2d[NSPIN][LNNR]

    //-----------------
	// wavefunctions' pointer
    //-----------------
    Local_Orbital_wfc* LOWF;
    //-----------------
	// Parallel Variables' pointer
    //-----------------
    const Parallel_Orbitals* ParaV;

    //temporary set it to public for ElecStateLCAO class, would be refactor later
    void cal_dk_k(const Grid_Technique &gt, const ModuleBase::matrix& wg_in);

private:

	// whether the DM array has been allocated
	bool init_DM;
	// whether the DM(R) array has been allocated
	bool init_DM_R;

	void cal_dk_gamma(void);

	// mohan add 2010-09-06
	int lgd_last;// sub-FFT-mesh orbitals number in previous step.
	int lgd_now;// sub-FFT-mesh orbitals number in this step.

	int nnrg_last = 0;// sub-FFT-mesh orbtials number in previous step, with k.
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

#ifdef __MPI
    int setAlltoallvParameter(MPI_Comm comm_2D, int blacs_ctxt, int nblk);
#endif
    void cal_dk_gamma_from_2D(void);
};

#endif
