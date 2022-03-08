#ifndef LOCAL_ORBITAL_WFC
#define LOCAL_ORBITAL_WFC

#include "grid_technique.h"
#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_orbital/ORB_control.h" // mohan add 2021-05-24

class Local_Orbital_wfc
{
public:

	Local_Orbital_wfc();
	~Local_Orbital_wfc();

    ///=========================================
    /// grid wfc 
    /// used to generate density matrix: LOC.DM_R,
	/// which is used to calculate the charge density. 
	/// which is got after the diagonalization of 
    /// std::complex Hamiltonian matrix.
    ///=========================================
    //( Old Name: WFC_K)
    std::complex<double>*** wfc_k_grid; // [NK, GlobalV::NBANDS, GlobalV::NLOCAL]	
    //( Old Name: WFC_K_POOL)
    std::complex<double>* wfc_k_grid2; // [NK*GlobalV::NBANDS*GlobalV::NLOCAL]

    ///=========================================
    /// 2d wfc
    /// directly output from elpa interface
    /// used to calculate density matrix LOC.dm_gamma and LOC.dm_k
    ///=========================================
    std::vector<ModuleBase::matrix> wfc_gamma;			// dm_gamma[is](iw1,iw2);
    std::vector<ModuleBase::ComplexMatrix> wfc_k;		// dm_k[ik](iw1,iw2);
    std::vector<ModuleBase::ComplexMatrix> wfc_k_laststep; // store the 2d wavefunction of the last step in TDDFT

    const Parallel_Orbitals *ParaV;


    void allocate_k(const Grid_Technique& gt,
        Local_Orbital_wfc &lowf);

    //=========================================
    // Init Cij, make it satisfy 2 conditions:
    // (1) Unit
    // (2) Orthogonal <i|S|j>= \delta{ij}
    //=========================================
	// void init_Cij(const bool change_c = 1);


    ///=========================================
    ///Parallel: map of index in 2D distribution: global<->local
    ///=========================================
    static int globalIndex(int localindex, int nblk, int nprocs, int myproc);
    static int localIndex(int globalindex, int nblk, int nprocs, int& myproc);
    
    ///=========================================
    ///Parallel: convert the distribution of wavefunction from 2D to grid
    ///=========================================
    //name will be changed
    
    ///for gamma_only, with output
    int q2CTOT(
        int myid,
        int naroc[2],
        int nb,
        int dim0,
        int dim1,
        int iprow,
        int ipcol,
        int loc_size,
        double* work,
        double** CTOT);
    
    ///for multi-k, no output
    int q2WFC_complex(
        int naroc[2],
        int nb,
        int dim0,
        int dim1,
        int iprow,
        int ipcol,
        std::complex<double>* work,
        std::complex<double>** WFC);
    
    ///for multi-k, with output
    int q2WFC_CTOT_complex(
        int myid,
        int naroc[2],
        int nb,
        int dim0,
        int dim1,
        int iprow,
        int ipcol,
        std::complex<double>* work,
        std::complex<double>** WFC,
        std::complex<double>** CTOT);
        
private:

	bool wfck_flag; 
	bool complex_flag;
	bool allocate_flag;

};

#endif
