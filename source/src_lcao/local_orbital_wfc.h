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

#ifdef __MPI
    ///=========================================
    ///Parallel: convert the distribution of wavefunction from 2D to grid
    ///=========================================
    /// For gamma_only, T = double; 
    /// For multi-k, T = complex<double>;
    /// Set myid and ctot when output is needed;
    /// Set wfc as nullptr when 2d-to-grid convertion is not needed.
    
    // Notice: here I reload this function rather than use template
    // (in which the implementation should be put in header file )
    // because sub-function `write_lowf_complex`contains GlobalC declared in `global.h`
    // which will cause lots of "not defined" if included in a header file.
    void wfc_2d_to_grid(int out_wfc_lcao, double* wfc_2d, double** wfc_grid);
    void wfc_2d_to_grid(int out_wfc_lcao, std::complex<double>* wfc_2d, std::complex<double>** wfc_grid, int ik);
#endif

private:
    template <typename T>
    int set_wfc_grid(
        int naroc[2],
        int nb,
        int dim0,
        int dim1,
        int iprow,
        int ipcol,
        T* work,
        T** wfc,
        int myid = -1,
        T** ctot = nullptr);

	bool wfck_flag; 
	bool complex_flag;
	bool allocate_flag;

};

#ifdef __MPI
template <typename T>
int Local_Orbital_wfc::set_wfc_grid(
    int naroc[2],
    int nb,
    int dim0,
    int dim1,
    int iprow,
    int ipcol,
    T* work,
    T** wfc,
    int myid,
    T** ctot)
{
    ModuleBase::TITLE(" Local_Orbital_wfc", "set_wfc_grid");
    if (!wfc && !ctot) return 0;
    for (int j = 0; j < naroc[1]; ++j)
    {
        int igcol=globalIndex(j, nb, dim1, ipcol);
        if(igcol>=GlobalV::NBANDS) continue;
        for(int i=0; i<naroc[0]; ++i)
        {
            int igrow=globalIndex(i, nb, dim0, iprow);
	        int mu_local=GlobalC::GridT.trace_lo[igrow];
            if (wfc && mu_local >= 0)
            {
                wfc[igcol][mu_local]=work[j*naroc[0]+i];
            }
            if (ctot && myid == 0)
                ctot[igcol][igrow] = work[j * naroc[0] + i];
        }
    }
    return 0;
}
#endif

#endif
