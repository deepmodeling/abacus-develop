#ifndef LOCAL_ORBITAL_CHARGE
#define LOCAL_ORBITAL_CHARGE

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/matrix.h"
#include "module_base/complexmatrix.h"
#include "module_base/parallel_common.h"
#include "module_hamilt_lcao/module_gint/grid_technique.h"
#include "module_hamilt_lcao/hamilt_lcaodft/record_adj.h"
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_wfc.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_hamilt.h"
#include "module_psi/psi.h"
#include "module_elecstate/elecstate.h"
#include "module_hamilt_lcao/hamilt_lcaodft/DM_gamma_2d_to_grid.h"
class Local_Orbital_Charge
{
	public:

	Local_Orbital_Charge();
	~Local_Orbital_Charge();

	// mohan added 2021-02-08
    void allocate_dm_wfc(const Grid_Technique& gt,
        elecstate::ElecState* pelec,
        Local_Orbital_wfc &lowf,
        psi::Psi<double>* psi,
        const K_Vectors& kv);
    void allocate_dm_wfc(const Grid_Technique& gt,
        elecstate::ElecState* pelec,
        Local_Orbital_wfc& lowf,
        psi::Psi<std::complex<double>>* psi,
        const K_Vectors& kv);
	//-----------------
	// in DM_gamma.cpp
	//-----------------
	void allocate_gamma(const int &lgd, psi::Psi<double>* psid, elecstate::ElecState* pelec, const int& nks);

    void gamma_file(psi::Psi<double>* psid, Local_Orbital_wfc &lowf, elecstate::ElecState* pelec);
    void cal_dk_gamma_from_2D_pub(void);

	//-----------------
	// in DM_k.cpp
	//-----------------
	void allocate_DM_k(const int& nks, const int& nnrg);

	// liaochen modify on 2010-3-23 
	// change its state from private to public
	double*** DM;	
	double** DM_R;

	// whether to printout density matrix
    static int out_dm; // output density matrix or not.
    static int out_dm1;
    
    //-----------------
	// in dm_2d.cpp
    //-----------------
    // dm stands for density matrix
    std::vector<ModuleBase::matrix> dm_gamma;			// dm_gamma[is](iw1,iw2);
    std::vector<ModuleBase::ComplexMatrix> dm_k;		// dm_k[ik](iw1,iw2);

    // use the original formula (Hamiltonian matrix) to calculate energy density matrix	
    std::vector<ModuleBase::ComplexMatrix> edm_k_tddft;

    void init_dm_2d(const int& nks);
    
    // dm(R) = wfc.T * wg * wfc.conj()*kphase, only used in multi-k 
    void cal_dm_R(
        std::vector<ModuleBase::ComplexMatrix>& dm_k,
        Record_adj& ra,
        double** dm2d,
        const K_Vectors& kv);     //output, dm2d[NSPIN][LNNR]

    //-----------------
	// wavefunctions' pointer
    //-----------------
    Local_Orbital_wfc* LOWF;
    //-----------------
	// Parallel Variables' pointer
    //-----------------
    const Parallel_Orbitals* ParaV;

    //temporary set it to public for ElecStateLCAO class, would be refactor later
    void cal_dk_k(const Grid_Technique &gt, const ModuleBase::matrix& wg_in, const K_Vectors& kv);

    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> DMR_sparse;

    void set_dm_k(int ik, std::complex<double>* dm_k_in); // set dm_k from a pointer
    void set_dm_gamma(int is, double* dm_gamma_in); // set dm_gamma from a pointer

private:

	// whether the DM array has been allocated
	bool init_DM;
	// whether the DM(R) array has been allocated
	bool init_DM_R;

	// mohan add 2010-09-06
	int lgd_last;// sub-FFT-mesh orbitals number in previous step.
	int lgd_now;// sub-FFT-mesh orbitals number in this step.

	int nnrg_last = 0;// sub-FFT-mesh orbtials number in previous step, with k.
	int nnrg_now; // sub-FFT-mesh orbitals number in this step, with k.

	// add by yshen on 9/22/2014
	// these variables are memory pool for DM series matrixes, 
	// so that these matrixes will be storaged continuously in the memory.
    double** DM_pool;

    DMgamma_2dtoGrid dm2g;

};

#endif
