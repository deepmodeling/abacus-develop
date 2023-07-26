//=========================================================
//AUTHOR : mohan
//DATE : 2009-09-16
//REFACTOR : Peize Lin, 2021.06.28
//=========================================================
#ifndef GINT_GAMMA_H
#define GINT_GAMMA_H
#include "gint.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "grid_technique.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_matrix.h"
#ifdef _OPENMP
#include <omp.h>
#endif

//=========================================================
// ModuleBase::Integral On 3D Grids, different from Grid_Integral
// Feature : Matrix Elements Of Local Potential For 
// Numerical Orbitals
//=========================================================

class Gint_Gamma : public Gint
{
	public:

	Gint_Gamma();
	~Gint_Gamma();

    //------------------------------------------------------
    // in gint_gamma_vl.cpp 
    //------------------------------------------------------
    // there is an additional step in calculating vlocal for gamma point
    // namely the redistribution of Hamiltonian from grid to 2D block format
    // hence we have an additional layer outside the unified interface
	void cal_vlocal(Gint_inout *inout, const bool new_e_iteration);

    //------------------------------------------------------
    // in gint_gamma_env.cpp 
    //------------------------------------------------------
	// calcualte the envelope function
	void cal_env(const double* wfc, double* rho);

private:

    double***  DM;   //pointer to LOC.DM

    ///------------------------------------------------------
    /// in gint_gamma_vl.cpp 
    ///------------------------------------------------------
    /// method for redistributing the Hamiltonian
    /// from grid to 2D format
    /// pass a setter function to customize row/col major and outputs
    void vl_grid_to_2D(const double* vl_grid,
        const Parallel_2D& p2d,
        const int loc_grid_dim,
        const bool new_e_iteration,
        double* vl_2d,
        std::function<void(const int&, const int&, const double&, double*)> setfunc);

    ///===============================
    /// Use MPI_Alltoallv to convert a grid distributed matrix
    /// to 2D - block cyclic distributed matrix.
    ///===============================
    int sender_index_size;
    int *sender_local_index;
    int sender_size;
    int *sender_size_process;
    int *sender_displacement_process;
    double* sender_buffer;

    int receiver_index_size;
    int *receiver_global_index;
    int receiver_size;
    int *receiver_size_process;
    int *receiver_displacement_process;
    double* receiver_buffer;

};

#endif
