//Wenfei Li, November 2022
//A new implementation of CG relaxation
#ifndef RELAX1_H
#define RELAX1_H

#include "module_base/matrix.h"
#include "module_base/matrix3.h"
#include "line_search.h"

class Relax
{
    public:

    Relax(){};
    ~Relax(){};

    //prepare for relaxation
    void init_relax(const int nat_in);
    //perform a single relaxation step
    bool relax_step(const ModuleBase::matrix& force, const ModuleBase::matrix &stress, const double etot_in);

    private:

    int istep; //count ionic step

    //setup gradient based on force and stress
    //constraints are considered here
    //also check if relaxation has converged
    //based on threshold in force & stress
    bool setup_gradient(const ModuleBase::matrix& force, const ModuleBase::matrix &stress);

    //check whether previous line search is done
    bool check_line_search();

    //if line search not done : perform line search
    void perform_line_search();

    //if line search done: find new search direction and make a trial move
    void new_direction();

    //move ions and lattice vectors
    void move_cell_ions(const bool is_new_dir);

    int nat; // number of atoms
    bool ltrial; // if last step is trial step

    double step_size;

    // Gradients; _p means previous step
    ModuleBase::matrix grad_ion;
    ModuleBase::matrix grad_cell;
    ModuleBase::matrix grad_ion_p;
    ModuleBase::matrix grad_cell_p;

    // Search directions; _p means previous step
    ModuleBase::matrix search_dr_ion;
    ModuleBase::matrix search_dr_cell;
    ModuleBase::matrix search_dr_ion_p;
    ModuleBase::matrix search_dr_cell_p;

    // Used for applyting constraints
    bool if_cell_moves;

    //Keeps track of how many CG trial steps have been performed,
    //namely the number of CG directions followed
    //Note : this should not be confused with number of ionic steps
    //which includes both trial and line search steps
    int cg_step;

    //in CG, search_dr = search_dr_p + grad * gamma
    double gamma;
    void calculate_gamma();

    //Intermediate variables
    //I put them here because they are used across different subroutines
    double sr_sr, srp_srp; //inner/cross products between search directions
    double gr_gr, gr_grp, grp_grp; //inner/cross products between gradients
    double gr_sr; //cross product between search direction and gradient
    double e1ord1, e1ord2, e2ord, e2ord2;
    double dmove,dmovel,dmoveh;
    double etot, etot_p;
    double force_thr_eva;

    bool brent_done; //if brent line search is finished

    double fac_force;
    double fac_stress;

    ModuleBase::Matrix3 latvec_save;
    Line_Search ls;
};

#endif