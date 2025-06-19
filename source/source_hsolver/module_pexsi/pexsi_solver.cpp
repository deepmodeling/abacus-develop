#include "module_base/parallel_global.h"
#ifdef __PEXSI
#include "module_parameter/parameter.h"
#include "pexsi_solver.h"

#include <mpi.h>
#include <cstring>
#include <vector>

#include "module_base/global_variable.h"
#include "simple_pexsi.h"

extern MPI_Comm DIAG_WORLD;
extern MPI_Comm GRID_WORLD;
namespace pexsi
{

int PEXSI_Solver::pexsi_npole = 0;
bool PEXSI_Solver::pexsi_inertia = 0;
int PEXSI_Solver::pexsi_nmax = 0;
// int PEXSI_Solver::pexsi_symbolic = 0;
bool PEXSI_Solver::pexsi_comm = 0;
bool PEXSI_Solver::pexsi_storage = 0;
int PEXSI_Solver::pexsi_ordering = 0;
int PEXSI_Solver::pexsi_row_ordering = 0;
int PEXSI_Solver::pexsi_nproc = 0;
bool PEXSI_Solver::pexsi_symm = 0;
bool PEXSI_Solver::pexsi_trans = 0;
int PEXSI_Solver::pexsi_method = 0;
int PEXSI_Solver::pexsi_nproc_pole = 0;
// double PEXSI_Solver::pexsi_spin = 2;
double PEXSI_Solver::pexsi_temp = 0.0;
double PEXSI_Solver::pexsi_gap = 0.0;
double PEXSI_Solver::pexsi_delta_e = 0.0;
double PEXSI_Solver::pexsi_mu_lower = 0.0;
double PEXSI_Solver::pexsi_mu_upper = 0.0;
double PEXSI_Solver::pexsi_mu = 0.0;
double PEXSI_Solver::pexsi_mu_thr = 0.0;
double PEXSI_Solver::pexsi_mu_expand = 0.0;
double PEXSI_Solver::pexsi_mu_guard = 0.0;
double PEXSI_Solver::pexsi_elec_thr = 0.0;
double PEXSI_Solver::pexsi_zero_thr = 0.0;

void PEXSI_Solver::prepare(const int blacs_text,
                           const int nb,
                           const int nrow,
                           const int ncol,
                           const double* h,
                           const double* s,
                           double*& _DM,
                           double*& _EDM)
{
    this->blacs_text = blacs_text;
    this->nb = nb;
    this->nrow = nrow;
    this->ncol = ncol;
    this->h = const_cast<double*>(h);
    this->s = const_cast<double*>(s);
    this->DM = _DM;
    this->EDM = _EDM;
    this->totalEnergyH = 0.0;
    this->totalEnergyS = 0.0;
    this->totalFreeEnergy = 0.0;
}

int PEXSI_Solver::solve(double mu0)
{
    MPI_Group grid_group;
    int myid, grid_np;
    MPI_Group world_group;
    MPI_Comm_rank(DIAG_WORLD, &myid);
    MPI_Comm_size(DIAG_WORLD, &grid_np);
    MPI_Comm_group(DIAG_WORLD, &world_group);

    int grid_proc_range[3]={0, (GlobalV::NPROC/grid_np)*grid_np-1, GlobalV::NPROC/grid_np};
    MPI_Group_range_incl(world_group, 1, &grid_proc_range, &grid_group);

    simplePEXSI(DIAG_WORLD,
                DIAG_WORLD,
                grid_group,
                this->blacs_text,
                PARAM.globalv.nlocal,
                this->nb,
                this->nrow,
                this->ncol,
                'c',
                this->h,
                this->s,
                PARAM.inp.nelec,
                "PEXSIOPTION",
                this->DM,
                this->EDM,
                this->totalEnergyH,
                this->totalEnergyS,
                this->totalFreeEnergy,
                mu,
                mu0);
    return 0;
}

const double PEXSI_Solver::get_totalFreeEnergy() const
{
    return totalFreeEnergy;
}

const double PEXSI_Solver::get_totalEnergyH() const
{
    return totalEnergyH;
}

const double PEXSI_Solver::get_totalEnergyS() const
{
    return totalEnergyS;
}

const double PEXSI_Solver::get_mu() const
{
    return mu;
}

} // namespace pexsi
#endif