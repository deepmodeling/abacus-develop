#include "module_base/parallel_global.h"
#ifdef __PEXSI
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
PEXSI_Solver::PEXSI_Solver(const int blacs_text,
                           const int nb,
                           const int nrow,
                           const int ncol,
                           const double* h,
                           const double* s,
                           double& totalEnergyH,
                           double& totalEnergyS,
                           double& totalFreeEnergy)
{
    this->blacs_text = blacs_text;
    this->nb = nb;
    this->nrow = nrow;
    this->ncol = ncol;
    this->h = new double[nrow * ncol];
    this->s = new double[nrow * ncol];
    std::memcpy(this->h, h, nrow * ncol * sizeof(double));
    std::memcpy(this->s, s, nrow * ncol * sizeof(double));
    this->DM = new double[nrow * ncol];
    this->EDM = new double[nrow * ncol];
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
                GlobalV::NLOCAL,
                this->nb,
                this->nrow,
                this->ncol,
                'C',
                this->h,
                this->s,
                GlobalV::nelec,
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

double* PEXSI_Solver::get_DM() const
{
    return DM;
}

double* PEXSI_Solver::get_EDM() const
{
    return EDM;
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