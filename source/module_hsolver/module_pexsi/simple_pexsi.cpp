// use PEXSI to solve a Kohn-Sham equation
// the H and S matrices are given by 2D block cyclic distribution
// the Density Matrix and Energy Density Matrix calculated by PEXSI are transformed to 2D block cyclic distribution
// #include "mpi.h"
#ifdef __PEXSI
#include <mpi.h>

#include <cfloat>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>

#include "c_pexsi_interface.h"
#include "dist_bcd_matrix.h"
#include "dist_ccs_matrix.h"
#include "dist_matrix_transformer.h"
#include "module_base/lapack_connector.h"
#include "module_base/timer.h"
#include "module_base/tool_quit.h"
#include "module_base/global_variable.h"
#include "module_hsolver/diago_pexsi.h"

namespace pexsi
{
inline void strtolower(char* sa, char* sb)
{
    char c;
    int len = strlen(sa);
    for (int i = 0; i < len; i++)
    {
        c = sa[i];
        sb[i] = tolower(c);
    }
    sb[len] = '\0';
}

inline void setDefaultOption(int* int_para, double* double_para)
{
    double_para[0] = 2;
    double_para[2] = 0;
    double_para[11] = DBL_MIN;
    int_para[3] = 0;
    int_para[6] = 1;
    int_para[8] = 0;
    int_para[9] = 0;
    int_para[11] = 0;
    int_para[12] = 0;
    int_para[14] = 2;
    int_para[15] = 1;
}

int loadPEXSIOption(MPI_Comm comm,
                    const std::string PexsiOptionFile,
                    PPEXSIOptions& options,
                    int& numProcessPerPole,
                    double& ZERO_Limit)
{

    // temp variable arrays read from conf file and will be bcast to all processors

    // parameter array of type int,
    //  0: numPole
    //  1: isInertiaCount
    //  2: maxPEXSIIter
    //  3: matrixType
    //  4: isSymbolicFactorize
    //  5: isConstructCommPattern
    //  6: solver
    //  7: symmetricStorage
    //  8: ordering
    //  9: rowOrdering
    // 10: npSymbFact
    // 11: symmetric
    // 12: transpose
    // 13: method
    // 14: nPoints
    // 15: verbosity
    // 16: numProcessPerPole
    int int_para[17];

    // parameter array of type double
    //  0: spin
    //  1: temperature
    //  2: gap
    //  3: deltaE
    //  4: muMin0
    //  5: muMax0
    //  6: mu0
    //  7: muInertiaTolerance
    //  8: muInertiaExpansion
    //  9: muPEXSISafeGuard
    // 10: numElectronPEXSITolerance
    // 11: ZERO_Limit
    double double_para[12];

    // read in PEXSI options from GlobalV
    int_para[0] = hsolver::DiagoPexsi<double>::pexsi_npole;
    int_para[1] = hsolver::DiagoPexsi<double>::pexsi_inertia;
    int_para[2] = hsolver::DiagoPexsi<double>::pexsi_nmax;
    int_para[3] = 0;
    int_para[4] = 1; // hsolver::DiagoPexsi<double>::pexsi_symbolic;
    int_para[5] = hsolver::DiagoPexsi<double>::pexsi_comm;
    int_para[6] = 0;
    int_para[7] = hsolver::DiagoPexsi<double>::pexsi_storage;
    int_para[8] = hsolver::DiagoPexsi<double>::pexsi_ordering;
    int_para[9] = hsolver::DiagoPexsi<double>::pexsi_row_ordering;
    int_para[10] = hsolver::DiagoPexsi<double>::pexsi_nproc;
    int_para[11] = hsolver::DiagoPexsi<double>::pexsi_symm;
    int_para[12] = hsolver::DiagoPexsi<double>::pexsi_trans;
    int_para[13] = hsolver::DiagoPexsi<double>::pexsi_method;
    int_para[14] = 2;
    int_para[15] = 0;
    int_para[16] = hsolver::DiagoPexsi<double>::pexsi_nproc_pole;

    double_para[0] = GlobalV::NSPIN; // hsolver::DiagoPexsi<double>::pexsi_spin;
    double_para[1] = hsolver::DiagoPexsi<double>::pexsi_temp;
    double_para[2] = hsolver::DiagoPexsi<double>::pexsi_gap;
    double_para[3] = hsolver::DiagoPexsi<double>::pexsi_delta_e;
    double_para[4] = hsolver::DiagoPexsi<double>::pexsi_mu_lower;
    double_para[5] = hsolver::DiagoPexsi<double>::pexsi_mu_upper;
    double_para[6] = hsolver::DiagoPexsi<double>::pexsi_mu;
    double_para[7] = hsolver::DiagoPexsi<double>::pexsi_mu_thr;
    double_para[8] = hsolver::DiagoPexsi<double>::pexsi_mu_expand;
    double_para[9] = hsolver::DiagoPexsi<double>::pexsi_mu_guard;
    double_para[10] = hsolver::DiagoPexsi<double>::pexsi_elec_thr;
    double_para[11] = hsolver::DiagoPexsi<double>::pexsi_zero_thr;

    options.numPole = int_para[0];
    options.isInertiaCount = int_para[1];
    options.maxPEXSIIter = int_para[2];
    options.matrixType = int_para[3];
    options.isSymbolicFactorize = int_para[4];
    options.isConstructCommPattern = int_para[5];
    options.solver = int_para[6];
    options.symmetricStorage = int_para[7];
    options.ordering = int_para[8];
    options.rowOrdering = int_para[9];
    options.npSymbFact = int_para[10];
    options.symmetric = int_para[11];
    options.transpose = int_para[12];
    options.method = int_para[13];
    options.nPoints = int_para[14];
    options.verbosity = int_para[15];
    numProcessPerPole = int_para[16];

    options.spin = double_para[0];
    options.temperature = double_para[1];
    options.gap = double_para[2];
    options.deltaE = double_para[3];
    options.muMin0 = double_para[4];
    options.muMax0 = double_para[5];
    options.mu0 = double_para[6];
    options.muInertiaTolerance = double_para[7];
    options.muInertiaExpansion = double_para[8];
    options.muPEXSISafeGuard = double_para[9];
    options.numElectronPEXSITolerance = double_para[10];
    ZERO_Limit = double_para[11];

    return 0;
}

void splitNProc2NProwNPcol(const int NPROC, int& nprow, int& npcol)
{
    int integral_part = (int)sqrt(NPROC);
    if (NPROC % integral_part == 0)
    {
        nprow = integral_part;
        npcol = NPROC / integral_part;
    }
    else
    {
        int flag;
        int i;
        int low = pow(integral_part, 2);
        int high = pow(integral_part + 1, 2);
        if ((NPROC - low) >= (high - NPROC))
        {
            flag = integral_part + 1;
        }
        else
        {
            flag = integral_part;
        }
        for (i = flag; i > 0; ++i)
        {
            if (NPROC % i == 0)
                break;
        }
        nprow = i;
        npcol = NPROC / i;
    }
}

int simplePEXSI(MPI_Comm comm_PEXSI,
                MPI_Comm comm_2D,
                MPI_Group group_2D,
                const int blacs_ctxt, // communicator parameters
                const int size,
                const int nblk,
                const int nrow,
                const int ncol,
                char layout, // matrix parameters
                double* H,
                double* S, // input matrices
                const double numElectronExact,
                const std::string PexsiOptionFile, // pexsi parameters file
                double*& DM,
                double*& EDM, // output matrices
                double& totalEnergyH,
                double& totalEnergyS,
                double& totalFreeEnergy) // output energy
{

    if (comm_2D == MPI_COMM_NULL && comm_PEXSI == MPI_COMM_NULL)
        return 0;
    int myid;
    std::ofstream f_log;
    if (comm_PEXSI != MPI_COMM_NULL)
    {
        MPI_Comm_rank(comm_PEXSI, &myid);
    }

    //  set up PEXSI parameter
    PPEXSIOptions options;
    PPEXSISetDefaultOptions(&options);
    int numProcessPerPole;
    double ZERO_Limit;
    loadPEXSIOption(comm_PEXSI, PexsiOptionFile, options, numProcessPerPole, ZERO_Limit);

    ModuleBase::timer::tick("Diago_LCAO_Matrix", "setup_PEXSI_plan");
    PPEXSIPlan plan;
    int info;
    int outputFileIndex;
    int pexsi_prow, pexsi_pcol;
    ModuleBase::timer::tick("Diago_LCAO_Matrix", "splitNProc2NProwNPcol");
    splitNProc2NProwNPcol(numProcessPerPole, pexsi_prow, pexsi_pcol);
    ModuleBase::timer::tick("Diago_LCAO_Matrix", "splitNProc2NProwNPcol");

    outputFileIndex = -1;
    ModuleBase::timer::tick("Diago_LCAO_Matrix", "PEXSIPlanInit");
    if (comm_PEXSI != MPI_COMM_NULL)
    {
        plan = PPEXSIPlanInitialize(comm_PEXSI, pexsi_prow, pexsi_pcol, outputFileIndex, &info);
    }
    ModuleBase::timer::tick("Diago_LCAO_Matrix", "PEXSIPlanInit");
    
    ModuleBase::timer::tick("Diago_LCAO_Matrix", "setup_PEXSI_plan");

    // create compressed column storage distribution matrix parameter
    // LiuXh modify 2021-03-30, add DONE(ofs_running,"xx") for test
    // DONE(ofs_running,"create compressed column storage distribution matrix parameter, begin");
    DistCCSMatrix DST_Matrix(comm_PEXSI, numProcessPerPole, size);
    // LiuXh modify 2021-03-30, add DONE(ofs_running,"xx") for test
    // DONE(ofs_running,"create compressed column storage distribution matrix parameter, finish");


    // create block cyclic distribution matrix parameter
    DistBCDMatrix SRC_Matrix(comm_2D, group_2D, blacs_ctxt, size, nblk, nrow, ncol, layout);
    // LiuXh modify 2021-03-30, add DONE(ofs_running,"xx") for test
    // DONE(ofs_running,"create block cyclic distribution matrix parameter, finish");
    double* HnzvalLocal = nullptr;
    double* SnzvalLocal = nullptr;
    double* DMnzvalLocal = nullptr;
    double* EDMnzvalLocal = nullptr;
    double* FDMnzvalLocal = nullptr;
    // transform H and S from 2D block cyclic distribution to compressed column sparse matrix
    // LiuXh modify 2021-03-30, add DONE(ofs_running,"xx") for test
    DistMatrixTransformer::transformBCDtoCCS(SRC_Matrix, H, S, ZERO_Limit, DST_Matrix, HnzvalLocal, SnzvalLocal);
    // MPI_Barrier(MPI_COMM_WORLD);
    // LiuXh modify 2021-03-30, add DONE(ofs_running,"xx") for test
    if (comm_PEXSI != MPI_COMM_NULL)
    {

        // Load H and S to PEXSI
        int isSIdentity = 0;
        PPEXSILoadRealHSMatrix(plan,
                               options,
                               size,
                               DST_Matrix.get_nnz(),
                               DST_Matrix.get_nnzlocal(),
                               DST_Matrix.get_numcol_local(),
                               DST_Matrix.get_colptr_local(),
                               DST_Matrix.get_rowind_local(),
                               HnzvalLocal,
                               isSIdentity,
                               SnzvalLocal,
                               &info);

        double mu;
        double nelec;
        double muMinInertia;
        double muMaxInertia;
        int numTotalPEXSIIter;
        int numTotalInertiaIter; // Number of total inertia[out]
        // LiuXh modify 2021-04-29, add DONE(ofs_running,"xx") for test
        ModuleBase::timer::tick("Diago_LCAO_Matrix", "PEXSIDFT");
        PPEXSIDFTDriver(plan,                 // PEXSI plan[in]
                        options,              // PEXSI Options[in]
                        numElectronExact,     // exact electron number[in]
                        &mu,                  // chemical potential[out]
                        &nelec,               // number of electrons[out]
                        &muMinInertia,        // Lower bound for mu after the last inertia[out]
                        &muMaxInertia,        // Upper bound for mu after the last inertia[out]
                        &numTotalInertiaIter, // Number of total inertia[out]
                        &numTotalPEXSIIter,   // number of total pexsi evaluation procedure[out]
                        &info);               // 0: successful; otherwise: unsuccessful
        // LiuXh modify 2021-04-29, add DONE(ofs_running,"xx") for test
        ModuleBase::timer::tick("Diago_LCAO_Matrix", "PEXSIDFT");

        // retrieve the results from the plan
        if (DMnzvalLocal != nullptr)
            delete[] DMnzvalLocal;
        if (EDMnzvalLocal != nullptr)
            delete[] EDMnzvalLocal;
        if (FDMnzvalLocal != nullptr)
            delete[] FDMnzvalLocal;
        DMnzvalLocal = new double[DST_Matrix.get_nnzlocal()];
        EDMnzvalLocal = new double[DST_Matrix.get_nnzlocal()];
        FDMnzvalLocal = new double[DST_Matrix.get_nnzlocal()];
        if (myid < numProcessPerPole)
        {
            PPEXSIRetrieveRealDFTMatrix(plan,
                                        DMnzvalLocal,
                                        EDMnzvalLocal,
                                        FDMnzvalLocal,
                                        &totalEnergyH,
                                        &totalEnergyS,
                                        &totalFreeEnergy,
                                        &info);
        }
        // clean PEXSI
        PPEXSIPlanFinalize(plan, &info);
    }

    // transform Density Matrix and Energy Density Matrix from compressed column sparse matrix
    // back to 2D block cyclic distribution if neccessary
    if (comm_2D != MPI_COMM_NULL)
    {
        delete[] DM;
        delete[] EDM;
        DM = new double[SRC_Matrix.get_nrow() * SRC_Matrix.get_ncol()];
        EDM = new double[SRC_Matrix.get_nrow() * SRC_Matrix.get_ncol()];
    }
    // LiuXh modify 2021-04-29, add DONE(ofs_running,"xx") for test
    ModuleBase::timer::tick("Diago_LCAO_Matrix", "TransMAT22D");
    DistMatrixTransformer::transformCCStoBCD(DST_Matrix, DMnzvalLocal, EDMnzvalLocal, SRC_Matrix, DM, EDM);
    ModuleBase::timer::tick("Diago_LCAO_Matrix", "TransMAT22D");
    // LiuXh modify 2021-04-29, add DONE(ofs_running,"xx") for test

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    delete[] DMnzvalLocal;
    delete[] EDMnzvalLocal;
    delete[] FDMnzvalLocal;
    delete[] HnzvalLocal;
    delete[] SnzvalLocal;
    MPI_Barrier(MPI_COMM_WORLD);
    return 0;
}
} // namespace pexsi
#endif