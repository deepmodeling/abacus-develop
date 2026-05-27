// use PEXSI to solve a Kohn-Sham equation
// the H and S matrices are given by 2D block cyclic distribution
#include "simple_pexsi.h"
#include "source_io/module_parameter/parameter.h"
// the Density Matrix and Energy Density Matrix calculated by PEXSI are transformed to 2D block cyclic distribution
// #include "mpi.h"
#ifdef __PEXSI
#include <mpi.h>

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>

#include "c_pexsi_interface.h"
#include "dist_bcd_matrix.h"
#include "dist_ccs_matrix.h"
#include "dist_matrix_transformer.h"
#include "source_base/module_external/lapack_connector.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#include "source_base/global_variable.h"
#include "pexsi_solver.h"

namespace pexsi
{
namespace
{
bool trace_pexsi_enabled(const char* category)
{
    const char* trace_env = std::getenv("ABACUS_PEXSI_TRACE");
    if (trace_env == nullptr)
    {
        return false;
    }
    const std::string trace_value(trace_env);
    return trace_value == "1" || trace_value == "all" || trace_value.find(category) != std::string::npos;
}

void trace_pexsi_stage(MPI_Comm comm, const char* path, const char* stage)
{
    if (!trace_pexsi_enabled("stage") || comm == MPI_COMM_NULL)
    {
        return;
    }
    int comm_rank = 0;
    int world_rank = 0;
    MPI_Comm_rank(comm, &comm_rank);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    if (comm_rank == 0)
    {
        std::cout << "PEXSI_TRACE world_rank=" << world_rank << " path=" << path << " stage=" << stage
                  << " time=" << MPI_Wtime() << std::endl;
    }
}

std::uint64_t fnv1a_append_int(std::uint64_t hash, const int value)
{
    std::uint32_t data = static_cast<std::uint32_t>(value);
    for (int byte = 0; byte < 4; ++byte)
    {
        hash ^= static_cast<unsigned char>((data >> (byte * 8)) & 0xffU);
        hash *= 1099511628211ULL;
    }
    return hash;
}

std::uint64_t local_pattern_hash(const int numColLocal,
                                 const int nnzLocal,
                                 const int* colptrLocal,
                                 const int* rowindLocal)
{
    std::uint64_t hash = 1469598103934665603ULL;
    hash = fnv1a_append_int(hash, numColLocal);
    hash = fnv1a_append_int(hash, nnzLocal);
    for (int i = 0; i < numColLocal + 1; ++i)
    {
        hash = fnv1a_append_int(hash, colptrLocal[i]);
    }
    for (int i = 0; i < nnzLocal; ++i)
    {
        hash = fnv1a_append_int(hash, rowindLocal[i]);
    }
    return hash;
}

void print_effective_pexsi_options_once(const PPEXSIOptions& options,
                                        const int numProcessPerPole,
                                        const double zeroLimit)
{
    static bool printed = false;
    if (printed || GlobalV::MY_RANK != 0)
    {
        return;
    }
    printed = true;
    GlobalV::ofs_running << "\n PEXSI effective options:"
                         << " numPole=" << options.numPole
                         << " temperature=" << options.temperature
                         << " numProcessPerPole=" << numProcessPerPole
                         << " muGuard=" << options.muPEXSISafeGuard
                         << " electronTolerance=" << options.numElectronPEXSITolerance
                         << " zeroLimit=" << zeroLimit << std::endl;
}
} // namespace

struct PexsiComplexReuseContextImpl
{
    struct Signature
    {
        int size = 0;
        int numProcessPerPole = 0;
        int pexsi_prow = 0;
        int pexsi_pcol = 0;
        int nnz = 0;
        int nnzLocal = 0;
        int numColLocal = 0;
        int matrixType = 0;
        int symmetric = 0;
        int symmetricStorage = 0;
        int ordering = 0;
        int rowOrdering = 0;
        int solver = 0;
        int npSymbFact = 0;
        int transpose = 0;
        std::uint64_t localPatternHash = 0;
        std::uint64_t globalPatternHash = 0;
    };

    PPEXSIPlan plan = 0;
    MPI_Comm comm = MPI_COMM_NULL;
    Signature signature;
    bool has_plan = false;

    ~PexsiComplexReuseContextImpl()
    {
        clear();
    }

    void clear()
    {
        int initialized = 0;
        int finalized = 0;
        MPI_Initialized(&initialized);
        if (initialized)
        {
            MPI_Finalized(&finalized);
        }
        if (initialized && !finalized)
        {
            if (has_plan)
            {
                int info = 0;
                trace_pexsi_stage(comm, "complex", "PPEXSIPlanFinalize.cached.begin");
                PPEXSIPlanFinalize(plan, &info);
                trace_pexsi_stage(comm, "complex", "PPEXSIPlanFinalize.cached.end");
            }
            if (comm != MPI_COMM_NULL)
            {
                MPI_Comm_free(&comm);
            }
        }
        plan = 0;
        comm = MPI_COMM_NULL;
        signature = Signature();
        has_plan = false;
    }

    static Signature make_signature(const PPEXSIOptions& options,
                                    const int size,
                                    const int numProcessPerPole,
                                    const int pexsi_prow,
                                    const int pexsi_pcol,
                                    const DistCCSMatrix& matrix,
                                    const std::uint64_t localPatternHash,
                                    const std::uint64_t globalPatternHash)
    {
        Signature sig;
        sig.size = size;
        sig.numProcessPerPole = numProcessPerPole;
        sig.pexsi_prow = pexsi_prow;
        sig.pexsi_pcol = pexsi_pcol;
        sig.nnz = matrix.get_nnz();
        sig.nnzLocal = matrix.get_nnzlocal();
        sig.numColLocal = matrix.get_numcol_local();
        sig.matrixType = options.matrixType;
        sig.symmetric = options.symmetric;
        sig.symmetricStorage = options.symmetricStorage;
        sig.ordering = options.ordering;
        sig.rowOrdering = options.rowOrdering;
        sig.solver = options.solver;
        sig.npSymbFact = options.npSymbFact;
        sig.transpose = options.transpose;
        sig.localPatternHash = localPatternHash;
        sig.globalPatternHash = globalPatternHash;
        return sig;
    }

    static bool same_signature(const Signature& lhs, const Signature& rhs)
    {
        return lhs.size == rhs.size && lhs.numProcessPerPole == rhs.numProcessPerPole
               && lhs.pexsi_prow == rhs.pexsi_prow && lhs.pexsi_pcol == rhs.pexsi_pcol
               && lhs.nnz == rhs.nnz && lhs.nnzLocal == rhs.nnzLocal
               && lhs.numColLocal == rhs.numColLocal && lhs.matrixType == rhs.matrixType
               && lhs.symmetric == rhs.symmetric && lhs.symmetricStorage == rhs.symmetricStorage
               && lhs.ordering == rhs.ordering && lhs.rowOrdering == rhs.rowOrdering
               && lhs.solver == rhs.solver && lhs.npSymbFact == rhs.npSymbFact
               && lhs.transpose == rhs.transpose && lhs.localPatternHash == rhs.localPatternHash
               && lhs.globalPatternHash == rhs.globalPatternHash;
    }

    bool comm_matches(MPI_Comm new_comm) const
    {
        if (!has_plan || comm == MPI_COMM_NULL || new_comm == MPI_COMM_NULL)
        {
            return false;
        }
        int result = MPI_UNEQUAL;
        MPI_Comm_compare(comm, new_comm, &result);
        if (result == MPI_IDENT || result == MPI_CONGRUENT)
        {
            return true;
        }
        if (result == MPI_SIMILAR)
        {
            int old_rank = -1;
            int new_rank = -2;
            MPI_Comm_rank(comm, &old_rank);
            MPI_Comm_rank(new_comm, &new_rank);
            return old_rank == new_rank;
        }
        return false;
    }

    PPEXSIPlan get_or_create_plan(MPI_Comm comm_in,
                                  const int pexsi_prow,
                                  const int pexsi_pcol,
                                  const int outputFileIndex,
                                  const Signature& new_signature,
                                  bool& reused)
    {
        reused = false;
        if (comm_in == MPI_COMM_NULL)
        {
            return 0;
        }

        const bool local_match = has_plan && comm_matches(comm_in) && same_signature(signature, new_signature);
        int local_reuse = local_match ? 1 : 0;
        int global_reuse = 0;
        MPI_Allreduce(&local_reuse, &global_reuse, 1, MPI_INT, MPI_MIN, comm_in);
        reused = global_reuse == 1;
        if (reused)
        {
            return plan;
        }

        clear();
        MPI_Comm_dup(comm_in, &comm);
        int info = 0;
        trace_pexsi_stage(comm, "complex", "PPEXSIPlanInitialize.cached.begin");
        plan = PPEXSIPlanInitialize(comm, pexsi_prow, pexsi_pcol, outputFileIndex, &info);
        trace_pexsi_stage(comm, "complex", "PPEXSIPlanInitialize.cached.end");
        signature = new_signature;
        has_plan = true;
        return plan;
    }
};

PexsiComplexReuseContext::PexsiComplexReuseContext()
    : impl_(new PexsiComplexReuseContextImpl())
{
}

PexsiComplexReuseContext::~PexsiComplexReuseContext() = default;

PexsiComplexReuseContext::PexsiComplexReuseContext(PexsiComplexReuseContext&&) noexcept = default;

PexsiComplexReuseContext& PexsiComplexReuseContext::operator=(PexsiComplexReuseContext&&) noexcept = default;

void PexsiComplexReuseContext::clear()
{
    if (impl_)
    {
        impl_->clear();
    }
}

inline void strtolower(char* sa, char* sb)
{
    char c = '\0';
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
    int_para[0] = pexsi::PEXSI_Solver::pexsi_npole;
    int_para[1] = pexsi::PEXSI_Solver::pexsi_inertia;
    int_para[2] = pexsi::PEXSI_Solver::pexsi_nmax;
    int_para[3] = 0;
    int_para[4] = 1; // pexsi::PEXSI_Solver::pexsi_symbolic;
    int_para[5] = pexsi::PEXSI_Solver::pexsi_comm;
    int_para[6] = 0;
    int_para[7] = pexsi::PEXSI_Solver::pexsi_storage;
    int_para[8] = pexsi::PEXSI_Solver::pexsi_ordering;
    int_para[9] = pexsi::PEXSI_Solver::pexsi_row_ordering;
    int_para[10] = pexsi::PEXSI_Solver::pexsi_nproc;
    int_para[11] = pexsi::PEXSI_Solver::pexsi_symm;
    int_para[12] = pexsi::PEXSI_Solver::pexsi_trans;
    int_para[13] = pexsi::PEXSI_Solver::pexsi_method;
    int_para[14] = 1;
    int_para[15] = 0;
    int_para[16] = pexsi::PEXSI_Solver::pexsi_nproc_pole;

    double_para[0] = 2.0;
    double_para[1] = pexsi::PEXSI_Solver::pexsi_temp;
    double_para[2] = pexsi::PEXSI_Solver::pexsi_gap;
    double_para[3] = pexsi::PEXSI_Solver::pexsi_delta_e;
    double_para[4] = pexsi::PEXSI_Solver::pexsi_mu_lower;
    double_para[5] = pexsi::PEXSI_Solver::pexsi_mu_upper;
    double_para[6] = pexsi::PEXSI_Solver::pexsi_mu;
    double_para[7] = pexsi::PEXSI_Solver::pexsi_mu_thr;
    double_para[8] = pexsi::PEXSI_Solver::pexsi_mu_expand;
    double_para[9] = pexsi::PEXSI_Solver::pexsi_mu_guard;
    double_para[10] = pexsi::PEXSI_Solver::pexsi_elec_thr;
    double_para[11] = pexsi::PEXSI_Solver::pexsi_zero_thr;

    const bool default_pexsi_temperature = std::abs(double_para[1] - 0.015) < 1.0e-14;
    const bool fermi_dirac_smearing = PARAM.inp.smearing_method == "fd"
                                      || PARAM.inp.smearing_method == "fermi-dirac";
    if (default_pexsi_temperature && fermi_dirac_smearing)
    {
        double_para[1] = PARAM.inp.smearing_sigma;
    }

    int comm_size = 1;
    if (comm != MPI_COMM_NULL)
    {
        MPI_Comm_size(comm, &comm_size);
    }
    int_para[16] = std::max(1, std::min(int_para[16], comm_size));

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
        int flag = 0;
        int i = 0;
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
                double& totalFreeEnergy, // output energy
                double& mu,
                double mu0)
{

    if (comm_2D == MPI_COMM_NULL && comm_PEXSI == MPI_COMM_NULL)
        return 0;
    int myid = 0;
    std::ofstream f_log;
    if (comm_PEXSI != MPI_COMM_NULL)
    {
        MPI_Comm_rank(comm_PEXSI, &myid);
    }

    //  set up PEXSI parameter
    PPEXSIOptions options;
    PPEXSISetDefaultOptions(&options);
    int numProcessPerPole = 0;
    double ZERO_Limit = 0.0;
    loadPEXSIOption(comm_PEXSI, PexsiOptionFile, options, numProcessPerPole, ZERO_Limit);
    options.mu0 = mu0;
    print_effective_pexsi_options_once(options, numProcessPerPole, ZERO_Limit);
    const bool detail_timing = std::getenv("ABACUS_PEXSI_REAL_DETAIL_TIMING") != nullptr;

    ModuleBase::timer::start("Diago_LCAO_Matrix", "setup_PEXSI_plan");
    PPEXSIPlan plan;
    int info = 0;
    int outputFileIndex = 0;
    int pexsi_prow, pexsi_pcol;
    ModuleBase::timer::start("Diago_LCAO_Matrix", "splitNProc2NProwNPcol");
    splitNProc2NProwNPcol(numProcessPerPole, pexsi_prow, pexsi_pcol);
    ModuleBase::timer::end("Diago_LCAO_Matrix", "splitNProc2NProwNPcol");

    outputFileIndex = -1;
    ModuleBase::timer::start("Diago_LCAO_Matrix", "PEXSIPlanInit");
    if (comm_PEXSI != MPI_COMM_NULL)
    {
        plan = PPEXSIPlanInitialize(comm_PEXSI, pexsi_prow, pexsi_pcol, outputFileIndex, &info);
    }
    ModuleBase::timer::end("Diago_LCAO_Matrix", "PEXSIPlanInit");
    
    ModuleBase::timer::end("Diago_LCAO_Matrix", "setup_PEXSI_plan");

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
    if (detail_timing)
    {
        ModuleBase::timer::start("Diago_LCAO_Matrix", "TransMAT2CCS");
    }
    DistMatrixTransformer::transformBCDtoCCS(SRC_Matrix, H, S, ZERO_Limit, DST_Matrix, HnzvalLocal, SnzvalLocal);
    if (detail_timing)
    {
        ModuleBase::timer::end("Diago_LCAO_Matrix", "TransMAT2CCS");
    }
    // MPI_Barrier(MPI_COMM_WORLD);
    // LiuXh modify 2021-03-30, add DONE(ofs_running,"xx") for test
    if (comm_PEXSI != MPI_COMM_NULL)
    {

        // Load H and S to PEXSI
        int isSIdentity = 0;
        if (detail_timing)
        {
            ModuleBase::timer::start("Diago_LCAO_Matrix", "PEXSI_LoadHS_R");
        }
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
        if (detail_timing)
        {
            ModuleBase::timer::end("Diago_LCAO_Matrix", "PEXSI_LoadHS_R");
        }

        double nelec = 0.0;
        double muMinInertia = 0.0;
        double muMaxInertia = 0.0;
        int numTotalPEXSIIter = 0;
        int numTotalInertiaIter = 0; // Number of total inertia[out]
        // LiuXh modify 2021-04-29, add DONE(ofs_running,"xx") for test
        ModuleBase::timer::start("Diago_LCAO_Matrix", "PEXSIDFT");
        PPEXSIDFTDriver2(plan,                 // PEXSI plan[in]
                        &options,              // PEXSI Options[in]
                        numElectronExact,     // exact electron number[in]
                        &mu,                  // chemical potential[out]
                        &nelec,               // number of electrons[out]
                        // &muMinInertia,        // Lower bound for mu after the last inertia[out]
                        // &muMaxInertia,        // Upper bound for mu after the last inertia[out]
                        &numTotalInertiaIter, // Number of total inertia[out]
                        // &numTotalPEXSIIter,   // number of total pexsi evaluation procedure[out]
                        &info);               // 0: successful; otherwise: unsuccessful
        // LiuXh modify 2021-04-29, add DONE(ofs_running,"xx") for test
        ModuleBase::timer::end("Diago_LCAO_Matrix", "PEXSIDFT");

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
            if (detail_timing)
            {
                ModuleBase::timer::start("Diago_LCAO_Matrix", "PEXSI_Retrieve_R");
            }
            PPEXSIRetrieveRealDFTMatrix(plan,
                                        DMnzvalLocal,
                                        EDMnzvalLocal,
                                        FDMnzvalLocal,
                                        &totalEnergyH,
                                        &totalEnergyS,
                                        &totalFreeEnergy,
                                        &info);
            if (detail_timing)
            {
                ModuleBase::timer::end("Diago_LCAO_Matrix", "PEXSI_Retrieve_R");
            }
        }
        // clean PEXSI
        if (detail_timing)
        {
            ModuleBase::timer::start("Diago_LCAO_Matrix", "PEXSI_Finalize");
        }
        PPEXSIPlanFinalize(plan, &info);
        if (detail_timing)
        {
            ModuleBase::timer::end("Diago_LCAO_Matrix", "PEXSI_Finalize");
        }
    }

    // transform Density Matrix and Energy Density Matrix from compressed column sparse matrix
    // back to 2D block cyclic distribution if neccessary
    if (comm_2D != MPI_COMM_NULL)
    {
        // delete[] DM;
        // delete[] EDM;
        // DM = new double[SRC_Matrix.get_nrow() * SRC_Matrix.get_ncol()];
        // EDM = new double[SRC_Matrix.get_nrow() * SRC_Matrix.get_ncol()];
    }
    // LiuXh modify 2021-04-29, add DONE(ofs_running,"xx") for test
    const MPI_Comm barrier_comm = comm_2D != MPI_COMM_NULL ? comm_2D : comm_PEXSI;
    ModuleBase::timer::start("Diago_LCAO_Matrix", "TransMAT22D");
    DistMatrixTransformer::transformCCStoBCD(DST_Matrix, DMnzvalLocal, EDMnzvalLocal, SRC_Matrix, DM, EDM);
    ModuleBase::timer::end("Diago_LCAO_Matrix", "TransMAT22D");
    // LiuXh modify 2021-04-29, add DONE(ofs_running,"xx") for test

    if (barrier_comm != MPI_COMM_NULL)
    {
        MPI_Barrier(barrier_comm);
        MPI_Barrier(barrier_comm);
    }
    delete[] DMnzvalLocal;
    delete[] EDMnzvalLocal;
    delete[] FDMnzvalLocal;
    delete[] HnzvalLocal;
    delete[] SnzvalLocal;
    if (barrier_comm != MPI_COMM_NULL)
    {
        MPI_Barrier(barrier_comm);
    }
    return 0;
}

int simplePEXSIComplex(MPI_Comm comm_PEXSI,
                       MPI_Comm comm_2D,
                       MPI_Group group_2D,
                       const int blacs_ctxt,
                       const int size,
                       const int nblk,
                       const int nrow,
                       const int ncol,
                       char layout,
                       std::complex<double>* H,
                       std::complex<double>* S,
                       const double numElectronExact,
                       const std::string PexsiOptionFile,
                       std::complex<double>*& DM,
                       std::complex<double>*& EDM,
                       double& totalEnergyH,
                       double& totalEnergyS,
                       double& totalFreeEnergy,
                       double& mu,
                       double mu0,
                       double* numElectronPEXSI,
                       double* numElectronDrvMuPEXSI,
                       PexsiComplexReuseContext* reuse_context)
{
    if (comm_2D == MPI_COMM_NULL && comm_PEXSI == MPI_COMM_NULL)
    {
        return 0;
    }

    int myid = 0;
    if (comm_PEXSI != MPI_COMM_NULL)
    {
        MPI_Comm_rank(comm_PEXSI, &myid);
    }

    PPEXSIOptions options;
    PPEXSISetDefaultOptions(&options);
    int numProcessPerPole = 0;
    double ZERO_Limit = 0.0;
    loadPEXSIOption(comm_PEXSI, PexsiOptionFile, options, numProcessPerPole, ZERO_Limit);
    options.symmetric = 1;
    options.symmetricStorage = 0;
    options.rowOrdering = 0;
    options.spin = PARAM.inp.nspin == 1 ? 2.0 : 1.0;
    options.mu0 = mu0;
    print_effective_pexsi_options_once(options, numProcessPerPole, ZERO_Limit);

    ModuleBase::timer::start("Diago_LCAO_Matrix", "setup_PEXSI_plan");
    PPEXSIPlan plan = 0;
    int info = 0;
    int outputFileIndex = -1;
    int pexsi_prow = 0;
    int pexsi_pcol = 0;
    ModuleBase::timer::start("Diago_LCAO_Matrix", "splitNProc2NProwNPcol");
    splitNProc2NProwNPcol(numProcessPerPole, pexsi_prow, pexsi_pcol);
    ModuleBase::timer::end("Diago_LCAO_Matrix", "splitNProc2NProwNPcol");

    ModuleBase::timer::end("Diago_LCAO_Matrix", "setup_PEXSI_plan");

    DistCCSMatrix DST_Matrix(comm_PEXSI, numProcessPerPole, size);
    DistBCDMatrix SRC_Matrix(comm_2D, group_2D, blacs_ctxt, size, nblk, nrow, ncol, layout);

    std::complex<double>* HnzvalLocal = nullptr;
    std::complex<double>* SnzvalLocal = nullptr;
    std::complex<double>* DMnzvalLocal = nullptr;
    std::complex<double>* EDMnzvalLocal = nullptr;
    std::complex<double>* FDMnzvalLocal = nullptr;

    ModuleBase::timer::start("Diago_LCAO_Matrix", "TransMAT2CCS");
    trace_pexsi_stage(comm_PEXSI, "complex", "TransMAT2CCS.begin");
    DistMatrixTransformer::transformBCDtoCCS(SRC_Matrix, H, S, ZERO_Limit, DST_Matrix, HnzvalLocal, SnzvalLocal);
    trace_pexsi_stage(comm_PEXSI, "complex", "TransMAT2CCS.end");
    ModuleBase::timer::end("Diago_LCAO_Matrix", "TransMAT2CCS");

    if (comm_PEXSI != MPI_COMM_NULL)
    {
        bool reuse_plan = false;
        ModuleBase::timer::start("Diago_LCAO_Matrix", "PEXSIPlanInit");
        if (reuse_context != nullptr)
        {
            const std::uint64_t local_pattern_hash_value = local_pattern_hash(DST_Matrix.get_numcol_local(),
                                                                               DST_Matrix.get_nnzlocal(),
                                                                               DST_Matrix.get_colptr_local(),
                                                                               DST_Matrix.get_rowind_local());
            std::uint64_t global_pattern_hash_value = local_pattern_hash_value;
            MPI_Allreduce(&local_pattern_hash_value,
                          &global_pattern_hash_value,
                          1,
                          MPI_UINT64_T,
                          MPI_BXOR,
                          comm_PEXSI);
            const auto signature = PexsiComplexReuseContextImpl::make_signature(options,
                                                                                size,
                                                                                numProcessPerPole,
                                                                                pexsi_prow,
                                                                                pexsi_pcol,
                                                                                DST_Matrix,
                                                                                local_pattern_hash_value,
                                                                                global_pattern_hash_value);
            plan = reuse_context->impl_->get_or_create_plan(comm_PEXSI,
                                                            pexsi_prow,
                                                            pexsi_pcol,
                                                            outputFileIndex,
                                                            signature,
                                                            reuse_plan);
            if (reuse_plan)
            {
                trace_pexsi_stage(comm_PEXSI, "complex", "PPEXSIPlanInitialize.reuse");
            }
        }
        else
        {
            trace_pexsi_stage(comm_PEXSI, "complex", "PPEXSIPlanInitialize.begin");
            plan = PPEXSIPlanInitialize(comm_PEXSI, pexsi_prow, pexsi_pcol, outputFileIndex, &info);
            trace_pexsi_stage(comm_PEXSI, "complex", "PPEXSIPlanInitialize.end");
        }
        ModuleBase::timer::end("Diago_LCAO_Matrix", "PEXSIPlanInit");

        int isSIdentity = 0;
        PPEXSIOptions load_options = options;
        if (reuse_plan)
        {
            load_options.isSymbolicFactorize = 0;
            load_options.isConstructCommPattern = 0;
        }
        ModuleBase::timer::start("Diago_LCAO_Matrix", "PEXSI_LoadHS_C");
        trace_pexsi_stage(comm_PEXSI, "complex", "PPEXSILoadComplexHSMatrix.begin");
        PPEXSILoadComplexHSMatrix(plan,
                                  load_options,
                                  size,
                                  DST_Matrix.get_nnz(),
                                  DST_Matrix.get_nnzlocal(),
                                  DST_Matrix.get_numcol_local(),
                                  DST_Matrix.get_colptr_local(),
                                  DST_Matrix.get_rowind_local(),
                                  reinterpret_cast<double*>(HnzvalLocal),
                                  isSIdentity,
                                  reinterpret_cast<double*>(SnzvalLocal),
                                  &info);
        trace_pexsi_stage(comm_PEXSI, "complex", "PPEXSILoadComplexHSMatrix.end");
        ModuleBase::timer::end("Diago_LCAO_Matrix", "PEXSI_LoadHS_C");
        if (!reuse_plan)
        {
            ModuleBase::timer::start("Diago_LCAO_Matrix", "PEXSI_Symbolic_C");
            trace_pexsi_stage(comm_PEXSI, "complex", "PPEXSISymbolicFactorizeComplexUnsymmetricMatrix.begin");
            PPEXSISymbolicFactorizeComplexUnsymmetricMatrix(plan,
                                                            options,
                                                            reinterpret_cast<double*>(HnzvalLocal),
                                                            &info);
            trace_pexsi_stage(comm_PEXSI, "complex", "PPEXSISymbolicFactorizeComplexUnsymmetricMatrix.end");
            ModuleBase::timer::end("Diago_LCAO_Matrix", "PEXSI_Symbolic_C");
        }
        else
        {
            trace_pexsi_stage(comm_PEXSI, "complex", "PPEXSISymbolicFactorizeComplexUnsymmetricMatrix.reuse");
        }

        double nelec = 0.0;
        double d_nelec_d_mu = 0.0;
        mu = mu0;
        PPEXSIOptions fermi_options = options;
        if (reuse_plan)
        {
            fermi_options.isSymbolicFactorize = 0;
            fermi_options.isConstructCommPattern = 0;
        }
        ModuleBase::timer::start("Diago_LCAO_Matrix", "PEXSIDFT");
        ModuleBase::timer::start("Diago_LCAO_Matrix", "PEXSI_FermiOp_C");
        trace_pexsi_stage(comm_PEXSI, "complex", "PPEXSICalculateFermiOperatorComplex.begin");
        PPEXSICalculateFermiOperatorComplex(plan,
                                            fermi_options,
                                            mu,
                                            numElectronExact,
                                            &nelec,
                                            &d_nelec_d_mu,
                                            &info);
        trace_pexsi_stage(comm_PEXSI, "complex", "PPEXSICalculateFermiOperatorComplex.end");
        ModuleBase::timer::end("Diago_LCAO_Matrix", "PEXSI_FermiOp_C");
        if (numElectronPEXSI != nullptr)
        {
            *numElectronPEXSI = nelec;
        }
        if (numElectronDrvMuPEXSI != nullptr)
        {
            *numElectronDrvMuPEXSI = d_nelec_d_mu;
        }
        ModuleBase::timer::end("Diago_LCAO_Matrix", "PEXSIDFT");

        DMnzvalLocal = new std::complex<double>[DST_Matrix.get_nnzlocal()];
        EDMnzvalLocal = new std::complex<double>[DST_Matrix.get_nnzlocal()];
        FDMnzvalLocal = new std::complex<double>[DST_Matrix.get_nnzlocal()];
        if (myid < numProcessPerPole)
        {
            ModuleBase::timer::start("Diago_LCAO_Matrix", "PEXSI_Retrieve_C");
            trace_pexsi_stage(comm_PEXSI, "complex", "PPEXSIRetrieveComplexDFTMatrix.begin");
            PPEXSIRetrieveComplexDFTMatrix(plan,
                                           reinterpret_cast<double*>(DMnzvalLocal),
                                           reinterpret_cast<double*>(EDMnzvalLocal),
                                           reinterpret_cast<double*>(FDMnzvalLocal),
                                           &totalEnergyH,
                                           &totalEnergyS,
                                           &totalFreeEnergy,
                                           &info);
            trace_pexsi_stage(comm_PEXSI, "complex", "PPEXSIRetrieveComplexDFTMatrix.end");
            ModuleBase::timer::end("Diago_LCAO_Matrix", "PEXSI_Retrieve_C");
        }
        if (reuse_context == nullptr)
        {
            ModuleBase::timer::start("Diago_LCAO_Matrix", "PEXSI_Finalize");
            trace_pexsi_stage(comm_PEXSI, "complex", "PPEXSIPlanFinalize.begin");
            PPEXSIPlanFinalize(plan, &info);
            trace_pexsi_stage(comm_PEXSI, "complex", "PPEXSIPlanFinalize.end");
            ModuleBase::timer::end("Diago_LCAO_Matrix", "PEXSI_Finalize");
        }
    }

    const MPI_Comm barrier_comm = comm_2D != MPI_COMM_NULL ? comm_2D : comm_PEXSI;
    ModuleBase::timer::start("Diago_LCAO_Matrix", "TransMAT22D");
    trace_pexsi_stage(barrier_comm, "complex", "TransMAT22D.begin");
    DistMatrixTransformer::transformCCStoBCD(DST_Matrix, DMnzvalLocal, EDMnzvalLocal, SRC_Matrix, DM, EDM);
    trace_pexsi_stage(barrier_comm, "complex", "TransMAT22D.end");
    ModuleBase::timer::end("Diago_LCAO_Matrix", "TransMAT22D");

    if (barrier_comm != MPI_COMM_NULL)
    {
        ModuleBase::timer::start("Diago_LCAO_Matrix", "PEXSI_Barrier");
        MPI_Barrier(barrier_comm);
        MPI_Barrier(barrier_comm);
        ModuleBase::timer::end("Diago_LCAO_Matrix", "PEXSI_Barrier");
    }
    delete[] DMnzvalLocal;
    delete[] EDMnzvalLocal;
    delete[] FDMnzvalLocal;
    delete[] HnzvalLocal;
    delete[] SnzvalLocal;
    if (barrier_comm != MPI_COMM_NULL)
    {
        ModuleBase::timer::start("Diago_LCAO_Matrix", "PEXSI_Barrier");
        MPI_Barrier(barrier_comm);
        ModuleBase::timer::end("Diago_LCAO_Matrix", "PEXSI_Barrier");
    }
    return 0;
}
} // namespace pexsi
#endif
