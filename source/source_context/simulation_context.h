#ifndef ABACUS_SOURCE_CONTEXT_SIMULATION_CONTEXT_H
#define ABACUS_SOURCE_CONTEXT_SIMULATION_CONTEXT_H

#include "source_context/context_types.h"

namespace ModuleContext
{

struct SimulationContext
{
    RunControl run;
    FileSystemLayout files;
    ParallelTopology parallel;
    LogStreams logs;
    BasisInfo basis;
    GridInfo grid;
    ElectronicRuntimeInfo electronic_runtime;
    FeatureRuntimeInfo features;
    SpinConfig spin;
    SolverConfig solver;
    MatrixOutputConfig matrix_output;
    GeneralOutputConfig output;
    DftUConfig dftu;
    ExactExchangeConfig exact_exchange;

    std::shared_ptr<ExactExchangeStateService> exact_exchange_state;
    std::shared_ptr<RestartService> restart;

    // Full input-domain type skeleton. The root aggregate is restricted to
    // Driver/ESolver orchestration and must not be passed into leaf modules.
    CellConfig cell;
    FileInputConfig file_input;
    ParallelConfig parallel_input;
    DeviceConfig device;
    ElectronicStructureConfig electronic_structure;
    ScfMixingConfig scf_mixing;
    LcaoConfig lcao;
    RelaxationConfig relaxation;
    MolecularDynamicsConfig molecular_dynamics;
    OfdftConfig ofdft;
    StochasticDftConfig stochastic_dft;
    DeepksConfig deepks;
    RealTimeTddftConfig real_time_tddft;
    LinearResponseTddftConfig linear_response_tddft;
    ChargeWavefunctionOutputConfig charge_wavefunction_output;
    RestartOutputConfig restart_output;
    PostprocessConfig postprocess;
    ModelConfig model;
    VdwConfig vdw;
    SpinConstraintConfig spin_constraint;
    QuasiAtomicOrbitalConfig quasi_atomic_orbital;
    PexsiConfig pexsi;
    TestConfig test;
    RdmftConfig rdmft;
    ExternalXcConfig external_xc;
    TimeDependentOfdftConfig time_dependent_ofdft;
    UncommonHardwareConfig uncommon_hardware;
};

} // namespace ModuleContext

#endif
