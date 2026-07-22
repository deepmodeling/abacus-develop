#include "source_context/simulation_context_builder.h"

#include "source_base/global_variable.h"
#include "source_hamilt/module_xc/exx_info.h"
#include "source_io/module_parameter/input_parameter.h"
#include "source_io/module_parameter/system_parameter.h"
#include "source_io/module_restart/restart.h"

#include <gtest/gtest.h>

#include <set>
#include <stdexcept>
#include <string>

namespace GlobalC
{
Exx_Info exx_info;
Restart restart;
} // namespace GlobalC

namespace
{

class SimulationContextBuilderTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        GlobalV::NPROC = 1;
        GlobalV::KPAR = 1;
        GlobalV::MY_RANK = 0;
        GlobalV::MY_POOL = 0;
        GlobalV::MY_BNDGROUP = 0;
        GlobalV::NPROC_IN_POOL = 1;
        GlobalV::NPROC_IN_BNDGROUP = 1;
        GlobalV::RANK_IN_POOL = 0;
        GlobalV::RANK_IN_BPGROUP = 0;
        GlobalV::DRANK = 0;
        GlobalV::DSIZE = 1;
        GlobalV::DCOLOR = 0;
        GlobalV::GRANK = 0;
        GlobalV::GSIZE = 1;
        GlobalC::exx_info = Exx_Info();
        GlobalC::restart = Restart();
    }

    ModuleContext::SimulationContext freeze(Input_para& input, System_para& system)
    {
        ModuleContext::SimulationContextBuilder builder(input, system);
        builder.capture_runtime(system);
        return builder.finalize(system);
    }
};

TEST_F(SimulationContextBuilderTest, ProjectsEveryInputFieldExactlyOnce)
{
    Input_para input;
    input.calculation = "nscf";
    input.suffix = "projection-test";
    input.pseudo_dir = "pseudo";
    input.device = "cpu";
    input.nelec = 12.5;
    input.mixing_beta = 0.37;
    input.nspin = 4;
    input.ks_solver = "cg";
    input.lcao_dr = 0.02;
    input.relax_nmax = 9;
    input.ref_cell_factor = 1.5;
    input.of_method = "cg1";
    input.nbands_sto = 32;
    input.deepks_model = "model.pt";
    input.td_dt = 0.1;
    input.lr_nstates = 7;
    input.out_level = "m";
    input.out_mat_hs2 = {1, 12, 2};
    input.out_chg = {1, 6};
    input.restart_save = true;
    input.dos_sigma = 0.12;
    input.efield_amp = 0.3;
    input.vdw_method = "d4";
    input.exx_hybrid_step = 13;
    input.dft_plus_u = 1;
    input.sc_thr = 2.0e-6;
    input.qo_switch = true;
    input.pexsi_npole = 28;
    input.nurse = 3;
    input.rdmft_power_alpha = 0.7;
    input.xc_exch_ext = {101.0, 0.9};
    input.of_cd = true;
    input.dsp_count = 6;

    System_para system;
    const ModuleContext::SimulationContext context = freeze(input, system);

    std::set<std::string> mapped_names;
#define MODULE_CONTEXT_INPUT_FIELD(domain, field)                                                        \
    do                                                                                                   \
    {                                                                                                    \
        EXPECT_TRUE(mapped_names.insert(#field).second) << "duplicate INPUT field " << #field;          \
    } while (false);
#include "source_context/input_field_mapping.inc"
#undef MODULE_CONTEXT_INPUT_FIELD
    EXPECT_EQ(504u, mapped_names.size());

    EXPECT_EQ("projection-test", context.cell.value<std::string>("suffix"));
    EXPECT_EQ("pseudo", context.file_input.value<std::string>("pseudo_dir"));
    EXPECT_DOUBLE_EQ(0.37, context.scf_mixing.value<double>("mixing_beta"));
    EXPECT_EQ("model.pt", context.deepks.value<std::string>("deepks_model"));
    EXPECT_EQ("d4", context.vdw.value<std::string>("vdw_method"));
    EXPECT_EQ(6, context.uncommon_hardware.value<int>("dsp_count"));
    EXPECT_EQ(0u, context.run.mapped_values.count("calculation"));
    EXPECT_EQ(0u, context.spin.mapped_values.count("nspin"));
    EXPECT_EQ(0u, context.matrix_output.mapped_values.count("out_mat_hs2"));

    EXPECT_EQ("nscf", context.run.calculation);
    EXPECT_EQ(4, context.spin.nspin);
    EXPECT_TRUE(context.matrix_output.hs_r.enabled);
    EXPECT_EQ(12, context.matrix_output.hs_r.precision);
    EXPECT_EQ(2, context.matrix_output.hs_r.mode);
}

TEST_F(SimulationContextBuilderTest, RequiresRuntimeCaptureBeforeFreeze)
{
    Input_para input;
    System_para system;
    ModuleContext::SimulationContextBuilder builder(input, system);

    EXPECT_THROW(builder.finalize(system), std::logic_error);
    EXPECT_FALSE(builder.finalized());

    system.nlocal = 17;
    system.npol = 2;
    builder.capture_runtime(system);
    EXPECT_TRUE(builder.runtime_captured());
    const ModuleContext::SimulationContext context = builder.finalize(system);
    EXPECT_EQ(17, context.basis.nlocal);
    EXPECT_EQ(2, context.basis.npol);
    EXPECT_TRUE(builder.finalized());
    EXPECT_THROW(builder.finalize(system), std::logic_error);
    EXPECT_THROW(builder.capture_runtime(system), std::logic_error);
}

TEST_F(SimulationContextBuilderTest, SourceMappingManifestsAreCompleteAndUnique)
{
    System_para system;
    std::set<std::string> system_fields;
#define MODULE_CONTEXT_SYSTEM_FIELD(domain, field) \
    static_cast<void>(system.field);                \
    EXPECT_TRUE(system_fields.insert(#field).second) << "duplicate System_para field " << #field;
#include "source_context/system_field_mapping.inc"
#undef MODULE_CONTEXT_SYSTEM_FIELD
    EXPECT_EQ(45u, system_fields.size());

    std::set<std::string> globalv_fields;
#define MODULE_CONTEXT_GLOBALV_FIELD(domain, field) \
    static_cast<void>(GlobalV::field);              \
    EXPECT_TRUE(globalv_fields.insert(#field).second) << "duplicate GlobalV field " << #field;
#include "source_context/globalv_field_mapping.inc"
#undef MODULE_CONTEXT_GLOBALV_FIELD
    EXPECT_EQ(18u, globalv_fields.size());

    std::set<std::string> restart_fields;
#define MODULE_CONTEXT_RESTART_FIELD(field)                                                \
    static_cast<void>(GlobalC::restart.field);                                             \
    EXPECT_TRUE(restart_fields.insert(#field).second) << "duplicate restart field " << #field;
#include "source_context/restart_field_mapping.inc"
#undef MODULE_CONTEXT_RESTART_FIELD
    EXPECT_EQ(8u, restart_fields.size());
}

TEST_F(SimulationContextBuilderTest, RejectsInconsistentLegacyWorldTopology)
{
    Input_para input;
    System_para system;
    system.nproc = 1;
    GlobalV::NPROC = 2;

    ModuleContext::SimulationContextBuilder builder(input, system);
    builder.capture_runtime(system);
    EXPECT_THROW(builder.finalize(system), std::logic_error);
}

TEST_F(SimulationContextBuilderTest, UsesLegacyGlobalVAsTheRuntimePoolTopology)
{
    Input_para input;
    System_para system;
    system.npool = 2;
    system.nproc_in_pool = 1;
    GlobalV::KPAR = 3;
    GlobalV::MY_POOL = 1;
    GlobalV::NPROC_IN_POOL = 4;

    ModuleContext::SimulationContextBuilder builder(input, system);
    builder.capture_runtime(system);
    const ModuleContext::SimulationContext context = builder.finalize(system);

    EXPECT_EQ(3, context.parallel.kpoint_pool_count);
    EXPECT_EQ(1, context.parallel.pool_index);
    EXPECT_EQ(4, context.parallel.processes_in_pool);
}

TEST_F(SimulationContextBuilderTest, AdaptersUseTheSingleLegacyStateInstances)
{
    GlobalC::exx_info.info_global.cal_exx = true;
    GlobalC::exx_info.info_global.hybrid_alpha = 0.25;
    GlobalC::exx_info.info_global.hse_omega = 0.11;
    GlobalC::exx_info.info_ri.real_number = true;
    GlobalC::restart.info_save.save_charge = true;
    GlobalC::restart.info_load.load_H = true;
    GlobalC::restart.folder = "restart-a/";

    Input_para input;
    System_para system;
    const ModuleContext::SimulationContext context = freeze(input, system);

    const ModuleContext::ExactExchangeState exx = context.exact_exchange_state->snapshot();
    std::set<std::string> exx_names;
#define MODULE_CONTEXT_EXX_STATE_FIELD(section, field)                                      \
    do                                                                                       \
    {                                                                                        \
        const std::string name = #section "." #field;                                      \
        EXPECT_TRUE(exx_names.insert(name).second) << "duplicate EXX state field " << name; \
    } while (false);
#include "source_context/exx_state_field_mapping.inc"
#undef MODULE_CONTEXT_EXX_STATE_FIELD
    EXPECT_EQ(41u, exx_names.size());
    EXPECT_TRUE(exx.enabled);
    EXPECT_TRUE(exx.real_number);
    EXPECT_DOUBLE_EQ(0.25, exx.hybrid_alpha);
    EXPECT_DOUBLE_EQ(0.11, exx.hse_omega);
    EXPECT_EQ(0u, exx.mapped_values.count("info_global.hybrid_alpha"));
    EXPECT_EQ(1u, exx.mapped_values.count("info_ri.C_threshold"));

    ModuleContext::RestartState restart = context.restart->snapshot();
    EXPECT_TRUE(restart.save_charge);
    EXPECT_TRUE(restart.load_hamiltonian);
    EXPECT_EQ("restart-a/", restart.folder);

    restart.save_charge = false;
    restart.load_hamiltonian_finished = true;
    restart.restart_exact_exchange = true;
    restart.folder = "restart-b/";
    context.restart->update(restart);
    EXPECT_FALSE(GlobalC::restart.info_save.save_charge);
    EXPECT_TRUE(GlobalC::restart.info_load.load_H_finish);
    EXPECT_TRUE(GlobalC::restart.info_load.restart_exx);
    EXPECT_EQ("restart-b/", GlobalC::restart.folder);
}

} // namespace
