#include "source_context/simulation_context_builder.h"

#include "source_base/global_variable.h"
#include "source_hamilt/module_xc/exx_info.h"
#include "source_io/module_parameter/input_parameter.h"
#include "source_io/module_parameter/system_parameter.h"
#include "source_io/module_restart/restart.h"

#include <stdexcept>

namespace ModuleContext
{
namespace
{

MatrixOutputRequest request_from(const std::vector<int>& value)
{
    MatrixOutputRequest request;
    if (!value.empty())
    {
        request.enabled = value[0] != 0;
    }
    if (value.size() > 1)
    {
        request.precision = value[1];
    }
    if (value.size() > 2)
    {
        request.mode = value[2];
    }
    return request;
}

class LegacyExactExchangeStateService : public ExactExchangeStateService
{
  public:
    ExactExchangeState snapshot() const override
    {
        ExactExchangeState state;
#define MODULE_CONTEXT_EXX_STATE_FIELD(section, field) \
        state.project(#section "." #field, GlobalC::exx_info.section.field);
#include "source_context/exx_state_field_mapping.inc"
#undef MODULE_CONTEXT_EXX_STATE_FIELD
        state.enabled = GlobalC::exx_info.info_global.cal_exx;
        state.real_number = GlobalC::exx_info.info_ri.real_number;
        state.hybrid_alpha = GlobalC::exx_info.info_global.hybrid_alpha;
        state.hse_omega = GlobalC::exx_info.info_global.hse_omega;
        return state;
    }
};

class LegacyRestartService : public RestartService
{
  public:
    RestartState snapshot() const override
    {
        RestartState state;
        state.save_charge = GlobalC::restart.info_save.save_charge;
        state.save_hamiltonian = GlobalC::restart.info_save.save_H;
        state.load_charge = GlobalC::restart.info_load.load_charge;
        state.load_charge_finished = GlobalC::restart.info_load.load_charge_finish;
        state.load_hamiltonian = GlobalC::restart.info_load.load_H;
        state.load_hamiltonian_finished = GlobalC::restart.info_load.load_H_finish;
        state.restart_exact_exchange = GlobalC::restart.info_load.restart_exx;
        state.folder = GlobalC::restart.folder;
        return state;
    }

    void update(const RestartState& state) override
    {
        GlobalC::restart.info_save.save_charge = state.save_charge;
        GlobalC::restart.info_save.save_H = state.save_hamiltonian;
        GlobalC::restart.info_load.load_charge = state.load_charge;
        GlobalC::restart.info_load.load_charge_finish = state.load_charge_finished;
        GlobalC::restart.info_load.load_H = state.load_hamiltonian;
        GlobalC::restart.info_load.load_H_finish = state.load_hamiltonian_finished;
        GlobalC::restart.info_load.restart_exx = state.restart_exact_exchange;
        GlobalC::restart.folder = state.folder;
    }
};

void capture_input(SimulationContext& context, const Input_para& input)
{
#define MODULE_CONTEXT_INPUT_FIELD(domain, field) context.domain.project(#field, input.field);
#include "source_context/input_field_mapping.inc"
#undef MODULE_CONTEXT_INPUT_FIELD

    context.run.calculation = input.calculation;
    context.run.esolver_type = input.esolver_type;
    context.run.output_level = input.out_level;
    context.run.cal_force = input.cal_force;
    context.run.cal_stress = input.cal_stress;

    context.spin.nspin = input.nspin;
    context.spin.spin_orbit = input.lspinorb;
    context.spin.noncollinear = input.noncolin;

    context.solver.ks_solver = input.ks_solver;
    context.solver.basis_type = input.basis_type;

    context.matrix_output.hs_k = request_from(input.out_mat_hs);
    context.matrix_output.hs_r = request_from(input.out_mat_hs2);
    context.matrix_output.kinetic_k = request_from(input.out_mat_tk);
    context.matrix_output.kinetic_r = request_from(input.out_mat_t);
    context.matrix_output.position_r = request_from(input.out_mat_r);
    context.matrix_output.dh = request_from(input.out_mat_dh);
    context.matrix_output.ds = request_from(input.out_mat_ds);
    context.matrix_output.h_t = request_from(input.out_mat_h_t);
    context.matrix_output.h_vnl = request_from(input.out_mat_h_vnl);
    context.matrix_output.h_vl = request_from(input.out_mat_h_vl);
    context.matrix_output.h_vh = request_from(input.out_mat_h_vh);
    context.matrix_output.h_vxc = request_from(input.out_mat_h_vxc);
    context.matrix_output.h_exx = request_from(input.out_mat_h_exx);
    context.matrix_output.dh_t = request_from(input.out_mat_dh_t);
    context.matrix_output.dh_vnl = request_from(input.out_mat_dh_vnl);
    context.matrix_output.dh_vl = request_from(input.out_mat_dh_vl);
    context.matrix_output.dh_vh = request_from(input.out_mat_dh_vh);
    context.matrix_output.dh_vxc = request_from(input.out_mat_dh_vxc);
    context.matrix_output.dh_exx = request_from(input.out_mat_dh_exx);
    context.matrix_output.vxc_r = request_from(input.out_mat_xc2);
    context.matrix_output.vxc_k = input.out_mat_xc;
    context.matrix_output.band_energy_terms = input.out_eband_terms;
    context.matrix_output.append = input.out_app_flag;
    context.matrix_output.digits = input.out_ndigits;

    context.output.level = input.out_level;
    context.output.electronic_frequency = input.out_freq_elec;
    context.output.ionic_frequency = input.out_freq_ion;
    context.output.tddft_frequency = input.out_freq_td;
    context.output.all_logs = input.out_alllog;

    context.dftu.enabled = input.dft_plus_u != 0;
    context.dftu.dmft_enabled = input.dft_plus_dmft;

    context.exact_exchange.separate_loop = input.exx_separate_loop;
    context.exact_exchange.hybrid_step = input.exx_hybrid_step;
    context.exact_exchange.mixing_beta = input.exx_mixing_beta;
    context.exact_exchange.symmetry_realspace = input.exx_symmetry_realspace;
}

void capture_system(SimulationContext& context, const System_para& system)
{
    context.run.start_time = system.start_time;
    context.files.input_card = system.global_in_card;
    context.files.structure_file = system.global_in_stru;
    context.files.output_directory = system.global_out_dir;
    context.files.readin_directory = system.global_readin_dir;
    context.files.structure_directory = system.global_stru_dir;
    context.files.matrix_directory = system.global_matrix_dir;
    context.files.wavefunction_directory = system.global_wfc_dir;
    context.files.mlkedf_descriptor_directory = system.global_mlkedf_descriptor_dir;
    context.files.deepks_label_directory = system.global_deepks_label_elec_dir;
    context.files.log_file = system.log_file;

    context.parallel.world_rank = system.myrank;
    context.parallel.world_size = system.nproc;
    context.parallel.pool_index = system.mypool;
    context.parallel.kpoint_pool_count = system.npool;
    context.parallel.processes_in_pool = system.nproc_in_pool;
    context.parallel.threads_per_process = system.nthread_per_proc;

    context.basis.nlocal = system.nlocal;
    context.basis.npol = system.npol;
    context.basis.gamma_only_pw = system.gamma_only_pw;
    context.basis.gamma_only_local = system.gamma_only_local;

    context.spin.domag = system.domag;
    context.spin.domag_z = system.domag_z;

    context.grid.charge_nx = system.ncx;
    context.grid.charge_ny = system.ncy;
    context.grid.charge_nz = system.ncz;
    context.grid.double_grid = system.double_grid;
    context.grid.radial_dq = system.dq;
    context.grid.radial_nqx = system.nqx;
    context.grid.radial_nqxq = system.nqxq;

    context.electronic_runtime.two_fermi = system.two_fermi;
    context.electronic_runtime.use_uspp = system.use_uspp;
    context.electronic_runtime.dos_minimum_is_explicit = system.dos_setemin;
    context.electronic_runtime.dos_maximum_is_explicit = system.dos_setemax;
    context.electronic_runtime.local_band_count = system.nbands_l;
    context.electronic_runtime.ks_run = system.ks_run;
    context.electronic_runtime.all_ks_run = system.all_ks_run;
    context.electronic_runtime.has_double_data = system.has_double_data;
    context.electronic_runtime.has_float_data = system.has_float_data;

    context.features.deepks_setorb = system.deepks_setorb;
    context.features.search_periodic_boundaries = system.search_pbc;
    context.features.lcao_kpoint_pool_count = system.kpar_lcao;

    context.output.molecular_dynamics_control = system.out_md_control;
    context.dftu.ramping_ry = system.uramping;
    context.dftu.hubbard_u_ry = system.hubbard_u;
}

} // namespace

SimulationContextBuilder::SimulationContextBuilder(const Input_para& input, const System_para& system)
{
    capture_input(context_, input);
    capture_system(context_, system);

    context_.logs.running = &GlobalV::ofs_running;
    context_.logs.warning = &GlobalV::ofs_warning;
    context_.logs.information = &GlobalV::ofs_info;
    context_.logs.device = &GlobalV::ofs_device;
    context_.exact_exchange_state.reset(new LegacyExactExchangeStateService());
    context_.restart.reset(new LegacyRestartService());
}

void SimulationContextBuilder::capture_runtime(const System_para& system)
{
    if (finalized_)
    {
        throw std::logic_error("SimulationContextBuilder::capture_runtime called after finalize");
    }
    capture_system(context_, system);
    runtime_captured_ = true;
}

SimulationContext SimulationContextBuilder::finalize(const Input_para& input, const System_para& system)
{
    if (finalized_)
    {
        throw std::logic_error("SimulationContextBuilder::finalize called more than once");
    }
    if (!runtime_captured_)
    {
        throw std::logic_error("SimulationContextBuilder::finalize called before runtime initialization");
    }
    if (system.myrank != GlobalV::MY_RANK || system.nproc != GlobalV::NPROC)
    {
        throw std::logic_error("legacy PARAM and GlobalV process topology disagree");
    }

    capture_input(context_, input);
    capture_system(context_, system);
    context_.parallel.world_rank = GlobalV::MY_RANK;
    context_.parallel.world_size = GlobalV::NPROC;
    context_.parallel.kpoint_pool_count = GlobalV::KPAR;
    context_.parallel.pool_index = GlobalV::MY_POOL;
    context_.parallel.processes_in_pool = GlobalV::NPROC_IN_POOL;
    context_.parallel.rank_in_pool = GlobalV::RANK_IN_POOL;
    context_.parallel.band_group_index = GlobalV::MY_BNDGROUP;
    context_.parallel.processes_in_band_group = GlobalV::NPROC_IN_BNDGROUP;
    context_.parallel.rank_in_band_group = GlobalV::RANK_IN_BPGROUP;
    context_.parallel.diagonalization_rank = GlobalV::DRANK;
    context_.parallel.diagonalization_size = GlobalV::DSIZE;
    context_.parallel.diagonalization_color = GlobalV::DCOLOR;
    context_.parallel.grid_rank = GlobalV::GRANK;
    context_.parallel.grid_size = GlobalV::GSIZE;

    finalized_ = true;
    return context_;
}

} // namespace ModuleContext
