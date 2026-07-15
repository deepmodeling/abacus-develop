#ifndef ABACUS_SOURCE_CONTEXT_CONTEXT_TYPES_H
#define ABACUS_SOURCE_CONTEXT_CONTEXT_TYPES_H

#include <cstddef>
#include <ctime>
#include <iosfwd>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <utility>
#include <vector>

namespace ModuleContext
{

class InputFieldValue
{
  private:
    struct HolderBase
    {
        virtual ~HolderBase() {}
        virtual const std::type_info& type() const = 0;
    };

    template <typename T>
    struct Holder : HolderBase
    {
        explicit Holder(const T& source) : value(source) {}
        const std::type_info& type() const override { return typeid(T); }
        T value;
    };

  public:
    InputFieldValue() {}

    template <typename T>
    explicit InputFieldValue(const T& value) : holder_(new Holder<T>(value))
    {
    }

    bool empty() const { return !holder_; }
    const std::type_info& type() const { return holder_ ? holder_->type() : typeid(void); }

    template <typename T>
    const T& get() const
    {
        if (!holder_ || holder_->type() != typeid(T))
        {
            throw std::bad_cast();
        }
        return static_cast<const Holder<T>&>(*holder_).value;
    }

  private:
    std::shared_ptr<const HolderBase> holder_;
};

struct FieldProjection
{
    template <typename T>
    void project(const std::string& legacy_name, const T& value)
    {
        mapped_values[legacy_name] = InputFieldValue(value);
    }

    template <typename T>
    const T& value(const std::string& legacy_name) const
    {
        const std::map<std::string, InputFieldValue>::const_iterator it = mapped_values.find(legacy_name);
        if (it == mapped_values.end())
        {
            throw std::out_of_range("unknown projected INPUT field: " + legacy_name);
        }
        return it->second.get<T>();
    }

    std::size_t mapped_field_count() const { return mapped_values.size(); }

    std::map<std::string, InputFieldValue> mapped_values;
};

struct InputConfigDomain : FieldProjection
{
};

struct RunControl : InputConfigDomain
{
    std::time_t start_time = 0;
    std::string calculation;
    std::string esolver_type;
    std::string output_level;
    bool cal_force = false;
    bool cal_stress = false;
};

struct FileSystemLayout
{
    std::string input_card;
    std::string structure_file;
    std::string output_directory;
    std::string readin_directory;
    std::string structure_directory;
    std::string matrix_directory;
    std::string wavefunction_directory;
    std::string mlkedf_descriptor_directory;
    std::string deepks_label_directory;
    std::string log_file;
};

struct ParallelTopology
{
    int world_size = 1;
    int world_rank = 0;
    int kpoint_pool_count = 1;
    int pool_index = 0;
    int processes_in_pool = 1;
    int rank_in_pool = 0;
    int band_group_index = 0;
    int processes_in_band_group = 1;
    int rank_in_band_group = 0;
    int diagonalization_rank = 0;
    int diagonalization_size = 1;
    int diagonalization_color = 0;
    int grid_rank = 0;
    int grid_size = 1;
    int threads_per_process = 1;
};

struct LogStreams
{
    std::ostream* running = nullptr;
    std::ostream* warning = nullptr;
    std::ostream* information = nullptr;
    std::ostream* device = nullptr;
};

struct BasisInfo
{
    int nlocal = 0;
    int npol = 1;
    bool gamma_only_pw = false;
    bool gamma_only_local = false;
};

struct GridInfo
{
    int charge_nx = 0;
    int charge_ny = 0;
    int charge_nz = 0;
    bool double_grid = false;
    double radial_dq = 0.0;
    int radial_nqx = 0;
    int radial_nqxq = 0;
};

struct ElectronicRuntimeInfo
{
    bool two_fermi = false;
    bool use_uspp = false;
    bool dos_minimum_is_explicit = false;
    bool dos_maximum_is_explicit = false;
    int local_band_count = 0;
    bool ks_run = false;
    bool all_ks_run = true;
    bool has_double_data = true;
    bool has_float_data = false;
};

struct FeatureRuntimeInfo
{
    bool deepks_setorb = false;
    bool search_periodic_boundaries = true;
    int lcao_kpoint_pool_count = 1;
};

struct SpinConfig : InputConfigDomain
{
    int nspin = 1;
    bool spin_orbit = false;
    bool noncollinear = false;
    bool domag = false;
    bool domag_z = false;
};

struct SolverConfig : InputConfigDomain
{
    std::string ks_solver;
    std::string basis_type;
};

struct MatrixOutputRequest
{
    bool enabled = false;
    int precision = 8;
    int mode = 0;
};

struct MatrixOutputConfig : InputConfigDomain
{
    MatrixOutputRequest hs_k;
    MatrixOutputRequest hs_r;
    MatrixOutputRequest kinetic_k;
    MatrixOutputRequest kinetic_r;
    MatrixOutputRequest position_r;
    MatrixOutputRequest dh;
    MatrixOutputRequest ds;
    MatrixOutputRequest h_t;
    MatrixOutputRequest h_vnl;
    MatrixOutputRequest h_vl;
    MatrixOutputRequest h_vh;
    MatrixOutputRequest h_vxc;
    MatrixOutputRequest h_exx;
    MatrixOutputRequest dh_t;
    MatrixOutputRequest dh_vnl;
    MatrixOutputRequest dh_vl;
    MatrixOutputRequest dh_vh;
    MatrixOutputRequest dh_vxc;
    MatrixOutputRequest dh_exx;
    MatrixOutputRequest vxc_r;
    bool vxc_k = false;
    bool band_energy_terms = false;
    bool append = true;
    int digits = 8;
    double sparse_threshold = 1.0e-10;
    bool binary = false;
};

struct GeneralOutputConfig : InputConfigDomain
{
    std::string level;
    int electronic_frequency = 0;
    int ionic_frequency = 0;
    int tddft_frequency = 0;
    bool all_logs = false;
    bool molecular_dynamics_control = false;
};

struct DftUConfig : InputConfigDomain
{
    bool enabled = false;
    bool dmft_enabled = false;
    double ramping_ry = -10.0 / 13.6;
    std::vector<double> hubbard_u_ry;
};

struct ExactExchangeConfig : InputConfigDomain
{
    bool separate_loop = true;
    int hybrid_step = 0;
    double mixing_beta = 1.0;
    bool symmetry_realspace = true;
};

struct ExactExchangeState : FieldProjection
{
    bool enabled = false;
    bool real_number = false;
    double hybrid_alpha = 0.0;
    double hse_omega = 0.0;
};

struct RestartState
{
    bool save_charge = false;
    bool save_hamiltonian = false;
    bool load_charge = false;
    bool load_charge_finished = false;
    bool load_hamiltonian = false;
    bool load_hamiltonian_finished = false;
    bool restart_exact_exchange = false;
    std::string folder;
};

class ExactExchangeStateService
{
  public:
    virtual ~ExactExchangeStateService() {}
    virtual ExactExchangeState snapshot() const = 0;
};

class RestartService
{
  public:
    virtual ~RestartService() {}
    virtual RestartState snapshot() const = 0;
    virtual void update(const RestartState& state) = 0;
};

// The legacy parser remains the only source of defaults and validation. The
// complete one-to-one field projection is retained in these domain containers;
// stable named members are added only when a migrated consumer needs them.
struct CellConfig : InputConfigDomain {};
struct FileInputConfig : InputConfigDomain {};
struct ParallelConfig : InputConfigDomain {};
struct DeviceConfig : InputConfigDomain {};
struct ElectronicStructureConfig : InputConfigDomain {};
struct ScfMixingConfig : InputConfigDomain {};
struct LcaoConfig : InputConfigDomain {};
struct RelaxationConfig : InputConfigDomain {};
struct MolecularDynamicsConfig : InputConfigDomain {};
struct OfdftConfig : InputConfigDomain {};
struct StochasticDftConfig : InputConfigDomain {};
struct DeepksConfig : InputConfigDomain {};
struct RealTimeTddftConfig : InputConfigDomain {};
struct LinearResponseTddftConfig : InputConfigDomain {};
struct ChargeWavefunctionOutputConfig : InputConfigDomain {};
struct RestartOutputConfig : InputConfigDomain {};
struct PostprocessConfig : InputConfigDomain {};
struct ModelConfig : InputConfigDomain {};
struct VdwConfig : InputConfigDomain {};
struct SpinConstraintConfig : InputConfigDomain {};
struct QuasiAtomicOrbitalConfig : InputConfigDomain {};
struct PexsiConfig : InputConfigDomain {};
struct TestConfig : InputConfigDomain {};
struct RdmftConfig : InputConfigDomain {};
struct ExternalXcConfig : InputConfigDomain {};
struct TimeDependentOfdftConfig : InputConfigDomain {};
struct UncommonHardwareConfig : InputConfigDomain {};

} // namespace ModuleContext

#endif
