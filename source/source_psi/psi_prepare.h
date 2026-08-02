#ifndef PSI_PREPARE_H
#define PSI_PREPARE_H
#include "source_hamilt/hamilt.h"
#include "source_psi/psi_base.h"
#include "source_psi/psi_prepare_base.h"

class UnitCell;
class Structure_Factor;
namespace ModulePW { class PW_Basis_K; }

namespace psi
{

template <typename T, typename Device = base_device::DEVICE_CPU>
class PSIPrepare : public PSIPrepareBase
{
  public:
    PSIPrepare(const std::string& init_wfc_in,
            const std::string& ks_solver_in,
            const std::string& basis_type_in,
            const int& rank,
            const UnitCell& ucell,
            const Structure_Factor& sf,
            const std::vector<int>& ik2iktot,
            const int& nkstot,
            const int& lmaxkb,
            const ModulePW::PW_Basis_K& pw_wfc);

    ~PSIPrepare(){};

    ///@brief prepare the wavefunction initialization
    ///@param random_seed seed for random initialization
    ///@param istep current ion/relax step; informational warnings are only
    ///       printed on the first step to avoid spamming relax output
    void prepare_init(const int& random_seed, const int istep);

    void initialize_psi(Psi<std::complex<double>>* psi,
                        psi::Psi<T, Device>* kspw_psi,
                        hamilt::Hamilt<T, Device>* p_hamilt,
                        std::ofstream& ofs_running);

    void initialize_lcao_in_pw(Psi<T>* psi_local, std::ofstream& ofs_running);

    std::unique_ptr<psi_base<T>> psi_initer;

  private:

    std::string init_wfc = "none";

    std::string ks_solver = "none";

    std::string basis_type = "none";

    const ModulePW::PW_Basis_K& pw_wfc;

    const std::vector<int>& ik2iktot_;
    const int nkstot_;

    const UnitCell& ucell;

    const Structure_Factor& sf;

    const int lmaxkb;

    Device* ctx = {};
    base_device::DEVICE_CPU* cpu_ctx = {};
    const int rank;

    using syncmem_complex_op = base_device::memory::synchronize_memory_op<T, Device, Device>;
    using syncmem_h2d_op = base_device::memory::synchronize_memory_op<T, Device, base_device::DEVICE_CPU>;
};

void allocate_psi(Psi<std::complex<double>>*& psi, const int& nks, const std::vector<int>& ngk, const int& nbands, const int& npwx);

} // namespace psi
#endif
