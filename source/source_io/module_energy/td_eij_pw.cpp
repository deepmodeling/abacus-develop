#include "source_io/module_energy/td_eij_pw.h"

#include "source_base/parallel_device.h"
#include "source_base/timer.h"
#include "source_io/module_energy/td_eij_io.h"

#include <algorithm>
#include <limits>
#include <vector>

namespace ModuleIO
{
template <typename T, typename Device>
void print_eij_pw(hsolver::HSOperator<T, Device>& op,
                  hsolver::TDOrbOps<T, Device>& orbitals,
                  psi::Psi<T, Device>& current,
                  const K_Vectors& kv,
                  int nspin,
                  int step,
                  int iter,
                  double threshold,
                  std::ostream& log)
{
    ModuleBase::timer::start("ModuleIO", "print_eij_pw");

    const int original_k = current.get_current_k();
    const int nband = current.get_nbands();
    const size_t elements = static_cast<size_t>(nband) * nband;
    const int nk = kv.get_nkstot() / (nspin == 2 ? 2 : 1);
    const Parallel_Kpoints& pools = kv.para_k;
    const bool main_rank = pools.my_pool == 0 && pools.rank_in_pool == 0;

    for (int global_k = 0; global_k < kv.get_nkstot(); ++global_k)
    {
        const int physical_k = global_k % nk;
        const int spin = global_k / nk;
        int owner = 0;
        int local_k = global_k;
#ifdef __MPI
        owner = pools.whichpool[physical_k];
        local_k = physical_k - pools.startk_pool[owner] + spin * pools.nks_pool[owner];
#endif

        std::vector<std::complex<double>> eij;
        if (pools.my_pool == owner)
        {
            eij = orbitals.cal_eij(op, current, local_k);
        }

#ifdef __MPI
        // Each pool already has the full matrix; transfer only one copy to the main log.
        if (owner != 0 && pools.rank_in_pool == 0 && (main_rank || pools.my_pool == owner))
        {
            if (main_rank)
            {
                eij.resize(elements);
            }
            MPI_Comm world = MPI_COMM_WORLD;
            const int source = pools.get_startpro_pool(owner);
            const size_t max_count = std::numeric_limits<int>::max();
            for (size_t offset = 0; offset < elements; offset += max_count)
            {
                const int count = static_cast<int>(std::min(max_count, elements - offset));
                if (main_rank)
                {
                    Parallel_Common::recv_data(eij.data() + offset, count, source, 0, world, MPI_STATUS_IGNORE);
                }
                else
                {
                    Parallel_Common::send_data(eij.data() + offset, count, 0, 0, world);
                }
            }
        }
#endif

        if (main_rank)
        {
            const EijOutputInfo output{threshold, step + 1, iter, physical_k + 1, nspin, spin + 1};
            print_eij(eij.data(), nband, nband, output, log);
        }
    }

    current.fix_k(original_k);
    op.update_k(original_k);

    ModuleBase::timer::end("ModuleIO", "print_eij_pw");
}

template void print_eij_pw(hsolver::HSOperator<std::complex<float>, base_device::DEVICE_CPU>&,
                           hsolver::TDOrbOps<std::complex<float>, base_device::DEVICE_CPU>&,
                           psi::Psi<std::complex<float>, base_device::DEVICE_CPU>&,
                           const K_Vectors&,
                           int,
                           int,
                           int,
                           double,
                           std::ostream&);
template void print_eij_pw(hsolver::HSOperator<std::complex<double>, base_device::DEVICE_CPU>&,
                           hsolver::TDOrbOps<std::complex<double>, base_device::DEVICE_CPU>&,
                           psi::Psi<std::complex<double>, base_device::DEVICE_CPU>&,
                           const K_Vectors&,
                           int,
                           int,
                           int,
                           double,
                           std::ostream&);
#if defined(__CUDA) || defined(__ROCM)
template void print_eij_pw(hsolver::HSOperator<std::complex<float>, base_device::DEVICE_GPU>&,
                           hsolver::TDOrbOps<std::complex<float>, base_device::DEVICE_GPU>&,
                           psi::Psi<std::complex<float>, base_device::DEVICE_GPU>&,
                           const K_Vectors&,
                           int,
                           int,
                           int,
                           double,
                           std::ostream&);
template void print_eij_pw(hsolver::HSOperator<std::complex<double>, base_device::DEVICE_GPU>&,
                           hsolver::TDOrbOps<std::complex<double>, base_device::DEVICE_GPU>&,
                           psi::Psi<std::complex<double>, base_device::DEVICE_GPU>&,
                           const K_Vectors&,
                           int,
                           int,
                           int,
                           double,
                           std::ostream&);
#endif
} // namespace ModuleIO
