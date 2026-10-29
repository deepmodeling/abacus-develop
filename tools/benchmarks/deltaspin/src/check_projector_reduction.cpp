#include "module_hamilt_pw/hamilt_pwdft/fs_nonlocal_tools.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_hamilt_pw/hamilt_pwdft/onsite_projector.h"
#include "module_base/parallel_comm.h"

#include <mpi.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <unistd.h>

using Complex = std::complex<double>;
using CPU = base_device::DEVICE_CPU;
using GPU = base_device::DEVICE_GPU;

void require(bool condition, const char* message)
{
    if (!condition) { std::cerr << message << '\n'; MPI_Abort(MPI_COMM_WORLD, 1); }
}

template <typename Device>
class Buffer
{
public:
    Complex* data = nullptr;
    const size_t size;
    explicit Buffer(size_t count) : size(count)
    {
        base_device::memory::resize_memory_op<Complex, Device>()(nullptr, data, size);
    }
    ~Buffer() { base_device::memory::delete_memory_op<Complex, Device>()(nullptr, data); }
    Buffer(const Buffer&) = delete;
    Buffer& operator=(const Buffer&) = delete;
    void upload(const std::vector<Complex>& values)
    {
        require(values.size() == size, "Input shape mismatch");
        base_device::memory::synchronize_memory_op<Complex, Device, CPU>()(
            nullptr, nullptr, data, values.data(), size);
    }
    std::vector<Complex> download() const
    {
        std::vector<Complex> values(size);
        base_device::memory::synchronize_memory_op<Complex, CPU, Device>()(
            nullptr, nullptr, values.data(), data, size);
        return values;
    }
};

template <typename Device>
void run(const char* device, int polarizations, int rank, int ranks)
{
    GlobalV::device_flag = device;
    GlobalV::NPOL = polarizations;
    auto& cell = GlobalC::ucell;
    cell.set_iat2iwt(polarizations);
    ModulePW::PW_Basis_K basis(device, "double");
    basis.initmpi(ranks, rank, MPI_COMM_WORLD);
    basis.initgrids(cell.lat0, cell.latvec, 16, 18, 20);
    ModuleBase::Vector3<double> points[2] = {{0.13, -0.21, 0.17}, {-0.2, 0.1, -0.07}};
    basis.initparameters(false, 15.0, 2, points, 1, true);
    basis.setuptransform();
    basis.collect_local_pw();
    Structure_Factor structure;
    structure.set(&basis, 0);
    structure.setup_structure_factor(&cell, &basis);
    constexpr int maximum_bands = 5;
    const int leading = basis.npwk_max + 7;
    psi::Psi<Complex, Device> waves(2, maximum_bands, polarizations * basis.npwk_max, basis.npwk);
    K_Vectors kpoints;
    ModuleBase::matrix weights(2, maximum_bands), energies(2, maximum_bands), angular(1, 1);
    ModuleBase::realArray radial(1, 1, 2048);
    for (int i = 0; i < radial.getSize(); ++i) { radial.ptr[i] = std::exp(-0.0001 * i * i); }
    Buffer<Device> projectors(basis.npwk_max);
    Buffer<Device> input(size_t(leading) * maximum_bands * polarizations);
    Buffer<Device> output(maximum_bands * polarizations + 3);
    hamilt::FS_Nonlocal_tools<double, Device> tools({1}, {0}, radial, angular, projectors.data,
        &cell, &waves, &kpoints, &basis, &structure, weights, energies);
    const Complex guard(17.25, -8.75);
    double error = 0;
    for (int ik : {0, 1, 0})
    {
        const int rows = basis.npwk[ik];
        for (int bands : {1, maximum_bands, 2})
        {
            const int fields = bands * polarizations;
            std::vector<Complex> values(input.size, guard), initialized(output.size, guard);
            for (int field = 0; field < fields; ++field)
            for (int row = 0; row < rows; ++row)
            {
                const double index = 1 + row + 13 * field + 7 * rank + 3 * ik;
                values[size_t(field) * leading + row] = {0.03 * std::sin(index), 0.02 * std::cos(0.7 * index)};
            }
            input.upload(values);
            output.upload(initialized);
            tools.cal_becp(ik, bands, output.data, input.data, leading);
            const auto matrix = projectors.download(), actual = output.download();
            std::vector<Complex> expected(fields);
            for (int field = 0; field < fields; ++field)
            {
                std::complex<long double> sum = 0;
                for (int row = 0; row < rows; ++row)
                {
                    sum += std::conj(std::complex<long double>(matrix[row]))
                           * std::complex<long double>(values[size_t(field) * leading + row]);
                }
                expected[field] = Complex(sum);
            }
            MPI_Allreduce(MPI_IN_PLACE, expected.data(), fields, MPI_C_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
            for (int field = 0; field < fields; ++field)
            {
                const double difference = std::abs(actual[field] - expected[field]);
                require(std::isfinite(difference), "Nonfinite projection");
                error = std::max(error, difference);
            }
            for (size_t i = fields; i < output.size; ++i) { require(actual[i] == guard, "Output guard changed"); }
            require(input.download() == values, "Wavefunction input changed");
        }
    }
    require(error < 2e-12, "Distributed projection differs from direct contraction");
    if (rank == 0)
    {
        std::cout << "device=" << device << " ranks=" << ranks << " polarizations=" << polarizations
                  << " cases=9 maximum_error=" << error << std::endl;
    }
}

void check_host_views(const char* asset_directory, int rank, int ranks)
{
    require(chdir(asset_directory) == 0, "Cannot enter the projector asset directory");
    GlobalV::device_flag = "gpu";
    GlobalV::NPOL = 2;
    auto& cell = GlobalC::ucell;
    cell.set_iat2iwt(2);
    cell.atoms[0].nwl = 1;
    cell.orbital_fn = new std::string[1]{"Fe_gga_8au_200.0Ry_4s2p2d1f.orb"};
    ModulePW::PW_Basis_K basis("gpu", "double");
    basis.initmpi(ranks, rank, MPI_COMM_WORLD);
    basis.initgrids(cell.lat0, cell.latvec, 16, 18, 20);
    ModuleBase::Vector3<double> points[2] = {{0.13, -0.21, 0.17}, {-0.2, 0.1, -0.07}};
    basis.initparameters(false, 15.0, 2, points, 1, true);
    basis.setuptransform();
    basis.collect_local_pw();
    Structure_Factor structure;
    structure.set(&basis, 0);
    structure.setup_structure_factor(&cell, &basis);
    psi::Psi<Complex, GPU> waves(2, 5, 2 * basis.npwk_max, basis.npwk);
    K_Vectors kpoints;
    ModuleBase::matrix weights(2, 5), energies(2, 5);
    auto* onsite = projectors::OnsiteProjector<double, GPU>::get_instance();
    require(onsite->get_h_becp() == nullptr, "Uninitialized host view is not null");
    onsite->init("./", &cell, waves, kpoints, basis, structure, 3.0, 2048, GlobalV::DQ, weights, energies);
    const int leading = basis.npwk_max + 7;
    Buffer<GPU> input(size_t(leading) * 10);
    std::vector<Complex> values(input.size);
    int generation = 0;
    for (int ik : {0, 1, 0})
    {
        onsite->tabulate_atomic(ik);
        for (int bands : {1, 5, 2})
        {
            ++generation;
            for (size_t i = 0; i < values.size(); ++i)
            {
                values[i] = {0.01 * std::sin(0.13 * (i + generation)),
                             0.02 * std::cos(0.17 * (i + 7 * generation + rank))};
            }
            input.upload(values);
            onsite->overlap_proj_psi(2 * bands, input.data, leading);
            const int active = 2 * bands * onsite->get_tot_nproj();
            std::vector<Complex> expected(active);
            base_device::memory::synchronize_memory_op<Complex, CPU, GPU>()(
                nullptr, nullptr, expected.data(), onsite->get_becp(), active);
            for (int access = 0; access < 2; ++access)
            {
                const auto* host = onsite->get_h_becp();
                for (int i = 0; i < active; ++i)
                {
                    require(std::isfinite(host[i].real()) && std::isfinite(host[i].imag()), "Nonfinite host view");
                    require(host[i] == expected[i], "Host view differs from current device overlaps");
                }
            }
            require(input.download() == values, "Host acquisition changed wavefunctions");
        }
    }
    if (rank == 0) { std::cout << "host_views ranks=" << ranks << " generations=9 accesses=18 bitwise=pass\n"; }
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    require(argc == 2, "Usage: check_projector_reduction ASSET_DIRECTORY");
    int rank = 0, ranks = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    POOL_WORLD = MPI_COMM_WORLD;
    GlobalV::NPROC_IN_POOL = ranks;
    GlobalV::RANK_IN_POOL = rank;
    GlobalV::MY_RANK = rank;
    GlobalV::NPROC = ranks;
    GlobalV::precision_flag = "double";
    auto& cell = GlobalC::ucell;
    cell.ntype = 1;
    cell.nat = 1;
    cell.iat2it = new int[1]{0};
    cell.iat2ia = new int[1]{0};
    cell.atoms = new Atom[1];
    cell.set_atom_flag = true;
    cell.atoms[0].na = 1;
    cell.atoms[0].nw = 1;
    cell.atoms[0].tau[0] = {0.13, 0.27, 0.38};
    cell.atoms[0].taud[0] = cell.atoms[0].tau[0];
    cell.lat0 = 6;
    cell.latvec = ModuleBase::Matrix3(1, 0, 0, 0, 1, 0, 0, 0, 1);
    cell.G = cell.latvec;
    cell.tpiba = 2 * std::acos(-1.0) / cell.lat0;
    cell.tpiba2 = cell.tpiba * cell.tpiba;
    cell.omega = cell.lat0 * cell.lat0 * cell.lat0;
    hsolver::createGpuBlasHandle();
    for (int polarizations : {1, 2})
    {
        run<CPU>("cpu", polarizations, rank, ranks);
        run<GPU>("gpu", polarizations, rank, ranks);
    }
    check_host_views(argv[1], rank, ranks);
    hsolver::destoryBLAShandle();
    MPI_Finalize();
}
