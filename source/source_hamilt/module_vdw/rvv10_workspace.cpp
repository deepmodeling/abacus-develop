#include "rvv10_workspace.h"

#include "source_base/matrix3.h"
#include "source_base/memory_recorder.h"
#include "source_base/parallel_device.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_hamilt/module_xc/xc_functional.h"

// vdwxc.h includes <complex.h> inside an extern "C" block. Include the C++ header first so that
// libc++ does not parse its templates with C linkage, then drop the I macro.
#include <complex.h>
#include <complex>
#include <vdwxc.h>
#ifdef __MPI
#include <vdwxc_mpi.h>
#endif
#undef I

#include <algorithm>

namespace vdw
{
namespace
{
#ifndef __LIBVDWXC_LOCAL_GRID
/// z planes [begin, end) of the libvdwxc slab of process rank: FFTW-MPI blocks of ceil(nz / nproc) planes.
void libvdwxc_slab(const int nz, const int nproc, const int rank, int& begin, int& end)
{
    const int block = (nz + nproc - 1) / nproc;
    const int first_plane = rank * block;
    begin = std::min(first_plane, nz);
    const int block_end = begin + block;
    end = std::min(block_end, nz);
}
#endif

/// The instance of the run. It is created with new and deleted in release() (not a static object): a static
/// object would be destroyed after MPI_Finalize, too late for the FFTW-MPI plans that the destructor frees.
Rvv10Workspace* workspace_instance = nullptr;
} // namespace

Rvv10Workspace& Rvv10Workspace::instance()
{
    if (workspace_instance == nullptr)
    {
        workspace_instance = new Rvv10Workspace;
    }
    return *workspace_instance;
}

void Rvv10Workspace::release()
{
    delete workspace_instance;
    workspace_instance = nullptr;
}

Rvv10Workspace::~Rvv10Workspace()
{
    if (data_ != nullptr)
    {
        vdwxc_finalize(&data_);
    }
}

bool Rvv10Workspace::same_grid(const ModulePW::PW_Basis& pw) const
{
    bool same = data_ != nullptr && pw.nx == nx_ && pw.ny == ny_ && pw.nz == nz_ && pw.poolnproc == nproc_;
#ifdef __MPI
    same = same && (pw.poolnproc == 1 || pw.pool_world == comm_);
#endif
    return same;
}

/// The cell is set on every call: a cell relaxation keeps the FFT grid but changes the cell. libvdwxc
/// integrates with the signed cell volume, and inverting a left-handed cell leaves rVV10 unchanged.
void Rvv10Workspace::set_cell(const ModulePW::PW_Basis& pw)
{
    ModuleBase::Matrix3 a = pw.latvec * pw.lat0;
    if (a.Det() < 0.0)
    {
        a = a * -1.0;
    }
    vdwxc_set_unit_cell(data_, pw.nz, pw.nx, pw.ny, a.e31, a.e32, a.e33, a.e11, a.e12, a.e13, a.e21, a.e22, a.e23);
}

void Rvv10Workspace::rebuild(const ModulePW::PW_Basis& pw)
{
    if (data_ != nullptr)
    {
        vdwxc_finalize(&data_);
    }
    data_ = vdwxc_new(VDWXC_RVV10);
    nx_ = pw.nx;
    ny_ = pw.ny;
    nz_ = pw.nz;
    nproc_ = pw.poolnproc;
    // the cell must be known before the FFT plans are made
    this->set_cell(pw);
    bool distributed = false;
#ifdef __MPI
    comm_ = pw.pool_world;
    distributed = pw.poolnproc > 1;
    if (distributed)
    {
        vdwxc_init_mpi(data_, comm_);
    }
#endif
    if (!distributed)
    {
        vdwxc_init_serial(data_);
    }
    this->query_slabs(pw);
}

void Rvv10Workspace::query_slabs(const ModulePW::PW_Basis& pw)
{
    const int nproc = pw.poolnproc;
    int slab[2] = {0, pw.nz}; // first z plane and end of the slab of this process
#ifdef __LIBVDWXC_LOCAL_GRID
    int shape[3] = {0, 0, 0};
    int offset[3] = {0, 0, 0};
    vdwxc_get_local_grid(data_, shape, offset);
    // the first libvdwxc axis is z (see set_cell), and whole z planes are exchanged
    if (offset[1] != 0 || offset[2] != 0 || shape[1] != pw.nx || shape[2] != pw.ny)
    {
        ModuleBase::WARNING_QUIT("Rvv10Workspace::query_slabs", "libvdwxc distributes the x or y axis of the grid");
    }
    slab[0] = offset[0];
    slab[1] = offset[0] + shape[0];
#else
    libvdwxc_slab(pw.nz, nproc, pw.poolrank, slab[0], slab[1]);
#endif
    std::vector<int> slabs = {slab[0], slab[1]};
#ifdef __MPI
    if (nproc > 1)
    {
        const int slab_values = 2 * nproc;
        slabs.resize(slab_values);
        Parallel_Common::allgather_data(slab, 2, slabs.data(), comm_);
    }
#endif
    // the slabs must follow each other in rank order and cover the z axis; an empty slab, wherever libvdwxc
    // places it, is put at the end of the previous one
    const int bound_count = nproc + 1;
    slab_bound_.assign(bound_count, 0);
    int covered = 0;
    for (int ip = 0; ip < nproc; ++ip)
    {
        const int begin = slabs[2 * ip];
        const int end = slabs[2 * ip + 1];
        if (end > begin)
        {
            if (begin != covered)
            {
                ModuleBase::WARNING_QUIT("Rvv10Workspace::query_slabs", "the libvdwxc slabs are not in rank order");
            }
            covered = end;
        }
        slab_bound_[ip + 1] = covered;
    }
    if (covered != pw.nz)
    {
        ModuleBase::WARNING_QUIT("Rvv10Workspace::query_slabs", "the libvdwxc slabs do not cover the z axis");
    }
}

void Rvv10Workspace::prepare(const ModulePW::PW_Basis& pw, const XCFunctionalParameters& parameters)
{
    if (this->same_grid(pw))
    {
        this->set_cell(pw);
    }
    else
    {
        this->rebuild(pw);
    }
    // the libvdwxc slabs are fixed with the grid and the pool, but the ABACUS planes depend also on the basis,
    // so the exchange is planned on every call (counts and resizes only)
    this->plan_exchange(pw);
    if (vdwxc_set_rvv10_parameters(data_, parameters.rvv10_b, parameters.rvv10_c) != 0)
    {
        ModuleBase::WARNING_QUIT("Rvv10Workspace::prepare", "libvdwxc rejected vdw_rvv10_b or vdw_rvv10_c");
    }
}

void Rvv10Workspace::plan_exchange(const ModulePW::PW_Basis& pw)
{
    const int nproc = pw.poolnproc;
    const int slab_start = slab_bound_[pw.poolrank];
    const int slab_end = slab_bound_[pw.poolrank + 1];
    const int plane_end = pw.startz_current + pw.nplane;
    send_count_.assign(nproc, 0);
    send_offset_.assign(nproc, 0);
    recv_count_.assign(nproc, 0);
    recv_offset_.assign(nproc, 0);
    for (int ip = 0; nproc > 1 && ip < nproc; ++ip)
    {
        // ABACUS planes of this process that lie in the slab of ip
        const int slab_start_of_ip = slab_bound_[ip];
        const int slab_end_of_ip = slab_bound_[ip + 1];
        int low = std::max(pw.startz_current, slab_start_of_ip);
        int high = std::min(plane_end, slab_end_of_ip);
        if (high > low)
        {
            send_count_[ip] = (high - low) * pw.nxy;
            send_offset_[ip] = (low - pw.startz_current) * pw.nxy;
        }
        // ABACUS planes of ip that lie in the slab of this process
        const int plane_end_of_ip = pw.startz[ip] + pw.numz[ip];
        low = std::max(pw.startz[ip], slab_start);
        high = std::min(plane_end_of_ip, slab_end);
        if (high > low)
        {
            recv_count_[ip] = (high - low) * pw.nxy;
            recv_offset_[ip] = (low - slab_start) * pw.nxy;
        }
    }
    // libvdwxc expects valid pointers also for an empty slab
    const int slab_points = (slab_end - slab_start) * pw.nxy;
    const int slab_size = std::max(1, slab_points);
    // the transposed ABACUS planes are only needed for the exchange between processes
    const int planes_size = (nproc > 1) ? std::max(1, pw.nrxx) : 0;
    planes_.resize(planes_size);
    rho_.resize(slab_size);
    sigma_.resize(slab_size);
    dedn_.resize(slab_size);
    dedsigma_.resize(slab_size);
    const long buffer_bytes = static_cast<long>(sizeof(double)) * (planes_size + 4L * slab_size);
    ModuleBase::Memory::record("Rvv10Workspace::buffers", buffer_bytes);
}

void Rvv10Workspace::to_libvdwxc(const ModulePW::PW_Basis& pw,
                                 const std::vector<double>& field,
                                 std::vector<double>& slab)
{
    std::vector<double>& planes = (pw.poolnproc > 1) ? planes_ : slab;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int ixy = 0; ixy < pw.nxy; ++ixy)
    {
        for (int iz = 0; iz < pw.nplane; ++iz)
        {
            planes[iz * pw.nxy + ixy] = field[ixy * pw.nplane + iz];
        }
    }
#ifdef __MPI
    if (pw.poolnproc > 1)
    {
        Parallel_Common::alltoallv_data(planes_.data(),
                                        send_count_.data(),
                                        send_offset_.data(),
                                        slab.data(),
                                        recv_count_.data(),
                                        recv_offset_.data(),
                                        pw.pool_world);
    }
#endif
}

void Rvv10Workspace::from_libvdwxc(const ModulePW::PW_Basis& pw,
                                   const std::vector<double>& slab,
                                   std::vector<double>& field)
{
#ifdef __MPI
    if (pw.poolnproc > 1)
    {
        Parallel_Common::alltoallv_data(slab.data(),
                                        recv_count_.data(),
                                        recv_offset_.data(),
                                        planes_.data(),
                                        send_count_.data(),
                                        send_offset_.data(),
                                        pw.pool_world);
    }
#endif
    const std::vector<double>& planes = (pw.poolnproc > 1) ? planes_ : slab;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int ixy = 0; ixy < pw.nxy; ++ixy)
    {
        for (int iz = 0; iz < pw.nplane; ++iz)
        {
            field[ixy * pw.nplane + iz] = planes[iz * pw.nxy + ixy];
        }
    }
}

double Rvv10Workspace::evaluate(const ModulePW::PW_Basis& pw,
                                const XCFunctionalParameters& parameters,
                                const std::vector<double>& rho,
                                const std::vector<double>& sigma,
                                std::vector<double>* dedn,
                                std::vector<double>& dedsigma,
                                double* stress)
{
    ModuleBase::timer::start("Rvv10Workspace", "evaluate");
    this->prepare(pw, parameters);
    this->to_libvdwxc(pw, rho, rho_);
    this->to_libvdwxc(pw, sigma, sigma_);
    double energy = 0.0;
    if (stress == nullptr)
    {
        energy = vdwxc_calculate(data_, rho_.data(), sigma_.data(), dedn_.data(), dedsigma_.data());
    }
    else
    {
        energy = vdwxc_stress(data_, stress, rho_.data(), sigma_.data(), dedn_.data(), dedsigma_.data());
    }
    // the exchange is collective, and every process of the pool asks for dE/dn or none does
    if (dedn != nullptr)
    {
        this->from_libvdwxc(pw, dedn_, *dedn);
    }
    this->from_libvdwxc(pw, dedsigma_, dedsigma);
#ifdef __MPI
    // the distributed backend returns the contribution of this process
    if (pw.poolnproc > 1)
    {
        Parallel_Common::reduce_data(&energy, 1, pw.pool_world);
        if (stress != nullptr)
        {
            Parallel_Common::reduce_data(stress, rvv10_stress_components, pw.pool_world);
        }
    }
#endif
    ModuleBase::timer::end("Rvv10Workspace", "evaluate");
    return energy;
}

} // namespace vdw
