#ifndef RVV10_WORKSPACE_H
#define RVV10_WORKSPACE_H

#ifdef __MPI
#include <mpi.h>
#endif
#include <vector>

struct XCFunctionalParameters;
struct vdwxc_data_obj; // libvdwxc's opaque handle; vdwxc_data in vdwxc.h is a pointer to it
namespace ModulePW
{
class PW_Basis;
}

namespace vdw
{

/// Cartesian components of the stress tensor returned by vdwxc_stress, stored as [3 * a + b].
const int rvv10_stress_components = 9;

/**
 * @brief The libvdwxc object of a run: one FFT grid and pool, kept between calls, rebuilt when either changes.
 *
 * libvdwxc gets the lattice vectors in the cyclic order (a3, a1, a2), so that the first axis of its grid is the
 * z axis that ABACUS distributes. libvdwxc splits that axis into slabs. When it provides vdwxc_get_local_grid
 * (__LIBVDWXC_LOCAL_GRID, set by CMake), every process asks for its slab and the slabs are gathered over the
 * pool and checked once per grid; otherwise the FFTW-MPI default blocks of ceil(nz/P) planes are assumed. The
 * slab of a process is stored as [z][x][y], the transpose of the ABACUS layout [x][y][z], and whole z planes
 * are moved between the two distributions by all-to-all exchanges.
 *
 * The single instance is created on first use and released by release(), which vdw::finalize_rvv10 calls at
 * the end of the run; it must happen before MPI_Finalize, because the FFTW-MPI plans are destroyed with it.
 */
class Rvv10Workspace
{
  public:
    /// The instance of the run, created on first use.
    static Rvv10Workspace& instance();
    /// Finalize the libvdwxc object and free the buffers; the next evaluate() builds them again.
    static void release();

    /**
     * @brief libvdwxc energy and derivatives of the ABACUS-distributed density.
     * @param pw          grid and pool of the density (chr.rhopw)
     * @param parameters  rVV10 switch and the parameters b and C
     * @param rho         total density on the nrxx points of this process
     * @param sigma       |grad rho|^2 on the same points
     * @param dedn        dE/dn on the same points, or nullptr when it is not needed
     * @param dedsigma    dE/dsigma on the same points
     * @param stress      with stress != nullptr also the explicit strain derivative, 9 Cartesian components
     *                    (Ha, not divided by the volume)
     * @return the energy in Ha, reduced over the pool
     */
    double evaluate(const ModulePW::PW_Basis& pw,
                    const XCFunctionalParameters& parameters,
                    const std::vector<double>& rho,
                    const std::vector<double>& sigma,
                    std::vector<double>* dedn,
                    std::vector<double>& dedsigma,
                    double* stress);

  private:
    Rvv10Workspace() = default;
    ~Rvv10Workspace();
    Rvv10Workspace(const Rvv10Workspace&) = delete;
    Rvv10Workspace& operator=(const Rvv10Workspace&) = delete;

    bool same_grid(const ModulePW::PW_Basis& pw) const;
    void set_cell(const ModulePW::PW_Basis& pw);
    void rebuild(const ModulePW::PW_Basis& pw);
    void query_slabs(const ModulePW::PW_Basis& pw);
    void prepare(const ModulePW::PW_Basis& pw, const XCFunctionalParameters& parameters);
    void plan_exchange(const ModulePW::PW_Basis& pw);
    void to_libvdwxc(const ModulePW::PW_Basis& pw, const std::vector<double>& field, std::vector<double>& slab);
    void from_libvdwxc(const ModulePW::PW_Basis& pw, const std::vector<double>& slab, std::vector<double>& field);


    vdwxc_data_obj* data_ = nullptr;
    int nx_ = 0;
    int ny_ = 0;
    int nz_ = 0;
    int nproc_ = 0;
#ifdef __MPI
    MPI_Comm comm_ = MPI_COMM_NULL;
#endif
    // the libvdwxc slab of process ip is the z planes [slab_bound_[ip], slab_bound_[ip + 1])
    std::vector<int> slab_bound_;
    // planes moved to (send) and from (recv) each process, in doubles
    std::vector<int> send_count_;
    std::vector<int> send_offset_;
    std::vector<int> recv_count_;
    std::vector<int> recv_offset_;
    // the ABACUS planes of this process as [z][xy], and the libvdwxc slab
    std::vector<double> planes_;
    std::vector<double> rho_;
    std::vector<double> sigma_;
    std::vector<double> dedn_;
    std::vector<double> dedsigma_;
};

} // namespace vdw

#endif // RVV10_WORKSPACE_H
