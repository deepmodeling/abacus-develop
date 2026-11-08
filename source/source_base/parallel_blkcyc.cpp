#include "source_base/parallel_blkcyc.h"

#include "source_base/module_external/blacs_connector.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_base/parallel_common.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>
#include <utility>

BlockCyclicDistribution::~BlockCyclicDistribution()
{
#ifdef __MPI
    release_grid();
#endif
}

BlockCyclicDistribution::BlockCyclicDistribution(
    BlockCyclicDistribution&& rhs) noexcept
{
    *this = std::move(rhs);
}

BlockCyclicDistribution& BlockCyclicDistribution::operator=(
    BlockCyclicDistribution&& rhs) noexcept
{
    if (this != &rhs)
    {
#ifdef __MPI
        release_grid();
        owns_ictxt_ = rhs.owns_ictxt_;
        rhs.owns_ictxt_ = false;
        ictxt_ = rhs.ictxt_;
        rhs.ictxt_ = UNOWNED_BLACS_CONTEXT;
#endif
        nrloc_ = rhs.nrloc_;
        ncloc_ = rhs.ncloc_;
        size_ = rhs.size_;
        nb_ = rhs.nb_;
        nrproc_ = rhs.nrproc_;
        ncproc_ = rhs.ncproc_;
        irproc_ = rhs.irproc_;
        icproc_ = rhs.icproc_;
        is_serial_ = rhs.is_serial_;
        irglb2loc_ = std::move(rhs.irglb2loc_);
        icglb2loc_ = std::move(rhs.icglb2loc_);
        irloc2glb_ = std::move(rhs.irloc2glb_);
        icloc2glb_ = std::move(rhs.icloc2glb_);
#ifdef __MPI
        std::copy(rhs.desc_.begin(), rhs.desc_.end(), desc_.begin());
#endif
    }
    return *this;
}

bool BlockCyclicDistribution::in_this_processor(const int irglb, const int icglb) const
{
    return irloc(irglb) != -1 && icloc(icglb) != -1;
}

int BlockCyclicDistribution::owner_processor(const int irglb, const int icglb) const
{
    assert(irglb >= 0 && irglb < nrglb() && icglb >= 0 && icglb < ncglb());
    if (is_serial_)
    {
        return 0;
    }
    const int irproc = (irglb / nb_) % nrproc_;
    const int icproc = (icglb / nb_) % ncproc_;
    return irproc * ncproc_ + icproc;
}

int BlockCyclicDistribution::nrglb() const
{
    if (!is_serial_)
    {
#ifdef __MPI
        return desc_[DESC_M];
#endif
    }
    return nrloc_;
}

int BlockCyclicDistribution::ncglb() const
{
    if (!is_serial_)
    {
#ifdef __MPI
        return desc_[DESC_N];
#endif
    }
    return ncloc_;
}

#ifdef __MPI
MPI_Comm BlockCyclicDistribution::comm() const
{
    // it is an error to call blacs_get with an invalid BLACS context
    if (ictxt_ < 0)
    {
        return MPI_COMM_NULL;
    }

    int sys_ctxt = 0;
    Cblacs_get(ictxt_, BLACS_WHAT_SYSTEM_CONTEXT, &sys_ctxt);
    // blacs_get with "what" = 10 takes a BLACS context and returns the index
    // of the associated system context (MPI communicator) that can be used by
    // blacs2sys_handle to get the MPI communicator.
    return Cblacs2sys_handle(sys_ctxt);
}

void BlockCyclicDistribution::init_grid(const MPI_Comm comm, const bool swap_grid)
{
    // determine the number of rows and columns of the process grid
    // by factorizing n = p * q such that p, q are closest and p <= q
    const int num_proc = Parallel_Common::comm_size(comm);
    nrproc_ = static_cast<int>(std::sqrt(num_proc + 0.5));
    while (ncproc_ = num_proc / nrproc_, nrproc_ * ncproc_ != num_proc)
    {
        --nrproc_;
    }

    if (swap_grid)
    {
        std::swap(nrproc_, ncproc_);
    }

    // initialize the BLACS grid accordingly
    ictxt_ = Csys2blacs_handle(comm);
    char order = 'R'; // row-major
    Cblacs_gridinit(&ictxt_, &order, nrproc_, ncproc_);
    Cblacs_gridinfo(ictxt_, &nrproc_, &ncproc_, &irproc_, &icproc_);
}

void BlockCyclicDistribution::release_grid()
{
    if (owns_ictxt_ && ictxt_ >= 0)
    {
        Cblacs_gridexit(ictxt_);
        ictxt_ = UNOWNED_BLACS_CONTEXT;
        owns_ictxt_ = false;
    }
}

void BlockCyclicDistribution::init_dist(const int nrglb,
                                                   const int ncglb,
                                                   const int nb)
{
    nb_ = nb;

    // number of local rows and columns
    const int zero = 0;
    int nrglb_arg = nrglb;
    int ncglb_arg = ncglb;
    int nb_arg = nb;
    nrloc_ = numroc_(&nrglb_arg, &nb_arg, &irproc_, &zero, &nrproc_);
    ncloc_ = numroc_(&ncglb_arg, &nb_arg, &icproc_, &zero, &ncproc_);
    size_ = static_cast<std::size_t>(nrloc_) * static_cast<std::size_t>(ncloc_);

    // initialize the ScaLAPACK descriptor
    int info = 0;
    int lld = std::max(nrloc_, 1);
    int desc_nb = nb_;
    int desc_ictxt = ictxt_;
    descinit_(desc_.data(),
              &nrglb_arg,
              &ncglb_arg,
              &desc_nb,
              &desc_nb,
              &zero,
              &zero,
              &desc_ictxt,
              &lld,
              &info);

    // generate the global-to-local and local-to-global index maps
    irloc2glb_.resize(nrloc_);
    irglb2loc_ = std::vector<int>(nrglb, -1);
    for (int i = 0; i < nrloc_; ++i)
    {
        irloc2glb_[i] = (i / nb * nrproc_ + irproc_) * nb + i % nb;
        irglb2loc_[irloc2glb_[i]] = i;
    }

    icloc2glb_.resize(ncloc_);
    icglb2loc_ = std::vector<int>(ncglb, -1);
    for (int j = 0; j < ncloc_; ++j)
    {
        icloc2glb_[j] = (j / nb * ncproc_ + icproc_) * nb + j % nb;
        icglb2loc_[icloc2glb_[j]] = j;
    }
}

int BlockCyclicDistribution::bind_descriptor(const int nrglb,
                                                        const int ncglb,
                                                        const int nb,
                                                        const int ictxt)
{
    // Reusing the context owned by this object must not destroy the grid or
    // drop ownership; otherwise the grid about to be reused is released.
    if (ictxt == ictxt_ && owns_ictxt_)
    {
        Cblacs_gridinfo(ictxt, &nrproc_, &ncproc_, &irproc_, &icproc_);
        init_dist(nrglb, ncglb, nb);
        return nrloc_ == 0 || ncloc_ == 0;
    }
    release_grid();
    ictxt_ = ictxt;
    owns_ictxt_ = false;
    Cblacs_gridinfo(ictxt, &nrproc_, &ncproc_, &irproc_, &icproc_);
    init_dist(nrglb, ncglb, nb);
    return nrloc_ == 0 || ncloc_ == 0;
}

int BlockCyclicDistribution::init(const int nrglb,
                                             const int ncglb,
                                             const int nb,
                                             const MPI_Comm comm,
                                             const bool swap_grid)
{
    release_grid();
    init_grid(comm, swap_grid);
    owns_ictxt_ = true;
    init_dist(nrglb, ncglb, nb);
    return nrloc_ == 0 || ncloc_ == 0;
}

int BlockCyclicDistribution::set(const int nrglb,
                                            const int ncglb,
                                            const int nb,
                                            const int ictxt)
{
    return bind_descriptor(nrglb, ncglb, nb, ictxt);
}
#endif

void BlockCyclicDistribution::set_serial(const int nrglb, const int ncglb)
{
    assert(nrglb > 0 && ncglb > 0);

    nb_ = 1;
    nrproc_ = 1;
    ncproc_ = 1;
    irproc_ = 0;
    icproc_ = 0;
    nrloc_ = nrglb;
    ncloc_ = ncglb;
    size_ = static_cast<std::size_t>(nrloc_) * static_cast<std::size_t>(ncloc_);
    irloc2glb_.resize(nrloc_);
    icloc2glb_.resize(ncloc_);
    std::iota(irloc2glb_.begin(), irloc2glb_.end(), 0);
    std::iota(icloc2glb_.begin(), icloc2glb_.end(), 0);
    irglb2loc_ = irloc2glb_;
    icglb2loc_ = icloc2glb_;
    is_serial_ = true;
#ifdef __MPI
    release_grid();
    // A serial layout must not reference any BLACS grid: for borrowers
    // release_grid() is a no-op, so clear the handle explicitly;
    // otherwise the object keeps a reference to the old grid, which also
    // dangles once the owner destroys it.
    ictxt_ = UNOWNED_BLACS_CONTEXT;
#endif
}
