#ifndef ABACUS_SOURCE_SOURCE_BASE_BLOCK_CYCLIC_DISTRIBUTION_INFORMATION_H
#define ABACUS_SOURCE_SOURCE_BASE_BLOCK_CYCLIC_DISTRIBUTION_INFORMATION_H

#include <array>
#include <cstddef>
#include <vector>

#ifdef __MPI
#include <mpi.h>
#endif

// Responsibilities:
// 1. Describe a block-cyclic layout: process grid, block size, local shape,
//    global/local index maps, and owner queries.
// 2. Own or borrow the BLACS context and maintain the ScaLAPACK descriptor.
// It does not own physical matrices, atom ownership, or trace information.
class BlockCyclicDistributionInformation
{
  public:
    BlockCyclicDistributionInformation() = default;
    ~BlockCyclicDistributionInformation();

    BlockCyclicDistributionInformation(const BlockCyclicDistributionInformation&) = delete;
    BlockCyclicDistributionInformation& operator=(const BlockCyclicDistributionInformation&) = delete;
    BlockCyclicDistributionInformation(BlockCyclicDistributionInformation&& rhs) noexcept;
    BlockCyclicDistributionInformation& operator=(BlockCyclicDistributionInformation&& rhs) noexcept;

    int nrloc() const
    {
        return nrloc_;
    }

    int ncloc() const
    {
        return ncloc_;
    }

    int nrglb() const;
    int ncglb() const;

    std::size_t size() const
    {
        return size_;
    }

    int nb() const
    {
        return nb_;
    }

    int nrproc() const
    {
        return nrproc_;
    }

    int ncproc() const
    {
        return ncproc_;
    }

    int irproc() const
    {
        return irproc_;
    }

    int icproc() const
    {
        return icproc_;
    }

    int irloc(const int irglb) const
    {
        return irglb2loc_[irglb];
    }

    int icloc(const int icglb) const
    {
        return icglb2loc_[icglb];
    }

    int irglb(const int irloc) const
    {
        return irloc2glb_[irloc];
    }

    int icglb(const int icloc) const
    {
        return icloc2glb_[icloc];
    }

    bool in_this_processor(const int irglb, const int icglb) const;
    int owner_processor(const int irglb, const int icglb) const;

    bool blacs_in_this_processor(const int irproc, const int icproc) const
    {
        return irproc == irproc_ && icproc == icproc_;
    }

    void set_iproc(const int irproc, const int icproc)
    {
        irproc_ = irproc;
        icproc_ = icproc;
    }

    bool is_serial() const
    {
        return is_serial_;
    }

    void set_serial(const int nrglb, const int ncglb);

#ifdef __MPI
    int bind_descriptor(const int nrglb,
                        const int ncglb,
                        const int nb,
                        const int ictxt);

    const int* desc() const
    {
        return desc_.data();
    }

    int ictxt() const
    {
        return ictxt_;
    }

    MPI_Comm comm() const;

    int* desc()
    {
        return desc_.data();
    }

    int init(const int nrglb,
             const int ncglb,
             const int nb,
             const MPI_Comm comm,
             const bool swap_grid);
    int set(const int nrglb,
            const int ncglb,
            const int nb,
            const int ictxt);
#endif

  private:
    static constexpr int kDescriptorEntryCount = 9;
    static constexpr int kDescriptorBlacsContextIndex = 1;
    static constexpr int kDescriptorGlobalRowsIndex = 2;
    static constexpr int kDescriptorGlobalColsIndex = 3;
    static constexpr int kDescriptorLeadingDimIndex = 8;
    static constexpr int kUnownedBlacsContext = -1;
    static constexpr int kBlacsWhatSystemContext = 10;

    int nrloc_ = 0;
    int ncloc_ = 0;
    std::size_t size_ = 0;
    int nb_ = 1;
    int nrproc_ = 0;
    int ncproc_ = 0;

#ifdef __MPI
    int ictxt_ = kUnownedBlacsContext;
    std::array<int, kDescriptorEntryCount> desc_ = {};
#endif

    std::vector<int> irglb2loc_;
    std::vector<int> icglb2loc_;
    std::vector<int> irloc2glb_;
    std::vector<int> icloc2glb_;
    int irproc_ = kUnownedBlacsContext;
    int icproc_ = kUnownedBlacsContext;
    bool is_serial_ = false;

#ifdef __MPI
    bool owns_ictxt_ = false;

    void init_grid(const MPI_Comm comm, const bool swap_grid);
    void init_dist(const int nrglb, const int ncglb, const int nb);
    void release_grid();
#endif
};

#endif
