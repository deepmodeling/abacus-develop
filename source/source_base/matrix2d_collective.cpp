#include "source_base/matrix2d_collective.h"

#ifdef __MPI
#include "source_base/matrix2d_mpi.h"
#include "source_base/parallel_common.h"

#include <cstring>
#include <stdexcept>
#include <type_traits>
#include <vector>

/**
 * @internal
 * @brief Communication adapters for Matrix2d distributed storage.
 *
 * ScaLAPACK has no drop-in gather/scatter operation for arbitrary frozen
 * Matrix2d layouts. Block ownership and BLACS descriptors remain defined by
 * `Parallel_Blkcyclic` through `Cblacs_gridinit`, `numroc_`, and `descinit_`.
 * The code below only serializes values by the existing owner maps and moves
 * the resulting contiguous rank segments with wrapper-managed collective
 * operations.
 */
namespace
{
namespace mpi = matrix2d_impl::mpi;

constexpr int LOCAL_METADATA_SIZE = 6;
// Width in bytes of a native int in the homogeneous metadata wire format.
constexpr int W_INT = sizeof(int);

/**
 * @brief Pack rank-local axis metadata for Matrix2d collective transfers.
 *
 * The metadata wire format is
 * `[nrloc, ncloc, nrglb, ncglb, iprow, ipcol,
 *   irglb[nrloc], icglb[ncloc]]`.
 * `make_local_buffer()` appends `nrloc * ncloc` logical row-major `T` values
 * with no BLACS leading-dimension padding. Every integer is a native `int`.
 * This byte-level protocol is valid only among homogeneous MPI processes that
 * use the same C++ ABI, scalar representation, and native `int` layout.
 */
template <typename T>
std::vector<int> make_local_metadata(const Matrix2d<T>& matrix)
{
    static_assert(std::is_trivially_copyable<T>::value,
                  "Matrix2d collective buffers require trivially copyable T");
    const std::shared_ptr<const Parallel_Blkcyclic> distrib
        = matrix.distribution();
    const std::size_t nmeta = static_cast<std::size_t>(LOCAL_METADATA_SIZE)
        + static_cast<std::size_t>(matrix.nrloc())
        + static_cast<std::size_t>(matrix.ncloc());
    std::vector<int> metadata(nmeta);
    std::size_t offset = 0;
    metadata[offset++] = matrix.nrloc();
    metadata[offset++] = matrix.ncloc();
    metadata[offset++] = matrix.nrglb();
    metadata[offset++] = matrix.ncglb();
    metadata[offset++] = distrib->irproc();
    metadata[offset++] = distrib->icproc();
    for (int irloc = 0; irloc < matrix.nrloc(); ++irloc)
    {
        metadata[offset++] = distrib->irglb(irloc);
    }
    for (int icloc = 0; icloc < matrix.ncloc(); ++icloc)
    {
        metadata[offset++] = distrib->icglb(icloc);
    }
    return metadata;
}

template <typename T>
std::vector<unsigned char> make_local_buffer(const Matrix2d<T>& matrix)
{
    const std::vector<int> metadata = make_local_metadata(matrix);
    const std::size_t nmeta = metadata.size();
    const std::size_t nbytes_meta = nmeta * W_INT;
    // Number of stored T values, not bytes.
    const std::size_t ndata = static_cast<std::size_t>(matrix.nrloc())
        * static_cast<std::size_t>(matrix.ncloc());
    const std::size_t nbytes_data = ndata * sizeof(T);
    const std::size_t nbytbuf = nbytes_meta + nbytes_data;
    std::vector<unsigned char> buffer(nbytbuf);

    std::memcpy(buffer.data(), metadata.data(), nbytes_meta);
    for (int irloc = 0; irloc < matrix.nrloc(); ++irloc)
    {
        for (int icloc = 0; icloc < matrix.ncloc(); ++icloc)
        {
            // The u prefix denotes a checked, unsigned local index.
            const std::size_t uirloc = static_cast<std::size_t>(irloc);
            const std::size_t uicloc = static_cast<std::size_t>(icloc);
            const std::size_t idx
                = uirloc * static_cast<std::size_t>(matrix.ncloc()) + uicloc;
            const T& value = matrix.atloc(irloc, icloc);
            const std::size_t offset = nbytes_meta + idx * sizeof(T);
            unsigned char* const destination
                = buffer.data() + offset;
            std::memcpy(destination, &value, sizeof(T));
        }
    }
    return buffer;
}

struct RankBufferView
{
    const unsigned char* data = nullptr;
    std::size_t size = 0;
};

struct LocalBufferView
{
    int nrloc = 0;
    int ncloc = 0;
    int nrglb = 0;
    int ncglb = 0;
    int iprow = 0;
    int ipcol = 0;
    std::vector<int> irglb;
    std::vector<int> icglb;
    const unsigned char* data = nullptr;
};

/**
 * @brief Pack values owned by one rank into a contiguous scatter segment.
 *
 * The block-cyclic owner map is defined by `Parallel_Blkcyclic` through the
 * BLACS/ScaLAPACK grid setup (`Cblacs_gridinit`, `numroc_`, and `descinit_`).
 * This helper does not reproduce that distribution algorithm; it only moves
 * values according to the already-computed local-to-global axis maps.
 */
template <typename T>
void pack_rank_values(const Matrix2d<T>& matglb,
                      const LocalBufferView& view,
                      unsigned char* destination)
{
    for (int irloc = 0; irloc < view.nrloc; ++irloc)
    {
        const int irglb = view.irglb[irloc];
        for (int icloc = 0; icloc < view.ncloc; ++icloc)
        {
            const int icglb = view.icglb[icloc];
            const std::size_t idx
                = static_cast<std::size_t>(irloc)
                    * static_cast<std::size_t>(view.ncloc)
                + static_cast<std::size_t>(icloc);
            const T& value = matglb.atloc(irglb, icglb);
            const std::size_t destination_offset = idx * sizeof(T);
            unsigned char* const destination_offset_ptr
                = destination + destination_offset;
            std::memcpy(destination_offset_ptr, &value, sizeof(T));
        }
    }
}

template <typename T>
LocalBufferView view_local_buffer(
    const RankBufferView& buffer)
{
    if (buffer.size < LOCAL_METADATA_SIZE * W_INT)
    {
        throw std::invalid_argument(
            "Matrix2d collective matrix buffer is truncated");
    }

    int header[LOCAL_METADATA_SIZE] = {};
    std::memcpy(header, buffer.data, sizeof(header));
    LocalBufferView view;
    view.nrloc = header[0];
    view.ncloc = header[1];
    view.nrglb = header[2];
    view.ncglb = header[3];
    view.iprow = header[4];
    view.ipcol = header[5];

    if (view.nrloc < 0 || view.ncloc < 0 || view.nrglb <= 0
        || view.ncglb <= 0)
    {
        throw std::invalid_argument(
            "Matrix2d collective matrix buffer has invalid dimensions");
    }

    const std::size_t nmeta
        = static_cast<std::size_t>(LOCAL_METADATA_SIZE) + view.nrloc + view.ncloc;
    const std::size_t nbytes_meta = nmeta * W_INT;
    const std::size_t ndata
        = static_cast<std::size_t>(view.nrloc) * static_cast<std::size_t>(view.ncloc);
    const std::size_t nbytes_expected = nbytes_meta + ndata * sizeof(T);
    if (buffer.size != nbytes_expected)
    {
        throw std::invalid_argument(
            "Matrix2d collective matrix buffer has an invalid size");
    }

    std::vector<int> metadata(nmeta);
    std::memcpy(metadata.data(), buffer.data, nbytes_meta);
    view.irglb.assign(metadata.begin() + LOCAL_METADATA_SIZE,
                      metadata.begin() + LOCAL_METADATA_SIZE + view.nrloc);
    view.icglb.assign(metadata.begin() + LOCAL_METADATA_SIZE + view.nrloc,
                      metadata.end());
    view.data = buffer.data + nbytes_meta;
    return view;
}

std::vector<RankBufferView> make_rank_buffer_views(
    const std::vector<unsigned char>& recvbuf,
    const std::vector<int>& nbytes_recv,
    const std::vector<int>& rdispls)
{
    std::vector<RankBufferView> buffers(nbytes_recv.size());
    for (std::size_t ip = 0; ip < nbytes_recv.size(); ++ip)
    {
        const std::size_t rank_begin
            = static_cast<std::size_t>(rdispls[ip]);
        const std::size_t rank_size
            = static_cast<std::size_t>(nbytes_recv[ip]);
        if (rank_begin > recvbuf.size()
            || rank_size > recvbuf.size() - rank_begin)
        {
            throw std::invalid_argument(
                "Matrix2d collective rank segment is outside the receive buffer");
        }
        buffers[ip].data = recvbuf.data() + rank_begin;
        buffers[ip].size = rank_size;
    }
    return buffers;
}

template <typename T>
Matrix2d<T> rebuild_global(
    const std::vector<RankBufferView>& buffers)
{
    std::vector<LocalBufferView> views;
    views.reserve(buffers.size());
    int nrglb = 0;
    int ncglb = 0;
    for (std::size_t ip = 0; ip < buffers.size(); ++ip)
    {
        views.push_back(view_local_buffer<T>(buffers[ip]));
        const LocalBufferView& view = views.back();
        if (ip == 0)
        {
            nrglb = view.nrglb;
            ncglb = view.ncglb;
        }
        else if (view.nrglb != nrglb || view.ncglb != ncglb)
        {
            throw std::invalid_argument("Matrix2d ranks disagree on global shape");
        }
    }

    Matrix2d<T> result(nrglb, ncglb);
    for (std::size_t ip = 0; ip < views.size(); ++ip)
    {
        const LocalBufferView& view = views[ip];
        for (int irloc = 0; irloc < view.nrloc; ++irloc)
        {
            const int irglb = view.irglb[irloc];
            if (irglb < 0 || irglb >= nrglb)
            {
                throw std::invalid_argument("Matrix2d row axis is outside global shape");
            }
            for (int icloc = 0; icloc < view.ncloc; ++icloc)
            {
                const int icglb = view.icglb[icloc];
                if (icglb < 0 || icglb >= ncglb)
                {
                    throw std::invalid_argument("Matrix2d column axis is outside global shape");
                }
                const std::size_t idx
                    = static_cast<std::size_t>(irloc)
                        * static_cast<std::size_t>(view.ncloc)
                    + static_cast<std::size_t>(icloc);
                const std::size_t source_offset = idx * sizeof(T);
                const unsigned char* const source
                    = view.data + source_offset;
                T value;
                std::memcpy(&value, source, sizeof(T));
                result.atloc(irglb, icglb) = value;
            }
        }
    }
    return result;
}

} // namespace

namespace matrix2d_impl
{
template <typename T>
Matrix2d<T> allgather(const Matrix2d<T>& matrix, const MPI_Comm comm)
{
    const std::vector<unsigned char> locbuf = make_local_buffer(matrix);
    const int nbytes_local_mpi = mpi::checked_mpi_count(locbuf.size());
    const std::vector<int> nbytes_recv
        = mpi::allgather_nbytes(nbytes_local_mpi, comm);
    const std::vector<int> rdispls = mpi::make_byte_offsets(nbytes_recv);
    const int nbytes_total = rdispls.back() + nbytes_recv.back();
    std::vector<unsigned char> recvbuf(
        static_cast<std::size_t>(nbytes_total));
    Parallel_Common::allgatherv_bytes(locbuf.data(),
                                      nbytes_local_mpi,
                                      recvbuf.data(),
                                      nbytes_recv.data(),
                                      rdispls.data(),
                                      comm);
    const std::vector<RankBufferView> buffers
        = make_rank_buffer_views(recvbuf, nbytes_recv, rdispls);
    return rebuild_global<T>(buffers);
}

template <typename T>
std::unique_ptr<Matrix2d<T>> gather_on_root(const Matrix2d<T>& matrix,
                                            const int iproot,
                                            const MPI_Comm comm)
{
    const int nproc = Parallel_Common::comm_size(comm);
    const int ip = Parallel_Common::comm_rank(comm);
    if (iproot < 0 || iproot >= nproc)
    {
        throw std::invalid_argument("gather root rank is outside the communicator");
    }

    const std::vector<unsigned char> locbuf = make_local_buffer(matrix);
    const int nbytes_local_mpi = mpi::checked_mpi_count(locbuf.size());
    const std::vector<int> nbytes_recv
        = mpi::allgather_nbytes(nbytes_local_mpi, comm);
    const std::vector<int> rdispls = mpi::make_byte_offsets(nbytes_recv);
    const int nbytes_total = rdispls.back() + nbytes_recv.back();
    std::vector<unsigned char> rxbuf;
    if (ip == iproot)
    {
        rxbuf.resize(nbytes_total);
    }

    unsigned char* const recvbuf = ip == iproot ? rxbuf.data() : nullptr;
    Parallel_Common::gatherv_bytes(locbuf.data(),
                                   nbytes_local_mpi,
                                   recvbuf,
                                   nbytes_recv.data(),
                                   rdispls.data(),
                                   iproot,
                                   comm);
    if (ip != iproot)
    {
        return std::unique_ptr<Matrix2d<T>>();
    }

    const std::vector<RankBufferView> buffers
        = make_rank_buffer_views(rxbuf, nbytes_recv, rdispls);
    return std::unique_ptr<Matrix2d<T>>(
        new Matrix2d<T>(rebuild_global<T>(buffers)));
}

template <typename T>
void scatter_from_host(const Matrix2d<T>& matglb, Matrix2d<T>& matloc)
{
    if (matglb.has_distribution())
    {
        throw std::invalid_argument("scatter requires a host global matrix");
    }

    const std::shared_ptr<const Parallel_Blkcyclic> distrib
        = matloc.distribution();
    const MPI_Comm comm = distrib->comm();
    if (comm == MPI_COMM_NULL)
    {
        throw std::invalid_argument("scatter requires a valid MPI comm");
    }

    const std::vector<int> metadata = make_local_metadata(matloc);
    const std::size_t nmeta = metadata.size();
    const std::size_t nbytes_meta = nmeta * W_INT;
    const int nmeta_mpi = mpi::checked_mpi_count(nbytes_meta);
    const std::vector<int> metacounts
        = mpi::allgather_nbytes(nmeta_mpi, comm);
    const std::vector<int> metadispls
        = mpi::make_byte_offsets(metacounts);
    const int metadata_total_bytes
        = metadispls.back() + metacounts.back();
    if (metadata_total_bytes % W_INT != 0)
    {
        throw std::invalid_argument("Matrix2d scatter metadata has an invalid byte count");
    }
    const int nmeta_all = metadata_total_bytes / W_INT;
    std::vector<int> all_metadata(nmeta_all);

    Parallel_Common::allgatherv_bytes(metadata.data(),
                                      nmeta_mpi,
                                      all_metadata.data(),
                                      metacounts.data(),
                                      metadispls.data(),
                                      comm);

    const int nproc = Parallel_Common::comm_size(comm);
    const int ip = Parallel_Common::comm_rank(comm);
    std::vector<int> sendcounts(nproc, 0);
    std::vector<LocalBufferView> views;
    views.reserve(nproc);
    std::size_t metadata_offset = 0;
    const int nrglb = matloc.nrglb();
    const int ncglb = matloc.ncglb();
    for (int ipsrc = 0; ipsrc < nproc; ++ipsrc)
    {
        if (metadata_offset + LOCAL_METADATA_SIZE > all_metadata.size())
        {
            throw std::invalid_argument("Matrix2d scatter metadata is truncated");
        }
        LocalBufferView view;
        view.nrloc = all_metadata[metadata_offset];
        view.ncloc = all_metadata[metadata_offset + 1];
        view.nrglb = all_metadata[metadata_offset + 2];
        view.ncglb = all_metadata[metadata_offset + 3];
        view.iprow = all_metadata[metadata_offset + 4];
        view.ipcol = all_metadata[metadata_offset + 5];
        metadata_offset += LOCAL_METADATA_SIZE;
        const std::size_t naxes
            = static_cast<std::size_t>(view.nrloc + view.ncloc);
        if (metadata_offset + naxes > all_metadata.size())
        {
            throw std::invalid_argument("Matrix2d scatter axes are truncated");
        }
        view.irglb.assign(all_metadata.begin() + metadata_offset,
                          all_metadata.begin() + metadata_offset + view.nrloc);
        view.icglb.assign(all_metadata.begin() + metadata_offset + view.nrloc,
                          all_metadata.begin() + metadata_offset + view.nrloc
                              + view.ncloc);
        metadata_offset += naxes;

        if (view.nrloc < 0 || view.ncloc < 0
            || view.nrglb != nrglb || view.ncglb != ncglb)
        {
            throw std::invalid_argument("Matrix2d ranks disagree on scatter global shape");
        }
        const std::size_t nelem_local
            = static_cast<std::size_t>(view.nrloc)
            * static_cast<std::size_t>(view.ncloc);
        const std::size_t nbytes_local = nelem_local * sizeof(T);
        sendcounts[ipsrc] = mpi::checked_mpi_count(nbytes_local);
        views.push_back(view);
    }
    if (metadata_offset != all_metadata.size())
    {
        throw std::invalid_argument("Matrix2d scatter metadata has trailing entries");
    }
    if (matglb.nrglb() != nrglb || matglb.ncglb() != ncglb)
    {
        throw std::invalid_argument("scatter global matrix shape does not match distribution");
    }

    const std::vector<int> sdispls = mpi::make_byte_offsets(sendcounts);
    std::vector<unsigned char> sendbuf;
    if (ip == 0)
    {
        const int total_bytes = sdispls.back() + sendcounts.back();
        sendbuf.resize(total_bytes);
        for (int ipdst = 0; ipdst < nproc; ++ipdst)
        {
            const LocalBufferView& view = views[ipdst];
            unsigned char* const destination
                = sendbuf.data() + static_cast<std::size_t>(sdispls[ipdst]);
            pack_rank_values(matglb, view, destination);
        }
    }

    const std::size_t nelem_recv
        = static_cast<std::size_t>(matloc.nrloc())
        * static_cast<std::size_t>(matloc.ncloc());
    const std::size_t nbytes_recv = nelem_recv * sizeof(T);
    std::vector<unsigned char> recvbuf(nbytes_recv);
    const int recv = mpi::checked_mpi_count(nbytes_recv);
    const unsigned char* receive_source = ip == 0 ? sendbuf.data() : nullptr;
    Parallel_Common::scatterv_bytes(receive_source,
                                    sendcounts.data(),
                                    sdispls.data(),
                                    recvbuf.data(),
                                    recv,
                                    0,
                                    comm);

    for (int irloc = 0; irloc < matloc.nrloc(); ++irloc)
    {
        for (int icloc = 0; icloc < matloc.ncloc(); ++icloc)
        {
            const std::size_t idx
                = static_cast<std::size_t>(irloc)
                    * static_cast<std::size_t>(matloc.ncloc())
                + static_cast<std::size_t>(icloc);
            const std::size_t receive_offset = idx * sizeof(T);
            const unsigned char* const source
                = recvbuf.data() + receive_offset;
            T value;
            std::memcpy(&value, source, sizeof(T));
            matloc.atloc(irloc, icloc) = value;
        }
    }
}

template Matrix2d<int> allgather<int>(const Matrix2d<int>& matrix,
                                       const MPI_Comm comm);
template Matrix2d<float> allgather<float>(const Matrix2d<float>& matrix,
                                           const MPI_Comm comm);
template Matrix2d<double> allgather<double>(const Matrix2d<double>& matrix,
                                             const MPI_Comm comm);
template Matrix2d<std::complex<float>> allgather<std::complex<float>>(
    const Matrix2d<std::complex<float>>& matrix,
    const MPI_Comm comm);
template Matrix2d<std::complex<double>> allgather<std::complex<double>>(
    const Matrix2d<std::complex<double>>& matrix,
    const MPI_Comm comm);

template std::unique_ptr<Matrix2d<int>> gather_on_root<int>(
    const Matrix2d<int>& matrix,
    const int iproot,
    const MPI_Comm comm);
template std::unique_ptr<Matrix2d<float>> gather_on_root<float>(
    const Matrix2d<float>& matrix,
    const int iproot,
    const MPI_Comm comm);
template std::unique_ptr<Matrix2d<double>> gather_on_root<double>(
    const Matrix2d<double>& matrix,
    const int iproot,
    const MPI_Comm comm);
template std::unique_ptr<Matrix2d<std::complex<float>>> gather_on_root<std::complex<float>>(
    const Matrix2d<std::complex<float>>& matrix,
    const int iproot,
    const MPI_Comm comm);
template std::unique_ptr<Matrix2d<std::complex<double>>> gather_on_root<std::complex<double>>(
    const Matrix2d<std::complex<double>>& matrix,
    const int iproot,
    const MPI_Comm comm);

template void scatter_from_host<int>(const Matrix2d<int>& matglb,
                                     Matrix2d<int>& matloc);
template void scatter_from_host<float>(const Matrix2d<float>& matglb,
                                       Matrix2d<float>& matloc);
template void scatter_from_host<double>(const Matrix2d<double>& matglb,
                                        Matrix2d<double>& matloc);
template void scatter_from_host<std::complex<float>>(
    const Matrix2d<std::complex<float>>& matglb,
    Matrix2d<std::complex<float>>& matloc);
template void scatter_from_host<std::complex<double>>(
    const Matrix2d<std::complex<double>>& matglb,
    Matrix2d<std::complex<double>>& matloc);
} // namespace matrix2d_impl
#endif
