#include "source_base/matrix2d_collective.h"

#ifdef __MPI
#include "source_base/parallel_common.h"

#include <cstring>
#include <limits>
#include <stdexcept>
#include <vector>

namespace
{
constexpr int LOCAL_METADATA_SIZE = 6;

void assert_valid_mpi_count(const std::size_t nbytes)
{
    if (nbytes > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    {
        throw std::invalid_argument(
            "Matrix2d collective matrix buffer exceeds MPI count range");
    }
}

int mpi_count(const std::size_t nbytes)
{
    return static_cast<int>(nbytes);
}

std::vector<int> make_byte_displs(const std::vector<int>& nbytes)
{
    std::vector<int> displs(nbytes.size(), 0);
    int offset = 0;
    for (std::size_t ip = 0; ip < nbytes.size(); ++ip)
    {
        displs[ip] = offset;
        if (nbytes[ip] > std::numeric_limits<int>::max() - offset)
        {
            throw std::invalid_argument(
                "Matrix2d collective matrix buffer size overflow");
        }
        offset += nbytes[ip];
    }
    return displs;
}

std::vector<int> allgather_nbytes(const int nbytes_local, const MPI_Comm comm)
{
    const int nproc = Parallel_Common::comm_size(comm);
    const int nbytes_int = static_cast<int>(sizeof(int));
    std::vector<int> recvcounts(static_cast<std::size_t>(nproc), nbytes_int);
    const std::vector<int> rdispls
        = make_byte_displs(recvcounts);
    std::vector<int> nbytes(static_cast<std::size_t>(nproc), 0);

    Parallel_Common::allgatherv_bytes(&nbytes_local,
                                      nbytes_int,
                                      nbytes.data(),
                                      recvcounts.data(),
                                      rdispls.data(),
                                      comm);
    return nbytes;
}

template <typename T>
std::vector<int> make_local_metadata(const Matrix2d<T>& matrix)
{
    const std::shared_ptr<const Parallel_Blkcyclic> distrib
        = matrix.distribution();
    std::vector<int> metadata;
    metadata.reserve(LOCAL_METADATA_SIZE + matrix.nrloc()
                     + matrix.ncloc());
    metadata.push_back(matrix.nrloc());
    metadata.push_back(matrix.ncloc());
    metadata.push_back(matrix.nrglb());
    metadata.push_back(matrix.ncglb());
    metadata.push_back(distrib->irproc());
    metadata.push_back(distrib->icproc());
    for (int irloc = 0; irloc < matrix.nrloc(); ++irloc)
    {
        metadata.push_back(distrib->irglb(irloc));
    }
    for (int icloc = 0; icloc < matrix.ncloc(); ++icloc)
    {
        metadata.push_back(distrib->icglb(icloc));
    }
    return metadata;
}

template <typename T>
std::vector<unsigned char> make_local_matrix_buffer(const Matrix2d<T>& matrix)
{
    const std::vector<int> metadata = make_local_metadata(matrix);
    const std::size_t nmeta = metadata.size();
    const std::size_t nmeta_bytes = nmeta * sizeof(int);
    const std::size_t data_elements = static_cast<std::size_t>(matrix.nrloc())
        * static_cast<std::size_t>(matrix.ncloc());
    const std::size_t data_bytes = data_elements * sizeof(T);
    const std::size_t nbytes_buffer = nmeta_bytes + data_bytes;
    std::vector<unsigned char> buffer(nbytes_buffer);

    std::memcpy(buffer.data(), metadata.data(), nmeta_bytes);
    for (int irloc = 0; irloc < matrix.nrloc(); ++irloc)
    {
        for (int icloc = 0; icloc < matrix.ncloc(); ++icloc)
        {
            const std::size_t uirloc = static_cast<std::size_t>(irloc);
            const std::size_t uicloc = static_cast<std::size_t>(icloc);
            const std::size_t idx
                = uirloc * static_cast<std::size_t>(matrix.ncloc()) + uicloc;
            const T& value = matrix.atloc(irloc, icloc);
            const std::size_t offset = nmeta_bytes + idx * sizeof(T);
            unsigned char* const destination
                = buffer.data() + offset;
            std::memcpy(destination, &value, sizeof(T));
        }
    }
    return buffer;
}

struct LocalMatrixBufferView
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

template <typename T>
LocalMatrixBufferView view_local_matrix_buffer(
    const std::vector<unsigned char>& buffer)
{
    if (buffer.size() < LOCAL_METADATA_SIZE * sizeof(int))
    {
        throw std::invalid_argument(
            "Matrix2d collective matrix buffer is truncated");
    }

    int header[LOCAL_METADATA_SIZE] = {};
    std::memcpy(header, buffer.data(), sizeof(header));
    LocalMatrixBufferView view;
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
    const std::size_t nmeta_bytes = nmeta * sizeof(int);
    const std::size_t data_elements
        = static_cast<std::size_t>(view.nrloc) * static_cast<std::size_t>(view.ncloc);
    const std::size_t expected_bytes = nmeta_bytes + data_elements * sizeof(T);
    if (buffer.size() != expected_bytes)
    {
        throw std::invalid_argument(
            "Matrix2d collective matrix buffer has an invalid size");
    }

    std::vector<int> metadata(nmeta);
    std::memcpy(metadata.data(), buffer.data(), nmeta_bytes);
    view.irglb.assign(metadata.begin() + LOCAL_METADATA_SIZE,
                      metadata.begin() + LOCAL_METADATA_SIZE + view.nrloc);
    view.icglb.assign(metadata.begin() + LOCAL_METADATA_SIZE + view.nrloc,
                      metadata.end());
    view.data = buffer.data() + nmeta_bytes;
    return view;
}

template <typename T>
Matrix2d<T> rebuild_global_matrix(
    const std::vector<std::vector<unsigned char>>& buffers)
{
    std::vector<LocalMatrixBufferView> views;
    views.reserve(buffers.size());
    int nrglb = 0;
    int ncglb = 0;
    for (std::size_t ip = 0; ip < buffers.size(); ++ip)
    {
        views.push_back(view_local_matrix_buffer<T>(buffers[ip]));
        const LocalMatrixBufferView& view = views.back();
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
        const LocalMatrixBufferView& view = views[ip];
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
                T value;
                std::memcpy(&value, view.data + idx * sizeof(T), sizeof(T));
                result.atloc(irglb, icglb) = value;
            }
        }
    }
    return result;
}

template <typename T>
std::vector<std::vector<unsigned char>> allgather_local_matrix_buffers(
    const Matrix2d<T>& matrix,
    const MPI_Comm comm)
{
    const std::vector<unsigned char> local_buffer
        = make_local_matrix_buffer(matrix);
    assert_valid_mpi_count(local_buffer.size());
    const int nbytes_local_mpi = mpi_count(local_buffer.size());
    const std::vector<int> recvcounts
        = allgather_nbytes(nbytes_local_mpi, comm);
    const std::vector<int> rdispls
        = make_byte_displs(recvcounts);
    const int total_bytes = rdispls.back() + recvcounts.back();
    std::vector<unsigned char> received(total_bytes);

    Parallel_Common::allgatherv_bytes(local_buffer.data(),
                                      nbytes_local_mpi,
                                      received.data(),
                                      recvcounts.data(),
                                      rdispls.data(),
                                      comm);

    std::vector<std::vector<unsigned char>> buffers(recvcounts.size());
    for (std::size_t ip = 0; ip < recvcounts.size(); ++ip)
    {
        const std::size_t begin = static_cast<std::size_t>(rdispls[ip]);
        const std::size_t end
            = begin + static_cast<std::size_t>(recvcounts[ip]);
        buffers[ip].assign(received.begin() + begin, received.begin() + end);
    }
    return buffers;
}
} // namespace

namespace matrix2d_impl
{
template <typename T>
Matrix2d<T> allgather(const Matrix2d<T>& matrix, const MPI_Comm comm)
{
    const std::vector<std::vector<unsigned char>> buffers
        = allgather_local_matrix_buffers(matrix, comm);
    return rebuild_global_matrix<T>(buffers);
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

    const std::vector<unsigned char> local_buffer
        = make_local_matrix_buffer(matrix);
    assert_valid_mpi_count(local_buffer.size());
    const int nbytes_local_mpi = mpi_count(local_buffer.size());
    const std::vector<int> recvcounts
        = allgather_nbytes(nbytes_local_mpi, comm);
    const std::vector<int> rdispls
        = make_byte_displs(recvcounts);
    const int total_bytes = rdispls.back() + recvcounts.back();
    std::vector<unsigned char> received;
    if (ip == iproot)
    {
        received.resize(total_bytes);
    }

    unsigned char* const recvbuf = ip == iproot ? received.data() : nullptr;
    Parallel_Common::gatherv_bytes(local_buffer.data(),
                                   nbytes_local_mpi,
                                   recvbuf,
                                   recvcounts.data(),
                                   rdispls.data(),
                                   iproot,
                                   comm);
    if (ip != iproot)
    {
        return std::unique_ptr<Matrix2d<T>>();
    }

    std::vector<std::vector<unsigned char>> buffers(recvcounts.size());
    for (std::size_t ipsrc = 0; ipsrc < recvcounts.size(); ++ipsrc)
    {
        const std::size_t begin
            = static_cast<std::size_t>(rdispls[ipsrc]);
        const std::size_t end
            = begin + static_cast<std::size_t>(recvcounts[ipsrc]);
        buffers[ipsrc].assign(received.begin() + begin, received.begin() + end);
    }
    return std::unique_ptr<Matrix2d<T>>(
        new Matrix2d<T>(rebuild_global_matrix<T>(buffers)));
}

template <typename T>
void scatter_from_host(const Matrix2d<T>& matglb, Matrix2d<T>& matloc)
{
    if (matglb.is_distributed())
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
    const std::size_t nmeta_bytes = nmeta * sizeof(int);
    assert_valid_mpi_count(nmeta_bytes);
    const int nmeta_mpi = mpi_count(nmeta_bytes);
    const std::vector<int> metacounts = allgather_nbytes(nmeta_mpi, comm);
    const std::vector<int> metadispls
        = make_byte_displs(metacounts);
    const int metadata_total_bytes
        = metadispls.back() + metacounts.back();
    const int nbytes_int = static_cast<int>(sizeof(int));
    if (metadata_total_bytes % nbytes_int != 0)
    {
        throw std::invalid_argument("Matrix2d scatter metadata has an invalid byte count");
    }
    const int nmeta_all = metadata_total_bytes / nbytes_int;
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
    std::vector<LocalMatrixBufferView> views;
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
        LocalMatrixBufferView view;
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
        assert_valid_mpi_count(nbytes_local);
        sendcounts[ipsrc] = mpi_count(nbytes_local);
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

    const std::vector<int> sdispls = make_byte_displs(sendcounts);
    std::vector<unsigned char> sendbuf;
    if (ip == 0)
    {
        const int total_bytes = sdispls.back() + sendcounts.back();
        sendbuf.resize(total_bytes);
        for (int ipdst = 0; ipdst < nproc; ++ipdst)
        {
            const LocalMatrixBufferView& view = views[ipdst];
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
                    const std::size_t send_offset
                        = static_cast<std::size_t>(sdispls[ipdst])
                            + idx * sizeof(T);
                    unsigned char* const destination
                        = sendbuf.data() + send_offset;
                    std::memcpy(destination, &value, sizeof(T));
                }
            }
        }
    }

    const std::size_t nelem_recv
        = static_cast<std::size_t>(matloc.nrloc())
        * static_cast<std::size_t>(matloc.ncloc());
    const std::size_t nbytes_recv = nelem_recv * sizeof(T);
    std::vector<unsigned char> recvbuf(nbytes_recv);
    assert_valid_mpi_count(nbytes_recv);
    const int recv = mpi_count(nbytes_recv);
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
