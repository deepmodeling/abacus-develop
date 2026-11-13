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

int checked_mpi_count(const std::size_t count)
{
    if (count > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    {
        throw std::invalid_argument("Matrix2d collective payload exceeds MPI count range");
    }
    return static_cast<int>(count);
}

std::vector<int> make_byte_displacements(const std::vector<int>& counts)
{
    std::vector<int> displacements(counts.size(), 0);
    int displacement = 0;
    for (std::size_t ip = 0; ip < counts.size(); ++ip)
    {
        displacements[ip] = displacement;
        if (counts[ip] > std::numeric_limits<int>::max() - displacement)
        {
            throw std::invalid_argument("Matrix2d collective payload size overflow");
        }
        displacement += counts[ip];
    }
    return displacements;
}

std::vector<int> allgather_byte_counts(const int local_count, const MPI_Comm comm)
{
    const int nproc = Parallel_Common::comm_size(comm);
    std::vector<int> receive_counts(nproc, checked_mpi_count(sizeof(int)));
    const std::vector<int> receive_displacements
        = make_byte_displacements(receive_counts);
    std::vector<int> counts(nproc, 0);

    Parallel_Common::allgatherv_bytes(&local_count,
                                      checked_mpi_count(sizeof(local_count)),
                                      counts.data(),
                                      receive_counts.data(),
                                      receive_displacements.data(),
                                      comm);
    return counts;
}

template <typename T>
std::vector<int> make_local_metadata(const Matrix2d<T>& matrix)
{
    const std::shared_ptr<const Parallel_Blkcyclic> dist
        = matrix.distribution();
    std::vector<int> metadata;
    metadata.reserve(LOCAL_METADATA_SIZE + matrix.nrloc()
                     + matrix.ncloc());
    metadata.push_back(matrix.nrloc());
    metadata.push_back(matrix.ncloc());
    metadata.push_back(matrix.nrglb());
    metadata.push_back(matrix.ncglb());
    metadata.push_back(dist->irproc());
    metadata.push_back(dist->icproc());
    for (int irloc = 0; irloc < matrix.nrloc(); ++irloc)
    {
        metadata.push_back(dist->irglb(irloc));
    }
    for (int icloc = 0; icloc < matrix.ncloc(); ++icloc)
    {
        metadata.push_back(dist->icglb(icloc));
    }
    return metadata;
}

template <typename T>
std::vector<unsigned char> make_local_payload(const Matrix2d<T>& matrix)
{
    const std::vector<int> metadata = make_local_metadata(matrix);
    const std::size_t metadata_bytes = metadata.size() * sizeof(int);
    const std::size_t data_elements = static_cast<std::size_t>(matrix.nrloc())
        * static_cast<std::size_t>(matrix.ncloc());
    const std::size_t data_bytes = data_elements * sizeof(T);
    std::vector<unsigned char> payload(metadata_bytes + data_bytes);

    std::memcpy(payload.data(), metadata.data(), metadata_bytes);
    for (int irloc = 0; irloc < matrix.nrloc(); ++irloc)
    {
        for (int icloc = 0; icloc < matrix.ncloc(); ++icloc)
        {
            const std::size_t irloc_index = static_cast<std::size_t>(irloc);
            const std::size_t icloc_index = static_cast<std::size_t>(icloc);
            const std::size_t element_index
                = irloc_index * static_cast<std::size_t>(matrix.ncloc())
                + icloc_index;
            const T& value = matrix.at_local(irloc, icloc);
            std::memcpy(payload.data() + metadata_bytes + element_index * sizeof(T),
                        &value,
                        sizeof(T));
        }
    }
    return payload;
}

struct LocalPayloadView
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
LocalPayloadView view_local_payload(const std::vector<unsigned char>& payload)
{
    if (payload.size() < LOCAL_METADATA_SIZE * sizeof(int))
    {
        throw std::invalid_argument("Matrix2d collective payload is truncated");
    }

    int header[LOCAL_METADATA_SIZE] = {};
    std::memcpy(header, payload.data(), sizeof(header));
    LocalPayloadView view;
    view.nrloc = header[0];
    view.ncloc = header[1];
    view.nrglb = header[2];
    view.ncglb = header[3];
    view.iprow = header[4];
    view.ipcol = header[5];

    if (view.nrloc < 0 || view.ncloc < 0 || view.nrglb <= 0
        || view.ncglb <= 0)
    {
        throw std::invalid_argument("Matrix2d collective payload has invalid dimensions");
    }

    const std::size_t metadata_count
        = static_cast<std::size_t>(LOCAL_METADATA_SIZE) + view.nrloc + view.ncloc;
    const std::size_t metadata_bytes = metadata_count * sizeof(int);
    const std::size_t data_elements
        = static_cast<std::size_t>(view.nrloc) * static_cast<std::size_t>(view.ncloc);
    const std::size_t expected_bytes = metadata_bytes + data_elements * sizeof(T);
    if (payload.size() != expected_bytes)
    {
        throw std::invalid_argument("Matrix2d collective payload has an invalid size");
    }

    std::vector<int> metadata(metadata_count);
    std::memcpy(metadata.data(), payload.data(), metadata_bytes);
    view.irglb.assign(metadata.begin() + LOCAL_METADATA_SIZE,
                      metadata.begin() + LOCAL_METADATA_SIZE + view.nrloc);
    view.icglb.assign(metadata.begin() + LOCAL_METADATA_SIZE + view.nrloc,
                      metadata.end());
    view.data = payload.data() + metadata_bytes;
    return view;
}

template <typename T>
Matrix2d<T> rebuild_global_matrix(
    const std::vector<std::vector<unsigned char>>& payloads)
{
    std::vector<LocalPayloadView> views;
    views.reserve(payloads.size());
    int nrglb = 0;
    int ncglb = 0;
    for (std::size_t ip = 0; ip < payloads.size(); ++ip)
    {
        views.push_back(view_local_payload<T>(payloads[ip]));
        const LocalPayloadView& view = views.back();
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
        const LocalPayloadView& view = views[ip];
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
                const std::size_t element_index
                    = static_cast<std::size_t>(irloc)
                        * static_cast<std::size_t>(view.ncloc)
                    + static_cast<std::size_t>(icloc);
                T value;
                std::memcpy(&value, view.data + element_index * sizeof(T), sizeof(T));
                result.at_local(irglb, icglb) = value;
            }
        }
    }
    return result;
}

template <typename T>
std::vector<std::vector<unsigned char>> allgather_payloads(
    const Matrix2d<T>& matrix,
    const MPI_Comm comm)
{
    const std::vector<unsigned char> local_payload = make_local_payload(matrix);
    const int local_count = checked_mpi_count(local_payload.size());
    std::vector<int> payload_counts = allgather_byte_counts(local_count, comm);
    const std::vector<int> payload_displacements
        = make_byte_displacements(payload_counts);
    const int total_bytes = payload_displacements.back() + payload_counts.back();
    std::vector<unsigned char> received(total_bytes);

    Parallel_Common::allgatherv_bytes(local_payload.data(),
                                      local_count,
                                      received.data(),
                                      payload_counts.data(),
                                      payload_displacements.data(),
                                      comm);

    std::vector<std::vector<unsigned char>> payloads(payload_counts.size());
    for (std::size_t ip = 0; ip < payload_counts.size(); ++ip)
    {
        const std::size_t begin = static_cast<std::size_t>(payload_displacements[ip]);
        const std::size_t end
            = begin + static_cast<std::size_t>(payload_counts[ip]);
        payloads[ip].assign(received.begin() + begin, received.begin() + end);
    }
    return payloads;
}
} // namespace

namespace matrix2d_impl
{
template <typename T>
Matrix2d<T> allgather(const Matrix2d<T>& matrix, const MPI_Comm comm)
{
    const std::vector<std::vector<unsigned char>> payloads
        = allgather_payloads(matrix, comm);
    return rebuild_global_matrix<T>(payloads);
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

    const std::vector<unsigned char> local_payload = make_local_payload(matrix);
    const int local_count = checked_mpi_count(local_payload.size());
    std::vector<int> payload_counts = allgather_byte_counts(local_count, comm);
    const std::vector<int> payload_displacements
        = make_byte_displacements(payload_counts);
    const int total_bytes = payload_displacements.back() + payload_counts.back();
    std::vector<unsigned char> received;
    if (ip == iproot)
    {
        received.resize(total_bytes);
    }

    unsigned char* const receive_buffer = ip == iproot ? received.data() : nullptr;
    Parallel_Common::gatherv_bytes(local_payload.data(),
                                   local_count,
                                   receive_buffer,
                                   payload_counts.data(),
                                   payload_displacements.data(),
                                   iproot,
                                   comm);
    if (ip != iproot)
    {
        return std::unique_ptr<Matrix2d<T>>();
    }

    std::vector<std::vector<unsigned char>> payloads(payload_counts.size());
    for (std::size_t ipsrc = 0; ipsrc < payload_counts.size(); ++ipsrc)
    {
        const std::size_t begin
            = static_cast<std::size_t>(payload_displacements[ipsrc]);
        const std::size_t end
            = begin + static_cast<std::size_t>(payload_counts[ipsrc]);
        payloads[ipsrc].assign(received.begin() + begin, received.begin() + end);
    }
    return std::unique_ptr<Matrix2d<T>>(
        new Matrix2d<T>(rebuild_global_matrix<T>(payloads)));
}

template <typename T>
void scatter_from_host(const Matrix2d<T>& matglb, Matrix2d<T>& matloc)
{
    if (matglb.is_distributed())
    {
        throw std::invalid_argument("scatter requires a host global matrix");
    }

    const std::shared_ptr<const Parallel_Blkcyclic> dist
        = matloc.distribution();
    const MPI_Comm comm = dist->comm();
    if (comm == MPI_COMM_NULL)
    {
        throw std::invalid_argument("scatter requires a valid MPI comm");
    }

    const std::vector<int> local_meta = make_local_metadata(matloc);
    const int local_metadata_bytes
        = checked_mpi_count(local_meta.size() * sizeof(int));
    std::vector<int> metadata_counts
        = allgather_byte_counts(local_metadata_bytes, comm);
    const std::vector<int> metadata_displacements
        = make_byte_displacements(metadata_counts);
    const int metadata_total_bytes
        = metadata_displacements.back() + metadata_counts.back();
    if (metadata_total_bytes % checked_mpi_count(sizeof(int)) != 0)
    {
        throw std::invalid_argument("Matrix2d scatter metadata has an invalid byte count");
    }
    std::vector<int> all_metadata(metadata_total_bytes / sizeof(int));

    Parallel_Common::allgatherv_bytes(local_meta.data(),
                                      local_metadata_bytes,
                                      all_metadata.data(),
                                      metadata_counts.data(),
                                      metadata_displacements.data(),
                                      comm);

    const int nproc = Parallel_Common::comm_size(comm);
    const int ip = Parallel_Common::comm_rank(comm);
    std::vector<int> send_counts(nproc, 0);
    std::vector<LocalPayloadView> views;
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
        LocalPayloadView view;
        view.nrloc = all_metadata[metadata_offset];
        view.ncloc = all_metadata[metadata_offset + 1];
        view.nrglb = all_metadata[metadata_offset + 2];
        view.ncglb = all_metadata[metadata_offset + 3];
        view.iprow = all_metadata[metadata_offset + 4];
        view.ipcol = all_metadata[metadata_offset + 5];
        metadata_offset += LOCAL_METADATA_SIZE;
        const std::size_t axes_count
            = static_cast<std::size_t>(view.nrloc + view.ncloc);
        if (metadata_offset + axes_count > all_metadata.size())
        {
            throw std::invalid_argument("Matrix2d scatter axes are truncated");
        }
        view.irglb.assign(all_metadata.begin() + metadata_offset,
                          all_metadata.begin() + metadata_offset + view.nrloc);
        view.icglb.assign(all_metadata.begin() + metadata_offset + view.nrloc,
                          all_metadata.begin() + metadata_offset + view.nrloc
                              + view.ncloc);
        metadata_offset += axes_count;

        if (view.nrloc < 0 || view.ncloc < 0
            || view.nrglb != nrglb || view.ncglb != ncglb)
        {
            throw std::invalid_argument("Matrix2d ranks disagree on scatter global shape");
        }
        const std::size_t local_elements
            = static_cast<std::size_t>(view.nrloc)
            * static_cast<std::size_t>(view.ncloc);
        send_counts[ipsrc] = checked_mpi_count(local_elements * sizeof(T));
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

    const std::vector<int> send_displacements = make_byte_displacements(send_counts);
    std::vector<unsigned char> send_buffer;
    if (ip == 0)
    {
        const int total_bytes = send_displacements.back() + send_counts.back();
        send_buffer.resize(total_bytes);
        for (int ipdst = 0; ipdst < nproc; ++ipdst)
        {
            const LocalPayloadView& view = views[ipdst];
            for (int irloc = 0; irloc < view.nrloc; ++irloc)
            {
                const int irglb = view.irglb[irloc];
                for (int icloc = 0; icloc < view.ncloc; ++icloc)
                {
                    const int icglb = view.icglb[icloc];
                    const std::size_t element_index
                        = static_cast<std::size_t>(irloc)
                            * static_cast<std::size_t>(view.ncloc)
                        + static_cast<std::size_t>(icloc);
                    const T& value = matglb.at_local(irglb, icglb);
                    std::memcpy(send_buffer.data()
                                    + send_displacements[ipdst]
                                    + element_index * sizeof(T),
                                &value,
                                sizeof(T));
                }
            }
        }
    }

    const std::size_t receive_elements
        = static_cast<std::size_t>(matloc.nrloc())
        * static_cast<std::size_t>(matloc.ncloc());
    std::vector<unsigned char> receive_buffer(receive_elements * sizeof(T));
    const unsigned char* receive_source = ip == 0 ? send_buffer.data() : nullptr;
    Parallel_Common::scatterv_bytes(receive_source,
                                    send_counts.data(),
                                    send_displacements.data(),
                                    receive_buffer.data(),
                                    checked_mpi_count(receive_buffer.size()),
                                    0,
                                    comm);

    for (int irloc = 0; irloc < matloc.nrloc(); ++irloc)
    {
        for (int icloc = 0; icloc < matloc.ncloc(); ++icloc)
        {
            const std::size_t element_index
                = static_cast<std::size_t>(irloc)
                    * static_cast<std::size_t>(matloc.ncloc())
                + static_cast<std::size_t>(icloc);
            T value;
            std::memcpy(&value,
                        receive_buffer.data() + element_index * sizeof(T),
                        sizeof(T));
            matloc.at_local(irloc, icloc) = value;
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
