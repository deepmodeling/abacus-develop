#include "gpu_neighbor_list.h"

#include "source_base/timer.h"

#include <cstddef>

int gpu_build_neighbor_list(int nall,
                                double cutoff,
                                const double* position,
                                int* neighbor_count,
                                int* neighbor_indices,
                                int* max_neighbors,
                                char* message,
                                int message_size);

int gpu_filter_neighbor_list(int nall,
                                 double cutoff,
                                 const double* position,
                                 const int* candidate_count,
                                 const int* candidate_indices,
                                 int candidate_max_neighbors,
                                 int* neighbor_count,
                                 int* neighbor_indices,
                                 int* max_neighbors,
                                 char* message,
                                 int message_size);

int gpu_build_neighbor_list_device(int nall,
                                       double cutoff,
                                       const double* position,
                                       void** device_neighbor_count,
                                       void** device_neighbor_indices,
                                       int* max_neighbors,
                                       char* message,
                                       int message_size);

int gpu_filter_neighbor_list_device(int nall,
                                        double cutoff,
                                        const double* position,
                                        int candidate_max_neighbors,
                                        const void* device_candidate_count,
                                        const void* device_candidate_indices,
                                        void** device_neighbor_count,
                                        void** device_neighbor_indices,
                                        int* max_neighbors,
                                        char* message,
                                        int message_size);

int gpu_upload_neighbor_list(int nall,
                                 int max_neighbors,
                                 const int* neighbor_count,
                                 const int* neighbor_indices,
                                 void** device_neighbor_count,
                                 void** device_neighbor_indices,
                                 char* message,
                                 int message_size);

void gpu_release_neighbor_list(void* device_neighbor_count,
                                   void* device_neighbor_indices);

namespace ModuleNeighList
{

GPU_NeighborList::~GPU_NeighborList()
{
    const bool aliases_candidate = device_neighbor_count_ == candidate_device_count_
                                   && device_neighbor_indices_ == candidate_device_indices_;
    const bool aliases_filtered = device_neighbor_count_ == filtered_device_count_
                                  && device_neighbor_indices_ == filtered_device_indices_;
    if (!aliases_candidate && !aliases_filtered)
    {
        gpu_release_neighbor_list(device_neighbor_count_, device_neighbor_indices_);
    }
    gpu_release_neighbor_list(candidate_device_count_, candidate_device_indices_);
    gpu_release_neighbor_list(filtered_device_count_, filtered_device_indices_);
}

bool GPU_NeighborList::build_device(int nall,
                                        double cutoff,
                                        const std::vector<double>& position,
                                        int& max_neighbors,
                                        std::string& error)
{
    ModuleBase::timer::start("GPU_NeighborList", "build_device");
    if (nall <= 0 || position.size() != static_cast<std::size_t>(3 * nall))
    {
        error = "invalid GPU device neighbor-list input";
        ModuleBase::timer::end("GPU_NeighborList", "build_device");
        return false;
    }

    char message[512] = {0};
    void* new_count = nullptr;
    void* new_indices = nullptr;
    const int status = gpu_build_neighbor_list_device(nall,
                                                          cutoff,
                                                          position.data(),
                                                          &new_count,
                                                          &new_indices,
                                                          &max_neighbors,
                                                          message,
                                                          sizeof(message));
    if (status != 0)
    {
        error = message;
        ModuleBase::timer::end("GPU_NeighborList", "build_device");
        return false;
    }
    gpu_release_neighbor_list(candidate_device_count_, candidate_device_indices_);
    gpu_release_neighbor_list(filtered_device_count_, filtered_device_indices_);
    candidate_device_count_ = new_count;
    candidate_device_indices_ = new_indices;
    candidate_device_nall_ = nall;
    candidate_device_max_neighbors_ = max_neighbors;
    filtered_device_count_ = nullptr;
    filtered_device_indices_ = nullptr;
    filtered_device_nall_ = 0;
    filtered_device_max_neighbors_ = 0;
    device_neighbor_count_ = candidate_device_count_;
    device_neighbor_indices_ = candidate_device_indices_;
    device_nall_ = nall;
    device_max_neighbors_ = max_neighbors;
    ModuleBase::timer::end("GPU_NeighborList", "build_device");
    return true;
}

bool GPU_NeighborList::filter_device(int nall,
                                         double cutoff,
                                         const std::vector<double>& position,
                                         int candidate_max_neighbors,
                                         int& max_neighbors,
                                         std::string& error)
{
    ModuleBase::timer::start("GPU_NeighborList", "filter_device");
    if (nall <= 0 || position.size() != static_cast<std::size_t>(3 * nall)
        || candidate_device_nall_ != nall
        || candidate_device_count_ == nullptr
        || (candidate_max_neighbors > 0 && candidate_device_indices_ == nullptr))
    {
        error = "invalid GPU device neighbor-list filter input";
        ModuleBase::timer::end("GPU_NeighborList", "filter_device");
        return false;
    }

    char message[512] = {0};
    void* new_count = nullptr;
    void* new_indices = nullptr;
    const int status = gpu_filter_neighbor_list_device(
        nall,
        cutoff,
        position.data(),
        candidate_max_neighbors,
        candidate_device_count_,
        candidate_device_indices_,
        &new_count,
        &new_indices,
        &max_neighbors,
        message,
        sizeof(message));
    if (status != 0)
    {
        error = message;
        ModuleBase::timer::end("GPU_NeighborList", "filter_device");
        return false;
    }
    gpu_release_neighbor_list(filtered_device_count_, filtered_device_indices_);
    filtered_device_count_ = new_count;
    filtered_device_indices_ = new_indices;
    filtered_device_nall_ = nall;
    filtered_device_max_neighbors_ = max_neighbors;
    device_neighbor_count_ = filtered_device_count_;
    device_neighbor_indices_ = filtered_device_indices_;
    device_nall_ = nall;
    device_max_neighbors_ = max_neighbors;
    ModuleBase::timer::end("GPU_NeighborList", "filter_device");
    return true;
}

void GPU_NeighborList::use_candidate_device()
{
    gpu_release_neighbor_list(filtered_device_count_, filtered_device_indices_);
    filtered_device_count_ = nullptr;
    filtered_device_indices_ = nullptr;
    filtered_device_nall_ = 0;
    filtered_device_max_neighbors_ = 0;
    device_neighbor_count_ = candidate_device_count_;
    device_neighbor_indices_ = candidate_device_indices_;
    device_nall_ = candidate_device_nall_;
    device_max_neighbors_ = candidate_device_max_neighbors_;
}

bool GPU_NeighborList::upload(int nall,
                                  int max_neighbors,
                                  const std::vector<int>& neighbor_count,
                                  const std::vector<int>& neighbor_indices,
                                  std::string& error)
{
    ModuleBase::timer::start("GPU_NeighborList", "upload");
    if (nall <= 0 || max_neighbors < 0
        || neighbor_count.size() != static_cast<std::size_t>(nall)
        || neighbor_indices.size() != static_cast<std::size_t>(nall)
                                             * static_cast<std::size_t>(max_neighbors))
    {
        error = "invalid GPU neighbor-list upload input";
        ModuleBase::timer::end("GPU_NeighborList", "upload");
        return false;
    }

    char message[512] = {0};
    void* new_device_count = nullptr;
    void* new_device_indices = nullptr;
    const int status = gpu_upload_neighbor_list(nall,
                                                    max_neighbors,
                                                    neighbor_count.data(),
                                                    neighbor_indices.empty() ? NULL : neighbor_indices.data(),
                                                    &new_device_count,
                                                    &new_device_indices,
                                                    message,
                                                    sizeof(message));
    if (status != 0)
    {
        error = message;
        ModuleBase::timer::end("GPU_NeighborList", "upload");
        return false;
    }
    gpu_release_neighbor_list(device_neighbor_count_, device_neighbor_indices_);
    device_neighbor_count_ = new_device_count;
    device_neighbor_indices_ = new_device_indices;
    device_nall_ = nall;
    device_max_neighbors_ = max_neighbors;
    ModuleBase::timer::end("GPU_NeighborList", "upload");
    return true;
}

const void* GPU_NeighborList::device_neighbor_count() const
{
    return device_neighbor_count_;
}

const void* GPU_NeighborList::device_neighbor_indices() const
{
    return device_neighbor_indices_;
}

int GPU_NeighborList::device_nall() const
{
    return device_nall_;
}

int GPU_NeighborList::device_max_neighbors() const
{
    return device_max_neighbors_;
}

bool GPU_NeighborList::build(int nall,
                                 double cutoff,
                                 const std::vector<double>& position,
                                 std::vector<int>& neighbor_count,
                                 std::vector<int>& neighbor_indices,
                                 int& max_neighbors,
                                 std::string& error) const
{
    ModuleBase::timer::start("GPU_NeighborList", "build");
    if (nall <= 0 || cutoff <= 0.0
        || position.size() != static_cast<std::size_t>(3 * nall))
    {
        error = "invalid GPU neighbor-list input";
        ModuleBase::timer::end("GPU_NeighborList", "build");
        return false;
    }

    neighbor_count.assign(static_cast<std::size_t>(nall), 0);
    neighbor_indices.clear();
    max_neighbors = 0;
    char message[512] = {0};
    const int status = gpu_build_neighbor_list(nall,
                                                   cutoff,
                                                   position.data(),
                                                   neighbor_count.data(),
                                                   NULL,
                                                   &max_neighbors,
                                                   message,
                                                   sizeof(message));
    if (status != 0)
    {
        error = message;
        ModuleBase::timer::end("GPU_NeighborList", "build");
        return false;
    }

    neighbor_indices.assign(static_cast<std::size_t>(nall)
                                * static_cast<std::size_t>(max_neighbors),
                            0);
    if (max_neighbors > 0)
    {
        const int fill_status = gpu_build_neighbor_list(nall,
                                                            cutoff,
                                                            position.data(),
                                                            neighbor_count.data(),
                                                            neighbor_indices.data(),
                                                            &max_neighbors,
                                                            message,
                                                            sizeof(message));
        if (fill_status != 0)
        {
            error = message;
            ModuleBase::timer::end("GPU_NeighborList", "build");
            return false;
        }
    }
    ModuleBase::timer::end("GPU_NeighborList", "build");
    return true;
}

bool GPU_NeighborList::filter(int nall,
                                  double cutoff,
                                  const std::vector<double>& position,
                                  const std::vector<int>& candidate_count,
                                  const std::vector<int>& candidate_indices,
                                  int candidate_max_neighbors,
                                  std::vector<int>& neighbor_count,
                                  std::vector<int>& neighbor_indices,
                                  int& max_neighbors,
                                  std::string& error) const
{
    ModuleBase::timer::start("GPU_NeighborList", "filter");
    if (nall <= 0 || cutoff <= 0.0
        || position.size() != static_cast<std::size_t>(3 * nall)
        || candidate_count.size() != static_cast<std::size_t>(nall)
        || candidate_max_neighbors < 0
        || candidate_indices.size() != static_cast<std::size_t>(nall)
                                         * static_cast<std::size_t>(candidate_max_neighbors))
    {
        error = "invalid GPU neighbor-list filter input";
        ModuleBase::timer::end("GPU_NeighborList", "filter");
        return false;
    }

    neighbor_count.assign(static_cast<std::size_t>(nall), 0);
    neighbor_indices.assign(candidate_indices.size(), 0);
    max_neighbors = 0;
    char message[512] = {0};
    const int status = gpu_filter_neighbor_list(nall,
                                                    cutoff,
                                                    position.data(),
                                                    candidate_count.data(),
                                                    candidate_indices.data(),
                                                    candidate_max_neighbors,
                                                    neighbor_count.data(),
                                                    neighbor_indices.empty() ? NULL : neighbor_indices.data(),
                                                    &max_neighbors,
                                                    message,
                                                    sizeof(message));
    if (status != 0)
    {
        error = message;
        ModuleBase::timer::end("GPU_NeighborList", "filter");
        return false;
    }
    neighbor_indices.resize(static_cast<std::size_t>(nall)
                            * static_cast<std::size_t>(max_neighbors));
    ModuleBase::timer::end("GPU_NeighborList", "filter");
    return true;
}

} // namespace ModuleNeighList
