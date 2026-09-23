#include "nep_gpu_backend.h"

#include "source_cell/module_neighlist/gpu_neighbor_list.h"
#include "nep.h"
#include "nep_gpu_kernels.h"

#include <algorithm>
#include <cmath>
#include <memory>

namespace ModuleESolver
{

class NEP_GPU_Backend::Impl
{
  public:
    Impl() : gpu_model_(NULL) {}

    ~Impl()
    {
        nep_gpu_destroy(gpu_model_);
    }

    NEP_GPU_Model* gpu_model_;
    GPU_NeighborList neighbor_list_;
    bool candidate_valid_ = false;
    int candidate_nlocal_ = 0;
    int candidate_nall_ = 0;
    double candidate_cutoff_ = 0.0;
    double candidate_skin_ = 0.0;
    std::vector<double> candidate_reference_position_;
    std::vector<int> candidate_count_;
    std::vector<int> candidate_indices_;
    int candidate_max_neighbors_ = 0;
};

NEP_GPU_Backend::NEP_GPU_Backend() : impl_(new Impl()) {}

NEP_GPU_Backend::~NEP_GPU_Backend()
{
    delete impl_;
}

bool NEP_GPU_Backend::initialize(const NEP& nep, int local_rank, std::string& error)
{
    if (nep.paramb.model_type != 0 || nep.zbl.enabled)
    {
        error = "the NEP GPU backend supports potential-only, non-ZBL NEP models";
        return false;
    }

    char message[512] = {0};
    int device_id = -1;
    if (nep_gpu_select_device(local_rank, &device_id, message, sizeof(message)) != 0)
    {
        error = message;
        return false;
    }

    impl_->gpu_model_ = nep_gpu_create(&nep, message, sizeof(message));
    if (impl_->gpu_model_ == NULL)
    {
        error = message;
        return false;
    }
    return true;
}

bool NEP_GPU_Backend::compute(int nlocal,
                              int nall,
                              double cutoff,
                              double skin,
                              const std::vector<int>& type,
                              const std::vector<double>& position,
                              std::vector<double>& energy,
                              std::vector<double>& force,
                              std::vector<double>& virial,
                              std::string& error)
{
    const double effective_skin = std::max(0.0, skin);
    bool rebuild_candidate = !impl_->candidate_valid_
                             || impl_->candidate_nlocal_ != nlocal
                             || impl_->candidate_nall_ != nall
                             || impl_->candidate_cutoff_ != cutoff
                             || impl_->candidate_skin_ != effective_skin
                             || effective_skin <= 0.0;
    if (!rebuild_candidate)
    {
        double max_displacement2 = 0.0;
        for (int dim = 0; dim < 3; ++dim)
        {
            for (int i = 0; i < nlocal; ++i)
            {
                const std::size_t index = static_cast<std::size_t>(dim * nall + i);
                const double displacement = position[index]
                                            - impl_->candidate_reference_position_[index];
                max_displacement2 = std::max(max_displacement2, displacement * displacement);
            }
        }
        rebuild_candidate = max_displacement2 >= 0.25 * effective_skin * effective_skin;
    }

    if (rebuild_candidate)
    {
        const double candidate_cutoff = cutoff + effective_skin;
        if (!impl_->neighbor_list_.build_device(nall,
                                                candidate_cutoff,
                                                position,
                                                impl_->candidate_max_neighbors_,
                                                error))
        {
            return false;
        }
        impl_->candidate_valid_ = true;
        impl_->candidate_nlocal_ = nlocal;
        impl_->candidate_nall_ = nall;
        impl_->candidate_cutoff_ = cutoff;
        impl_->candidate_skin_ = effective_skin;
        impl_->candidate_reference_position_ = position;
    }

    int max_count = 0;
    if (effective_skin > 0.0
        && (rebuild_candidate || impl_->candidate_max_neighbors_ > 0))
    {
        if (!impl_->neighbor_list_.filter_device(nall,
                                                 cutoff,
                                                 position,
                                                 impl_->candidate_max_neighbors_,
                                                 max_count,
                                                 error))
        {
            return false;
        }
    }
    else
    {
        impl_->neighbor_list_.use_candidate_device();
        max_count = impl_->candidate_max_neighbors_;
    }

    std::vector<double> x(static_cast<std::size_t>(nall), 0.0);
    std::vector<double> y(static_cast<std::size_t>(nall), 0.0);
    std::vector<double> z(static_cast<std::size_t>(nall), 0.0);
    for (int i = 0; i < nall; ++i)
    {
        x[static_cast<std::size_t>(i)] = position[static_cast<std::size_t>(i)];
        y[static_cast<std::size_t>(i)] = position[static_cast<std::size_t>(nall + i)];
        z[static_cast<std::size_t>(i)] = position[static_cast<std::size_t>(2 * nall + i)];
    }

    energy.assign(static_cast<std::size_t>(nlocal), 0.0);
    force.assign(static_cast<std::size_t>(3 * nlocal), 0.0);
    virial.assign(static_cast<std::size_t>(9 * nlocal), 0.0);
    NEP_GPU_Result result;
    result.potential = energy.data();
    result.force = force.data();
    result.virial = virial.data();
    char message[512] = {0};
    const int status = nep_gpu_compute_device_neighbors(impl_->gpu_model_,
                                                        nlocal,
                                                        nall,
                                                        max_count,
                                                        max_count,
                                                        type.data(),
                                                        x.data(),
                                                        y.data(),
                                                        z.data(),
                                                        impl_->neighbor_list_.device_neighbor_count(),
                                                        impl_->neighbor_list_.device_neighbor_indices(),
                                                        impl_->neighbor_list_.device_neighbor_count(),
                                                        impl_->neighbor_list_.device_neighbor_indices(),
                                                        result,
                                                        message,
                                                        sizeof(message));
    if (status != 0)
    {
        error = message;
        return false;
    }
    return true;
}

} // namespace ModuleESolver
