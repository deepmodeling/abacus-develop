#ifndef HYBRID_IMAGE_DOMAIN_H
#define HYBRID_IMAGE_DOMAIN_H

#include <RI/global/Tensor.h>
#include <array>
#include <complex>
#include <map>
#include <vector>

namespace RI_2D_Comm
{
using ImageCell = std::array<int, 3>;

ImageCell centered_image(const ImageCell& cell, const ImageCell& period);

// Use the existing LibRI max-norm threshold, with no additional approximation.
template <typename T>
ImageCell retained_image_extent(std::map<int, std::map<std::pair<int, ImageCell>, RI::Tensor<T>>>& blocks,
                                const double threshold);

ImageCell hybrid_image_period(const ImageCell& c_extent, const ImageCell& v_extent);

template <typename T>
void unfold_hybrid_images(std::vector<std::map<int, std::map<std::pair<int, ImageCell>, RI::Tensor<T>>>>& density,
                          const ImageCell& physical_period, const ImageCell& image_period,
                          const std::array<double, 3>& mesh_origin,
                          const std::array<double, 3>& input_phase_per_cell);

template <typename T>
void decode_hybrid_images(std::map<int, std::map<std::pair<int, ImageCell>, RI::Tensor<T>>>& blocks,
                          const ImageCell& period);
} // namespace RI_2D_Comm
#endif
