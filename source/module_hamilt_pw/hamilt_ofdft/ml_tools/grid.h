#ifndef GRID_H
#define GRID_H

#include <torch/torch.h>

class Grid
{
  public:
    void initGrid(const int fftdim,
                  const int ndata,
                  const std::string *cell,
                  const double *a,
                  const torch::Device device,
                  double *volume);

    // fft grid
    std::vector<std::vector<torch::Tensor>> fft_grid; // ntrain*3*fftdim*fftdim*fftdim
    std::vector<torch::Tensor> fft_gg;

  private:
    void initGrid_(const int fftdim,
                   const int ndata,
                   const std::string *cell,
                   const double *a,
                   const torch::Device device,
                   double *volume,
                   std::vector<std::vector<torch::Tensor>> &grid,
                   std::vector<torch::Tensor> &gg);
    void initScRecipGrid(const int fftdim,
                         const double a,
                         const int index,
                         const torch::Device device,
                         double *volume,
                         std::vector<std::vector<torch::Tensor>> &grid,
                         std::vector<torch::Tensor> &gg);
    void initFccRecipGrid(const int fftdim,
                          const double a,
                          const int index,
                          const torch::Device device,
                          double *volume,
                          std::vector<std::vector<torch::Tensor>> &grid,
                          std::vector<torch::Tensor> &gg);
    void initBccRecipGrid(const int fftdim,
                          const double a,
                          const int index,
                          const torch::Device device,
                          double *volume,
                          std::vector<std::vector<torch::Tensor>> &grid,
                          std::vector<torch::Tensor> &gg);
};
#endif