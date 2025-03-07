#include "./grid.h"

void Grid::initGrid(const int fftdim,
                    const int ndata,
                    const std::string *cell,
                    const double *a,
                    const torch::Device device,
                    double *volume)
{
    this->initGrid_(fftdim, ndata, cell, a, device, volume, this->fft_grid, this->fft_gg);
    std::cout << "Init grid done" << std::endl;
}

void Grid::initGrid_(const int fftdim,
                     const int ndata,
                     const std::string *cell,
                     const double *a,
                     const torch::Device device,
                     double *volume,
                     std::vector<std::vector<torch::Tensor>> &grid,
                     std::vector<torch::Tensor> &gg)
{
    this->fft_grid = std::vector<std::vector<torch::Tensor>>(ndata);
    this->fft_gg = std::vector<torch::Tensor>(ndata);
    for (int i = 0; i < ndata; ++i)
    {
        this->fft_grid[i] = std::vector<torch::Tensor>(3);
    }

    for (int it = 0; it < ndata; ++it)
    {
        if (cell[it] == "sc"){
            this->initScRecipGrid(fftdim, a[it], it, device, volume, grid, gg);
        }
        else if (cell[it] == "fcc"){
            this->initFccRecipGrid(fftdim, a[it], it, device, volume, grid, gg);
        }
        else if (cell[it] == "bcc"){
            this->initBccRecipGrid(fftdim, a[it], it, device, volume, grid, gg);
        }
    }
}

void Grid::initScRecipGrid(const int fftdim,
                           const double a,
                           const int index,
                           const torch::Device device,
                           double *volume,
                           std::vector<std::vector<torch::Tensor>> &grid,
                           std::vector<torch::Tensor> &gg)
{
    volume[index] = std::pow(a, 3);
    torch::Tensor fre = torch::fft::fftfreq(fftdim, a / fftdim).to(device) * 2. * M_PI;
    grid[index] = torch::meshgrid({fre, fre, fre});
    gg[index] = grid[index][0] * grid[index][0] + grid[index][1] * grid[index][1] + grid[index][2] * grid[index][2];
}

void Grid::initFccRecipGrid(const int fftdim,
                            const double a,
                            const int index,
                            const torch::Device device,
                            double *volume,
                            std::vector<std::vector<torch::Tensor>> &grid,
                            std::vector<torch::Tensor> &gg)
{
    // std::cout << "init grid" << std::endl;
    volume[index] = std::pow(a, 3) / 4.;
    double coef = 1. / sqrt(2.);
    // std::cout << "fftfreq" << std::endl;
    torch::Tensor fre = torch::fft::fftfreq(fftdim, a * coef / fftdim).to(device) * 2. * M_PI;
    auto originalGrid = torch::meshgrid({fre, fre, fre});
    grid[index][0] = coef * (-originalGrid[0] + originalGrid[1] + originalGrid[2]);
    grid[index][1] = coef * (originalGrid[0] - originalGrid[1] + originalGrid[2]);
    grid[index][2] = coef * (originalGrid[0] + originalGrid[1] - originalGrid[2]);
    // std::cout << "gg" << std::endl;
    gg[index] = grid[index][0] * grid[index][0] + grid[index][1] * grid[index][1] + grid[index][2] * grid[index][2];
}

void Grid::initBccRecipGrid(const int fftdim,
                            const double a,
                            const int index,
                            const torch::Device device,
                            double *volume,
                            std::vector<std::vector<torch::Tensor>> &grid,
                            std::vector<torch::Tensor> &gg)
{
    volume[index] = std::pow(a, 3) / 2.;
    double coef = sqrt(3.) / 2.;
    torch::Tensor fre = torch::fft::fftfreq(fftdim, a * coef / fftdim).to(device) * 2. * M_PI;
    auto originalGrid = torch::meshgrid({fre, fre, fre});
    grid[index][0] = coef * (originalGrid[1] + originalGrid[2]);
    grid[index][1] = coef * (originalGrid[0] + originalGrid[2]);
    grid[index][2] = coef * (originalGrid[0] + originalGrid[1]);
    gg[index] = grid[index][0] * grid[index][0] + grid[index][1] * grid[index][1] + grid[index][2] * grid[index][2];
}