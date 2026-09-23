#ifndef ESOLVER_NEP_GPU_BACKEND_H
#define ESOLVER_NEP_GPU_BACKEND_H

#include <string>
#include <vector>

class NEP;

namespace ModuleESolver
{

class NEP_GPU_Backend
{
  public:
    NEP_GPU_Backend();
    ~NEP_GPU_Backend();

    bool initialize(const NEP& nep, int local_rank, std::string& error);
    bool compute(int nlocal,
                 int nall,
                 double cutoff,
                 double skin,
                 const std::vector<int>& type,
                 const std::vector<double>& position,
                 std::vector<double>& energy,
                 std::vector<double>& force,
                 std::vector<double>& virial,
                 std::string& error);

  private:
    class Impl;
    Impl* impl_;
};

} // namespace ModuleESolver

#endif
