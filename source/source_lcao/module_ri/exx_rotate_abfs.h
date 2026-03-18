
#ifndef EXX_ROTATE_ABFS_H
#define EXX_ROTATE_ABFS_H

#include "LRI_CV.h"
// #include "module_xc/exx_info.h"
// #include "module_basis/module_ao/ORB_atomic_lm.h"
#include "source_base/matrix.h"
#include "Exx_LRI.h"
// #include "module_ri/Exx_LRI.h"
// #include <RI/physics/Exx.h>
#include <RI/ri/RI_Tools.h>
#include <array>
#include <map>
#include <mpi.h>
#include <vector>

template <typename Tdata>
class Moment_abfs
{
  private:
    using TA = int;
    using Tcell = int;
    static constexpr std::size_t Ndim = 3;
    using TC = std::array<Tcell, Ndim>;
    using Tq = std::array<double, Ndim>;
    using TAC = std::pair<TA, TC>;
    using TAq = std::pair<TA, Tq>;

  public:
    Moment_abfs(Exx_Info::Exx_Info_RI& info_in) : info(info_in) {};
    ~Moment_abfs() {};
    void cal_VR(
        const UnitCell& ucell,
        const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& orb_in,
        const std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA, std::array<Tcell, Ndim>>>>>& list_r,
        const std::vector<double>& orb_cutoff,
        const double Rc,
        LRI_CV<Tdata>& cv,
        std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Vs_cut);
    void discard0_VR(
        const UnitCell& ucell,
        const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& orb_in,
        const std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA, std::array<Tcell, Ndim>>>>>& list_r,
        const std::vector<double>& orb_cutoff,
        const double Rc,
        LRI_CV<Tdata>& cv,
        std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Vs_cut);
    void cal_multipole(const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& orb_in);
    void rotate_abfs(std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& orb_in);
    double sum_triple_Y_YLM_real(int l1,
                                 int m1, // real m1, not index
                                 int l2,
                                 int m2,                         // real m2, not index
                                 const std::vector<double>& rly, // real Y_LM(R)
                                 const ORB_gaunt_table& MGT,
                                 const double distance);
    double cal_cl1l2(int l1, int l2) const;
    /// double factorial
    double dfact(const int& l) const;
    int factorial(const int& n) const;
    double ln_factorial(int n) const;

    void out_pure_ri_tensor(const std::string fn, RI::Tensor<std::complex<double>>& olp, const double threshold);
    void out_pure_ri_tensor(const std::string fn, RI::Tensor<double>& olp, const double threshold);

    std::vector<std::vector<std::vector<double>>> multipole;

  private:
    // std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> VR;
    Exx_Info::Exx_Info_RI& info;
};
#include "exx_rotate_abfs.hpp"

#endif
