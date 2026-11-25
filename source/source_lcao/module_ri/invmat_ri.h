//=======================
// AUTHOR : Peize Lin
// DATE :   2022-08-17
//=======================

#ifndef ABACUS_SOURCE_LCAO_MODULE_RI_INVMAT_RI_H
#define ABACUS_SOURCE_LCAO_MODULE_RI_INVMAT_RI_H



#include <RI/global/Tensor.h>
#include <vector>
#include <array>

template <typename Tdata>
class Inverse_Matrix
{
  public:
    enum class Method
    {
        potrf,
        syev
    };
    void cal_inverse(const Method& method, const double& threshold_condition_number = 0.);

    /// Invert with a fixed-rank spectral cutoff and differentiate the retained subspace.
    std::array<RI::Tensor<Tdata>, 3> cal_inverse_with_derivatives(
        const std::array<RI::Tensor<Tdata>, 3>& derivatives,
        const double threshold_condition_number);

    void input(const RI::Tensor<Tdata>& m);
    void input(const std::vector<std::vector<RI::Tensor<Tdata>>>& ms);
    RI::Tensor<Tdata> output() const;
    std::vector<std::vector<RI::Tensor<Tdata>>> output(const std::vector<size_t>& n0,
                                                       const std::vector<size_t>& n1) const;

  private:
    void using_potrf();
    void using_syev(const double& threshold_condition_number);
    void copy_down_triangle();
    RI::Tensor<Tdata> A;
};


#endif // ABACUS_SOURCE_LCAO_MODULE_RI_INVMAT_RI_H
