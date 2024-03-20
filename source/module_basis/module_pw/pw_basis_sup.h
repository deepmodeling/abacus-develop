#ifndef PWBASIS_SUP_H
#define PWBASIS_SUP_H

#include "module_base/complexmatrix.h"

namespace ModulePW
{

/**
 * @brief Special pw_basis class for sup girds, which is constrcuted in order to be consistent with the smooth grids
 * in terms of sticks. Easy for conversion between smooth and sup grids in reciprocal space.
 * @author liuyu on 2023-10-12
 * @details
 * Math:
 * plane waves: <r|g>=1/sqrt(V) * exp(igr)
 * f(r) = 1/sqrt(V) * \sum_g{c(g)*exp(igr)}
 * c(g) = \int f(r)*exp(-igr) dr
 *
 * USAGE:
 * Similar to PW_Basis, but we need to set up the smooth grids first.
 */
class PW_Basis_Sup : public PW_Basis
{

  public:
    PW_Basis_Sup()
    {
    }
    PW_Basis_Sup(std::string device_, std::string precision_) : PW_Basis(device_, precision_)
    {
        classname = "PW_Basis_Sup";
    }
    ~PW_Basis_Sup();

    // distribute plane waves and grids and set up fft according to the smooth grids
    void setuptransform(const ModulePW::PW_Basis* pw_rho);

  protected:
    // distribute plane waves to different processors according to the smooth grids
    void distribute_g(const ModulePW::PW_Basis* pw_rho);

    // method 3: ONLY for dense grids in uspp
    // consider the consistence of sticks between dense and smooth grids
    void distribution_method3(const ModulePW::PW_Basis* pw_rho);

    // Distribute sticks to cores in method 3.
    void divide_sticks_3(const int* st_length2D, // st_length2D[ixy], number of planewaves in stick on (x, y).
                         const int* st_i,        // x or x + fftnx (if x < 0) of stick.
                         const int* st_j,        // y or y + fftny (if y < 0) of stick.
                         const int* st_length,   // the stick on (x, y) consists of st_length[x*fftny+y] planewaves.
                         const int* fftixy2ip_s, // fftixy2ip of smooth grids
                         const int& nx_s,        // nx of smooth grids
                         const int& ny_s         // ny of smooth grids
    );

    void get_ig2isz_is2fftixy(int* st_bottom2D, // minimum z of stick, stored in 1d array with this->nstot elements.
                              int* st_length2D, // the stick on (x, y) consists of st_length[x*fftny+y] planewaves.
                              const ModulePW::PW_Basis* pw_rho);
}; // class PW_Basis_Sup

} // namespace ModulePW
#endif // PWBASIS_SUP_H