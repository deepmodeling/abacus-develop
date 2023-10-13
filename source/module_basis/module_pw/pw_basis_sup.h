#ifndef PWBASIS_SUP_H
#define PWBASIS_SUP_H

#include "module_base/complexmatrix.h"
#include "pw_basis.h"

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

    // distribute plane waves and grids and set up fft
    void setuptransform(const int* fftixy2ip_s, // fftixy2ip of smooth grids
                        const int& nx_s,        // nx of smooth grids
                        const int& ny_s         // ny of smooth grids
    );

    // get igs2igd and igd2igs
    void link_igs_igd(const ModuleBase::Vector3<double>* gcar_s, // G vectors  of smooth grids
                      const int& npw_s                           // npw of smooth grids
    );

    // recip of smooth grids to recip of dense grids
    void recip_gs2gd(const ModuleBase::ComplexMatrix& data_s, // data of smooth grids
                     ModuleBase::ComplexMatrix& data_d        // data of dense grids
    ) const;

    // recip of dense grids to recip of smooth grids
    void recip_gd2gs(const ModuleBase::ComplexMatrix& data_d, // data of dense grids
                     ModuleBase::ComplexMatrix& data_s        // data of smooth grids
    ) const;

  protected:
    // distribute plane waves to different processors
    void distribute_g(const int* fftixy2ip_s, // fftixy2ip of smooth grids
                      const int& nx_s,        // nx of smooth grids
                      const int& ny_s         // ny of smooth grids
    );

    // method 3: ONLY for dense grids in uspp
    // consider the consistence of sticks between dense and smooth grids
    void distribution_method3(const int* fftixy2ip_s, // fftixy2ip of smooth grids
                              const int& nx_s,        // nx of smooth grids
                              const int& ny_s         // ny of smooth grids
    );

    // Distribute sticks to cores in method 3.
    void divide_sticks_3(const int* st_length2D, // st_length2D[ixy], number of planewaves in stick on (x, y).
                         const int* st_i,        // x or x + fftnx (if x < 0) of stick.
                         const int* st_j,        // y or y + fftny (if y < 0) of stick.
                         const int* st_length,   // the stick on (x, y) consists of st_length[x*fftny+y] planewaves.
                         const int* fftixy2ip_s, // fftixy2ip of smooth grids
                         const int& nx_s,        // nx of smooth grids
                         const int& ny_s         // ny of smooth grids
    );

    int* igs2igd = nullptr; // ig of smooth grids to ig of dense grids
    int* igd2igs = nullptr; // ig of dense grids to ig of smooth grids

}; // class PW_Basis_Sup

} // namespace ModulePW
#endif // PWBASIS_SUP_H