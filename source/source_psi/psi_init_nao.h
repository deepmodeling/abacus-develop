#ifndef PSI_INIT_NAO_H
#define PSI_INIT_NAO_H
#include "source_base/cubic_spline.h"
#include "source_base/realarray.h"
#include "source_base/spherical_bessel_transformer.h"
#include "psi_base.h"

#include <memory>
/*
Psi (planewave based wavefunction) initializer: numerical atomic orbital method
*/
template <typename T>
class psi_init_nao : public psi_base<T>
{
  private:
    using Real = typename GetTypeReal<T>::type;

  public:
    psi_init_nao()
    {
        this->method_ = "nao";
    };
    ~psi_init_nao(){};

    virtual void init_psig(T* psig, const int& ik) override;

    /// @brief initialize the psi_base with external data and methods
    virtual void initialize(const Structure_Factor* sf,             //< structure factor
                            const ModulePW::PW_Basis_K* pw_wfc,         //< planewave basis
                            const UnitCell* p_ucell,                     //< unit cell
                            const std::vector<int>& ik2iktot,             //< ik2iktot: local->global k-point mapping
                            const int& nkstot,                          //< nkstot: total number of k-points
                            const int& random_seed,                      //< random seed
                            const int& lmaxkb,                          //< lmaxkb: max angular momentum for non-local projectors
                            const int& rank,                            //< MPI rank
                            const int& npol,                            //< npol
                            const int& nbands) override;                //< nbands

    void read_external_orbs(const std::string* orbital_files, const int& rank);

    virtual void tabulate() override;

    std::vector<std::string> external_orbs() const
    {
        return orbital_files_;
    }

    std::vector<std::vector<int>> nr() const
    {
        return nr_;
    }

    std::vector<int> nr(const int& itype) const
    {
        return nr_[itype];
    }

    int nr(const int& itype, const int& ichi) const
    {
        return nr_[itype][ichi];
    }

    std::vector<std::vector<std::vector<double>>> chi() const
    {
        return chi_;
    }

    std::vector<std::vector<double>> chi(const int& itype) const
    {
        return chi_[itype];
    }

    std::vector<double> chi(const int& itype, const int& ichi) const
    {
        return chi_[itype][ichi];
    }

    double chi(const int& itype, const int& ichi, const int& ir) const
    {
        return chi_[itype][ichi][ir];
    }

    std::vector<std::vector<std::vector<double>>> rgrid() const
    {
        return rgrid_;
    }

    std::vector<std::vector<double>> rgrid(const int& itype) const
    {
        return rgrid_[itype];
    }

    std::vector<double> rgrid(const int& itype, const int& ichi) const
    {
        return rgrid_[itype][ichi];
    }

    double rgrid(const int& itype, const int& ichi, const int& ir) const
    {
        return rgrid_[itype][ichi][ir];
    }

  protected:
    /// @brief allocate memory for overlap table
    void allocate_ao_table();
    std::vector<std::string> orbital_files_;
    /// @brief cubic spline for interpolation
    std::unique_ptr<ModuleBase::CubicSpline> cubspl_;
    /// @brief radial map, [itype][l][izeta] -> i
    ModuleBase::realArray projmap_;
    /// @brief number of realspace grids per type per chi, [itype][ichi]
    std::vector<std::vector<int>> nr_;
    /// @brief data of numerical atomic orbital per type per chi per position, [itype][ichi][ir]
    std::vector<std::vector<std::vector<double>>> chi_;
    /// @brief r of numerical atomic orbital per type per chi per position, [itype][ichi][ir]
    std::vector<std::vector<std::vector<double>>> rgrid_;
    /// @brief useful for atomic-like methods
    ModuleBase::SphericalBesselTransformer sbt;
};
#endif
