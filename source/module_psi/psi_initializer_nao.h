#ifndef PSI_INITIALIZER_NAO_H
#define PSI_INITIALIZER_NAO_H
#include "module_base/cubic_spline.h"
#include "module_base/realarray.h"
#include "psi_initializer.h"

#include <memory>
/*
Psi (planewave based wavefunction) initializer: numerical atomic orbital method
*/
template <typename T>
class psi_initializer_nao : public psi_initializer<T>
{
  private:
    using Real = typename GetTypeReal<T>::type;

  public:
    psi_initializer_nao()
    {
        this->method_ = "nao";
    };
    ~psi_initializer_nao(){};

    virtual void init_psig(T* psig,  const int& ik) override;

    /// @brief initialize the psi_initializer with external data and methods
    virtual void initialize(const Structure_Factor*,             //< structure factor
                            const ModulePW::PW_Basis_K*,         //< planewave basis
                            const UnitCell*,                     //< unit cell
                            const Parallel_Kpoints*,             //< parallel kpoints
                            const int& = 1,                //< random seed
                            const pseudopot_cell_vnl* = nullptr, //< nonlocal pseudopotential
                            const int& = 0) override;      //< MPI rank

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