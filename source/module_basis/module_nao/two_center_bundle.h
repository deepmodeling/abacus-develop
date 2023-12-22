#ifndef TWO_CENTER_BUNDLE_H
#define TWO_CENTER_BUNDLE_H

#include "module_basis/module_nao/two_center_integrator.h"
#include "module_basis/module_ao/ORB_read.h"

#include <memory>
#include <string>

class TwoCenterBundle
{
  public:
    TwoCenterBundle() = default;
    ~TwoCenterBundle() = default;

    //void build(const int nfile_orb,
    //           const std::string* file_orb,
    //           const int nfile_pp,
    //           const std::string* file_pp,
    //           const int nfile_desc = 0,
    //           const std::string* file_desc = nullptr);

    void build(const int ntype,
               const std::string* file_orb,
               Numerical_Nonlocal* const nl,
               const int nfile_desc = 0,
               const std::string* file_desc = nullptr);

    /**
     * @brief Overwrites the content of a LCAO_Orbitals object (e.g. GlobalC::ORB)
     * with the current object.
     *
     * This function provides an interface to the corresponding object in the old module_ao.
     */
    void to_LCAO_Orbitals(LCAO_Orbitals&) const;

    std::unique_ptr<TwoCenterIntegrator> kinetic_orb;
    std::unique_ptr<TwoCenterIntegrator> overlap_orb;
    std::unique_ptr<TwoCenterIntegrator> overlap_orb_beta;
    std::unique_ptr<TwoCenterIntegrator> overlap_orb_alpha;

  private:
    std::unique_ptr<RadialCollection> orb_;
    std::unique_ptr<RadialCollection> beta_;
    std::unique_ptr<RadialCollection> alpha_;
};

#endif
