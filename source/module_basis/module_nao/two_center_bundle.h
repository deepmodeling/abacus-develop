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

    // NOTE: some variables might be set only on RANK-0
    void build_orb(int ntype, const std::string* file_orb0);
    void build_beta(int ntype, Numerical_Nonlocal* nl);
    void build_alpha(int ndesc = 0, std::string* file_desc0 = nullptr);

    void tabulate();

    /**
     * @brief Overwrites the content of a LCAO_Orbitals object (e.g. GlobalC::ORB)
     * with the current object.
     *
     * This function provides an interface to the corresponding object in the old module_ao.
     */
    void to_LCAO_Orbitals(LCAO_Orbitals& orb,
                          const double lcao_ecut,
                          const double lcao_dk,
                          const double lcao_dr,
                          const double lcao_rmax
                          ) const;

    std::unique_ptr<TwoCenterIntegrator> kinetic_orb;
    std::unique_ptr<TwoCenterIntegrator> overlap_orb;
    std::unique_ptr<TwoCenterIntegrator> overlap_orb_beta;
    std::unique_ptr<TwoCenterIntegrator> overlap_orb_alpha;

    std::unique_ptr<RadialCollection> orb_;
    std::unique_ptr<RadialCollection> beta_;
    std::unique_ptr<RadialCollection> alpha_;
};

#endif
