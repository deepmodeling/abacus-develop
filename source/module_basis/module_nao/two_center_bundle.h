#ifndef W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_BASIS_MODULE_NAO_TWO_CENTER_BUNDLE_H
#define W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_BASIS_MODULE_NAO_TWO_CENTER_BUNDLE_H

#include "module_basis/module_ao/ORB_read.h"
#include "module_basis/module_nao/two_center_integrator.h"

#include <memory>
#include <string>

class TwoCenterBundle
{
  public:
    TwoCenterBundle() = default;
    ~TwoCenterBundle() = default;
    TwoCenterBundle& operator=(TwoCenterBundle&&) = default;

    // NOTE: some variables might be set only on RANK-0
    void build_orb(int ntype, const std::string* file_orb0);
    void build_beta(int ntype, Numerical_Nonlocal* nl);
    void build_alpha(int ndesc = 0, std::string* file_desc0 = nullptr);
    void build_orb_onsite(const double& radius);

    void tabulate();

    // Unlike the tabulate() above, this overload function computes
    // two-center integration table by direct integration with Simpson's
    // rule, which was the algorithm used prior to v3.3.4.
    void tabulate(const double lcao_ecut, const double lcao_dk, const double lcao_dr, const double lcao_rmax);

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
                          const double lcao_rmax) const;

    std::unique_ptr<TwoCenterIntegrator> kinetic_orb;
    std::unique_ptr<TwoCenterIntegrator> overlap_orb;
    std::unique_ptr<TwoCenterIntegrator> overlap_orb_beta;
    std::unique_ptr<TwoCenterIntegrator> overlap_orb_alpha;
    std::unique_ptr<TwoCenterIntegrator> overlap_orb_onsite;

    std::unique_ptr<RadialCollection> orb_;
    std::unique_ptr<RadialCollection> beta_;
    std::unique_ptr<RadialCollection> alpha_;
    std::unique_ptr<RadialCollection> orb_onsite_;
};

#endif
