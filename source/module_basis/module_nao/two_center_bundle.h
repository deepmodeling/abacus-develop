#ifndef TWO_CENTER_BUNDLE_H
#define TWO_CENTER_BUNDLE_H

#include "module_basis/module_nao/two_center_integrator.h"

#include <memory>
#include <string>

class TwoCenterBundle
{
  public:
    TwoCenterBundle() {}
    ~TwoCenterBundle();

    void build(const int nfile_orb,
               const std::string* file_orb,
               const int nfile_pp,
               const std::string* file_pp,
               const int nfile_desc = 0,
               const std::string* file_desc = nullptr);

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
