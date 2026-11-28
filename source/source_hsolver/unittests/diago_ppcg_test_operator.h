#ifndef ABACUS_PPCG_TEST_OPERATOR_H
#define ABACUS_PPCG_TEST_OPERATOR_H

#include "../diago_ppcg.h"
#include "../hs_operator.h"

#include <algorithm>
#include <functional>
#include <iostream>

namespace ppcg_test
{
// Adapt analytical test matrices to the same operator interface used by PW.
// Identity overlap is a fixture choice; production supplies its real S operator.
template <typename Scalar>
class CallbackHSOperator : public hsolver::HSOperator<Scalar, base_device::DEVICE_CPU>
{
  public:
    using Function = std::function<void(const Scalar *, Scalar *, int, int)>;

    CallbackHSOperator(const Function &h, const Function &s) : h_(h), s_(s) {}

    void update_k(const int) override {}

    void hpsi(
        const Scalar *in,
        Scalar *out,
        const int ld,
        const int nvec) const override
    {
        h_(in, out, ld, nvec);
    }

    void spsi(
        const Scalar *in,
        Scalar *out,
        const int ld,
        const int nvec) const override
    {
        if (s_)
        {
            s_(in, out, ld, nvec);
        }
        else
        {
            const int count = ld * nvec;
            std::copy(in, in + count, out);
        }
    }

  private:
    Function h_;
    Function s_;
};

template <typename Scalar>
double diagonalize(
    hsolver::DiagoPPCG<Scalar, base_device::DEVICE_CPU> &solver,
    const typename CallbackHSOperator<Scalar>::Function &h,
    const typename CallbackHSOperator<Scalar>::Function &s,
    int ld,
    int bands,
    int dimension,
    Scalar *psi,
    typename hsolver::DiagoPPCG<Scalar, base_device::DEVICE_CPU>::Real *eigenvalues,
    const std::vector<double> &thresholds,
    const typename hsolver::DiagoPPCG<Scalar, base_device::DEVICE_CPU>::Real *preconditioner,
    std::ostream &log = std::clog)
{
    const CallbackHSOperator<Scalar> op(h, s);
    const hsolver::PPCGTraceContext context = {0, 0};
    return solver.diag(op, ld, bands, dimension, psi, eigenvalues, thresholds, preconditioner, log, context);
}
} // namespace ppcg_test

#endif
