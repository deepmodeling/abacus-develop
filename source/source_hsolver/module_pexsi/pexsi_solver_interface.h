#ifndef PEXSI_SOLVER_INTERFACE_H
#define PEXSI_SOLVER_INTERFACE_H

namespace pexsi
{
class IPexsiSolver
{
  public:
    virtual ~IPexsiSolver() = default;

    virtual void prepare(const int blacs_text,
                         const int nb,
                         const int nrow,
                         const int ncol,
                         const double* h,
                         const double* s,
                         double*& DM,
                         double*& EDM)
        = 0;
    virtual int solve(double mu0) = 0;
    virtual double get_totalFreeEnergy() const = 0;
    virtual double get_totalEnergyH() const = 0;
    virtual double get_totalEnergyS() const = 0;
    virtual double get_mu() const = 0;
};
} // namespace pexsi

#endif // PEXSI_SOLVER_INTERFACE_H
