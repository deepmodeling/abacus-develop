#ifndef STERNHEIMER_WEAK_AUGMENTED_H
#define STERNHEIMER_WEAK_AUGMENTED_H

#include "source_lcao/module_ri/sternheimer_rpa.h"

#include <cstddef>
#include <functional>
#include <memory>
#include <vector>

namespace ModuleRI
{

// Low-level algebra only, with no production response/FFT/weak-form assembly.
// U=[O,F] is fine-metric orthonormal; E is normalized so E*E=I. All supplied
// coefficients and inner products are Euclidean: never insert dVc or dVf here.
// This helper cannot verify the basis normalization or full-Pc positivity.
class SternheimerWeakAugmented
{
  public:
    using Complex = SternheimerRPA::Complex;
    using Vector = SternheimerRPA::Vector;
    using Apply = std::function<void(const Vector&, Vector&)>;

    struct Data
    {
        int nocc = 0;
        int nvirtual = 0;
        int ncoarse = 0;
        // All matrices are packed COLUMN-MAJOR, index row + leading_dim*col.
        // hu: nu x nu; l and k: nu x ncoarse; nu=nocc+nvirtual.
        // hu=h(U,U), l=U*E, k=h(U,E). Occupied rows precede virtual rows.
        Vector hu;
        Vector l;
        Vector k;
    };

    struct Vertices
    {
        Vector f;
        Vector w;
    };

    struct Expansion
    {
        // Response U*u + E*e = F*y + (E-U*L)*x, not a fine-grid vector.
        Vector u;
        Vector e;
    };

    // Owns the supplied matrices, with no mutable shared application scratch.
    // Budget covers retained matrix capacities plus metric eigensolver scratch;
    // excludes caller-owned input copies, object/allocator and LAPACK overhead.
    // Diagnoses min eig(S_W)=1-lambda_max(L L*) in nu space, never dense coarse
    // space. Throws domain_error on singular/ill-conditioned/indefinite metric.
    // metric_tolerance is a rejection threshold, NOT a regularization or shift.
    explicit SternheimerWeakAugmented(Data data,
                                      std::size_t max_payload_bytes = 512ULL * 1024 * 1024,
                                      double metric_tolerance = 1.0e-12,
                                      double hermitian_tolerance = 1.0e-12);
    SternheimerWeakAugmented(const SternheimerWeakAugmented&) = delete;
    SternheimerWeakAugmented& operator=(const SternheimerWeakAugmented&) = delete;

    int nocc() const;
    int nvirtual() const;
    int ncoarse() const;
    std::size_t storage_bytes() const;
    // Returns 1 for the empty complement; it is not a Hamiltonian spectrum test.
    double minimum_metric_eigenvalue() const;

    // Augmented vectors always concatenate [y;x]. Input/output may alias.
    void apply_metric(const Vector& coefficients, Vector& output) const;
    // Exact positive-metric coordinate normalization T=(I-L*L)^(-1/2).
    // Low-rank application; no coarse square matrix, shift or mode removal.
    void apply_complement_inverse_sqrt(const Vector& x, Vector& output) const;
    // gF=(gU)_F; gW=gE-L* gU. Inputs are already projected fine-grid vertices.
    Vertices project_vertices(const Vector& g_u, const Vector& g_e) const;
    Expansion expand_coordinates(const Vector& coefficients) const;

    class Worker
    {
      public:
        struct SolveResult
        {
            Vector coefficients;
            SternheimerRPA::SolverResult schur;
            bool converged = false;
            // Recomputed residual of BOTH original augmented equations,
            // normalized by ||[gF;gW]|| (absolute residual if that norm is zero).
            double absolute_residual = 0.0;
            double relative_residual = 0.0;
        };

        // CE applies h(E,E), must be fixed, linear and Hermitian, and must return
        // exactly ncoarse finite entries. Hermiticity/linearity of an arbitrary
        // callback is a caller contract, not established by finite probes here.
        // A worker owns its callback and LU(A-zI). Concurrent workers may share
        // Blocks; any shared callback state must itself be thread-safe. Callbacks
        // receive distinct input/output vectors, with output initially zeroed.
        // z=eps_source-i*omega; signed omega supports the adjoint branch. No
        // energy shift, eigenmode truncation or pseudoinverse is performed.
        // Rejects singular LU and a ZGECON reciprocal condition estimate at or
        // below double epsilon, reporting an unresolved resolvent, not fixing it.
        Worker(std::shared_ptr<const SternheimerWeakAugmented> blocks,
               Apply ce,
               double eps_source,
               double omega,
               std::size_t max_workspace_bytes = 512ULL * 1024 * 1024);
        Worker(const Worker&) = delete;
        Worker& operator=(const Worker&) = delete;

        void apply_hamiltonian(const Vector& coefficients, Vector& output) const;
        void apply_shifted(const Vector& coefficients, Vector& output) const;
        void apply_schur(const Vector& x, Vector& output) const;
        Vector schur_rhs(const Vertices& g) const;
        Vector reconstruct(const Vertices& g, const Vector& x) const;

        // Starts at zero and uses GMRES on T Q T p = T b, x=T p, where T is
        // the exact inverse square root of the complement overlap. Residuals
        // are still checked in the original untransformed equations.
        // No kinetic preconditioner is used. Nonconvergence is explicit; numerical
        // breakdown, invalid/nonfinite data and singular LU raise exceptions.
        SolveResult solve(const Vertices& g,
                          const SternheimerRPA::SolverOptions& options,
                          int restart_dimension = 50) const;

        // Conservative numerical-payload bound including LU, temporary vectors
        // and GMRES/least-squares scratch. Excludes shared blocks, caller-owned
        // inputs/outputs, callback storage/scratch, allocator/runtime overhead.
        // This is not an RSS bound. Checked before LU and before GMRES allocation.
        std::size_t workspace_bytes_required(int restart_dimension) const;

      private:
        Vector solve_virtual(const Vector& rhs) const;
        void apply_complement(const Vector& x, Vector& output) const;
        void validate_vertices(const Vertices& g) const;

        std::shared_ptr<const SternheimerWeakAugmented> blocks_;
        Apply ce_;
        Complex z_;
        std::size_t max_workspace_bytes_;
        Vector lu_;
        std::vector<int> pivots_;
    };

  private:
    void validate_coefficients(const Vector& coefficients) const;
    void apply_complement_metric(const Vector& x, Vector& output) const;
    void apply_b(const Vector& x, Vector& output) const;
    void apply_b_adjoint(const Vector& y, Vector& output) const;

    Data data_;
    Vector inverse_sqrt_kernel_;
    int nu_ = 0;
    double minimum_metric_eigenvalue_ = 1.0;
};

} // namespace ModuleRI

#endif
