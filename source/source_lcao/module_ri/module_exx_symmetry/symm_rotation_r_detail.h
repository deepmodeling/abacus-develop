#ifndef ABACUS_SOURCE_LCAO_MODULE_RI_MODULE_EXX_SYMMETRY_SYMM_ROTATION_R_DETAIL_H
#define ABACUS_SOURCE_LCAO_MODULE_RI_MODULE_EXX_SYMMETRY_SYMM_ROTATION_R_DETAIL_H

// Small helpers that the real-space rotation templates need to see. They are
// non-template, so they live in a plain header instead of the implementation
// header `symm_rotation_r.h`, which holds only what must be visible at the
// instantiation site.

#include "source_base/timer.h"

#include <complex>

namespace ModuleSymmetry
{
    /// Balance ModuleBase::timer for a scope that may exit by throwing.
    ///
    /// The restore helpers validate their input and throw on a malformed
    /// tensor; without this the matching `timer::end` would be skipped and the
    /// process-global entry would stay open, so the next `timer::start` on the
    /// same name aborts.
    class ScopedTimer
    {
      public:
        ScopedTimer(const char* class_name, const char* name)
            : class_name_(class_name), name_(name)
        {
            ModuleBase::timer::start(class_name_, name_);
        }
        ~ScopedTimer()
        {
            ModuleBase::timer::end(class_name_, name_);
        }
        ScopedTimer(const ScopedTimer&) = delete;
        ScopedTimer& operator=(const ScopedTimer&) = delete;

      private:
        const char* class_name_;
        const char* name_;
    };

    /// Elementwise complex conjugation used by the time-reversal branch of
    /// restore_HR_nspin4. Overloaded rather than specialized, so a real Tdata
    /// compiles to the identity.
    inline float conj_elem(const float v) { return v; }
    inline double conj_elem(const double v) { return v; }
    inline std::complex<float> conj_elem(const std::complex<float>& v) { return std::conj(v); }
    inline std::complex<double> conj_elem(const std::complex<double>& v) { return std::conj(v); }
} // namespace ModuleSymmetry

#endif
