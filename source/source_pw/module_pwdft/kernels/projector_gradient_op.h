#ifndef PW_PROJECTOR_GRADIENT_OP_H
#define PW_PROJECTOR_GRADIENT_OP_H
#include "source_base/module_device/types.h"
#include <complex>
#include <cmath>
namespace hamilt
{
/** @brief Evaluate a normalized real solid harmonic and its Cartesian derivatives. */
template <typename Real>
#if defined(__CUDACC__) || defined(__HIPCC__)
__host__ __device__
#endif
void gradient_harmonic(const int l, const int lm, const Real* q, Real* value)
{
    const int index = lm - l * l;
    const int m = (index + 1) / 2;
    const bool imaginary = index > 0 && index % 2 == 0;
    const Real radius = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]);
    for (int d = 0; d < 4; ++d) { value[d] = Real(0); }
    if (radius < Real(1e-9))
    {
        if (l == 0) { value[0] = Real(0.28209479177387814347); }
        return;
    }
    const Real x = q[0]/radius;
    const Real y = q[1]/radius;
    const Real z = q[2]/radius;
    Real re[4] = {Real(1), Real(0), Real(0), Real(0)};
    Real im[4] = {Real(0), Real(0), Real(0), Real(0)};
    // Cartesian solid harmonics avoid polar-axis singularities.
    for (int k = 1; k <= m; ++k)
    {
        const Real factor = Real(-(2*k-1));
        const Real old_re = re[0];
        const Real old_im = im[0];
        for (int d = 0; d < 4; ++d)
        {
            const Real r = re[d];
            re[d] = factor * (x*r-y*im[d] + (d==1 ? old_re : Real(0)) - (d==2 ? old_im : Real(0)));
            im[d] = factor * (x*im[d]+y*r + (d==1 ? old_im : Real(0)) + (d==2 ? old_re : Real(0)));
        }
    }
    Real previous[4] = {};
    Real current[4];
    for (int d = 0; d < 4; ++d) { current[d] = imaginary ? im[d] : re[d]; }
    for (int n = m+1; n <= l; ++n)
    {
        Real next[4];
        for (int d = 0; d < 4; ++d)
        {
            const Real coordinate = d==1 ? x : (d==2 ? y : z);
            next[d] = (Real(2*n-1)*(z*current[d] + (d==3 ? current[0] : Real(0)))
                       - Real(n+m-1)*(previous[d] + (d>0 ? Real(2)*coordinate*previous[0] : Real(0)))) / Real(n-m);
        }
        for (int d = 0; d < 4; ++d) { previous[d]=current[d]; current[d]=next[d]; }
    }
    Real normalization = Real(2*l+1) / Real(12.56637061435917295385);
    for (int k=l-m+1; k<=l+m; ++k) { normalization /= Real(k); }
    normalization = sqrt(normalization * (m==0 ? Real(1) : Real(2)));
    value[0] = normalization * current[0];
    for (int d=1; d<4; ++d)
    {
        value[d] = normalization * (current[d] - Real(l)*current[0]*q[d-1]/radius) / radius;
    }
}

/** @brief Device projector-gradient assembly; metadata is [type,beta,l,lm,atom]. */
template <typename Real, typename Device>
class projector_gradient_op
{
  public:
    void operator()(int npw, int ld, int nkb, int nbeta, int nq, Real dq, Real tpiba,
                    const int* metadata, const Real* q, const Real* tab, const Real* derivative,
                    const std::complex<Real>* sk, std::complex<Real>* out) const;
};

#if defined(__CUDA) || defined(__ROCM)
template <typename Real>
class projector_gradient_op<Real, base_device::DEVICE_GPU>
{
  public:
    void operator()(int npw, int ld, int nkb, int nbeta, int nq, Real dq, Real tpiba,
                    const int* metadata, const Real* q, const Real* tab, const Real* derivative,
                    const std::complex<Real>* sk, std::complex<Real>* out) const;
};
#endif

template <typename Real>
#if defined(__CUDACC__) || defined(__HIPCC__)
__host__ __device__
#endif
Real gradient_interpolate(const Real* table, const Real position)
{
    const int index = static_cast<int>(position);
    const Real x = position - Real(index);
    return table[index]*(Real(1)-x)*(Real(2)-x)*(Real(3)-x)/Real(6)
           + table[index+1]*x*(Real(2)-x)*(Real(3)-x)/Real(2)
           - table[index+2]*x*(Real(1)-x)*(Real(3)-x)/Real(2)
           + table[index+3]*x*(Real(1)-x)*(Real(2)-x)/Real(6);
}

template <typename Real, typename Complex>
#if defined(__CUDACC__) || defined(__HIPCC__)
__host__ __device__
#endif
void gradient_element(const int ig, const int projector, const int npw, const int ld, const int nkb,
                      const int nbeta, const int nq, const Real dq, const Real tpiba,
                      const int* metadata, const Real* q, const Real* tab, const Real* derivative,
                      const Complex* sk, Complex* out)
{
    const int* meta = metadata + 5*projector;
    const int l = meta[2];
    const int lm = meta[3];
    const Real* momentum = q+3*ig;
    const Real radius = sqrt(momentum[0]*momentum[0]+momentum[1]*momentum[1]+momentum[2]*momentum[2]);
    const Real position = radius*tpiba/dq;
    Real harmonic[4];
    gradient_harmonic(l, lm, momentum, harmonic);
    const Real* radial = tab + (meta[0]*nbeta+meta[1])*nq;
    const Real* radial_d = derivative + (meta[0]*nbeta+meta[1])*nq;
    const Real v = gradient_interpolate(radial, position);
    const Real dv = gradient_interpolate(radial_d, position);
    const int phase_index = l % 4;
    const Complex phase = phase_index==0 ? Complex(1,0) : (phase_index==1 ? Complex(0,-1)
                                         : (phase_index==2 ? Complex(-1,0) : Complex(0,1)));
    for (int d=0; d<3; ++d)
    {
        Real gradient = Real(0);
        if (radius < Real(1e-8))
        {
            if (lm==1 && d==2) { gradient = dv*Real(0.48860251190291992159); }
            if ((lm==2 && d==0) || (lm==3 && d==1)) { gradient = -dv*Real(0.48860251190291992159); }
        }
        else
        {
            gradient = harmonic[0]*dv*momentum[d]/radius + harmonic[d+1]*v/tpiba;
        }
        out[(d*nkb+projector)*ld+ig] = gradient * sk[meta[4]*npw+ig] * phase;
    }
}
} // namespace hamilt
#endif
