#ifndef SPHERICAL_BESSEL_TRANSFORMER_H_ 
#define SPHERICAL_BESSEL_TRANSFORMER_H_

#include <fftw3.h>

namespace ModuleBase {

    /* 
     * A class to perform spherical Bessel transforms.
     *
     *
     * The spherical Bessel transform of a function F(x) is defined as
     *
     *                               / +inf     2           
     *          G(y) = sqrt(2/pi) *  |      dx x  F(x) j (x)
     *                               /  0               l
     *
     * where j (x) is the l-th order spherical Bessel funciton of the first kind. 
     *        l
     *
     * This class interprets the input array as
     *                   p
     *          in[i] = x [i] F(x[i])   (p is an input argument)
     *
     * and, on finish, fills the output array with 
     *
     *          out[j] = G(y[j])
     *
     * Currently the supported grids must follow
     *
     *          x[i] = i*dx             i = 0, 1, 2, ..., N
     *          y[j] = j*pi/(N*dx)      j = 0, 1, 2, ..., N
     *
     * That is, the input grid must be uniform and starts from 0, and there is no 
     * freedom to choose the output grid once the input grid is specified.
     *
     *
     * Usage:
     *
     * SphericalBesselTransformer sbt;
     *
     * // This flag will optimize the transforms at the cost of introducing large 
     * // overhead during planning the FFTs. Alternatively, one may use FFTW_ESTIMATE,
     * // which leads to less optimized transforms and much less overhead (or simply 
     * // skip this step as FFTW_ESTIMATE is used by default).
     * sbt.set_plan_flag(FFTW_MEASURE)
     *
     * // FFTW plan is created first and reused for consecutive same-sized transforms
     * sbt.radrfft(0, 1000, ...);
     * sbt.radrfft(1, 1000, ...);
     * sbt.radrfft(2, 1000, ...);
     *
     * // FFTW plan has to be re-created for a new size
     * sbt.radrfft(0, 2000, ...);
     *                                                                              */
    class SphericalBesselTransformer {
    
    public:

        SphericalBesselTransformer() {};
        ~SphericalBesselTransformer() { rfft_clean(); }

        SphericalBesselTransformer(SphericalBesselTransformer const&) = delete;
        SphericalBesselTransformer& operator=(SphericalBesselTransformer const&) = delete;
    
        /* 
         * Performs an l-th order spherical Bessel transform via real-input fast Fourier transforms. 
         *
         * This function computes the spherical Bessel transform F(x) -> G(y) with input values
         *
         *                   p
         *          in[i] = x [i] F(x[i]) 
         *
         * where p is an arbitrary integer, and
         *
         *                     cutoff
         *          x[i] = i * -------          i = 0, 1, 2,..., ngrid-1.
         *                     ngrid-1
         * 
         * On exit, out[j] = G(y[j]) where
         *
         *                      pi
         *          y[j] = j * ------           j = 0, 1, 2,..., ngrid-1.
         *                     cutoff
         *
         *
         * Caveats:
         *
         * 1. This function does not allocate memory for output array; it must be pre-allocated.
         * 2. F(x) is supposed to be exactly zero at and after cutoff. Results would not make
         *    sense if input is truncated at a place where F(x) is still significantly non-zero.
         *                                                                                      */
        void radrfft(
                int l,
                int ngrid,
                double cutoff,
                double* in,
                double* out,
                int p = 0
        );


        /* 
         * Sets the FFTW planner flag.
         *
         * Saved fftw_plan will be destroyed if it was created with a flag other than new_flag.
         *                                                                                      */
        void set_plan_flag(unsigned new_flag);
    

    private:
    
        /* Buffers for in-place real-input FFT (interpreted as double* on input) */
        fftw_complex* f_ = nullptr;
    
        /* FFTW plan saved for reuse */
        fftw_plan rfft_plan_ = nullptr;

        /* Size of the planned rFFT */
        int sz_planned_ = -1;

        /* Planner flag used to create rfft_plan_ */
        unsigned fftw_plan_flag_ = FFTW_ESTIMATE;
    
        /* Applies an in-place 1-d real-input discrete Fourier transform to f_ */
        void rfft_in_place();
    
        /* Buffer (f_) allocation and plan creation for an N-element real-input FFT */
        void rfft_prepare(int N);

        /* Buffer deallocation and plan destruction */
        void rfft_clean();
    
        /* 
         * Returns the polynomial coefficient of the n-th power term in the sin+cos
         * representation of the l-th order spherical Bessel functions of the first kind.
         *
         * The l-th order spherical Bessel function of the first kind can be expressed as
         *
         *                  sin(x)*P(x) + cos(x)*Q(x)
         *          j (x) = --------------------------
         *           l               l+1
         *                          x
         *
         * where P(x) and Q(x) are polynomials of order no more than l. This function returns
         * the coefficients within those polynomials.
         *
         * Caveats:
         * 1. Coefficients grow very quickly as l increases. Currently l is capped at 10 since
         * Some coefficients exceed INT_MAX (2^31-1) for l >= 11.
         *                                                                                      */
        int spherical_bessel_sincos_polycoef(bool get_sine, int l, int n);
    
    };
}

#endif
