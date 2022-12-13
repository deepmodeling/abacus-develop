#ifndef ELPA_H
#define ELPA_H

#ifdef __ELPA

#include <complex>
#include <elpa/elpa_version.h>
#include <limits.h>

struct elpa_struct;
typedef struct elpa_struct *elpa_t;

struct elpa_autotune_struct;
typedef struct elpa_autotune_struct *elpa_autotune_t;

#include <elpa/elpa_constants.h>
#include <elpa/elpa_generated_c_api.h>
// ELPA only provides a C interface header, causing inconsistence of complex
// between C99 (e.g. double complex) and C++11 (std::complex).
// Thus, we have to define a wrapper of complex over the c api
// for compatiability.
#define complex _Complex
#include <elpa/elpa_generated.h>
#undef complex
#include <elpa/elpa_generic.h>

#define ELPA_2STAGE_REAL_GPU    ELPA_2STAGE_REAL_NVIDIA_GPU
#define ELPA_2STAGE_COMPLEX_GPU ELPA_2STAGE_COMPLEX_NVIDIA_GPU

const char *elpa_strerr(int elpa_error);

#endif
#endif
