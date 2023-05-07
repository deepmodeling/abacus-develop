#pragma once

#include <elpa/elpa_version.h>

#if ELPA_API_VERSION >= 20221101
#include <elpa/elpa.h>
#else
extern "C"
{
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
    // #include <elpa/elpa_generic.h>
#undef complex
    const char *elpa_strerr(int elpa_error);
}

#include "elpa_generic.hpp" // This is a wrapper for `elpa/elpa_generic.h`.
#endif
