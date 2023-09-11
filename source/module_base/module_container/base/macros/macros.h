#ifndef BASE_MACROS_MACROS_H_
#define BASE_MACROS_MACROS_H_

#include <stdint.h>

#if __CUDA
#include <base/macros/cuda.h> 
#endif 

#include <ATen/core/tensor_types.h>

#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
  TypeName(const TypeName&) = delete;      \
  void operator=(const TypeName&) = delete

#define DISALLOW_MOVE_AND_ASSIGN(TypeName) \
  TypeName(TypeName&&) = delete;           \
  void operator=(TypeName&&) = delete

#define DISALLOW_COPY_MOVE_AND_ASSIGN(TypeName) \
  DISALLOW_COPY_AND_ASSIGN(TypeName);           \
  DISALLOW_MOVE_AND_ASSIGN(TypeName)

#define DISALLOW_IMPLICIT_CONSTRUCTORS(TypeName) \
  TypeName() = delete;                           \
  DISALLOW_COPY_MOVE_AND_ASSIGN(TypeName)

#define MAX_SIZE_T UINT64_MAX


// The macro TEMPLATE_1() expands to a switch statement conditioned on
// TYPE_ENUM. Each case expands the STMTS after a typedef for T.
#define SINGLE_ARG(...) __VA_ARGS__

#define CASE_2(TYPE, DEVICE, STMTS)               \
  case (int(DataTypeToEnum<TYPE>::value) * 10 +   \
        int(DeviceTypeToEnum<DEVICE>::value)): {  \
    typedef TYPE T_;                              \
    typedef DEVICE DEVICE_;                       \
    STMTS;                                        \
    break;                                        \
  }

#define CASE_3(TYPE, DEVICE_OUT, DEVICE_IN, STMTS)      \
  case (int(DataTypeToEnum<TYPE>::value) * 100 +        \
        int(DeviceTypeToEnum<DEVICE_OUT>::value) * 10 + \
        int(DeviceTypeToEnum<DEVICE_IN>::value)): {     \
    typedef TYPE T_;                                    \
    typedef DEVICE_OUT DEVICE_OUT_;                     \
    typedef DEVICE_IN DEVICE_IN_;                       \
    STMTS;                                              \
    break;                                              \
  }

#define CASES_ALL_WITH_DEFAULT_2_CPU(TYPE_ENUM, DEVICE_ENUM, STMTS, DEFAULT) \
  switch (int(TYPE_ENUM) * 10 + int(DEVICE_ENUM)) {                      \
    CASE_2(float, DEVICE_CPU, SINGLE_ARG(STMTS))                         \
    CASE_2(double, DEVICE_CPU, SINGLE_ARG(STMTS))                        \
    CASE_2(int, DEVICE_CPU, SINGLE_ARG(STMTS))                           \
    CASE_2(int64_t, DEVICE_CPU, SINGLE_ARG(STMTS))                       \
    CASE_2(std::complex<float>, DEVICE_CPU, SINGLE_ARG(STMTS))           \
    CASE_2(std::complex<double>, DEVICE_CPU, SINGLE_ARG(STMTS))          \
    default:                                                             \
      DEFAULT;                                                           \
      break;                                                             \
  }

#define CASES_BLAS_WITH_DEFAULT_2_CPU(TYPE_ENUM, DEVICE_ENUM, STMTS, DEFAULT) \
  switch (int(TYPE_ENUM) * 10 + int(DEVICE_ENUM)) {                      \
    CASE_2(float, DEVICE_CPU, SINGLE_ARG(STMTS))                         \
    CASE_2(double, DEVICE_CPU, SINGLE_ARG(STMTS))                        \
    CASE_2(std::complex<float>, DEVICE_CPU, SINGLE_ARG(STMTS))           \
    CASE_2(std::complex<double>, DEVICE_CPU, SINGLE_ARG(STMTS))          \
    default:                                                             \
      DEFAULT;                                                           \
      break;                                                             \
  }

#define CASES_ALL_WITH_DEFAULT_2(TYPE_ENUM, DEVICE_ENUM, STMTS, DEFAULT) \
  switch (int(TYPE_ENUM) * 10 + int(DEVICE_ENUM)) {                      \
    CASE_2(float, DEVICE_CPU, SINGLE_ARG(STMTS))                         \
    CASE_2(float, DEVICE_GPU, SINGLE_ARG(STMTS))                         \
    CASE_2(double, DEVICE_CPU, SINGLE_ARG(STMTS))                        \
    CASE_2(double, DEVICE_GPU, SINGLE_ARG(STMTS))                        \
    CASE_2(int, DEVICE_CPU, SINGLE_ARG(STMTS))                           \
    CASE_2(int, DEVICE_GPU, SINGLE_ARG(STMTS))                           \
    CASE_2(int64_t, DEVICE_CPU, SINGLE_ARG(STMTS))                       \
    CASE_2(int64_t, DEVICE_GPU, SINGLE_ARG(STMTS))                       \
    CASE_2(std::complex<float>, DEVICE_CPU, SINGLE_ARG(STMTS))           \
    CASE_2(std::complex<float>, DEVICE_GPU, SINGLE_ARG(STMTS))           \
    CASE_2(std::complex<double>, DEVICE_CPU, SINGLE_ARG(STMTS))          \
    CASE_2(std::complex<double>, DEVICE_GPU, SINGLE_ARG(STMTS))          \
    default:                                                             \
      DEFAULT;                                                           \
      break;                                                             \
  }

#define CASES_BLAS_WITH_DEFAULT_2(TYPE_ENUM, DEVICE_ENUM, STMTS, DEFAULT) \
  switch (int(TYPE_ENUM) * 10 + int(DEVICE_ENUM)) {                      \
    CASE_2(float, DEVICE_CPU, SINGLE_ARG(STMTS))                         \
    CASE_2(float, DEVICE_GPU, SINGLE_ARG(STMTS))                         \
    CASE_2(double, DEVICE_CPU, SINGLE_ARG(STMTS))                        \
    CASE_2(double, DEVICE_GPU, SINGLE_ARG(STMTS))                        \
    CASE_2(std::complex<float>, DEVICE_CPU, SINGLE_ARG(STMTS))           \
    CASE_2(std::complex<float>, DEVICE_GPU, SINGLE_ARG(STMTS))           \
    CASE_2(std::complex<double>, DEVICE_CPU, SINGLE_ARG(STMTS))          \
    CASE_2(std::complex<double>, DEVICE_GPU, SINGLE_ARG(STMTS))          \
    default:                                                             \
      DEFAULT;                                                           \
      break;                                                             \
  }

#define CASES_CZ_WITH_DEFAULT_2(TYPE_ENUM, DEVICE_ENUM, STMTS, DEFAULT)  \
  switch (int(TYPE_ENUM) * 10 + int(DEVICE_ENUM)) {                      \
    CASE_2(std::complex<float>, DEVICE_CPU, SINGLE_ARG(STMTS))           \
    CASE_2(std::complex<float>, DEVICE_GPU, SINGLE_ARG(STMTS))           \
    CASE_2(std::complex<double>, DEVICE_CPU, SINGLE_ARG(STMTS))          \
    CASE_2(std::complex<double>, DEVICE_GPU, SINGLE_ARG(STMTS))          \
    default:                                                             \
      DEFAULT;                                                           \
      break;                                                             \
  }

#define CASES_ALL_WITH_DEFAULT_3(TYPE_ENUM, DEVICE_ENUM_OUT, DEVICE_ENUM_IN, STMTS, DEFAULT) \
  switch (int(TYPE_ENUM) * 100 + int(DEVICE_ENUM_OUT) * 10 + int(DEVICE_ENUM_IN)) {          \
    CASE_3(float, DEVICE_CPU, DEVICE_CPU, SINGLE_ARG(STMTS))                                 \
    CASE_3(float, DEVICE_GPU, DEVICE_CPU, SINGLE_ARG(STMTS))                                 \
    CASE_3(double, DEVICE_CPU, DEVICE_CPU, SINGLE_ARG(STMTS))                                \
    CASE_3(double, DEVICE_GPU, DEVICE_CPU, SINGLE_ARG(STMTS))                                \
    CASE_3(int, DEVICE_CPU, DEVICE_CPU, SINGLE_ARG(STMTS))                                   \
    CASE_3(int, DEVICE_GPU, DEVICE_CPU, SINGLE_ARG(STMTS))                                   \
    CASE_3(int64_t, DEVICE_CPU, DEVICE_CPU, SINGLE_ARG(STMTS))                               \
    CASE_3(int64_t, DEVICE_GPU, DEVICE_CPU, SINGLE_ARG(STMTS))                               \
    CASE_3(std::complex<float>, DEVICE_CPU, DEVICE_CPU, SINGLE_ARG(STMTS))                   \
    CASE_3(std::complex<float>, DEVICE_CPU, DEVICE_GPU, SINGLE_ARG(STMTS))                   \
    CASE_3(std::complex<double>, DEVICE_CPU, DEVICE_CPU, SINGLE_ARG(STMTS))                  \
    CASE_3(std::complex<double>, DEVICE_CPU, DEVICE_GPU, SINGLE_ARG(STMTS))                  \
    CASE_3(float, DEVICE_CPU, DEVICE_GPU, SINGLE_ARG(STMTS))                                 \
    CASE_3(float, DEVICE_GPU, DEVICE_GPU, SINGLE_ARG(STMTS))                                 \
    CASE_3(double, DEVICE_CPU, DEVICE_GPU, SINGLE_ARG(STMTS))                                \
    CASE_3(double, DEVICE_GPU, DEVICE_GPU, SINGLE_ARG(STMTS))                                \
    CASE_3(int, DEVICE_CPU, DEVICE_GPU, SINGLE_ARG(STMTS))                                   \
    CASE_3(int, DEVICE_GPU, DEVICE_GPU, SINGLE_ARG(STMTS))                                   \
    CASE_3(int64_t, DEVICE_CPU, DEVICE_GPU, SINGLE_ARG(STMTS))                               \
    CASE_3(int64_t, DEVICE_GPU, DEVICE_GPU, SINGLE_ARG(STMTS))                               \
    CASE_3(std::complex<float>, DEVICE_GPU, DEVICE_CPU, SINGLE_ARG(STMTS))                   \
    CASE_3(std::complex<float>, DEVICE_GPU, DEVICE_GPU, SINGLE_ARG(STMTS))                   \
    CASE_3(std::complex<double>, DEVICE_GPU, DEVICE_CPU, SINGLE_ARG(STMTS))                  \
    CASE_3(std::complex<double>, DEVICE_GPU, DEVICE_GPU, SINGLE_ARG(STMTS))                  \
    default:                                                                                 \
      DEFAULT;                                                                               \
      break;                                                                                 \
  }


#if __CUDA || __ROCM
#define TEMPLATE_ALL_2(TYPE_ENUM, DEVICE_ENUM, ...)                        \
    CASES_ALL_WITH_DEFAULT_2(TYPE_ENUM, DEVICE_ENUM, (__VA_ARGS__),        \
                           std::cerr << "Unexpected type: " << TYPE_ENUM; exit(EXIT_FAILURE));
#define TEMPLATE_BLAS_2(TYPE_ENUM, DEVICE_ENUM, ...)                        \
    CASES_BLAS_WITH_DEFAULT_2(TYPE_ENUM, DEVICE_ENUM, (__VA_ARGS__),        \
                           std::cerr << "Unexpected type: " << TYPE_ENUM; exit(EXIT_FAILURE));
#else 
#define TEMPLATE_ALL_2(TYPE_ENUM, DEVICE_ENUM, ...)                    \
CASES_ALL_WITH_DEFAULT_2_CPU(TYPE_ENUM, DEVICE_ENUM, (__VA_ARGS__),    \
                       std::cerr << "Unexpected type: " << TYPE_ENUM; exit(EXIT_FAILURE));
#define TEMPLATE_BLAS_2(TYPE_ENUM, DEVICE_ENUM, ...)                        \
CASES_BLAS_WITH_DEFAULT_2_CPU(TYPE_ENUM, DEVICE_ENUM, (__VA_ARGS__),        \
                       std::cerr << "Unexpected type: " << TYPE_ENUM; exit(EXIT_FAILURE));
#endif 

#define TEMPLATE_CZ_2(TYPE_ENUM, DEVICE_ENUM, ...)                         \
    CASES_CZ_WITH_DEFAULT_2(TYPE_ENUM, DEVICE_ENUM, (__VA_ARGS__),         \
                           std::cerr << "Unexpected type: " << TYPE_ENUM; exit(EXIT_FAILURE));

#define TEMPLATE_ALL_3(TYPE_ENUM, DEVICE_ENUM_OUT, DEVICE_ENUM_IN, ...)                        \
    CASES_ALL_WITH_DEFAULT_3(TYPE_ENUM, DEVICE_ENUM_OUT, DEVICE_ENUM_IN, (__VA_ARGS__),        \
                           std::cerr << "Unexpected type: " << TYPE_ENUM; exit(EXIT_FAILURE));


#endif // BASE_MACROS_MACROS_H_