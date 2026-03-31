//
// Created by rhx on 25-6-26.
//

#ifndef VEC_MUL_VEC_OP_H
#define VEC_MUL_VEC_OP_H
namespace hamilt {

template <typename T, typename Device>
struct vec_mul_vec_complex_op
{
    // Multiply a vector with a complex vector
    void operator()(const T *vec1, const T *vec2, T *out, int n);
};
} // namespace hamilt
#endif //VEC_MUL_VEC_OP_H
