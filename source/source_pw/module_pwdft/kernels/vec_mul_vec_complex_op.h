//
// Created by rhx on 25-6-26.
//

#ifndef VEC_MUL_VEC_OP_H
#define VEC_MUL_VEC_OP_H
namespace hamilt {

template <typename T, typename Device>
struct vec_mul_vec_complex_op
{
    // Single element-wise operation: out[i] = vec1[i] * vec2[i]
    void operator()(const T *vec1, const T *vec2, T *out, int n);

    // Batched element-wise operation: process multiple vector pairs in parallel
    // vec1_batch, vec2_batch, out_batch: contiguous batches (batch_size × n elements)
    // n: elements per vector
    // batch_size: number of vector pairs to process
    void operator_batch(
        const T *vec1_batch,
        const T *vec2_batch,
        T *out_batch,
        int n,
        int batch_size);
};
} // namespace hamilt
#endif //VEC_MUL_VEC_OP_H
