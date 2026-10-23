#ifndef ABACUS_SOURCE_RELAX_SOCKET_FRAME_UTILS_H
#define ABACUS_SOURCE_RELAX_SOCKET_FRAME_UTILS_H

#include "source_relax/socket_frame.h"

namespace FrameUtils
{
constexpr int kMatrixDimension = 3;
extern const int kMaxJacobiSweeps;

bool is_finite_matrix(const SocketFrame::Matrix9& values);
double column_norm_squared(const SocketFrame::Matrix9& values, int column);
double column_dot(const SocketFrame::Matrix9& values, int first, int second);
bool columns_are_orthogonal(const SocketFrame::Matrix9& values);
void rotate_columns(SocketFrame::Matrix9& values, int first, int second, double cosine, double sine);
bool one_sided_jacobi(SocketFrame::Matrix9& columns, SocketFrame::Matrix9& right_vectors);
long double scaled_determinant(const SocketFrame::Matrix9& values);
double received_inverse_residual(const SocketFrame::Matrix9& cell,
                                 const SocketFrame::Matrix9& inverse,
                                 bool transpose_inverse);
} // namespace FrameUtils

#endif
