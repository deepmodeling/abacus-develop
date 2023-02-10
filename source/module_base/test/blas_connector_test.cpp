#include "../blas_connector.h"
#include "gtest/gtest.h"

#include <algorithm>
#include <array>
#include <complex>
#include <cstdlib>
TEST(blas_connector, sscal_) {
    typedef float T;
    const int size = 8;
    const T scale = 2;
    const int incx = 1;
    std::array<T, size> result, answer;
    std::generate(result.begin(), result.end(),
                  []() { return std::rand() / T(RAND_MAX); });
    for (int i = 0; i < size; i++)
        answer[i] = result[i] * scale;
    sscal_(&size, &scale, result.data(), &incx);
    for (int i = 0; i < size; i++)
        EXPECT_FLOAT_EQ(answer[i], result[i]);
}

TEST(blas_connector, dscal_) {
    typedef double T;
    const int size = 8;
    const T scale = 2;
    const int incx = 1;
    std::array<T, size> result, answer;
    std::generate(result.begin(), result.end(),
                  []() { return std::rand() / T(RAND_MAX); });
    for (int i = 0; i < size; i++)
        answer[i] = result[i] * scale;
    dscal_(&size, &scale, result.data(), &incx);
    for (int i = 0; i < size; i++)
        EXPECT_DOUBLE_EQ(answer[i], result[i]);
}

TEST(blas_connector, cscal_) {
    typedef std::complex<float> T;
    const int size = 8;
    const T scale = {2, 3};
    const int incx = 1;
    std::array<T, size> result, answer;
    std::generate(result.begin(), result.end(), []() {
        return T{static_cast<float>(std::rand() / float(RAND_MAX)),
                 static_cast<float>(std::rand() / float(RAND_MAX))};
    });
    for (int i = 0; i < size; i++)
        answer[i] = result[i] * scale;
    cscal_(&size, &scale, result.data(), &incx);
    for (int i = 0; i < size; i++) {
        EXPECT_FLOAT_EQ(answer[i].real(), result[i].real());
        EXPECT_FLOAT_EQ(answer[i].imag(), result[i].imag());
    }
}

TEST(blas_connector, zscal_) {
    typedef std::complex<double> T;
    const int size = 8;
    const T scale = {2, 3};
    const int incx = 1;
    std::array<T, size> result, answer;
    std::generate(result.begin(), result.end(), []() {
        return T{static_cast<double>(std::rand() / double(RAND_MAX)),
                 static_cast<double>(std::rand() / double(RAND_MAX))};
    });
    for (int i = 0; i < size; i++)
        answer[i] = result[i] * scale;
    zscal_(&size, &scale, result.data(), &incx);
    for (int i = 0; i < size; i++) {
        EXPECT_DOUBLE_EQ(answer[i].real(), result[i].real());
        EXPECT_DOUBLE_EQ(answer[i].imag(), result[i].imag());
    }
}

TEST(blas_connector, Scal) {
    const int size = 8;
    const std::complex<double> scale = {2, 3};
    const int incx = 1;
    std::complex<double> result[8], answer[8];
    for (int i=0; i< size; i++) {
        result[i] = std::complex<double>{static_cast<double>(std::rand() / double(RAND_MAX)),
                 static_cast<double>(std::rand() / double(RAND_MAX))};
    };
    for (int i = 0; i < size; i++)
        answer[i] = result[i] * scale;
    BlasConnector bs;
    bs.scal(size,scale,result,incx);
    // incx is the spacing between elements if result
    for (int i = 0; i < size; i++) {
        EXPECT_DOUBLE_EQ(answer[i].real(), result[i].real());
        EXPECT_DOUBLE_EQ(answer[i].imag(), result[i].imag());
    }
}


TEST(blas_connector, daxpy_) {
    typedef double T;
    const int size = 8;
    const T scale = 2;
    const int incx = 1;
    const int incy = 1;
    std::array<T, size> x_const, result, answer;
    std::generate(x_const.begin(), x_const.end(),
                  []() { return std::rand() / double(RAND_MAX); });
    std::generate(result.begin(), result.end(),
                  []() { return std::rand() / double(RAND_MAX); });
    for (int i = 0; i < size; i++)
        answer[i] = x_const[i] * scale + result[i];
    daxpy_(&size, &scale, x_const.data(), &incx, result.data(), &incy);
    for (int i = 0; i < size; i++)
        EXPECT_DOUBLE_EQ(answer[i], result[i]);
}

TEST(blas_connector, zaxpy_) {
    typedef std::complex<double> T;
    const int size = 8;
    const T scale = {2, 3};
    const int incx = 1;
    const int incy = 1;
    std::array<T, size> x_const, result, answer;
    std::generate(x_const.begin(), x_const.end(), []() {
        return T{static_cast<double>(std::rand() / double(RAND_MAX)),
                 static_cast<double>(std::rand() / double(RAND_MAX))};
    });
    std::generate(result.begin(), result.end(), []() {
        return T{static_cast<double>(std::rand() / double(RAND_MAX)),
                 static_cast<double>(std::rand() / double(RAND_MAX))};
    });
    for (int i = 0; i < size; i++)
        answer[i] = x_const[i] * scale + result[i];
    zaxpy_(&size, &scale, x_const.data(), &incx, result.data(), &incy);
    for (int i = 0; i < size; i++) {
        EXPECT_DOUBLE_EQ(answer[i].real(), result[i].real());
        EXPECT_DOUBLE_EQ(answer[i].imag(), result[i].imag());
    }
}

TEST(blas_connector, dcopy_) {
    typedef double T;
    long const size = 8;
    int const incx = 1;
    int const incy = 1;
    std::array<T, size> x_const, result, answer;
    std::generate(x_const.begin(), x_const.end(),
                  []() { return std::rand() / double(RAND_MAX); });
    for (int i = 0; i < size; i++)
        answer[i] = x_const[i];
    dcopy_(&size, x_const.data(), &incx, result.data(), &incy);
    for (int i = 0; i < size; i++)
        EXPECT_DOUBLE_EQ(answer[i], result[i]);
}

TEST(blas_connector, zcopy_) {
    typedef std::complex<double> T;
    long const size = 8;
    int const incx = 1;
    int const incy = 1;
    std::array<T, size> x_const, result, answer;
    std::generate(x_const.begin(), x_const.end(), []() {
        return T{static_cast<double>(std::rand() / double(RAND_MAX)),
                 static_cast<double>(std::rand() / double(RAND_MAX))};
    });
    for (int i = 0; i < size; i++)
        answer[i] = x_const[i];
    zcopy_(&size, x_const.data(), &incx, result.data(), &incy);
    for (int i = 0; i < size; i++) {
        EXPECT_DOUBLE_EQ(answer[i].real(), result[i].real());
        EXPECT_DOUBLE_EQ(answer[i].imag(), result[i].imag());
    }
}

TEST(blas_connector, copy) {
    long const size = 8;
    int const incx = 1;
    int const incy = 1;
    std::complex<double> result[8], answer[8];
    for (int i = 0; i < size; i++)
    {
	    answer[i] = std::complex<double>{static_cast<double>(std::rand() / double(RAND_MAX)),
		    static_cast<double>(std::rand() / double(RAND_MAX))};
    }
    BlasConnector bs;
    bs.copy(size, answer, incx, result, incy);
    for (int i = 0; i < size; i++) {
        EXPECT_DOUBLE_EQ(answer[i].real(), result[i].real());
        EXPECT_DOUBLE_EQ(answer[i].imag(), result[i].imag());
    }
}

TEST(blas_connector, dgemv_) {
    typedef double T;
    const char transa_m = 'N';
    const char transa_n = 'T';
    const int size_m = 3;
    const int size_n = 4;
    const T alpha_const = 2;
    const T beta_const = 3;
    const int lda = size_m;
    const int incx = 1;
    const int incy = 1;
    std::array<T, size_m> x_const_m, result_m, answer_m, c_dot_m{};
    std::array<T, size_n> x_const_n, result_n, answer_n, c_dot_n{};
    std::generate(x_const_n.begin(), x_const_n.end(),
                  []() { return std::rand() / double(RAND_MAX); });
    std::generate(result_n.begin(), result_n.end(),
                  []() { return std::rand() / double(RAND_MAX); });
    std::generate(x_const_m.begin(), x_const_m.end(),
                  []() { return std::rand() / double(RAND_MAX); });
    std::generate(result_m.begin(), result_m.end(),
                  []() { return std::rand() / double(RAND_MAX); });
    std::array<T, size_n * lda> a_const;
    std::generate(a_const.begin(), a_const.end(),
                  []() { return std::rand() / double(RAND_MAX); });
    for (int i = 0; i < size_m; i++) {
        for (int j = 0; j < size_n; j++) {
            c_dot_m[i] += a_const[i + j * lda] * x_const_n[j];
        }
        answer_m[i] = alpha_const * c_dot_m[i] + beta_const * result_m[i];
    }
    dgemv_(&transa_m, &size_m, &size_n, &alpha_const, a_const.data(), &lda,
           x_const_n.data(), &incx, &beta_const, result_m.data(), &incy);

    for (int j = 0; j < size_n; j++) {
        for (int i = 0; i < size_m; i++) {
            c_dot_n[j] += a_const[i + j * lda] * x_const_m[i];
        }
        answer_n[j] = alpha_const * c_dot_n[j] + beta_const * result_n[j];
    }
    dgemv_(&transa_n, &size_m, &size_n, &alpha_const, a_const.data(), &lda,
           x_const_m.data(), &incx, &beta_const, result_n.data(), &incy);

    for (int i = 0; i < size_m; i++)
        EXPECT_DOUBLE_EQ(answer_m[i], result_m[i]);
    for (int j = 0; j < size_n; j++)
        EXPECT_DOUBLE_EQ(answer_n[j], result_n[j]);
}

TEST(blas_connector, zgemv_) {
    typedef std::complex<double> T;
    const char transa_m = 'N';
    const char transa_n = 'T';
    const char transa_h = 'C';
    const int size_m = 3;
    const int size_n = 4;
    const T alpha_const = {2, 3};
    const T beta_const = {3, 4};
    const int lda = 5;
    const int incx = 1;
    const int incy = 1;
    std::array<T, size_m> x_const_m, x_const_c, result_m, answer_m, c_dot_m{};
    std::array<T, size_n> x_const_n, result_n, result_c, answer_n, answer_c,
        c_dot_n{}, c_dot_c{};
    std::generate(x_const_n.begin(), x_const_n.end(), []() {
        return T{static_cast<double>(std::rand() / double(RAND_MAX)),
                 static_cast<double>(std::rand() / double(RAND_MAX))};
    });
    std::generate(result_n.begin(), result_n.end(), []() {
        return T{static_cast<double>(std::rand() / double(RAND_MAX)),
                 static_cast<double>(std::rand() / double(RAND_MAX))};
    });
    std::generate(x_const_m.begin(), x_const_m.end(), []() {
        return T{static_cast<double>(std::rand() / double(RAND_MAX)),
                 static_cast<double>(std::rand() / double(RAND_MAX))};
    });
    std::generate(result_m.begin(), result_m.end(), []() {
        return T{static_cast<double>(std::rand() / double(RAND_MAX)),
                 static_cast<double>(std::rand() / double(RAND_MAX))};
    });
    std::array<T, size_n * lda> a_const;
    std::generate(a_const.begin(), a_const.end(), []() {
        return T{static_cast<double>(std::rand() / double(RAND_MAX)),
                 static_cast<double>(std::rand() / double(RAND_MAX))};
    });
    for (int i = 0; i < size_m; i++) {
        for (int j = 0; j < size_n; j++) {
            c_dot_m[i] += a_const[i + j * lda] * x_const_n[j];
        }
        answer_m[i] = alpha_const * c_dot_m[i] + beta_const * result_m[i];
    }
    zgemv_(&transa_m, &size_m, &size_n, &alpha_const, a_const.data(), &lda,
           x_const_n.data(), &incx, &beta_const, result_m.data(), &incy);

    for (int j = 0; j < size_n; j++) {
        for (int i = 0; i < size_m; i++) {
            c_dot_n[j] += a_const[i + j * lda] * x_const_m[i];
        }
        answer_n[j] = alpha_const * c_dot_n[j] + beta_const * result_n[j];
    }
    zgemv_(&transa_n, &size_m, &size_n, &alpha_const, a_const.data(), &lda,
           x_const_m.data(), &incx, &beta_const, result_n.data(), &incy);

    for (int j = 0; j < size_n; j++) {
        for (int i = 0; i < size_m; i++) {
            c_dot_c[j] += conj(a_const[i + j * lda]) * x_const_c[i];
        }
        answer_c[j] = alpha_const * c_dot_c[j] + beta_const * result_c[j];
    }
    zgemv_(&transa_h, &size_m, &size_n, &alpha_const, a_const.data(), &lda,
           x_const_c.data(), &incx, &beta_const, result_c.data(), &incy);

    for (int i = 0; i < size_m; i++) {
        EXPECT_DOUBLE_EQ(answer_m[i].real(), result_m[i].real());
        EXPECT_DOUBLE_EQ(answer_m[i].imag(), result_m[i].imag());
    }
    for (int j = 0; j < size_n; j++) {
        EXPECT_DOUBLE_EQ(answer_n[j].real(), result_n[j].real());
        EXPECT_DOUBLE_EQ(answer_n[j].imag(), result_n[j].imag());
    }
    for (int j = 0; j < size_n; j++) {
        EXPECT_DOUBLE_EQ(answer_c[j].real(), result_c[j].real());
        EXPECT_DOUBLE_EQ(answer_c[j].imag(), result_c[j].imag());
    }
}

TEST(blas_connector, dgemm_) {
    typedef double T;
    const char transa_m = 'N';
    const char transb_m = 'N';
    const int size_m = 3;
    const int size_n = 4;
    const int size_k = 5;
    const T alpha_const = 2;
    const T beta_const = 3;
    const int lda = 6;
    const int ldb = 5;
    const int ldc = 4;
    std::array<T, size_k * lda> a_const;
    std::array<T, size_n * ldb> b_const;
    std::array<T, size_n * ldc> c_dot{}, answer, result;
    std::generate(a_const.begin(), a_const.end(),
                  []() { return std::rand() / double(RAND_MAX); });
    std::generate(b_const.begin(), b_const.end(),
                  []() { return std::rand() / double(RAND_MAX); });
    std::generate(result.begin(), result.end(),
                  []() { return std::rand() / double(RAND_MAX); });
    for (int i = 0; i < size_m; i++) {
        for (int j = 0; j < size_n; j++) {
            for (int k = 0; k < size_k; k++) {
                c_dot[i + j * ldc] +=
                    a_const[i + k * lda] * b_const[k + j * ldb];
            }
            answer[i + j * ldc] = alpha_const * c_dot[i + j * ldc] +
                                  beta_const * result[i + j * ldc];
        }
    }
    dgemm_(&transa_m, &transb_m, &size_m, &size_n, &size_k, &alpha_const,
           a_const.data(), &lda, b_const.data(), &ldb, &beta_const,
           result.data(), &ldc);

    for (int i = 0; i < size_m; i++)
        for (int j = 0; j < size_n; j++) {
            EXPECT_DOUBLE_EQ(answer[i + j * ldc], result[i + j * ldc]);
        }
}

TEST(blas_connector, zgemm_) {
    typedef std::complex<double> T;
    const char transa_m = 'N';
    const char transb_m = 'N';
    const int size_m = 3;
    const int size_n = 4;
    const int size_k = 5;
    const T alpha_const = {2, 3};
    const T beta_const = {3, 4};
    const int lda = 6;
    const int ldb = 5;
    const int ldc = 4;
    std::array<T, size_k * lda> a_const;
    std::array<T, size_n * ldb> b_const;
    std::array<T, size_n * ldc> c_dot{}, answer, result;
    std::generate(a_const.begin(), a_const.end(), []() {
        return T{static_cast<double>(std::rand() / double(RAND_MAX)),
                 static_cast<double>(std::rand() / double(RAND_MAX))};
    });
    std::generate(b_const.begin(), b_const.end(), []() {
        return T{static_cast<double>(std::rand() / double(RAND_MAX)),
                 static_cast<double>(std::rand() / double(RAND_MAX))};
    });
    std::generate(result.begin(), result.end(), []() {
        return T{static_cast<double>(std::rand() / double(RAND_MAX)),
                 static_cast<double>(std::rand() / double(RAND_MAX))};
    });
    for (int i = 0; i < size_m; i++) {
        for (int j = 0; j < size_n; j++) {
            for (int k = 0; k < size_k; k++) {
                c_dot[i + j * ldc] +=
                    a_const[i + k * lda] * b_const[k + j * ldb];
            }
            answer[i + j * ldc] = alpha_const * c_dot[i + j * ldc] +
                                  beta_const * result[i + j * ldc];
        }
    }
    zgemm_(&transa_m, &transb_m, &size_m, &size_n, &size_k, &alpha_const,
           a_const.data(), &lda, b_const.data(), &ldb, &beta_const,
           result.data(), &ldc);

    for (int i = 0; i < size_m; i++)
        for (int j = 0; j < size_n; j++) {
            EXPECT_DOUBLE_EQ(answer[i + j * ldc].real(),
                             result[i + j * ldc].real());
            EXPECT_DOUBLE_EQ(answer[i + j * ldc].imag(),
                             result[i + j * ldc].imag());
        }
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
