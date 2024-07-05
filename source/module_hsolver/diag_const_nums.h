#ifndef DIAG_CONST_NUMS
#define DIAG_CONST_NUMS

template <typename T>
struct const_nums
{
    const_nums();
    T zero;
    T one;
    T neg_one;
};

#endif


// #ifndef DIAG_CONST_NUMS_H
// #define DIAG_CONST_NUMS_H

// #include <complex>

// template <typename T>
// struct const_nums
// {
//     T zero;
//     T one;
//     T neg_one;
//     const_nums();
// };

// // 特化模板以支持 double 类型
// template <>
// const_nums<double>::const_nums() : zero(0.0), one(1.0), neg_one(-1.0) {}

// // 特化模板以支持 std::complex<double> 类型
// template <>
// const_nums<std::complex<double>>::const_nums() : zero(std::complex<double>(0.0, 0.0)), one(std::complex<double>(1.0, 0.0)), neg_one(std::complex<double>(-1.0, 0.0)) {}

// // 特化模板以支持 std::complex<float> 类型
// template <>
// const_nums<std::complex<float>>::const_nums() : zero(std::complex<float>(0.0, 0.0)), one(std::complex<float>(1.0, 0.0)), neg_one(std::complex<float>(-1.0, 0.0)) {}

// #endif // DIAG_CONST_NUMS_H