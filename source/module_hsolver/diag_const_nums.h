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