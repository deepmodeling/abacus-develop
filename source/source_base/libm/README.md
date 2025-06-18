# Math fucntions source code

These math functions is ported form glibc-2.36/math (https://github.com/bminor/glibc/tree/release/2.36/master/math)

The new version of sin/cos/exp implementation is much faster than glibc-2.27.

In order to avoid the performance impact of different versions of glibc installed on systems and different versions of math functions provided by compilers, we choose to migrate these source codes.
