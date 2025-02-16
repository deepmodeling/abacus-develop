#include "NCLibxc.h"
#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <complex>
#include <iomanip>
#include <stdexcept>
#include <fstream>


///////////////////////////////////////////////////////////////////////////////////
//Matrix


// Pauli matrix initialization
const Matrix2x2 NCLibxc::sigma_x = {{{std::complex<double>(0.0, 0.0), std::complex<double>(1.0, 0.0)},
                            {std::complex<double>(1.0, 0.0), std::complex<double>(0.0, 0.0)}}};

const Matrix2x2 NCLibxc::sigma_y = {{{std::complex<double>(0.0, 0.0), std::complex<double>(0.0, -1.0)},
                            {std::complex<double>(0.0, 1.0), std::complex<double>(0.0, 0.0)}}};

const Matrix2x2 NCLibxc::sigma_z = {{{std::complex<double>(1.0, 0.0), std::complex<double>(0.0, 0.0)},
                            {std::complex<double>(0.0, 0.0), std::complex<double>(-1.0, 0.0)}}};

// matrix addition
Matrix2x2 NCLibxc::add(const Matrix2x2 &a, const Matrix2x2 &b)
{
    Matrix2x2 result;
    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            result[i][j] = a[i][j] + b[i][j];
        }
    }
    return result;
}

// matrix scalar multiplication
Matrix2x2 NCLibxc::scalar_multiply(const std::complex<double> &scalar, const Matrix2x2 &matrix)
{
    Matrix2x2 result;
    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            result[i][j] = scalar * matrix[i][j];
        }
    }
    return result;
}


// construct a Pauli matrix based on x, y, z
Matrix2x2 NCLibxc::construct_pauli_matrix(double x, double y, double z)
{
    Matrix2x2 pauli_matrix = add(add(scalar_multiply(x, sigma_x), scalar_multiply(y, sigma_y)), scalar_multiply(z, sigma_z));
    return pauli_matrix;
}

// identity matrix
Matrix2x2 NCLibxc::identity_matrix()
{
    Matrix2x2 result;
    result[0][0] = 1.0;
    result[0][1] = 0.0;
    result[1][0] = 0.0;
    result[1][1] = 1.0;
    return result;
}
///////////////////////////////////////////////////////////////////////////////////