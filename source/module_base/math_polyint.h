#ifndef MATH_POLYINT_H
#define MATH_POLYINT_H

#include "realarray.h"

namespace ModuleBase
{

// mohan add 2021-05-07
class PolyInt
{

	public:

	PolyInt();
	~PolyInt();

    //========================================================
    // Polynomial_Interpolation
    //========================================================

    /**
     * @brief Lagrange interpolation
     * 
     * @param table [in] three dimension matrix, the data in 3rd dimension is used to do prediction
     * @param dim1 [in] index of 1st dimension of table/y
     * @param dim2 [in] index of 2nd dimension of table/y
     * @param y [out] three dimension matrix to store the predicted value
     * @param dim_y [in] index of 3rd dimension of y to store predicted value
     * @param table_length [in] length of 3rd dimension of table
     * @param table_interval [in] interval of 3rd dimension of table
     * @param x [in] the position in 3rd dimension to be predicted
     */
    static void Polynomial_Interpolation
    (
        const ModuleBase::realArray &table,
        const int &dim1,
        const int &dim2,
        ModuleBase::realArray &y,
        const int &dim_y,
        const int &table_length,
        const double &table_interval,
        const double &x
    );

    /**
     * @brief Lagrange interpolation
     * 
     * @param table [in] three dimension matrix, the data in 3rd dimension is used to do prediction
     * @param dim1 [in] index of 1st dimension of table
     * @param dim2 [in] index of 2nd dimension of table
     * @param table_length [in] length of 3rd dimension of table
     * @param table_interval [in] interval of 3rd dimension of table
     * @param x [in] the position in 3rd dimension to be predicted
     * @return double the predicted value
     */
    static double Polynomial_Interpolation
    (
        const ModuleBase::realArray &table,
        const int &dim1,
        const int &dim2,
        const int &table_length,
        const double &table_interval,
        const double &x            
    );

    /**
     * @brief Lagrange interpolation
     * 
     * @param table [in] four dimension matrix, the data in 4th dimension is used to do prediction
     * @param dim1 [in] index of 1st dimension of table
     * @param dim2 [in] index of 2nd dimension of table 
     * @param dim3 [in] index of 3rd dimension of table 
     * @param table_length [in] length of 4th dimension of table
     * @param table_interval [in] interval of 4th dimension of table
     * @param x [in] the position in 4th dimension to be predicted
     * @return double the predicted value
     * @author pengfei Li
     * @date 2018-3-23
     */
    static double Polynomial_Interpolation            
    (
        const ModuleBase::realArray &table,
        const int &dim1,
        const int &dim2,
        const int &dim3,
        const int &table_length,
        const double &table_interval,
        const double &x            
    );

    /**
     * @brief  Lagrange interpolation
     * 
     * @param table [in] the data used to do prediction
     * @param table_length [in] length of table
     * @param table_interval [in] interval of table
     * @param x [in] the position to be predicted
     * @return double the predicted value
     */
	static double Polynomial_Interpolation
	(
        const double *table,
        const int &table_length,
        const double &table_interval,
        const double &x            
    );

    /**
     * @brief Lagrange interpolation
     * 
     * @param xpoint [in] array of postion
     * @param ypoint [in] array of data to do prediction
     * @param table_length [in] length of xpoint
     * @param x [in] position to be predicted
     * @return double predicted value
     */
    static double Polynomial_Interpolation_xy
    (
        const double *xpoint,
        const double *ypoint,
        const int table_length,
        const double &x            
    );

};
}
#endif
