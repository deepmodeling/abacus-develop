#include "basic_funcs.h"

#include <iostream>

double maxval_abs_2d(const std::vector<std::vector<double>>& array)
{
    double max = 0;
    for (const auto& row: array)
    {
        for (double value: row)
        {
            max = std::max(max, std::abs(value));
        }
    }
    return max;
}

void maxloc_abs_2d(const std::vector<std::vector<double>>& array, std::vector<int>& result)
{
    double max = 0;
    int size_1 = array.size();
    int size_2 = array[0].size();
    for (int i = 0; i < size_1; i++)
    {
        for (int j = 0; j < size_2; j++)
        {
            if ((max < abs(array[i][j])))
            {
                max = abs(array[i][j]);
                result[0] = i;
                result[1] = j;
            }
        }
    }
}

double sum_2d(const std::vector<std::vector<double>>& array)
{
    double sum = 0;
    for (const auto& row: array)
    {
        for (const auto& element: row)
        {
            sum += element;
        }
    }
    return sum;
}

void scalar_multiply_2d(const std::vector<std::vector<double>>& array,
                        double scalar,
                        std::vector<std::vector<double>>& result)
{
    int size_1 = array.size();
    int size_2 = array[0].size();
    result.resize(size_1);
    for (int i = 0; i < size_1; i++)
    {
        result[i].resize(size_2);
        for (int j = 0; j < size_2; j++)
        {
            result[i][j] = scalar * array[i][j];
        }
    }
}

void add_scalar_multiply_2d(const std::vector<std::vector<double>>& array_1,
                            const std::vector<std::vector<double>>& array_2,
                            double scalar,
                            std::vector<std::vector<double>>& result)
{
    int size_1 = array_1.size();
    int size_2 = array_1[0].size();
    result.resize(size_1);
    for (int i = 0; i < size_1; i++)
    {
        result[i].resize(size_2);
        for (int j = 0; j < size_2; j++)
        {
            result[i][j] = array_1[i][j] + scalar * array_2[i][j];
        }
    }
}

void subtract_2d(const std::vector<std::vector<double>>& array_1,
                 const std::vector<std::vector<double>>& array_2,
                 std::vector<std::vector<double>>& result)
{
    int size_1 = array_1.size();
    int size_2 = array_1[0].size();
    result.resize(size_1);
    for (int i = 0; i < size_1; i++)
    {
        result[i].resize(size_2);
        for (int j = 0; j < size_2; j++)
        {
            result[i][j] = array_1[i][j] - array_2[i][j];
        }
    }
}

void fill_scalar_2d(double scalar, std::vector<std::vector<double>>& result)
{
    for (auto& row: result)
    {
        std::fill(row.begin(), row.end(), scalar);
    }
}

void where_fill_scalar_2d(const std::vector<std::vector<int>>& array_mask,
                          int mask,
                          double scalar,
                          std::vector<std::vector<double>>& result)
{
    int size_1 = array_mask.size();
    int size_2 = array_mask[0].size();
    result.resize(size_1);
    for (int i = 0; i < size_1; i++)
    {
        result[i].resize(size_2);
        for (int j = 0; j < size_2; j++)
        {
            if (array_mask[i][j] == mask)
            {
                result[i][j] = scalar;
            }
        }
    }
}

void where_fill_scalar_else_2d(const std::vector<std::vector<int>>& array_mask,
                               int mask,
                               double scalar,
                               const std::vector<std::vector<double>>& rest,
                               std::vector<std::vector<double>>& result)
{
    int size_1 = array_mask.size();
    int size_2 = array_mask[0].size();
    result.resize(size_1);
    for (int i = 0; i < size_1; i++)
    {
        result[i].resize(size_2);
        for (int j = 0; j < size_2; j++)
        {
            result[i][j] = (array_mask[i][j] == mask) ? scalar : rest[i][j];
        }
    }
}

void print_2d(const std::vector<std::vector<double>>& array)
{
    for (const auto& row: array)
    {
        for (const auto& element: row)
        {
            std::cout << element << " ";
        }
        std::cout << std::endl;
    }
}

void optimize_lambda(const std::vector<std::vector<double>> &M_CONSTR,
                     const std::vector<std::vector<int>> &CONSTRL,
                     const int NIONS, const int NTYP, const std::vector<int> &NITYP,
                     const double INISC,
                     const double SCDIFF,
                     const std::vector<double> &SCCONV_GRAD,
                     const int NSC,
                     const int NSCMIN,
                     const double SCCUT,
                     const int N,
                     std::vector<std::vector<double>> &MW,
                     std::vector<std::vector<double>> &OUT_LAMBDA
                     //  TODO datatype
                     //<> CHTOT,
                     //<> CHTOTL,
                     //<> W)

{
    // TODO datatype
    //<> CHTOT_RESERVE;
    //<> CHTOTL_RESERVE;
    //<> W_RESERVE;

    double epsilon;
    double alpha_trial, alpha_opt, alpha_plus;
    double beta;
    double mean_error, mean_error_old, rms_error;
    double restrict_current;
    double boundary;
    double sum_k, sum_k2;
    double g;

    int num_component, num_step, i_step;

    //    std::vector<std::vector<double>> INITIAL_LAMBDA(3, std::vector<double>(NIONS));
    //    std::vector<std::vector<double>> DELTA_LAMBDA(3, std::vector<double>(NIONS));
    //    std::vector<std::vector<double>> nu(3, std::vector<double>(NIONS));
    //    std::vector<std::vector<double>> dnu(3, std::vector<double>(NIONS));
    //    std::vector<std::vector<double>> dnu_last_step(3, std::vector<double>(NIONS));
    //    std::vector<std::vector<double>> target_spin(3, std::vector<double>(NIONS));
    //    std::vector<std::vector<double>> spin(3, std::vector<double>(NIONS));
    //    std::vector<std::vector<double>> spin_plus(3, std::vector<double>(NIONS));
    //    std::vector<std::vector<double>> delta_spin(3, std::vector<double>(NIONS));
    //    std::vector<std::vector<double>> delta_spin_old(3, std::vector<double>(NIONS));
    //    std::vector<std::vector<double>> search(3, std::vector<double>(NIONS));
    //    std::vector<std::vector<double>> search_old(3, std::vector<double>(NIONS));
    //    std::vector<std::vector<double>> target_spin_mask(3, std::vector<double>(NIONS));
    //    std::vector<std::vector<double>> spin_mask(3, std::vector<double>(NIONS));
    //    std::vector<std::vector<double>> spin_plus_mask(3, std::vector<double>(NIONS));
    //    std::vector<std::vector<std::vector<std::vector<double>>>> spin_nu_gradient(
    //        3, std::vector<std::vector<std::vector<double>>>(
    //               NIONS, std::vector<std::vector<double>>(
    //                          3, std::vector<double>(NIONS))));
    //    std::vector<std::vector<double>> spin_nu_gradient_diag(3, std::vector<double>(NIONS));
    //    std::vector<std::vector<double>> spin_change(3, std::vector<double>(NIONS));
    //    std::vector<std::vector<double>> nu_change(3, std::vector<double>(NIONS));
    //    std::vector<std::vector<int>> max_gradient_index(2, std::vector<int>(NTYP));
    //    std::vector<double> max_gradient(NTYP);
    //    std::vector<double> bound_gradient(NTYP);
    //
    //    std::vector<std::vector<double>> temp_1(3, std::vector<double>(NIONS));
    //    std::vector<std::vector<double>> temp_2(3, std::vector<double>(NIONS));
    //    std::vector<int> temp_3(2);
    //
    //    target_spin = M_CONSTR;
    //    fill_scalar_2d(0.0, dnu);
    //    fill_scalar_2d(0.0, dnu_last_step);
    //    num_component = sum_2d(CONSTRL);
    //    epsilon = SCDIFF;
    //    bound_gradient = SCCONV_GRAD;
    //    alpha_trial = INISC;
    //    num_step = NSC;
    //    restrict_current = SCCUT;
    //
    //    // TODO
    //    CHTOT_RESERVE = CHTOT;
    //    CHTOTL_RESERVE = CHTOTL;
    //    W_RESERVE = W;
    //
    //    std::cout << "===============================================================================";
    //    << std::endl;
    //    std::cout << "Inner optimization for lambda begins ..." << std::endl;
    //    std::cout << "Covergence criterion for the iteration: " << epsilon << std::endl;
    //
    //    // lambda loop
    //    for (int i_step = 0; i_step < num_step; i_step++)
    //    {
    //        if (i_step == 0)
    //        {
    //            nu = OUT_LAMBDA;
    //            where_fill_scalar_else_2d(CONSTRL, 0, 0.0, OUT_LAMBDA, INITIAL_LAMBDA);
    //            spin = MW;
    //            std::cout << "initial lambda:" << std::endl;
    //            print_2d(INITIAL_LAMBDA);
    //            std::cout << "initial spin: " << std::endl;
    //            print_2d(spin);
    //            std::cout << "target spin: " << std::endl;
    //            print_2d(target_spin);
    //        }
    //        else
    //        {
    //            std::cout << "optimal delta lambda: " << std::endl;
    //            print_2d(DELTA_LAMBDA);
    //
    //            add_scalar_multiply_2d(INITIAL_LAMBDA, DELTA_LAMBDA, 1.0, temp_1);
    //            // TODO, also in-place change density CHTOT and orbital W, const 3 regular 3
    //
    //            // basically, use CHTOTL_RESERVE and temp_1(LAMBDA) recalculate V, then calculate H, then
    //            // diagonalize H (in the subspace spanned by W_RESERVE), then use obtained orbitals W
    //            // calculate density CHTOT and orbital mag MW.
    //
    //            // Note that using CHTOTL instead of CHTOT is to recreate the H that has W_RESERVE as eigenvectors
    //            calculate_MW_from_lambda(temp_1, CHTOTL_RESERVE, W_RESERVE, MW, CHTOT, W);
    //
    //            subtract_2d(MW, spin, spin_change);
    //            subtract_2d(DELTA_LAMBDA, dnu_last_step, nu_change);
    //            where_fill_scalar_2d(CONSTRL, 0, 0.0, spin_change);
    //            where_fill_scalar_2d(CONSTRL, 0, 0.0, nu_change);
    //
    //            for (int i = 0; i < 3; i++)
    //            {
    //                for (int j = 0; j < NIONS; j++)
    //                {
    //                    for (int k = 0; k < 3; k++)
    //                    {
    //                        for (int l = 0; l < NIONS; l++)
    //                        {
    //                            spin_nu_gradient[i][j][k][l] = spin_change[i][j] / nu_change[k][l]
    //                        }
    //                    }
    //                }
    //            }
    //            for (int i = 0; i < 3; i++)
    //            {
    //                for (int j = 0; j < NIONS; j++)
    //                {
    //                    spin_nu_gradient_diag[i][j] = spin_nu_gradient[i][j][i][j]
    //                }
    //            }
    //            int i = 0;
    //            for (int j = 0; j < NTYP; j++)
    //            {
    //                for (int k = 0; k < 3; k++)
    //                {
    //                    for (int l = 0; l < NIONS; l++)
    //                    {
    //                        if (i <= l && l < (i + NITYP[j]))
    //                        {
    //                            temp_1[k][l] = spin_nu_gradient_diag[k][l]
    //                        }
    //                        else
    //                        {
    //                            temp_1[k][l] = 0.0;
    //                        }
    //                    }
    //                }
    //                maxloc_abs_2d(temp_1, temp_3);
    //                for (int k = 0; k < 2; k++)
    //                {
    //                    for (int l = 0; l < NTYP; l++)
    //                    {
    //                        if (l == j)
    //                        {
    //                            max_gradient_index[k][l] = temp_3[k];
    //                        }
    //                    }
    //                }
    //                max_gradient[j] = maxval_abs_2d(temp_1);
    //                i = i + NITYP[j];
    //            }
    //
    //            std::cout << "diagonal gradient:" << std::endl;
    //            print_2d(spin_nu_gradient_diag);
    //            std::cout << "maximum gradient appears at: " << std::endl;
    //            for (int i = 0; i < NTYP; i++)
    //            {
    //                std::cout << "(" << max_gradient_index[0][i] << "," << max_gradient_index[1][i] << ")" <<
    //                std::endl;
    //            }
    //            std::cout << "maximum gradient: " << std::endl;
    //            for (int i = 0; i < NTYP; i++)
    //            {
    //                std::cout << max_gradient[i] << std::endl;
    //            }
    //
    //            for (int i = 0; i < NTYP; i++)
    //            {
    //                if (i_step >= NSCMIN && bound_gradient[i] > 0 && max_gradient[i] < bound_gradient[i])
    //                {
    //                    std::cout << "Reach limitation of current step ( maximum gradient < " << bound_gradient[i] <<
    //                    " in atom type " //<< i << " ), exit." << std::endl;
    //                    // roll back to the last step
    //                    // TODO
    //                    CHTOT = CHTOT_last_step;
    //                    add_scalar_multiply_2d(INITIAL_LAMBDA, dnu_last_step, 1.0, OUT_LAMBDA);
    //                    goto cg_stop;
    //                }
    //            }
    //
    //            spin = MW;
    //            std::cout << "current spin:" << std::endl;
    //            print_2d(spin);
    //            std::cout << "target spin: " << std::endl;
    //            print_2d(target_spin);
    //        }
    //
    //        subtract_2d(spin, target_spin, delta_spin);
    //        where_fill_scalar_2d(CONSTRL, 0, 0.0, delta_spin);
    //        search = delta_spin;
    //        for (int i = 0; i < 3; i++)
    //        {
    //            for (int j = 0; j < NIONS; j++)
    //            {
    //                temp_1[i][j] = pow(delta_spin[i][j], 2);
    //            }
    //        }
    //        mean_error = sum_2d(temp_1) / num_component;
    //        rms_error = sqrt(mean_error);
    //
    //        std::cout << "Step (Outer -- Inner) =  " << N << " -- " << i_step + 1 << "       RMS =", rms_error <<
    //        std::endl;
    //
    //        if (rms_error < epsilon || i_step == num_step - 1)
    //        {
    //            if (rms_error < epsilon)
    //            {
    //                std::cout << "Meet convergence criterion ( < " << epsilon << " ), exit." << std::endl;
    //            }
    //            else if (i_step == num_step - 1)
    //            {
    //                std::cout << "Reach maximum number of steps ( " << num_step << " ), exit." << std::endl;
    //            }
    //            else
    //            {
    //                ;
    //            }
    //            add_scalar_multiply_2d(INITIAL_LAMBDA, DELTA_LAMBDA, 1.0, OUT_LAMBDA);
    //            goto cg_stop;
    //        }
    //
    //        if (i_step >= 1)
    //        {
    //            beta = mean_error / mean_error_old;
    //            temp_1 = search;
    //            add_scalar_multiply_2d(temp_1, search_old, beta, search);
    //        }
    //
    //        boundary = abs(alpha_trial * maxval_abs_2d(search));
    //        std::cout << "restriction of this step = " << restrict_current << std::endl;
    //        std::cout << "alpha_trial before restrict = " << alpha_trial << std::endl;
    //        std::cout << "boundary before = ", boundary << std::endl;
    //        std::cout << "trial need restriction: false" << std::endl;
    //        std::cout << "delta delta lambda:" << std::endl;
    //        temp_1 = scalar_multiply_2d(search, alpha_trial, temp_1);
    //        print_2d(temp_1);
    //
    //        CHTOT_last_step = CHTOT;
    //        dnu_last_step = dnu;
    //        temp_1 = dnu;
    //        add_scalar_multiply_2d(temp_1, search, alpha_trial, dnu);
    //        DELTA_LAMBDA = dnu;
    //
    //        std::cout << "trial delta lambda:" << std::endl;
    //        print_2d(DELTA_LAMBDA);
    //
    //        if (LDESC == true)
    //        {
    //            std::cout << "(Debug) before-trial-step spin:" << std::endl;
    //            print_2d(MW);
    //            std::cout << "(Debug) target spin:" << std::endl;
    //            print_2d(M_CONSTR);
    //        }
    //
    //        add_scalar_multiply_2d(INITIAL_LAMBDA, DELTA_LAMBDA, 1.0, temp_1);
    //        // TODO
    //        calculate_MW_from_lambda(temp_1, CHTOTL_RESERVE, W_RESERVE, MW, CHTOT, W);
    //
    //        spin_plus = MW;
    //        std::cout << "current spin(trial):" << std::endl;
    //        print_2d(spin_plus);
    //
    //        where_fill_scalar_else_2d(CONSTRL, 0, 0.0, target_spin, target_spin_mask);
    //        where_fill_scalar_else_2d(CONSTRL, 0, 0.0, spin, spin_mask);
    //        where_fill_scalar_else_2d(CONSTRL, 0, 0.0, spin_plus, spin_plus_mask);
    //        for (int i = 0; i < 3; i++)
    //        {
    //            for (int j = 0; j < NIONS; j++)
    //            {
    //                temp_1[i][j] = (target_spin_mask[i][j] - spin_mask[i][j]) * (spin_plus_mask[i][j] -
    //                spin_mask[i][j]); temp_2[i][j] = pow(spin_mask[i][j] - spin_plus_mask[i][j], 2)
    //            }
    //        }
    //        sum_k = sum_2d(temp_1);
    //        sum_k2 = sum_2d(temp_2);
    //        alpha_opt = sum_k * alpha_trial / sum_k2;
    //
    //        boundary = abs(alpha_trial * maxval_abs_2d(search));
    //        std::cout << "alpha_opt before restrict = " << alpha_opt << std::endl;
    //        std::cout << "boundary before = ", boundary << std::endl;
    //
    //        if (SCCUT > 0 && boundary > restrict_current)
    //        {
    //            alpha_opt = copysign(1.0, alpha_opt) * restrict_current / maxval_abs_2d(search);
    //            boundary = abs(alpha_opt * maxval_abs_2d(search));
    //            std::cout << "restriction needed: true" << std::endl;
    //            std::cout << "alpha_opt after restrict = " << alpha_opt << std::endl;
    //            std::cout << "boundary after = ", boundary << std::endl;
    //        }
    //        else
    //        {
    //            std::cout << "restriction needed: false" << std::endl;
    //        }
    //
    //        alpha_plus = alpha_opt - alpha_trial;
    //        std::cout << "delta delta lambda:" << std::endl;
    //        scalar_multiply_2d(search, alpha_plus, temp_1);
    //        print_2d(temp_1);
    //
    //        temp_2 = dnu;
    //        add_scalar_multiply_2d(temp_2, temp_1, 1.0, dnu);
    //        DELTA_LAMBDA = dnu;
    //
    //        search_old = search;
    //        delta_spin_old = delta_spin;
    //        mean_error_old = mean_error;
    //
    //        g = 1.5 * abs(alpha_opt) / alpha_trial;
    //        if (g > 2.0)
    //        {
    //            g = 2;
    //        }
    //        else if (g < 0.5)
    //        {
    //            g = 0.5
    //        }
    //        else
    //        {
    //            ;
    //        }
    //        alpha_trial = alpha_trial * pow(g, 0.7);
    //    }
    //
    // cg_stop:;
    //
    //    // TODO
    //    CHTOTL = CHTOTL_RESERVE;
    //
    //    if (LDESC == true)
    //    {
    //        std::cout << "(Debug) after-optimization spin: (print in the inner loop):" << std::endl;
    //        print_2d(MW);
    //        std::cout << "(Debug) target spin:" << std::endl;
    //        print_2d(M_CONSTR);
    //    }
    //
    //    std::cout << "Inner optimization for lambda ends." << std::endl;
    //    std::cout << "===============================================================================" << std::endl;
}