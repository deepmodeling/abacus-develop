#include "../line_search.h"

#include <string>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
/************************************************
 *  unit test of class Line_Search
 ***********************************************/

/**
 * - Tested Functions:
 *   - Line_Search::line_search()
 *   - Line_Search::first_order()
 *   - Line_Search::third_order()
 *   - Line_Search::init_brent()
 *   - Line_Search::update_brent()
 *   - Line_Search::brent()
 */

// test of all cases of third_order()
class TestTO : public testing::Test
{
  protected:
    std::vector<double> xnew_arr;

    void SetUp()
    {
        Line_Search ls;

        double x, xnew, y, f;

        // 3rd order, harmonic, dmove > 0
        x = 0;
        y = 0;
        f = -2;
        ls.line_search(1, x, y, f, xnew, 1e-10);
        xnew_arr.push_back(xnew);

        x = 1;
        y = 0;
        f = 2;
        ls.line_search(0, x, y, f, xnew, 1e-10);
        xnew_arr.push_back(xnew);

        // 3rd order, harmonic, dmove < 0
        x = 0;
        y = 0;
        f = 1;
        ls.line_search(1, x, y, f, xnew, 1e-10);
        xnew_arr.push_back(xnew);

        x = 1;
        y = 2;
        f = 3;
        ls.line_search(0, x, y, f, xnew, 1e-10);
        xnew_arr.push_back(xnew);

        // 3rd order, harmonic, dmove > 4
        x = 0;
        y = 0;
        f = -20;
        ls.line_search(1, x, y, f, xnew, 1e-10);
        xnew_arr.push_back(xnew);

        x = 1;
        y = -9;
        f = -18;
        ls.line_search(0, x, y, f, xnew, 1e-10);
        xnew_arr.push_back(xnew);

        // 3rd order, anharmonic, use dmove1 & dmove
        x = 0;
        y = 0;
        f = -1;
        ls.line_search(1, x, y, f, xnew, 1e-10);
        xnew_arr.push_back(xnew);

        x = 1;
        y = 2;
        f = 3.5;
        ls.line_search(0, x, y, f, xnew, 1e-10);
        xnew_arr.push_back(xnew);

        // 3rd order, anharmonic, use dmove2 & dmoveh
        x = 0;
        y = 0;
        f = 4.5;
        ls.line_search(1, x, y, f, xnew, 1e-10);
        xnew_arr.push_back(xnew);

        x = 1;
        y = 2;
        f = -1;
        ls.line_search(0, x, y, f, xnew, 1e-10);
        xnew_arr.push_back(xnew);

        // 3rd order, anharmonic, dmoveh > 4
        x = 0;
        y = 0;
        f = -20;
        ls.line_search(1, x, y, f, xnew, 1e-10);
        xnew_arr.push_back(xnew);

        x = 1;
        y = -5;
        f = -18;
        ls.line_search(0, x, y, f, xnew, 1e-10);
        xnew_arr.push_back(xnew);
    }
};

TEST_F(TestTO, LineSearch)
{
    std::vector<double> xnew_arr_ref = {1, 0.5, 1, 4, 1, 4, 1, 0.1180828963, 1, 0.8181818182, 1, 4};

    for (int i = 0; i < xnew_arr.size(); i++)
    {
        EXPECT_NEAR(xnew_arr[i], xnew_arr_ref[i], 1e-8);
    }
}

class TestLS : public testing::Test
{
  protected:
    std::vector<double> xnew_arr;
    Line_Search ls;

    double x = 0, xnew = 0;

    // here are 150 random numbers as values of x,y,f
    // try to catch as many cases as possible
    std::vector<double> rand
        = {7.84824536,  9.64663789,  3.02069873,  11.26507127, 0.57183349,  13.78382791, 11.35073038, 13.93731114,
           17.79133148, 8.21541546,  14.93677637, 18.98757025, 3.9381696,   16.94282377, 15.3640563,  7.42831671,
           17.3805471,  8.80330788,  0.99259854,  6.1560873,   2.81682224,  8.53026114,  6.09436657,  11.7394269,
           12.90193598, 5.89789934,  8.64853882,  4.73933146,  3.36094493,  14.25637065, 14.60564235, 13.66375342,
           4.4094395,   13.24080597, 6.94100136,  6.74562183,  4.93984237,  3.28646424,  12.95371468, 2.28366581,
           14.10115138, 17.89334777, 17.85486276, 11.66796459, 3.83858509,  17.67856229, 14.09360901, 6.32802805,
           9.81252062,  13.00909712, 18.59789075, 2.20899305,  17.91715464, 17.65862859, 0.49912591,  15.47188694,
           3.74374207,  14.36382216, 3.95503455,  9.02487212,  4.93351326,  8.21473905,  7.10799664,  6.13204507,
           1.85529643,  5.98573821,  8.77475213,  4.84374422,  8.91737412,  12.26443836, 16.89717106, 19.53642954,
           4.69604119,  12.24565644, 14.7098488,  12.98241239, 8.19982363,  19.52610301, 6.18040644,  9.14778002,
           0.3568356,   7.09769476,  19.18409854, 15.39537538, 8.41229496,  12.21441438, 17.09138446, 10.10983535,
           10.93238563, 13.22155721, 8.86644395,  12.70139813, 10.73434842, 11.82537801, 2.14964817,  14.59896678,
           15.41267417, 2.81137988,  3.71229432,  1.83056217,  7.41129002,  1.37004508,  17.76889137, 1.89905894,
           14.16419654, 4.33564906,  9.08449031,  3.62527869,  12.79424987, 18.55812915, 3.74804462,  2.75107002,
           2.60441846,  7.82141378,  4.16003289,  9.68131555,  12.5484824,  9.68999296,  9.7243542,   18.39408735,
           19.00070811, 11.87098669, 3.78418022,  6.1277818,   10.78841329, 5.24128999,  16.13326588, 3.58595401,
           5.21059945,  19.06084549, 0.5703917,   9.66890009,  0.73642789,  2.74387619,  5.27688081,  14.07082298,
           9.87873532,  10.62883574, 10.89594589, 5.06857927,  9.08512489,  8.75285552,  12.71893686, 9.90408215,
           2.96173193,  0.8821268,   2.75896956,  18.05164767, 17.20957441, 8.86361389};
};

TEST_F(TestLS, LineSearch)
{
    int ind = 0;
    for (int i = 0; i < 4; i++)
    {
        bool restart = false;
        if (i % 10 == 0)
            restart = true;

        if (i < 3)
        {
            ls.line_search(restart, rand[ind], rand[ind + 1], rand[ind + 2], xnew, 1e-10);
            ind = ind + 3;
            xnew_arr.push_back(xnew);
        }
        else
        {
            testing::internal::CaptureStdout();
            EXPECT_EXIT(ls.line_search(restart, rand[ind], rand[ind + 1], rand[ind + 2], xnew, 1e-10),
                        ::testing::ExitedWithCode(0),
                        "");
            std::string output = testing::internal::GetCapturedStdout();
            EXPECT_THAT(output, testing::HasSubstr("something wrong with Brent line search!"));
        }
    }
    /*
    std::vector<double> xnew_arr_ref
        = {8.84824536,    11.84824536, 11.5220486,   1.94478562,   -4.61632212,   14.40861093,  -11.8788378,
           23.60558634,   21.64528566, -11.58587758, 15.60564235,  18.60564235,   -11.66208483, -3.02868731,
           22.1076108,    17.32596135, -5.91956272,  -12.99806209, 0.03909330126, 42.09321466,  5.93351326,
           8.93351326,    14.06016625, 19.24381082,  -10.44075315, 29.55515479,   -7.42360546,  8.9322714,
           11.04149536,   13.50491613, 9.86644395,   12.86644395,  16.63592206,   -25.33366183, 49.64554977,
           -0.2852472955, 29.71145149, -17.33528968, 6.97795863,   20.7499131,    20.00070811,  23.00070811,
           36.14423404,   24.91600471, -35.91240731, 40.73961316,  4.54619171,    4.46667478,   -8.62051525,
           48.23147915};
    for (int i = 0; i < xnew_arr.size(); i++)
    {
        EXPECT_NEAR(xnew_arr[i], xnew_arr_ref[i], 1e-8);
    }
    */
}
