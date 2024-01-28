#include"../math_ylmreal.h"
#include"../ylm.h"
#include"../vector3.h"
#include"../matrix.h"
#include"gtest/gtest.h"
#include<math.h>
#include "module_psi/psi.h"

#define doublethreshold 1e-12

/************************************************
*  unit test of class YlmReal and Ylm
***********************************************/

/**
 * For lmax <5 cases, the reference values are calculated by the formula from 
 * https://formulasearchengine.com/wiki/Table_of_spherical_harmonics. Note, these
 * formula lack of the Condon–Shortley phase (-1)^m, and in this unit test, item 
 * (-1)^m is multiplied.
 * For lmax >=5, the reference values are calculated by YlmReal::Ylm_Real.
 *
 * - Tested functions of class YlmReal
 *      - Ylm_Real
 *      - Ylm_Real2
 *      - rlylm 
 *      - YlmRealTemplate (double and float)
 *
 * - Tested functions of class Ylm
 *      - get_ylm_real
 *      - sph_harm
 *      - rl_sph_harm
 *      - grad_rl_sph_harm
 *      - equality_value_test: test the eqaulity of Ylm function between rl_sph_harm (spherical input) and  get_ylm_real (Cartesian input) 
 *      - equality_gradient_test:test the eqaulity of Ylm gradient function between grad_rl_sph_harm(spherical input) and  rlylm (Cartesian input)
 * 
 */



//mock functions of WARNING_QUIT and WARNING
namespace ModuleBase
{
    void WARNING_QUIT(const std::string &file,const std::string &description) {return ;}
    void WARNING(const std::string &file,const std::string &description) {return ;}
}


class YlmRealTest : public testing::Test
{
    protected:

    int lmax = 7;   //maximum angular quantum number
    int ng = 4;     //test the 4 selected points on the sphere
    int nylm = 64;      //total Ylm number;

    ModuleBase::matrix ylm;         //Ylm
    ModuleBase::matrix *dylm;       //dYlm/dx, dYlm/dy, dYlm/dz
    ModuleBase::Vector3<double> *g; //vectors of the 4 points
    double *ref;        //reference of Ylm
    double *rly;        //Ylm
    double (*rlgy)[3];  //the gradient of Ylm  
    std::vector<double> rlyvector; //Ylm
    std::vector<std::vector<double>> rlgyvector; //the gradient of Ylm

    //Ylm function
    inline double norm(const double &x, const double &y, const double &z) {return sqrt(x*x + y*y + z*z);}
    double y00(const double &x, const double &y, const double &z)  {return 1.0/2.0/sqrt(M_PI);}
    double y10(const double &x, const double &y, const double &z)  {double r=norm(x,y,z); return sqrt(3.0/(4.0*M_PI)) * z / r;}
    double y11(const double &x, const double &y, const double &z)  {double r=norm(x,y,z); return -1.0*sqrt(3.0/(4.*M_PI)) * x / r;}
    double y1m1(const double &x, const double &y, const double &z) {double r=norm(x,y,z); return -1.0*sqrt(3./(4.*M_PI)) * y / r;} // y1m1 means Y1,-1
    double y20(const double &x, const double &y, const double &z)  {double r=norm(x,y,z); return 1./4. * sqrt(5./M_PI) * (-1.*x*x - y*y + 2.*z*z) / (r*r);}
    double y21(const double &x, const double &y, const double &z)  {double r=norm(x,y,z); return -1.0*1./2. * sqrt(15./M_PI) * (z*x) / (r*r);}
    double y2m1(const double &x, const double &y, const double &z) {double r=norm(x,y,z); return -1.0*1./2. * sqrt(15./M_PI) * (z*y) / (r*r);}
    double y22(const double &x, const double &y, const double &z)  {double r=norm(x,y,z); return 1./4. * sqrt(15./M_PI) * (x*x - y*y) / (r*r);}
    double y2m2(const double &x, const double &y, const double &z) {double r=norm(x,y,z); return 1./2. * sqrt(15./M_PI) * (x*y) / (r*r);}
    double y30(const double &x, const double &y, const double &z)  {double r=norm(x,y,z); return 1./4. * sqrt(7./M_PI) * z*(2.*z*z-3.*x*x-3.*y*y) / (r*r*r);}
    double y31(const double &x, const double &y, const double &z)  {double r=norm(x,y,z); return -1.0*1./4. * sqrt(21./2./M_PI) * x*(4.*z*z-x*x-y*y) / (r*r*r);}
    double y3m1(const double &x, const double &y, const double &z) {double r=norm(x,y,z); return -1.0*1./4. * sqrt(21./2./M_PI) * y*(4.*z*z-x*x-y*y) / (r*r*r);}
    double y32(const double &x, const double &y, const double &z)  {double r=norm(x,y,z); return 1./4. * sqrt(105./M_PI) * (x*x - y*y)*z / (r*r*r);}
    double y3m2(const double &x, const double &y, const double &z) {double r=norm(x,y,z); return 1./2. * sqrt(105./M_PI) * x*y*z / (r*r*r);}
    double y33(const double &x, const double &y, const double &z)  {double r=norm(x,y,z); return -1.0*1./4. * sqrt(35./2./M_PI) * x*(x*x - 3.*y*y) / (r*r*r);}
    double y3m3(const double &x, const double &y, const double &z) {double r=norm(x,y,z); return -1.0*1./4. * sqrt(35./2./M_PI) * y*(3.*x*x - y*y) / (r*r*r);}
    double y40(const double &x, const double &y, const double &z)  {double r=norm(x,y,z); return 3./16.*sqrt(1./M_PI) * (35.*z*z*z*z - 30.*z*z*r*r + 3*r*r*r*r) / (r*r*r*r);}
    double y41(const double &x, const double &y, const double &z)  {double r=norm(x,y,z); return -1.0*3./4.*sqrt(5./2./M_PI) * x*z*(7.*z*z - 3*r*r) / (r*r*r*r);}
    double y4m1(const double &x, const double &y, const double &z) {double r=norm(x,y,z); return -1.0*3./4.*sqrt(5./2./M_PI) * y*z*(7.*z*z - 3.*r*r) / (r*r*r*r);}
    double y42(const double &x, const double &y, const double &z)  {double r=norm(x,y,z); return 3./8.*sqrt(5./M_PI) * (x*x-y*y)*(7.*z*z-r*r) / (r*r*r*r);}
    double y4m2(const double &x, const double &y, const double &z) {double r=norm(x,y,z); return 3./4.*sqrt(5./M_PI) * x*y*(7.*z*z - r*r) / (r*r*r*r);}
    double y43(const double &x, const double &y, const double &z)  {double r=norm(x,y,z); return -1.0*3./4.*sqrt(35./2./M_PI) * x*z*(x*x - 3.*y*y) / (r*r*r*r);}
    double y4m3(const double &x, const double &y, const double &z) {double r=norm(x,y,z); return -1.0*3./4.*sqrt(35./2./M_PI) * y*z*(3.*x*x - y*y) / (r*r*r*r);}
    double y44(const double &x, const double &y, const double &z)  {double r=norm(x,y,z); return 3./16.*sqrt(35./M_PI) * (x*x*(x*x - 3.*y*y) - y*y*(3.*x*x-y*y)) / (r*r*r*r);}
    double y4m4(const double &x, const double &y, const double &z) {double r=norm(x,y,z); return 3./4.*sqrt(35./M_PI) * x*y*(x*x - y*y) / (r*r*r*r);}

    //the reference values are calculated by ModuleBase::Ylm::grad_rl_sph_harm
    //1st dimension: example, 2nd dimension: Ylm, 3rd dimension: dx/dy/dz 
    double rlgyref[4][64][3] = {
        {   { 0.00000e+00,  0.00000e+00,  0.00000e+00}, { 0.00000e+00,  0.00000e+00,  4.88603e-01}, {-4.88603e-01,  0.00000e+00,  0.00000e+00}, 
            { 0.00000e+00, -4.88603e-01,  0.00000e+00}, {-6.30783e-01,  0.00000e+00,  0.00000e+00}, {-0.00000e+00,  0.00000e+00, -1.09255e+00}, 
            { 0.00000e+00, -0.00000e+00,  0.00000e+00}, { 1.09255e+00,  0.00000e+00,  0.00000e+00}, {-0.00000e+00,  1.09255e+00, -0.00000e+00}, 
            {-0.00000e+00,  0.00000e+00, -1.11953e+00}, { 1.37114e+00,  0.00000e+00, -0.00000e+00}, { 0.00000e+00,  4.57046e-01,  0.00000e+00}, 
            { 0.00000e+00,  0.00000e+00,  1.44531e+00}, {-0.00000e+00,  0.00000e+00,  0.00000e+00}, {-1.77013e+00,  0.00000e+00, -0.00000e+00}, 
            { 0.00000e+00, -1.77013e+00,  0.00000e+00}, { 1.26943e+00,  0.00000e+00, -0.00000e+00}, { 0.00000e+00,  0.00000e+00,  2.00714e+00}, 
            { 0.00000e+00,  0.00000e+00,  0.00000e+00}, {-1.89235e+00,  0.00000e+00,  0.00000e+00}, {-0.00000e+00, -9.46175e-01,  0.00000e+00}, 
            {-0.00000e+00,  0.00000e+00, -1.77013e+00}, { 0.00000e+00, -0.00000e+00,  0.00000e+00}, { 2.50334e+00,  0.00000e+00,  0.00000e+00}, 
            {-0.00000e+00,  2.50334e+00,  0.00000e+00}, { 0.00000e+00,  0.00000e+00,  1.75425e+00}, {-2.26473e+00,  0.00000e+00,  0.00000e+00}, 
            { 0.00000e+00, -4.52947e-01,  0.00000e+00}, {-0.00000e+00,  0.00000e+00, -2.39677e+00}, {-0.00000e+00, -0.00000e+00,  0.00000e+00}, 
            { 2.44619e+00,  0.00000e+00,  0.00000e+00}, { 0.00000e+00,  1.46771e+00,  0.00000e+00}, { 0.00000e+00,  0.00000e+00,  2.07566e+00}, 
            {-0.00000e+00,  0.00000e+00,  0.00000e+00}, {-3.28191e+00,  0.00000e+00,  0.00000e+00}, { 0.00000e+00, -3.28191e+00,  0.00000e+00}, 
            {-1.90708e+00,  0.00000e+00,  0.00000e+00}, {-0.00000e+00,  0.00000e+00, -2.91311e+00}, { 0.00000e+00, -0.00000e+00,  0.00000e+00}, 
            { 2.76362e+00,  0.00000e+00, -0.00000e+00}, {-0.00000e+00,  9.21205e-01,  0.00000e+00}, { 0.00000e+00,  0.00000e+00,  2.76362e+00}, 
            { 0.00000e+00,  0.00000e+00,  0.00000e+00}, {-3.02739e+00,  0.00000e+00,  0.00000e+00}, {-0.00000e+00, -2.01826e+00,  0.00000e+00}, 
            {-0.00000e+00,  0.00000e+00, -2.36662e+00}, { 0.00000e+00, -0.00000e+00,  0.00000e+00}, { 4.09910e+00,  0.00000e+00,  0.00000e+00}, 
            {-0.00000e+00,  4.09910e+00,  0.00000e+00}, {-0.00000e+00,  0.00000e+00, -2.38995e+00}, { 3.16161e+00,  0.00000e+00, -0.00000e+00}, 
            { 0.00000e+00,  4.51658e-01,  0.00000e+00}, { 0.00000e+00,  0.00000e+00,  3.31900e+00}, {-0.00000e+00,  0.00000e+00,  0.00000e+00}, 
            {-3.28564e+00,  0.00000e+00,  0.00000e+00}, { 0.00000e+00, -1.40813e+00,  0.00000e+00}, {-0.00000e+00,  0.00000e+00, -3.11349e+00}, 
            {-0.00000e+00, -0.00000e+00,  0.00000e+00}, { 3.63241e+00,  0.00000e+00, -0.00000e+00}, { 0.00000e+00,  2.59458e+00,  0.00000e+00}, 
            { 0.00000e+00,  0.00000e+00,  2.64596e+00}, { 0.00000e+00,  0.00000e+00,  0.00000e+00}, {-4.95014e+00,  0.00000e+00, -0.00000e+00}, 
            { 0.00000e+00, -4.95014e+00,  0.00000e+00}
        },
        {
            { 0.00000e+00,  0.00000e+00,  0.00000e+00}, { 0.00000e+00,  0.00000e+00,  4.88603e-01}, {-4.88603e-01,  0.00000e+00,  0.00000e+00}, 
            { 0.00000e+00, -4.88603e-01,  0.00000e+00}, { 0.00000e+00, -6.30783e-01,  0.00000e+00}, {-0.00000e+00,  0.00000e+00,  0.00000e+00}, 
            { 0.00000e+00, -0.00000e+00, -1.09255e+00}, { 0.00000e+00, -1.09255e+00,  0.00000e+00}, { 1.09255e+00,  0.00000e+00, -0.00000e+00}, 
            { 0.00000e+00, -0.00000e+00, -1.11953e+00}, { 4.57046e-01,  0.00000e+00,  0.00000e+00}, { 0.00000e+00,  1.37114e+00, -0.00000e+00}, 
            { 0.00000e+00, -0.00000e+00, -1.44531e+00}, { 0.00000e+00,  0.00000e+00,  0.00000e+00}, { 1.77013e+00,  0.00000e+00,  0.00000e+00}, 
            { 0.00000e+00,  1.77013e+00,  0.00000e+00}, { 0.00000e+00,  1.26943e+00, -0.00000e+00}, { 0.00000e+00,  0.00000e+00,  0.00000e+00}, 
            { 0.00000e+00,  0.00000e+00,  2.00714e+00}, { 0.00000e+00,  1.89235e+00, -0.00000e+00}, {-9.46175e-01,  0.00000e+00,  0.00000e+00}, 
            { 0.00000e+00,  0.00000e+00,  0.00000e+00}, { 0.00000e+00,  0.00000e+00,  1.77013e+00}, { 0.00000e+00,  2.50334e+00, -0.00000e+00}, 
            {-2.50334e+00,  0.00000e+00,  0.00000e+00}, { 0.00000e+00,  0.00000e+00,  1.75425e+00}, {-4.52947e-01,  0.00000e+00,  0.00000e+00}, 
            { 0.00000e+00, -2.26473e+00,  0.00000e+00}, { 0.00000e+00,  0.00000e+00,  2.39677e+00}, {-0.00000e+00,  0.00000e+00,  0.00000e+00}, 
            {-1.46771e+00,  0.00000e+00,  0.00000e+00}, { 0.00000e+00, -2.44619e+00,  0.00000e+00}, { 0.00000e+00,  0.00000e+00,  2.07566e+00}, 
            {-0.00000e+00,  0.00000e+00,  0.00000e+00}, {-3.28191e+00,  0.00000e+00,  0.00000e+00}, { 0.00000e+00, -3.28191e+00,  0.00000e+00}, 
            { 0.00000e+00, -1.90708e+00,  0.00000e+00}, {-0.00000e+00,  0.00000e+00,  0.00000e+00}, { 0.00000e+00, -0.00000e+00, -2.91311e+00}, 
            { 0.00000e+00, -2.76362e+00,  0.00000e+00}, { 9.21205e-01,  0.00000e+00,  0.00000e+00}, {-0.00000e+00,  0.00000e+00,  0.00000e+00}, 
            { 0.00000e+00, -0.00000e+00, -2.76362e+00}, { 0.00000e+00, -3.02739e+00,  0.00000e+00}, { 2.01826e+00,  0.00000e+00,  0.00000e+00}, 
            {-0.00000e+00,  0.00000e+00,  0.00000e+00}, { 0.00000e+00, -0.00000e+00, -2.36662e+00}, { 0.00000e+00, -4.09910e+00,  0.00000e+00}, 
            { 4.09910e+00,  0.00000e+00,  0.00000e+00}, { 0.00000e+00, -0.00000e+00, -2.38995e+00}, { 4.51658e-01,  0.00000e+00,  0.00000e+00}, 
            { 0.00000e+00,  3.16161e+00, -0.00000e+00}, { 0.00000e+00, -0.00000e+00, -3.31900e+00}, { 0.00000e+00,  0.00000e+00,  0.00000e+00}, 
            { 1.40813e+00,  0.00000e+00,  0.00000e+00}, { 0.00000e+00,  3.28564e+00, -0.00000e+00}, { 0.00000e+00, -0.00000e+00, -3.11349e+00}, 
            { 0.00000e+00,  0.00000e+00,  0.00000e+00}, { 2.59458e+00,  0.00000e+00,  0.00000e+00}, { 0.00000e+00,  3.63241e+00, -0.00000e+00}, 
            { 0.00000e+00,  0.00000e+00, -2.64596e+00}, { 0.00000e+00,  0.00000e+00,  0.00000e+00}, { 4.95014e+00,  0.00000e+00,  0.00000e+00}, 
            { 0.00000e+00,  4.95014e+00, -0.00000e+00}           
        },
        {
            { 0.00000e+00,  0.00000e+00,  0.00000e+00}, { 0.00000e+00,  0.00000e+00,  4.88603e-01}, {-4.88603e-01,  0.00000e+00,  0.00000e+00}, 
            { 0.00000e+00, -4.88603e-01,  0.00000e+00}, { 0.00000e+00,  0.00000e+00,  1.26157e+00}, {-1.09255e+00,  0.00000e+00,  0.00000e+00}, 
            { 0.00000e+00, -1.09255e+00,  0.00000e+00}, { 0.00000e+00,  0.00000e+00,  2.22045e-16}, {-0.00000e+00,  0.00000e+00, -0.00000e+00}, 
            { 0.00000e+00,  0.00000e+00,  2.23906e+00}, {-1.82818e+00,  0.00000e+00,  0.00000e+00}, { 0.00000e+00, -1.82818e+00,  0.00000e+00}, 
            { 0.00000e+00,  0.00000e+00,  8.81212e-16}, {-0.00000e+00,  0.00000e+00,  0.00000e+00}, {-1.84324e-16,  0.00000e+00,  0.00000e+00}, 
            { 0.00000e+00,  5.55112e-17,  0.00000e+00}, { 0.00000e+00,  0.00000e+00,  3.38514e+00}, {-2.67619e+00,  0.00000e+00,  0.00000e+00}, 
            { 0.00000e+00, -2.67619e+00,  0.00000e+00}, { 0.00000e+00,  0.00000e+00,  2.30756e-15}, {-0.00000e+00,  0.00000e+00,  0.00000e+00}, 
            {-5.52973e-16,  0.00000e+00,  0.00000e+00}, { 0.00000e+00,  1.66533e-16,  0.00000e+00}, { 0.00000e+00,  0.00000e+00,  0.00000e+00}, 
            {-0.00000e+00,  0.00000e+00,  0.00000e+00}, { 0.00000e+00,  0.00000e+00,  4.67801e+00}, {-3.62357e+00,  0.00000e+00,  0.00000e+00}, 
            { 0.00000e+00, -3.62357e+00,  0.00000e+00}, { 0.00000e+00,  0.00000e+00,  4.87108e-15}, {-0.00000e+00,  0.00000e+00,  0.00000e+00}, 
            {-1.22267e-15,  0.00000e+00,  0.00000e+00}, { 0.00000e+00,  3.68219e-16,  0.00000e+00}, { 0.00000e+00,  0.00000e+00,  0.00000e+00}, 
            {-0.00000e+00,  0.00000e+00,  0.00000e+00}, { 4.93038e-32,  0.00000e+00,  0.00000e+00}, { 0.00000e+00, -6.16298e-33,  0.00000e+00}, 
            { 0.00000e+00,  0.00000e+00,  6.10264e+00}, {-4.66097e+00,  0.00000e+00,  0.00000e+00}, { 0.00000e+00, -4.66097e+00,  0.00000e+00}, 
            { 0.00000e+00,  0.00000e+00,  8.98664e-15}, {-0.00000e+00,  0.00000e+00,  0.00000e+00}, {-2.30221e-15,  0.00000e+00,  0.00000e+00}, 
            { 0.00000e+00,  6.93334e-16,  0.00000e+00}, { 0.00000e+00,  0.00000e+00,  0.00000e+00}, {-0.00000e+00,  0.00000e+00,  0.00000e+00}, 
            { 1.77767e-31,  0.00000e+00,  0.00000e+00}, { 0.00000e+00, -2.22209e-32,  0.00000e+00}, { 0.00000e+00,  0.00000e+00,  0.00000e+00}, 
            {-0.00000e+00,  0.00000e+00,  0.00000e+00}, { 0.00000e+00,  0.00000e+00,  7.64784e+00}, {-5.78122e+00,  0.00000e+00,  0.00000e+00}, 
            { 0.00000e+00, -5.78122e+00,  0.00000e+00}, { 0.00000e+00,  0.00000e+00,  1.51096e-14}, {-0.00000e+00,  0.00000e+00,  0.00000e+00}, 
            {-3.91011e-15,  0.00000e+00,  0.00000e+00}, { 0.00000e+00,  1.17757e-15,  0.00000e+00}, { 0.00000e+00,  0.00000e+00,  0.00000e+00}, 
            {-0.00000e+00,  0.00000e+00,  0.00000e+00}, { 4.67737e-31,  0.00000e+00,  0.00000e+00}, { 0.00000e+00, -5.84671e-32,  0.00000e+00}, 
            { 0.00000e+00,  0.00000e+00,  0.00000e+00}, { 0.00000e+00,  0.00000e+00,  0.00000e+00}, { 1.13319e-47,  0.00000e+00,  0.00000e+00}, 
            { 0.00000e+00, -1.41649e-48,  0.00000e+00}
        },
        {
            { 0.00000e+00,  0.00000e+00,  0.00000e+00}, { 0.00000e+00,  0.00000e+00,  4.88603e-01}, {-4.88603e-01,  0.00000e+00,  0.00000e+00}, 
            { 0.00000e+00, -4.88603e-01,  0.00000e+00}, { 3.64183e-01,  3.64183e-01, -7.28366e-01}, { 6.30783e-01, -0.00000e+00,  6.30783e-01}, 
            {-0.00000e+00,  6.30783e-01,  6.30783e-01}, {-6.30783e-01,  6.30783e-01, -1.66533e-16}, {-6.30783e-01, -6.30783e-01,  0.00000e+00}, 
            {-7.46353e-01, -7.46353e-01,  0.00000e+00}, { 0.00000e+00,  3.04697e-01, -1.21879e+00}, { 3.04697e-01,  0.00000e+00, -1.21879e+00}, 
            { 9.63537e-01, -9.63537e-01,  4.01253e-16}, { 9.63537e-01,  9.63537e-01,  9.63537e-01}, {-4.44089e-16,  1.18009e+00, -2.22045e-16}, 
            {-1.18009e+00, -1.11022e-16,  0.00000e+00}, { 4.88603e-01,  4.88603e-01,  1.30294e+00}, {-1.03006e+00, -7.72548e-01,  7.72548e-01}, 
            {-7.72548e-01, -1.03006e+00,  7.72548e-01}, {-7.28366e-01,  7.28366e-01, -5.25363e-16}, {-3.64183e-01, -3.64183e-01, -2.18510e+00}, 
            { 7.69185e-16, -2.04397e+00, -6.81324e-01}, { 2.04397e+00,  1.92296e-16,  6.81324e-01}, { 9.63537e-01,  9.63537e-01, -1.44756e-16}, 
            {-9.63537e-01,  9.63537e-01, -5.55112e-17}, { 5.19779e-01,  5.19779e-01, -1.81923e+00}, { 1.40917e+00,  8.05238e-01,  8.05238e-01}, 
            { 8.05238e-01,  1.40917e+00,  8.05238e-01}, { 0.00000e+00, -4.44089e-16,  3.24739e-16}, {-1.06523e+00, -1.06523e+00,  2.13046e+00}, 
            {-2.17439e-01,  1.73951e+00,  1.73951e+00}, {-1.73951e+00,  2.17439e-01, -1.73951e+00}, {-1.84503e+00, -1.84503e+00, -9.22517e-01}, 
            { 1.84503e+00, -1.84503e+00,  6.58625e-16}, { 1.45863e+00,  1.11022e-15,  0.00000e+00}, {-8.88178e-16,  1.45863e+00,  0.00000e+00}, 
            {-1.46807e+00, -1.46807e+00,  5.87227e-01}, {-4.48502e-01, -3.36617e-16, -2.24251e+00}, {-3.36617e-16, -4.48502e-01, -2.24251e+00}, 
            { 7.09144e-01, -7.09144e-01,  1.87222e-16}, { 2.12743e+00,  2.12743e+00, -9.38779e-16}, { 7.09144e-01, -5.11006e-16, -2.12743e+00}, 
            { 1.02201e-15, -7.09144e-01,  2.12743e+00}, { 1.81260e+00,  1.81260e+00,  2.58943e+00}, {-2.07154e+00,  2.07154e+00, -1.66969e-15}, 
            {-3.03637e+00, -2.31111e-15, -6.07275e-01}, { 1.84889e-15, -3.03637e+00, -6.07275e-01}, { 1.05183e+00, -1.05183e+00,  5.77778e-17}, 
            { 1.05183e+00,  1.05183e+00,  4.03986e-17}, { 1.27464e+00,  1.27464e+00,  1.69952e+00}, {-1.28472e+00, -1.20442e+00,  1.92707e+00}, 
            {-1.20442e+00, -1.28472e+00,  1.92707e+00}, {-8.52285e-01,  8.52285e-01, -6.74704e-16}, {-1.50789e+00, -1.50789e+00, -2.95022e+00}, 
            {-1.11260e+00, -2.08612e+00,  9.27164e-01}, { 2.08612e+00,  1.11260e+00, -9.27164e-01}, {-3.07506e-01, -3.07506e-01, -3.69007e+00}, 
            { 1.23002e+00, -1.23002e+00,  2.28018e-15}, { 3.69007e+00, -1.53753e-01,  1.84503e+00}, {-1.53753e-01,  3.69007e+00,  1.84503e+00}, 
            {-2.35197e+00,  2.35197e+00, -8.00513e-16}, {-2.35197e+00, -2.35197e+00, -7.83988e-01}, { 1.37903e-15, -1.46671e+00,  9.77875e-17}, 
            { 1.46671e+00,  1.14919e-15,  1.34475e-16}
        }
    };

    void SetUp()
    {
        ylm.create(nylm,ng);
        dylm = new ModuleBase::matrix[3];
        for(int i = 0 ; i < 3 ; ++i)    dylm[i].create(nylm,ng);
        g = new ModuleBase::Vector3<double>[ng];
        g[0].set(1.0,0.0,0.0);
        g[1].set(0.0,1.0,0.0);
        g[2].set(0.0,0.0,1.0);
        g[3].set(-1.0,-1.0,-1.0);

        rly = new double[nylm];
        rlyvector.resize(nylm);
        rlgy = new double[nylm][3];
        rlgyvector.resize(nylm,std::vector<double>(3));
        ref = new double[64*4]{
            y00(g[0].x, g[0].y, g[0].z),  y00(g[1].x, g[1].y, g[1].z),  y00(g[2].x, g[2].y, g[2].z),  y00(g[3].x, g[3].y, g[3].z),  
            y10(g[0].x, g[0].y, g[0].z),  y10(g[1].x, g[1].y, g[1].z),  y10(g[2].x, g[2].y, g[2].z),  y10(g[3].x, g[3].y, g[3].z),  
            y11(g[0].x, g[0].y, g[0].z),  y11(g[1].x, g[1].y, g[1].z),  y11(g[2].x, g[2].y, g[2].z),  y11(g[3].x, g[3].y, g[3].z),  
            y1m1(g[0].x, g[0].y, g[0].z), y1m1(g[1].x, g[1].y, g[1].z), y1m1(g[2].x, g[2].y, g[2].z), y1m1(g[3].x, g[3].y, g[3].z), 
            y20(g[0].x, g[0].y, g[0].z),  y20(g[1].x, g[1].y, g[1].z),  y20(g[2].x, g[2].y, g[2].z),  y20(g[3].x, g[3].y, g[3].z),  
            y21(g[0].x, g[0].y, g[0].z),  y21(g[1].x, g[1].y, g[1].z),  y21(g[2].x, g[2].y, g[2].z),  y21(g[3].x, g[3].y, g[3].z),  
            y2m1(g[0].x, g[0].y, g[0].z), y2m1(g[1].x, g[1].y, g[1].z), y2m1(g[2].x, g[2].y, g[2].z), y2m1(g[3].x, g[3].y, g[3].z),                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
            y22(g[0].x, g[0].y, g[0].z),  y22(g[1].x, g[1].y, g[1].z),  y22(g[2].x, g[2].y, g[2].z),  y22(g[3].x, g[3].y, g[3].z),  
            y2m2(g[0].x, g[0].y, g[0].z), y2m2(g[1].x, g[1].y, g[1].z), y2m2(g[2].x, g[2].y, g[2].z), y2m2(g[3].x, g[3].y, g[3].z), 
            y30(g[0].x, g[0].y, g[0].z),  y30(g[1].x, g[1].y, g[1].z),  y30(g[2].x, g[2].y, g[2].z),  y30(g[3].x, g[3].y, g[3].z),  
            y31(g[0].x, g[0].y, g[0].z),  y31(g[1].x, g[1].y, g[1].z),  y31(g[2].x, g[2].y, g[2].z),  y31(g[3].x, g[3].y, g[3].z),  
            y3m1(g[0].x, g[0].y, g[0].z), y3m1(g[1].x, g[1].y, g[1].z), y3m1(g[2].x, g[2].y, g[2].z), y3m1(g[3].x, g[3].y, g[3].z), 
            y32(g[0].x, g[0].y, g[0].z),  y32(g[1].x, g[1].y, g[1].z),  y32(g[2].x, g[2].y, g[2].z),  y32(g[3].x, g[3].y, g[3].z),  
            y3m2(g[0].x, g[0].y, g[0].z), y3m2(g[1].x, g[1].y, g[1].z), y3m2(g[2].x, g[2].y, g[2].z), y3m2(g[3].x, g[3].y, g[3].z), 
            y33(g[0].x, g[0].y, g[0].z),  y33(g[1].x, g[1].y, g[1].z),  y33(g[2].x, g[2].y, g[2].z),  y33(g[3].x, g[3].y, g[3].z),  
            y3m3(g[0].x, g[0].y, g[0].z), y3m3(g[1].x, g[1].y, g[1].z), y3m3(g[2].x, g[2].y, g[2].z), y3m3(g[3].x, g[3].y, g[3].z), 
            y40(g[0].x, g[0].y, g[0].z),  y40(g[1].x, g[1].y, g[1].z),  y40(g[2].x, g[2].y, g[2].z),  y40(g[3].x, g[3].y, g[3].z),  
            y41(g[0].x, g[0].y, g[0].z),  y41(g[1].x, g[1].y, g[1].z),  y41(g[2].x, g[2].y, g[2].z),  y41(g[3].x, g[3].y, g[3].z),  
            y4m1(g[0].x, g[0].y, g[0].z), y4m1(g[1].x, g[1].y, g[1].z), y4m1(g[2].x, g[2].y, g[2].z), y4m1(g[3].x, g[3].y, g[3].z), 
            y42(g[0].x, g[0].y, g[0].z),  y42(g[1].x, g[1].y, g[1].z),  y42(g[2].x, g[2].y, g[2].z),  y42(g[3].x, g[3].y, g[3].z),  
            y4m2(g[0].x, g[0].y, g[0].z), y4m2(g[1].x, g[1].y, g[1].z), y4m2(g[2].x, g[2].y, g[2].z), y4m2(g[3].x, g[3].y, g[3].z), 
            y43(g[0].x, g[0].y, g[0].z),  y43(g[1].x, g[1].y, g[1].z),  y43(g[2].x, g[2].y, g[2].z),  y43(g[3].x, g[3].y, g[3].z),  
            y4m3(g[0].x, g[0].y, g[0].z), y4m3(g[1].x, g[1].y, g[1].z), y4m3(g[2].x, g[2].y, g[2].z), y4m3(g[3].x, g[3].y, g[3].z), 
            y44(g[0].x, g[0].y, g[0].z),  y44(g[1].x, g[1].y, g[1].z),  y44(g[2].x, g[2].y, g[2].z),  y44(g[3].x, g[3].y, g[3].z),  
            y4m4(g[0].x, g[0].y, g[0].z), y4m4(g[1].x, g[1].y, g[1].z), y4m4(g[2].x, g[2].y, g[2].z), y4m4(g[3].x, g[3].y, g[3].z), 
              0.000000000000000,    0.000000000000000,    0.935602579627389,    0.090028400200397, 
             -0.452946651195697,   -0.000000000000000,   -0.000000000000000,   -0.348678494661834, 
             -0.000000000000000,   -0.452946651195697,   -0.000000000000000,   -0.348678494661834, 
             -0.000000000000000,    0.000000000000000,    0.000000000000000,   -0.000000000000000, 
             -0.000000000000000,   -0.000000000000000,    0.000000000000000,   -0.000000000000000, 
              0.489238299435250,    0.000000000000000,   -0.000000000000000,   -0.376615818502422, 
              0.000000000000000,   -0.489238299435250,   -0.000000000000000,    0.376615818502422, 
              0.000000000000000,    0.000000000000000,    0.000000000000000,    0.532615198330370, 
              0.000000000000000,    0.000000000000000,    0.000000000000000,   -0.000000000000000, 
             -0.656382056840170,   -0.000000000000000,   -0.000000000000000,   -0.168427714314628, 
             -0.000000000000000,   -0.656382056840170,   -0.000000000000000,   -0.168427714314628, 
             -0.317846011338142,   -0.317846011338142,    1.017107236282055,    0.226023830284901, 
             -0.000000000000000,   -0.000000000000000,   -0.000000000000000,    0.258942827786103, 
             -0.000000000000000,   -0.000000000000000,   -0.000000000000000,    0.258942827786103, 
              0.460602629757462,   -0.460602629757462,    0.000000000000000,   -0.000000000000000, 
              0.000000000000000,    0.000000000000000,    0.000000000000000,   -0.409424559784410, 
             -0.000000000000000,   -0.000000000000000,   -0.000000000000000,    0.136474853261470, 
             -0.000000000000000,    0.000000000000000,   -0.000000000000000,   -0.136474853261470, 
             -0.504564900728724,   -0.504564900728724,    0.000000000000000,   -0.598002845308118, 
             -0.000000000000000,   -0.000000000000000,    0.000000000000000,    0.000000000000000, 
             -0.000000000000000,   -0.000000000000000,   -0.000000000000000,    0.350610246256556, 
             -0.000000000000000,   -0.000000000000000,   -0.000000000000000,    0.350610246256556, 
              0.683184105191914,   -0.683184105191914,    0.000000000000000,   -0.000000000000000, 
              0.000000000000000,    0.000000000000000,    0.000000000000000,   -0.202424920056864, 
              0.000000000000000,    0.000000000000000,    1.092548430592079,   -0.350435072502801, 
              0.451658037912587,    0.000000000000000,   -0.000000000000000,    0.046358202625865, 
              0.000000000000000,    0.451658037912587,   -0.000000000000000,    0.046358202625865, 
              0.000000000000000,   -0.000000000000000,    0.000000000000000,    0.000000000000000, 
              0.000000000000000,    0.000000000000000,    0.000000000000000,    0.492067081245654, 
             -0.469376801586882,   -0.000000000000000,   -0.000000000000000,    0.187354445356332, 
             -0.000000000000000,    0.469376801586882,   -0.000000000000000,   -0.187354445356332, 
              0.000000000000000,    0.000000000000000,    0.000000000000000,    0.355076798886913, 
              0.000000000000000,    0.000000000000000,    0.000000000000000,   -0.000000000000000, 
              0.518915578720260,    0.000000000000000,   -0.000000000000000,   -0.443845998608641, 
              0.000000000000000,    0.518915578720260,   -0.000000000000000,   -0.443845998608641, 
              0.000000000000000,   -0.000000000000000,    0.000000000000000,    0.000000000000000, 
              0.000000000000000,    0.000000000000000,    0.000000000000000,    0.452635881587108, 
             -0.707162732524596,    0.000000000000000,   -0.000000000000000,    0.120972027847095, 
             -0.000000000000000,    0.707162732524596,   -0.000000000000000,   -0.120972027847095  
         } ; 
    }

    void TearDown()
    {
        delete [] dylm;
        delete [] g;
        delete [] ref;
        delete [] rly;
        delete [] rlgy;
    }
};

TEST_F(YlmRealTest,Constructor)
{
    EXPECT_NO_THROW(ModuleBase::YlmReal YR);
}

TEST_F(YlmRealTest,YlmReal)
{
    ModuleBase::YlmReal::Ylm_Real(nylm,ng,g,ylm);
    for(int i=0;i<nylm;++i)
    {
        for(int j=0;j<ng;++j) 
        {
            EXPECT_NEAR(ylm(i,j),ref[i*ng+j],doublethreshold)  << "Ylm[" << i << "], example " << j << " not pass";
        }
    } 
}

TEST_F(YlmRealTest,YlmRealTemplate)
{
    psi::DEVICE_CPU * cpu_ctx = {};
    ModuleBase::YlmReal::Ylm_Real(cpu_ctx, nylm, ng, reinterpret_cast<double*>(g), ylm.c);
    for(int i=0;i<nylm;++i)
    {
        for(int j=0;j<ng;++j)
        {
            EXPECT_NEAR(ylm(i,j),ref[i*ng+j],doublethreshold)  << "Ylm[" << i << "], example " << j << " not pass";
        }
    }
}

TEST_F(YlmRealTest,gradYlmReal)
{
    ModuleBase::YlmReal::grad_Ylm_Real(nylm,ng,g,ylm,dylm[0],dylm[1],dylm[2]);
    for(int i=0;i<nylm;++i)
    {
        for(int j=0;j<ng;++j) 
        {
            EXPECT_NEAR(ylm(i,j),ref[i*ng+j],doublethreshold)  << "Ylm[" << i << "], example " << j << " not pass";
        }
    }
    ModuleBase::matrix ylmplus(nylm,1);
    ModuleBase::matrix ylmminus(nylm,1);
    double step = 1e-7;
    for(int id = 0 ; id < 3 ; ++id)
    {
        for(int j=0;j<ng;++j) 
        {
            ModuleBase::Vector3<double> gplus = g[j];
            ModuleBase::Vector3<double> gminus = g[j];
            gplus[id] += step/2;
            gminus[id] -= step/2;
            ModuleBase::YlmReal::Ylm_Real(nylm,1,&gplus,ylmplus);
            ModuleBase::YlmReal::Ylm_Real(nylm,1,&gminus,ylmminus);
            for(int i=0;i<nylm;++i)
            {
                if(std::abs(ylmplus(i,0)) < 1e-6 && std::abs(ylmminus(i,0)) < 1e-6) continue;
                double diff = (ylmplus(i,0) - ylmminus(i,0))/step;
                EXPECT_NEAR(diff,dylm[id](i,j),1e-6) << "dYlm[" << id << "][" << i << "], example " << j << " not pass";
            }
        }
    }
}


TEST_F(YlmRealTest,YlmReal2)
{
    ModuleBase::YlmReal::Ylm_Real2(nylm,ng,g,ylm);
    for(int i=0;i<nylm;++i)
    {
        for(int j=0;j<ng;++j) 
        {
            EXPECT_NEAR(ylm(i,j),ref[i*ng+j],doublethreshold) << "Ylm[" << i << "], example " << j << " not pass";
        }
    } 
}


TEST_F(YlmRealTest,YlmRealRlylm)
{    
    for(int j=0;j<ng;++j)
    {
        ModuleBase::YlmReal::rlylm(lmax,g[j].x,g[j].y,g[j].z,rly);
        for(int i=0;i<nylm;++i)
        {
            EXPECT_NEAR(rly[i],ref[i*ng+j],doublethreshold) << "Ylm[" << i << "], example " << j << " not pass";
        }
    }
}


TEST_F(YlmRealTest,YlmGetYlmReal)
{    
    for(int j=0;j<ng;++j)
    {
        ModuleBase::Ylm::get_ylm_real(lmax+1,g[j],rly);
        for(int i=0;i<nylm;++i)
        {
            EXPECT_NEAR(rly[i],ref[i*ng+j],1e-5)  << "Ylm[" << i << "], example " << j << " not pass";
        }
    }
}

TEST_F(YlmRealTest,YlmSphHarm)
{    
    ModuleBase::Ylm::set_coefficients ();
    for(int j=0;j<ng;++j)
    {
        double r = sqrt(g[j].x * g[j].x + g[j].y * g[j].y + g[j].z * g[j].z);
        ModuleBase::Ylm::sph_harm(lmax,g[j].x/r,g[j].y/r,g[j].z/r,rlyvector);
        for(int i=0;i<nylm;++i)
        {
            EXPECT_NEAR(rlyvector[i],ref[i*ng+j],doublethreshold)  << "Ylm[" << i << "], example " << j << " not pass";
            
        }
    }
}

TEST_F(YlmRealTest,YlmRlSphHarm)
{    
    ModuleBase::Ylm::set_coefficients ();
    for(int j=0;j<ng;++j)
    {
        double r = sqrt(g[j].x * g[j].x + g[j].y * g[j].y + g[j].z * g[j].z);
        ModuleBase::Ylm::rl_sph_harm(lmax,g[j].x/r,g[j].y/r,g[j].z/r,rlyvector);
        for(int i=0;i<nylm;++i)
        {
            EXPECT_NEAR(rlyvector[i],ref[i*ng+j],doublethreshold)  << "Ylm[" << i << "], example " << j << " not pass";
            
        }
    }
}
//used to be test1 in ylm.h
TEST_F(YlmRealTest,YlmGradRlSphHarm)
{    
    ModuleBase::Ylm::set_coefficients ();
    for(int j=0;j<ng;++j)
    {
        double r = sqrt(g[j].x * g[j].x + g[j].y * g[j].y + g[j].z * g[j].z);
        ModuleBase::Ylm::grad_rl_sph_harm(lmax,g[j].x/r,g[j].y/r,g[j].z/r,rlyvector,rlgyvector);
        for(int i=0;i<nylm;++i)
        {
            EXPECT_NEAR(rlyvector[i],ref[i*ng+j],doublethreshold)  << "Ylm[" << i << "], example " << j << " not pass";
            for(int k=0;k<3;++k) {EXPECT_NEAR(rlgyvector[i][k],rlgyref[j][i][k],1e-5);}
            
        }
    }
}

//used to be test1 in ylm.h
TEST_F(YlmRealTest, equality_value_test)
{

    
    ModuleBase::Vector3<double> R (20.0, 0.0, 0.0);
	const double xdr = R.x/R.norm();
	const double ydr = R.y/R.norm();
	const double zdr = R.z/R.norm();
	const int L = 9;
	const double rl = std::pow( R.norm(), L);
	//std::cout << " rl=" << rl << std::endl;
	ModuleBase::Ylm::set_coefficients();
	
	int nu = 100;
	
	// Peize Lin change rlya 2016-08-26
	std::vector<double> rlya;
	double rlyb[400];
	ModuleBase::Ylm::ZEROS( rlyb, 400);
	
	ModuleBase::Ylm::rl_sph_harm(L, xdr, ydr, zdr, rlya);
	ModuleBase::Ylm::get_ylm_real(L+1, R, rlyb);
	
	for (int i=0; i < nu; i++)
	{
		double diff = fabs(rlya[i]-rlyb[i]);
        EXPECT_LT(diff,1e-8);
	}

}

//used to be test2 in ylm.h
TEST_F(YlmRealTest, equality_gradient_test)
{

    
    ModuleBase::Vector3<double> R (0.1,-0.2,0.5);
	ModuleBase::Ylm::set_coefficients();
	
	//int nu = 100;

	std::vector<double> rlya;
	double rlyb[400];
	
	std::vector<std::vector<double>> grlya;
	double grlyb[400][3];
	
	ModuleBase::Ylm::grad_rl_sph_harm (9, R.x, R.y, R.z, rlya, grlya);
	ModuleBase::Ylm::rlylm (10, R.x, R.y, R.z, rlyb, grlyb);
	
	for (int i = 0; i < 100; i++)
	{
		double diffx = fabs(grlya[i][2]-grlyb[i][2]);
        EXPECT_LT(diffx,1e-8);
	}

}
TEST_F(YlmRealTest,YlmRealTemplatefloat)
{
    ModuleBase::Vector3<float> *gg;
    gg = new ModuleBase::Vector3<float>[ng];
    gg[0].set(1.0,0.0,0.0);
    gg[1].set(0.0,1.0,0.0);
    gg[2].set(0.0,0.0,1.0);
    gg[3].set(-1.0,-1.0,-1.0);
    float*ccc;
    ccc=new float[nylm*ng];
    psi::DEVICE_CPU * cpu_ctx = {};
    ModuleBase::YlmReal::Ylm_Real<float, psi::DEVICE_CPU>(cpu_ctx, nylm, ng, reinterpret_cast<float*>(gg), ccc);
    for(int i=0;i<nylm;++i)
    {
        for(int j=0;j<ng;++j)
        {
            EXPECT_NEAR(ccc[i*ng+j],ref[i*ng+j],2.45e-06)  << "Ylm[" << i << "], example " << j << " not pass";
        }
    }
    delete [] gg;
    delete [] ccc;
}
