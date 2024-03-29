#include"../math_integral.h"
#include"gtest/gtest.h"

#include<math.h>
#include <fstream>
#include <iostream>

#include "module_base/constants.h"

#define doublethreshold 1e-12



/************************************************
*  unit test of class Integral
***********************************************/

/**
 * Tested functions:
 *  - function Simpson_Integral with rab as input
 *  - function Simpson_Integral with dr as input
 *  - function Simpson_Integral_0toall
 *  - function Simpson_Integral_alltoinf
 *  - function Simpson_Integral with f(x) = sin(x)
 *  - function Simpson_Integral with f(x) = 1 / (1 + x^2)
 *  - function Simpson_Integral with f(x) = exp(x)
 */

// generate irregular grid with sinx
void sinspace(double start, double end, const int nums, double* xv, double* h){
    double astart = asin(start);
    double aend = asin(end);
    double step = (aend - astart) / (nums - 1);

    for(int i = 0; i < nums; ++i){
        h[i] = sin(astart + i * step);
    }
    // calculate the difference
    xv[0] = start;
    for(int i = 0; i< nums - 1; ++i){
        h[i] = h[i+1] - h[i];
        xv[i+1] = xv[i] + h[i]; 
    }
}

class SimpsonIntegralSinx : public testing::Test
{
    /**
     * test the integral of sinx between [0,PI],
     * devide to mesh-1 parts
     */
    
    protected:

    double*     func;
    double*     rab;
    int         mesh = 10001;
    double      dr   = M_PI/(mesh-1);
    double      asum;
    double*     asumlist;
    double      expectvalue = 2.0;

    void SetUp()
    {
        func        = new double[mesh];
        rab         = new double[mesh];
        asumlist    = new double[mesh];

        for (int i=0;i<mesh;i++)
        {
            func[i] = sin(M_PI*i/(mesh-1));
            rab[i] = dr;
        }
    }

    void TearDown()
    {
        delete[] func;
        delete[] rab;
        delete[] asumlist;
    }
};

TEST_F(SimpsonIntegralSinx,Constructor)
{
    EXPECT_NO_THROW(ModuleBase::Integral intgral);
}

TEST_F(SimpsonIntegralSinx,SimpsonIntegralRab)
{
    ModuleBase::Integral::Simpson_Integral(mesh,func,rab,asum);
    EXPECT_NEAR(asum,expectvalue,doublethreshold);
    EXPECT_DEATH(ModuleBase::Integral::Simpson_Integral(mesh-1,func,rab,asum),"");
}

TEST_F(SimpsonIntegralSinx,SimpsonIntegralDr)
{
    ModuleBase::Integral::Simpson_Integral(mesh,func,dr,asum);
    EXPECT_NEAR(asum,expectvalue,doublethreshold);
    EXPECT_DEATH(ModuleBase::Integral::Simpson_Integral(mesh-1,func,dr,asum),"");
}

TEST_F(SimpsonIntegralSinx,SimpsonIntegral0toall)
{
    int halfmesh = round(mesh/2);
    ModuleBase::Integral::Simpson_Integral_0toall(mesh,func,rab,asumlist);
    EXPECT_NEAR(asumlist[mesh-1],expectvalue,doublethreshold);
    EXPECT_NEAR(asumlist[mesh-2],1.9999999506519738901,doublethreshold);
    EXPECT_NEAR(asumlist[halfmesh],expectvalue/2.0,doublethreshold);
    EXPECT_NEAR(asumlist[halfmesh-1],0.99968584073722277505,doublethreshold);
}


TEST_F(SimpsonIntegralSinx,SimpsonIntegralAlltoinf)
{
    int halfmesh = round(mesh/2);
    ModuleBase::Integral::Simpson_Integral_alltoinf(mesh,func,rab,asumlist);
    EXPECT_NEAR(asumlist[0],expectvalue,doublethreshold);
    EXPECT_NEAR(asumlist[1],1.9999999506519754444,doublethreshold);
    EXPECT_NEAR(asumlist[halfmesh],expectvalue/2.0,doublethreshold);
    EXPECT_NEAR(asumlist[halfmesh-1],1.0003141592627740053,doublethreshold);
}

TEST_F(SimpsonIntegralSinx, UniformGridOdd)
{
    double start = 0.0;
    int ngrid_max = std::pow(3, 7);
    double* f = new double[ngrid_max];
    for (int ngrid = 3; ngrid <= ngrid_max; ngrid *= 3)
    {
        for (double end = 0.1; end < 10*ModuleBase::PI; end += 0.5)
        {
            double dx = (end-start) / (ngrid-1);
            for (int i = 0; i != ngrid; ++i)
            {
                f[i] = std::sin(start + i*dx);
            
            }
            // error bound of composite Simpson's rule (max|sin^(4)(x)| = 1)
            double tol = (end-start) * std::pow(dx, 4) / 180;
            EXPECT_NEAR(std::cos(start)-std::cos(end), ModuleBase::Integral::simpson(ngrid, f, dx), std::max(tol, doublethreshold));
        }
    }
    delete[] f;
}

TEST_F(SimpsonIntegralSinx, UniformGridEven)
{
    double start = 0.0;
    int ngrid_max = std::pow(2, 12);
    double* f = new double[ngrid_max];
    for (int ngrid = 4; ngrid <= ngrid_max; ngrid *= 2)
    {
        for (double end = 0.1; end < 10*ModuleBase::PI; end += 0.5)
        {
            double dx = (end-start) / (ngrid-1);
            for (int i = 0; i != ngrid; ++i)
            {
                f[i] = std::sin(start + i*dx);
            
            }
            // error bound of composite Simpson's rule (max|sin^(4)(x)| = 1)
            double tol = ((ngrid-4)*dx-start) * std::pow(dx, 4) / 180 + (end-(ngrid-4)*dx) * std::pow(dx, 4) / 80;
            EXPECT_NEAR(std::cos(start)-std::cos(end), ModuleBase::Integral::simpson(ngrid, f, dx), std::max(tol, doublethreshold));
        }
    }
    delete[] f;
}

TEST_F(SimpsonIntegralSinx, LogGridOdd)
{
    double start = 0.0;
    int ngrid_max = std::pow(3, 7);
    double* f = new double[ngrid_max];
    double* h = new double[ngrid_max];
    for (int ngrid = 3; ngrid <= ngrid_max; ngrid *= 3)
    {
        for (double w = 0.1; w < 10*ModuleBase::PI; w += 0.5)
        {
            double b = std::pow(1.0 + w, 1.0/(ngrid-1));
            for (int i = 0; i != ngrid; ++i)
            {
                double x = std::pow(b, i) - 1.0 + start;
                f[i] = std::sin(x);
                h[i] = std::pow(b, i+1) - std::pow(b, i); // the last one is not used
            }
            // FIXME error bound for irregularly-spaced Simpson's rule should be derived
            // the one below is a crude estimate
            double tol = w * std::pow(h[ngrid-2], 4);
            EXPECT_NEAR(std::cos(start)-std::cos(start+w), ModuleBase::Integral::simpson(ngrid, f, h), std::max(tol, doublethreshold));
            
            //printf("%6i %8.4f %20.12f\n", ngrid, w, ModuleBase::Integral::simpson(ngrid, f, h));
            //double tmp = 0.0;
            //ModuleBase::Integral::Simpson_Integral(ngrid, f, h, tmp);
            //EXPECT_NEAR(std::cos(start)-std::cos(end), tmp, std::max(tol, doublethreshold));
        }
    }
    delete[] f;
    delete[] h;
}

TEST_F(SimpsonIntegralSinx, LogGridEven)
{
    double start = 0.0;
    int ngrid_max = std::pow(2, 12);
    double* f = new double[ngrid_max];
    double* h = new double[ngrid_max];
    for (int ngrid = 4; ngrid <= ngrid_max; ngrid *= 2)
    {
        for (double w = 0.1; w < 10*ModuleBase::PI; w += 0.5)
        {
            double b = std::pow(1.0 + w, 1.0/(ngrid-1));
            for (int i = 0; i != ngrid; ++i)
            {
                double x = std::pow(b, i) - 1.0 + start;
                f[i] = std::sin(x);
                h[i] = std::pow(b, i+1) - std::pow(b, i); // the last one is not used
            }
            // FIXME error bound for irregularly-spaced Simpson's rule should be derived
            // the one below is a crude estimate
            double tol = w * std::pow(h[ngrid-2], 4);
            EXPECT_NEAR(std::cos(start)-std::cos(start+w), ModuleBase::Integral::simpson(ngrid, f, h), std::max(tol, doublethreshold));

            //double tmp = 0.0;
            //ModuleBase::Integral::Simpson_Integral(ngrid, f, h, tmp);
            //EXPECT_NEAR(std::cos(start)-std::cos(end), tmp, std::max(tol, doublethreshold));
        }
    }
    delete[] f;
    delete[] h;
}

TEST_F(SimpsonIntegralSinx, FourPoints)
{
    /*!
     * This test checks whether "simpson" yields the exact integral for a cubic polynomial 
     * sampled by four points. This case is handled by Simpson's 3/8 rule and should be exact.
     *                                                                                              */
    double A = 0.1; 
    double B = 0.25;
    double C = 0.40;
    double D = 0.77;
    double x0 = 0.52;
    
    auto f = [&](double x) { return A*std::pow(x-x0,3) + B*std::pow(x-x0,2) + C*(x-x0) + D;};
    auto I = [&](double x) { return A*std::pow(x-x0,4)/4 + B*std::pow(x-x0,3)/3 + C*std::pow(x-x0,2)/2 + D*(x-x0);};
    
    double x[4];
    double y[4];

    // unevenly-spaced 4 points
    x[0] = 0.11;
    x[1] = 0.22;
    x[2] = 0.44;
    x[3] = 0.88;
    
    for (int i = 0; i != 4; ++i)
    {
        y[i] = f(x[i]);
    }

    double ref = I(x[3]) - I(x[0]);
    double h[3] = {x[1]-x[0], x[2]-x[1], x[3]-x[2]};
    EXPECT_DOUBLE_EQ(ref, ModuleBase::Integral::simpson(4, y, h));

    // evenly-spaced 4 points
    x[0] = 0.1;
    x[1] = 0.2;
    x[2] = 0.3;
    x[3] = 0.4;
    for (int i = 0; i != 4; ++i)
    {
        y[i] = f(x[i]);
    }
    
    ref = I(x[3]) - I(x[0]);
    EXPECT_DOUBLE_EQ(ref, ModuleBase::Integral::simpson(4, y, x[1]-x[0]));
}

TEST_F(SimpsonIntegralSinx, ThreePoints)
{
    /*!
     * This test checks whether "simpson" yields the exact integral for a quadratic polynomial 
     * sampled by three points.
     *                                                                                              */
    double A = 0.17; 
    double B = 0.29;
    double C = 0.35;
    double x0 = 0.5;
    
    auto f = [&](double x) { return A*std::pow(x-x0,2) + B*(x-x0) + C;};
    auto I = [&](double x) { return A*std::pow(x-x0,3)/3 + B*std::pow(x-x0,2)/2 + C*(x-x0);};
    
    double x[3];

    // unevenly-spaced 3 points
    x[0] = 0.1;
    x[1] = 0.2;
    x[2] = 0.4;
    
    double y[3];
    y[0] = f(x[0]);
    y[1] = f(x[1]);
    y[2] = f(x[2]);

    double h[2] = {x[1]-x[0], x[2]-x[1]};

    double ref = I(x[2]) - I(x[0]);
    EXPECT_DOUBLE_EQ(ref, ModuleBase::Integral::simpson(3, y, h));

    // evenly-spaced 4 points
    x[0] = 0.1;
    x[1] = 0.2;
    x[2] = 0.3;
    y[0] = f(x[0]);
    y[1] = f(x[1]);
    y[2] = f(x[2]);
    
    ref = I(x[2]) - I(x[0]);
    EXPECT_DOUBLE_EQ(ref, ModuleBase::Integral::simpson(3, y, x[1]-x[0]));
}

class SimpsonIntegralITF : public testing::Test{

};

TEST_F(SimpsonIntegralITF, UniformGridOdd)
{
    double start = -1.0, end = 1.0;
    const int ngrid_max = 10000;
    double *f = new double[ngrid_max];
    int ind = 0;
    std::ofstream file_o("data/itf_uni_out.bin", std::ios::binary);
    for (int ngrid = 5; ngrid <= ngrid_max; ngrid += 2) {
        const double dx = (end - start) / (ngrid - 1);
        for (int i = 0; i < ngrid; ++i) {
            double x = start + i * dx;
            f[i] = 1.0 / (1.0 + x * x);
        }
        double tol = (end-start) * std::pow(dx, 4) * 24 / 180;
        EXPECT_NEAR(std::atan(end) - std::atan(start), ModuleBase::Integral::simpson(ngrid, f, dx), std::max(tol, doublethreshold));
        double err = std::abs(ModuleBase::Integral::simpson(ngrid, f, dx) - ModuleBase::Integral::simpson(ngrid - 2, f, dx)) / ModuleBase::Integral::simpson(ngrid, f, dx);
        file_o.write(reinterpret_cast<char*>(&err), sizeof(double));
    }
    delete[] f;
    file_o.close();
}

TEST_F(SimpsonIntegralITF, SinGridOdd)
{
    double start = -1.0, end = 1.0;
    const int ngrid_max = 10000;
    double *xv = new double[ngrid_max];
    double *h = new double[ngrid_max];
    double *f = new double[ngrid_max];
    std::ofstream file_o("data/itf_sin_out.bin", std::ios::binary);
    for (int ngrid = 5; ngrid <= ngrid_max; ngrid += 2) {
        sinspace(start, end, ngrid, xv, h);
        for (int i = 0; i < ngrid; ++i) {
            f[i] = 1.0 / (1.0 + xv[i] * xv[i]);
        }
        
        // crude estimate for irregular-grid error bound
        double dx = h[ngrid / 2];
        double tol = (end-start) * std::pow(dx, 4);
        EXPECT_NEAR(std::atan(end) - std::atan(start), ModuleBase::Integral::simpson(ngrid, f, h), std::max(tol, doublethreshold));
        double err = std::abs(ModuleBase::Integral::simpson(ngrid, f, h) - ModuleBase::Integral::simpson(ngrid - 2, f, h)) / ModuleBase::Integral::simpson(ngrid, f, h);
        file_o.write(reinterpret_cast<char*>(&err), sizeof(double));
    }
    
    delete[] xv;
    delete[] h;
    delete[] f;
    file_o.close();
}

class SimpsonIntegralExp : public testing::Test{

};

TEST_F(SimpsonIntegralExp, UniformGridOdd)
{
    double start = 0.0, end = 1.0;
    const int ngrid_max = 10000;
    double *f = new double[ngrid_max];
    std::ofstream file_o("data/exp_uni_out.bin", std::ios::binary);
    for (int ngrid = 5; ngrid <= ngrid_max; ngrid += 2) {
        const double dx = (end - start) / (ngrid - 1);
        for (int i = 0; i < ngrid; ++i) {
            double x = start + i * dx;
            f[i] = std::exp(x);
        }
        double tol = (end-start) * std::exp(1) * std::pow(dx, 4) / 180;
        EXPECT_NEAR(std::exp(end) - std::exp(start), ModuleBase::Integral::simpson(ngrid, f, dx), std::max(tol, doublethreshold));
         double err = std::abs(ModuleBase::Integral::simpson(ngrid, f, dx) - ModuleBase::Integral::simpson(ngrid - 2, f, dx)) / ModuleBase::Integral::simpson(ngrid, f, dx);
        file_o.write(reinterpret_cast<char*>(&err), sizeof(double));
    }
    delete[] f;
    file_o.close();
}

TEST_F(SimpsonIntegralExp, SinGridOdd)
{
    double start = 0.0, end = 1.0;
    const int ngrid_max = 10000;
    double *xv = new double[ngrid_max];
    double *h = new double[ngrid_max];
    double *f = new double[ngrid_max];
    std::ofstream file_o("data/exp_sin_out.bin", std::ios::binary);

   // skip ngrid = 3 since the errors exceeds the threshold
    for (int ngrid = 5; ngrid <= ngrid_max; ngrid += 2) {
        sinspace(start, end, ngrid, xv, h);
        for (int i = 0; i < ngrid; ++i) {
            f[i] = std::exp(xv[i]);
        }

        double dx = h[ngrid / 2];
        double tol = (end-start) * std::pow(dx, 4);
        EXPECT_NEAR(std::exp(end) - std::exp(start), ModuleBase::Integral::simpson(ngrid, f, h), std::max(tol, doublethreshold));
        double err = std::abs(ModuleBase::Integral::simpson(ngrid, f, h) - ModuleBase::Integral::simpson(ngrid - 2, f, h)) / ModuleBase::Integral::simpson(ngrid, f, h);
        file_o.write(reinterpret_cast<char*>(&err), sizeof(double));
    }
    
    delete[] xv;
    delete[] h;
    delete[] f;
    file_o.close();
}