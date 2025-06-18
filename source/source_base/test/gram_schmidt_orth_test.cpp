#include"../gram_schmidt_orth.h"
#include"../gram_schmidt_orth-inl.h"
#include"gtest/gtest.h"


#define DOUBLETHRESHOLD 1e-8


/************************************************
*  unit test of class Gram_Schmidt_Orth
***********************************************/

/**
 * Based on an linearly independent, but not orthonormal, 
 * set of functions x:{x1,x2,x3,...}, we can construct an 
 * orthonormal set X:{X1, X2, X3, ...} by using Gram-Schmidt 
 * orthogonalization.
 * The new set X should has below properties:
 * 1. X1 = x1/||x1||
 * 2. <Xi,Xj> = 1 if (i == j) else 0
 * 
 * Note:in this class, for coordinate of sphere, the inner product 
 * of two radial function f(r) and g(r) equals the integral of r^2*f(r)*g(r)
 * $$ (f(r),g(r)) = {\int}r^2f(r)g(r)dr $$
 * 
 */

class GramSchmidtOrth 
{
    public:
    int nbasis;
    int ndim;
    double dr;
    std::vector<double> r2;
    double norm0;
    ModuleBase::Gram_Schmidt_Orth<double,double>::Coordinate coordinate;
    std::vector<double> rab;
    std::vector<std::vector<double>> basis;

    GramSchmidtOrth(int nbasis, int ndim, double dr,
                        ModuleBase::Gram_Schmidt_Orth<double,double>::Coordinate coordinate):
                        nbasis(nbasis),ndim(ndim),dr(dr),coordinate(coordinate)
    {
        basis.resize(nbasis,std::vector<double>(ndim));
        rab.resize(ndim,dr);
        r2.resize(ndim,1.0);

        norm0 = sqrt(1.0/3.0 * pow(dr*(static_cast<double>(ndim-1)),3.0));
        if (ModuleBase::Gram_Schmidt_Orth<double,double>::Coordinate::Sphere == this->coordinate)
        {
            for(int i=0;i<ndim;++i) {r2[i] = dr*i*dr*i;}
            norm0 = sqrt(1.0/5.0 * pow(dr*(static_cast<double>(ndim-1)),5.0));
        }

        //build the function basis
        for(int i=0;i<nbasis;++i)
        {
            for(int j=0;j<ndim;++j)
            {
                //function: f_i(x) = x^(i+1)
                basis[i][j] = pow(static_cast<double>(j) * dr, static_cast<double>(i+1));
            }
        }
    }

    //calculate the inner product of two vector
    double inner_product(std::vector<double> a, std::vector<double> b)
    {
        double ip;
        std::vector<double> mul_func = ModuleBase::Mathzone::Pointwise_Product(a,b);
        std::vector<double> mul_func1 = ModuleBase::Mathzone::Pointwise_Product(mul_func,r2);
        ModuleBase::Integral::Simpson_Integral(mul_func1.size(),ModuleBase::GlobalFunc::VECTOR_TO_PTR(mul_func1),ModuleBase::GlobalFunc::VECTOR_TO_PTR(rab),ip);       
        return ip;
    }
    
};

class GramSchmidtOrthTest : public ::testing::TestWithParam<GramSchmidtOrth> {};


TEST_P(GramSchmidtOrthTest,CalOrth)
{
    GramSchmidtOrth gsot = GetParam();
    ModuleBase::Gram_Schmidt_Orth<double,double> gso_sphere(gsot.rab,gsot.coordinate);
    std::vector<std::vector<double>> old_basis = gsot.basis; 
    std::vector<std::vector<double>> new_basis = gso_sphere.cal_orth(old_basis);
    
    //==========================================================
    // VERIFY X0=x0/|x0|
    // the integral of old_basis[0] = {\int}_{0}^{dr*(ndim-1)} r^2*r*r dr
    // =1/5*r^5|_{0}^{dr*(ndim-1)}
    //==========================================================
    for(int i=0;i<gsot.ndim;i++)
    {
        EXPECT_NEAR(old_basis[0][i]/gsot.norm0,new_basis[0][i],DOUBLETHRESHOLD) << "the first basis is wrong";
    }

    //==========================================================
    // VERIFY <Xi,Xj> = 0 for i!=j
    //==========================================================
    int niter = 1;
    int maxiter = 1;
    bool pass = false;
    double maxip; 

    //do iteration.
    while (true)
    {    
        int nbasis = new_basis.size();
        maxip = std::abs(gsot.inner_product(new_basis[nbasis-1],new_basis[nbasis-2]));           
        for(int i=0;i<nbasis-1;++i)
        {
            for(int j=i+1;j<nbasis;++j)
            {
                double ip = gsot.inner_product(new_basis[i],new_basis[j]);
                //std::cout << "i=" << i << ", j=" << j << ": " << ip << std::endl;
                if(std::abs(ip) > maxip) {maxip = std::abs(ip);}
            }            
        }
        if (maxip < DOUBLETHRESHOLD) {pass = true; break;};
        if (niter >= maxiter) {break;}

        niter += 1; 
        old_basis = gso_sphere.cal_orth(new_basis); new_basis = old_basis;
    }

    //std::cout <<  "nbasis=" << gsot.nbasis << "niter=" << niter << " max_inner_product=" << std::setprecision(15) << maxip << std::endl;
    EXPECT_TRUE(pass) << "nbasis=" << gsot.nbasis << "niter=" << niter << " max_inner_product=" << std::setprecision(15) << maxip;
}

INSTANTIATE_TEST_SUITE_P(VerifyOrth,GramSchmidtOrthTest,::testing::Values(
        GramSchmidtOrth(10,101,0.1,ModuleBase::Gram_Schmidt_Orth<double,double>::Coordinate::Sphere),
        GramSchmidtOrth(20,1001,0.01,ModuleBase::Gram_Schmidt_Orth<double,double>::Coordinate::Sphere),
        GramSchmidtOrth(50,10001,0.001,ModuleBase::Gram_Schmidt_Orth<double,double>::Coordinate::Sphere),
        GramSchmidtOrth(10,10001,0.001,ModuleBase::Gram_Schmidt_Orth<double,double>::Coordinate::Cartesian),
        GramSchmidtOrth(20,1001,0.01,ModuleBase::Gram_Schmidt_Orth<double,double>::Coordinate::Cartesian),
        GramSchmidtOrth(50,101,0.1,ModuleBase::Gram_Schmidt_Orth<double,double>::Coordinate::Cartesian)
));

