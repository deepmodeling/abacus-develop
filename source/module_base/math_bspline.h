#ifndef MATH_BSPLINE_H
#define MATH_BSPLINE_H

namespace ModuleBase
{

/**
 * @brief A class to treat Cardinal B-spline interpolation.
 *
 * @author qianrui created 2021-09-14
 * @details see: J. Chem. Phys. 103, 8577 (1995).
 * Math:
 * Only uniform nodes are considered: xm-x[m-1]=Dx(>= 0) for control node: X={x0,x1,...,xm};
 * Any function p(x) can be written by
 * p(x)=\sum_i{ci*M_ik(x)} (k->infinity),
 * where M_ik is the i-th k-order Cardinal B-spline base function
 * and ci is undetermined coefficient.
 * M_i0 =  H(x-xi)-H(x-x[i+1]), H(x): step function
 *                  x-xi                    x[i+k+1]-x
 *      M_ik(x)=  ---------*M_i(k-1)(x)+ ----------------*M_[i+1][k-1](x)  ( xi <= x <= x[i+1] )
 *                x[i+k]-xi               x[i+k+1]-x[i+1]
 *   For uniform nodes: M_[i+1]k(x+Dx)=M_ik(x)
 *   If we define Bk[n] stores M_ik(x+n*Dx) for x in (xi,xi+Dx):
 *               x+n*Dx-xi               xi+(k-n+1)*Dx-x
 *      Bk[n] = -----------*B(k-1)[n] + -----------------*B(k-1)[n-1]
 *                 k*Dx                        k*Dx
 * USAGE: 
 *   ModuleBase::Bspline bp;
 *   bp.init(10,0.7,2); //Dx = 0.7, xi = 2
 *   bp.getbslpine(0.5); //x = 0.5
 *   cout<<bp.bezier_ele(3)<<endl; //print M_ik[xi+3*Dx+x]: M_i[10](4.6)
 *
 */
class Bspline
{
  private:
    int norder; // the order of bezier base; norder >= 0
    double Dx; // Dx: the interval of control node
    double xi; // xi: the starting point
    double *bezier; // bezier[n] = Bk[n]

  public:
    Bspline();
    ~Bspline();

    void init(int norderin, double Dxin, double xiin);

    // Get the result of i-th bezier base functions for different input x+xi+n*Dx.
    // x should be in [0,Dx]
    // n-th result is stored in bezier[n];
    void getbspline(double x);

    // get the element of bezier
    double bezier_ele(int n);
};
} // namespace ModuleBase
#endif
