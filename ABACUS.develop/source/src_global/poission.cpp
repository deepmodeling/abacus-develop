#include "poission.h"
#include "integral.h"
#include <cassert>
#include "constants.h"
using namespace std;

void Poission::SolPoissonEq
(
    const double* rho,
    const double* r,
    const int& mesh,
    double* pot
)
{
    assert(mesh > 0);

    double a, b, c;

    const double tiny = 1e-12;

    //auxilliary array
    double* rad_f1 = new double[mesh];
    double* rad_f2 = new double[mesh];

    //initialization
    for(int ir = 0; ir < mesh; ir++)
    {
        rad_f2[ir] = r[ir] * rho[ir];
        rad_f1[ir] = rad_f2[ir] * r[ir];
    }

    //value at the beginning
    a = r[0];
    b = r[mesh-1];
    pot[0] = Integral_G::Gauss_Legendre(a, b, rad_f2, r, mesh) * 4.0 * PI * e2;

    //value at the end
    assert(r[mesh-1] > tiny);
    pot[mesh-1] = Integral_G::Gauss_Legendre(a, b, rad_f1, r, mesh) * 4.0 * PI / r[mesh-1] * e2;
	
	//points in the interval
    for(int ir = 1; ir < mesh-1; ir++)
    {
        a = r[0];
        b = r[ir];
        c = r[mesh-1];

        //integrate inside
        const double inside = Integral_G::Gauss_Legendre(a, b, rad_f1, r, mesh) / r[ir];

        //integrate outside
        const double outside = Integral_G::Gauss_Legendre(b, c, rad_f2, r, mesh);

        //inside + outside
        pot[ir] = (inside + outside) * 4.0 * PI * e2;
    }

    delete[] rad_f1;
    delete[] rad_f2;

    return;
}
	
