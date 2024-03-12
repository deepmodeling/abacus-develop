#ifndef YLM_REAL_H
#define YLM_REAL_H
#include "../src_tools/vector3.h"
#include "../src_tools/realarray.h"
#include "../src_tools/matrix.h"
using namespace std;

void Ylm_Real
(
	const int lmax2, 			// lmax2 = (lmax+1)^2
	const int ng,				//
	const Vector3<double> *g, 	// g_cartesian_vec(x,y,z)
	matrix &ylm 				// output
);

long double Fact(const int n);
int Semi_Fact(const int n);

#endif
