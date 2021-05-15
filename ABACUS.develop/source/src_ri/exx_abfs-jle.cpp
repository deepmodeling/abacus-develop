#include "exx_abfs-jle.h"

#include "../src_pw/global.h"
#include "module_ORB/ORB_read.h"
#include "../src_global/global_function.h"
#include "../src_global/mathzone.h"
#include "../src_global/math_sphbes.h" // mohan add 2021-05-06

int Exx_Abfs::Jle::Lmax = 2;
double Exx_Abfs::Jle::Ecut_exx = 60;
double Exx_Abfs::Jle::tolerence = 1.0e-12;	

void Exx_Abfs::Jle::init_jle( const double kmesh_times )
{
	jle.resize( ucell.ntype );

	for (int T = 0;  T < ucell.ntype ; T++)
	{
		jle[T].resize( Lmax+1 );
		for (int L=0; L <= Lmax ; ++L)
		{
			const size_t ecut_number 
				= static_cast<size_t>( sqrt( Ecut_exx ) * ORB.Phi[T].getRcut() / PI ); // Rydberg Unit.

			jle[T][L].resize( ecut_number );

			vector<double> en(ecut_number, 0.0);
			Sphbes::Spherical_Bessel_Roots(ecut_number, L, tolerence, VECTOR_TO_PTR(en), ORB.Phi[T].getRcut());

			for(size_t E=0; E!=ecut_number; ++E)
			{
				vector<double> jle_r( ORB.Phi[T].PhiLN(0,0).getNr() );
				Sphbes::Spherical_Bessel(
					ORB.Phi[T].PhiLN(0,0).getNr(), 
					ORB.Phi[T].PhiLN(0,0).getRadial(), 
					en[E], 
					L, 
					VECTOR_TO_PTR(jle_r));
				jle[T][L][E].set_orbital_info(
					ORB.Phi[T].PhiLN(0,0).getLabel(),
					ORB.Phi[T].PhiLN(0,0).getType(),
					L,
					E,						// N?
					ORB.Phi[T].PhiLN(0,0).getNr(),
					ORB.Phi[T].PhiLN(0,0).getRab(),
					ORB.Phi[T].PhiLN(0,0).getRadial(),
					Numerical_Orbital_Lm::Psi_Type::Psi,
					VECTOR_TO_PTR(jle_r),
					static_cast<int>(ORB.Phi[T].PhiLN(0,0).getNk() * kmesh_times) | 1,
					ORB.Phi[T].PhiLN(0,0).getDk(),
					ORB.Phi[T].PhiLN(0,0).getDruniform(),
					false,
					true, FORCE);
			}
		}
	}
}
