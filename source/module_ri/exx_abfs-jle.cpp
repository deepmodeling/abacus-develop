#include "exx_abfs-jle.h"

#include "../module_hamilt_pw/hamilt_pwdft/global.h"
#include "../module_basis/module_ao/ORB_read.h"
#include "../module_base/global_function.h"
#include "../module_base/mathzone.h"
#include "../module_base/math_sphbes.h" // mohan add 2021-05-06

bool Exx_Abfs::Jle::generate_matrix = false;
int Exx_Abfs::Jle::Lmax = 2;
double Exx_Abfs::Jle::Ecut_exx = 60;
double Exx_Abfs::Jle::tolerence = 1.0e-12;	

void Exx_Abfs::Jle::init_jle( const double kmesh_times )
{
	jle.resize( GlobalC::ucell.ntype );

	for (int T = 0;  T < GlobalC::ucell.ntype ; T++)
	{
		jle[T].resize( Lmax+1 );
		for (int L=0; L <= Lmax ; ++L)
		{
			const size_t ecut_number 
				= static_cast<size_t>( sqrt( Ecut_exx ) * GlobalC::ORB.Phi[T].getRcut() / ModuleBase::PI ); // Rydberg Unit.

			jle[T][L].resize( ecut_number );

			std::vector<double> en(ecut_number, 0.0);
			ModuleBase::Sphbes::Spherical_Bessel_Roots(ecut_number, L, tolerence, ModuleBase::GlobalFunc::VECTOR_TO_PTR(en), GlobalC::ORB.Phi[T].getRcut());

			for(size_t E=0; E!=ecut_number; ++E)
			{
				std::vector<double> jle_r( GlobalC::ORB.Phi[T].PhiLN(0,0).getNr() );
				ModuleBase::Sphbes::Spherical_Bessel(
					GlobalC::ORB.Phi[T].PhiLN(0,0).getNr(), 
					GlobalC::ORB.Phi[T].PhiLN(0,0).getRadial(), 
					en[E], 
					L, 
					ModuleBase::GlobalFunc::VECTOR_TO_PTR(jle_r));
				jle[T][L][E].set_orbital_info(
					GlobalC::ORB.Phi[T].PhiLN(0,0).getLabel(),
					GlobalC::ORB.Phi[T].PhiLN(0,0).getType(),
					L,
					E,						// N?
					GlobalC::ORB.Phi[T].PhiLN(0,0).getNr(),
					GlobalC::ORB.Phi[T].PhiLN(0,0).getRab(),
					GlobalC::ORB.Phi[T].PhiLN(0,0).getRadial(),
					Numerical_Orbital_Lm::Psi_Type::Psi,
					ModuleBase::GlobalFunc::VECTOR_TO_PTR(jle_r),
					static_cast<int>(GlobalC::ORB.Phi[T].PhiLN(0,0).getNk() * kmesh_times) | 1,
					GlobalC::ORB.Phi[T].PhiLN(0,0).getDk(),
					GlobalC::ORB.Phi[T].PhiLN(0,0).getDruniform(),
					false,
					true, GlobalV::CAL_FORCE);
			}
		}
	}
}
