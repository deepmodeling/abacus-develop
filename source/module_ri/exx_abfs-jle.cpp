#include "exx_abfs-jle.h"

#include "module_parameter/parameter.h"
#include "../module_hamilt_pw/hamilt_pwdft/global.h"
#include "../module_basis/module_ao/ORB_read.h"
#include "../module_base/global_function.h"
#include "../module_base/mathzone.h"
#include "../module_base/math_sphbes.h" // mohan add 2021-05-06

bool Exx_Abfs::Jle::generate_matrix = false;
int Exx_Abfs::Jle::Lmax = 2;
double Exx_Abfs::Jle::Ecut_exx = 60;
double Exx_Abfs::Jle::tolerence = 1.0e-12;	

void Exx_Abfs::Jle::init_jle(const double kmesh_times, 
							 const UnitCell& ucell,
							 const LCAO_Orbitals& orb)
{
	jle.resize( ucell.ntype );

	for (int T = 0;  T < ucell.ntype ; T++)
	{
		jle[T].resize( Lmax+1 );
		for (int L=0; L <= Lmax ; ++L)
		{
			const size_t ecut_number 
				= static_cast<size_t>( sqrt( Ecut_exx ) * orb.Phi[T].getRcut() / ModuleBase::PI ); // Rydberg Unit.

			jle[T][L].resize( ecut_number );

			std::vector<double> en(ecut_number, 0.0);
			ModuleBase::Sphbes::Spherical_Bessel_Roots(ecut_number, L, tolerence, ModuleBase::GlobalFunc::VECTOR_TO_PTR(en), orb.Phi[T].getRcut());

			for(size_t E=0; E!=ecut_number; ++E)
			{
				std::vector<double> jle_r( orb.Phi[T].PhiLN(0,0).getNr() );
				ModuleBase::Sphbes::Spherical_Bessel(
					orb.Phi[T].PhiLN(0,0).getNr(), 
					orb.Phi[T].PhiLN(0,0).getRadial(), 
					en[E], 
					L, 
					ModuleBase::GlobalFunc::VECTOR_TO_PTR(jle_r));
				jle[T][L][E].set_orbital_info(
					orb.Phi[T].PhiLN(0,0).getLabel(),
					orb.Phi[T].PhiLN(0,0).getType(),
					L,
					E,						// N?
					orb.Phi[T].PhiLN(0,0).getNr(),
					orb.Phi[T].PhiLN(0,0).getRab(),
					orb.Phi[T].PhiLN(0,0).getRadial(),
					Numerical_Orbital_Lm::Psi_Type::Psi,
					ModuleBase::GlobalFunc::VECTOR_TO_PTR(jle_r),
					static_cast<int>(orb.Phi[T].PhiLN(0,0).getNk() * kmesh_times) | 1,
					orb.Phi[T].PhiLN(0,0).getDk(),
					orb.Phi[T].PhiLN(0,0).getDruniform(),
					false,
					true, PARAM.inp.cal_force);
			}
		}
	}
}
