#include "exx_abfs-jle.h"

#include "source_io/module_parameter/parameter.h"
#include "../../source_basis/module_ao/ORB_read.h"
#include "../../source_cell/unitcell.h"
#include "../../source_base/mathzone.h"
#include "../../source_base/math_sphbes.h" // mohan add 2021-05-06
#include "source_base/tool_title.h"

std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>
Exx_Abfs::Jle::init_jle(
	const Exx_Info::Exx_Info_Opt_ABFs &info,
	const double kmesh_times, 
	const UnitCell& ucell,
	const LCAO_Orbitals& orb)
{
	ModuleBase::TITLE("Exx_Abfs::Jle","init_jle");
	std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> jle( ucell.ntype );

	for(int T=0; T<ucell.ntype; ++T)
	{
		if(info.abfs_Lmax+1>0)
			{ jle[T].resize( info.abfs_Lmax+1 ); }
		for(int L=0; L<=info.abfs_Lmax; ++L)
		{
			const size_t ecut_number 
				= static_cast<size_t>( std::sqrt( info.ecut_exx ) * orb.Phi[T].getRcut() / ModuleBase::PI ); // Rydberg Unit.
			if(ecut_number<=0)
				{ continue; }

			jle[T][L].resize( ecut_number );

			std::vector<double> en(ecut_number, 0.0);
			ModuleBase::Sphbes::Spherical_Bessel_Roots(ecut_number, L, info.tolerence, en.data(), orb.Phi[T].getRcut());

			for(size_t E=0; E<ecut_number; ++E)
			{
				std::vector<double> jle_r( orb.Phi[T].PhiLN(0,0).getNr() );
				ModuleBase::Sphbes::Spherical_Bessel(
					orb.Phi[T].PhiLN(0,0).getNr(), 
					orb.Phi[T].PhiLN(0,0).getRadial(), 
					en[E], 
					L, 
					jle_r.data());
				jle[T][L][E].set_orbital_info(
					orb.Phi[T].PhiLN(0,0).getLabel(),
					orb.Phi[T].PhiLN(0,0).getType(),
					L,
					E,						// N?
					orb.Phi[T].PhiLN(0,0).getNr(),
					orb.Phi[T].PhiLN(0,0).getRab(),
					orb.Phi[T].PhiLN(0,0).getRadial(),
					Numerical_Orbital_Lm::Psi_Type::Psi,
					jle_r.data(),
					static_cast<int>(orb.Phi[T].PhiLN(0,0).getNk() * kmesh_times) | 1,
					orb.Phi[T].PhiLN(0,0).getDk(),
					orb.Phi[T].PhiLN(0,0).getDruniform(),
					false,
					true,
					PARAM.inp.cal_force);
			}
		}
	}

	return jle;
}
