#include "read_stru.h"
#include "source_base/timer.h"
#include "source_base/vector3.h"
#include "source_base/mathzone.h"

bool unitcell::check_tau(const Atom* atoms,
		const int& ntype,
		const double& lat0) 
{
	ModuleBase::TITLE("UnitCell","check_tau");
	ModuleBase::timer::tick("UnitCell","check_tau");

	ModuleBase::Vector3<double> diff = 0.0;
	double norm = 0.0;
	double tolerence_bohr = 1.0e-3;

	for(int T1=0; T1< ntype; T1++)
	{
		for(int I1=0; I1< atoms[T1].na; I1++)
		{    
			double shortest_norm = 10000.0; // a large number
			for(int T2=0; T2<ntype; T2++)
			{
				for(int I2=0; I2<atoms[T2].na; I2++)
				{
					if(T1==T2 && I1==I2)
					{
						shortest_norm = 0.0;
					}
					else
					{
						diff = atoms[T1].tau[I1] - atoms[T2].tau[I2];
						norm = diff.norm() * lat0;
						if( shortest_norm > norm )
						{
							shortest_norm = norm;
						}
						if( norm < tolerence_bohr ) // unit is Bohr
						{    
							GlobalV::ofs_warning << " two atoms are too close!" << std::endl;
							GlobalV::ofs_warning << " type:" << atoms[T1].label << " atom " << I1 + 1 << std::endl; 
							GlobalV::ofs_warning << " type:" << atoms[T2].label << " atom " << I2 + 1 << std::endl; 
							GlobalV::ofs_warning << " distance = " << norm << " Bohr" << std::endl;
							return false;
						}
					}
				}
			}
		}
	}
	ModuleBase::timer::tick("UnitCell","check_tau");
	return true;
}


void unitcell::check_dtau(Atom* atoms,
		const int& ntype,
		const double& lat0,
		ModuleBase::Matrix3& latvec)
{
	for(int it=0; it<ntype; it++)
	{
		Atom* atom1 = &atoms[it];
		for(int ia=0; ia<atoms[it].na; ia++)
		{
			// mohan add 2011-04-07            
			// fmod(x,1.0) set the result between the [0,1.0),
			// while the x may be the negtivate value,thus we add 10000.
			atom1->taud[ia].x=fmod(atom1->taud[ia].x + 10000,1.0);
			atom1->taud[ia].y=fmod(atom1->taud[ia].y + 10000,1.0);
			atom1->taud[ia].z=fmod(atom1->taud[ia].z + 10000,1.0);

			double cx2=0.0;
			double cy2=0.0;
			double cz2=0.0;

			ModuleBase::Mathzone::Direct_to_Cartesian(
					atom1->taud[ia].x, atom1->taud[ia].y, atom1->taud[ia].z,
					latvec.e11, latvec.e12, latvec.e13,
					latvec.e21, latvec.e22, latvec.e23,
					latvec.e31, latvec.e32, latvec.e33,
					cx2, cy2, cz2);

			atom1->tau[ia].x = cx2;
			atom1->tau[ia].y = cy2;
			atom1->tau[ia].z = cz2;

		}
	}
	return;
}
