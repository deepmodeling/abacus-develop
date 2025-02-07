#include "read_stru.h"
#include "module_base/timer.h"
#include "module_base/vector3.h"
#include "module_base/mathzone.h"


namespace unitcell
{
    bool check_tau(const Atom* atoms,
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

    void check_dtau(Atom* atoms,
                    const int& ntype,
                    const double& lat0,
                    ModuleBase::Matrix3& latvec)
    {
        for(int it=0; it<ntype; it++)
        {
            Atom* atom1 = &atoms[it];
            for(int ia=0; ia<atoms[it].na; ia++)
            {
                //fmod(x,1.0) set the result as [0,1),
                //if x=2.3,fmod(2.3,1.0)=0.3,fmod(2.3,1.0)+1=1.3,
                // fmod(fmod(2.3,1.0)+1,1.0)=0.3
                //if x=-0.7,fmod(-0.7,1.0)=-0.7 or 0.3,fmod(-0.7,1.0)+1=1.3,
                // fmod(fmod(-0.7,1.0)+1,1.0)=0.3
                double dx2 = fmod(fmod(atom1->taud[ia].x,1.0)+1.0,1.0);
                double dy2 = fmod(fmod(atom1->taud[ia].y,1.0)+1.0,1.0);
                double dz2 = fmod(fmod(atom1->taud[ia].z,1.0)+1.0,1.0);

                atom1->taud[ia].x = dx2;
                atom1->taud[ia].y = dy2;
                atom1->taud[ia].z = dz2;

                double cx2=0.0;
                double cy2=0.0;
                double cz2=0.0;

                ModuleBase::Mathzone::Direct_to_Cartesian(
                atom1->taud[ia].x, 
                atom1->taud[ia].y, 
                atom1->taud[ia].z,
                latvec.e11, 
                latvec.e12, 
                latvec.e13,
                latvec.e21, 
                latvec.e22, 
                latvec.e23,
                latvec.e31, 
                latvec.e32, 
                latvec.e33,
                cx2, 
                cy2, 
                cz2);

                atom1->tau[ia].x = cx2;
                atom1->tau[ia].y = cy2;
                atom1->tau[ia].z = cz2;
                
            }
        }

    }

}