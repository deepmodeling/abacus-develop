#include "symmetry.h"
using namespace ModuleSymmetry;

#include "module_parameter/parameter.h"

//---------------------------------------------------
// The lattice will be transformed to a 'standard
// cystallographic setting', the relation between
// 'origin' and 'transformed' lattice vectors will
// be givin in matrix form
//---------------------------------------------------
int Symmetry::standard_lat(
    ModuleBase::Vector3<double> &a,
    ModuleBase::Vector3<double> &b,
    ModuleBase::Vector3<double> &c,
    double *cel_const) const
{
    static bool first = true;
    // there are only 14 types of Bravais lattice.
    int type = 15;
	//----------------------------------------------------
    // used to calculte the volume to judge whether 
	// the lattice vectors corrispond the right-hand-sense
	//----------------------------------------------------
    double volume = 0;
    //the lattice vectors have not been changed

    const double aa = a * a;
    const double bb = b * b;
    const double cc = c * c;
    const double ab = a * b; //std::vector: a * b * cos(alpha)
    const double bc = b * c; //std::vector: b * c * cos(beta)
    const double ca = c * a; //std::vector: c * a * cos(gamma)
    double norm_a = a.norm();
    double norm_b = b.norm();
    double norm_c = c.norm();
    double gamma = ab /( norm_a * norm_b ); // cos(gamma)
    double alpha  = bc /( norm_b * norm_c ); // cos(alpha)
    double beta = ca /( norm_a * norm_c ); // cos(beta)
    double amb = sqrt( aa + bb - 2 * ab );	//amb = |a - b|
    double bmc = sqrt( bb + cc - 2 * bc );
    double cma = sqrt( cc + aa - 2 * ca );
    double apb = sqrt( aa + bb + 2 * ab );  //amb = |a + b|
    double bpc = sqrt( bb + cc + 2 * bc );
    double cpa = sqrt( cc + aa + 2 * ca );
    double apbmc = sqrt( aa + bb + cc + 2 * ab - 2 * bc - 2 * ca );	//apbmc = |a + b - c|
    double bpcma = sqrt( bb + cc + aa + 2 * bc - 2 * ca - 2 * ab );
    double cpamb = sqrt( cc + aa + bb + 2 * ca - 2 * ab - 2 * bc );
    double abc = ab + bc + ca;

	if (first)
    {
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"NORM_A",norm_a);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"NORM_B",norm_b);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"NORM_C",norm_c);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"ALPHA (DEGREE)", acos(alpha)/ModuleBase::PI*180.0 );
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"BETA  (DEGREE)" ,acos(beta)/ModuleBase::PI*180.0  );
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"GAMMA (DEGREE)" ,acos(gamma)/ModuleBase::PI*180.0 );
        first = false;
    }

    Symm_Other::right_hand_sense(a, b, c);
	ModuleBase::GlobalFunc::ZEROS(cel_const, 6);
	const double small = PARAM.inp.symmetry_prec;

	//---------------------------	
	// 1. alpha == beta == gamma 
	//---------------------------	
	if( equal(alpha, gamma) && equal(alpha, beta) )
	{
		//--------------
		// a == b == c 
		//--------------
		if( equal(norm_a, norm_b) && equal(norm_b, norm_c))
		{
			//---------------------------------------
			// alpha == beta == gamma == 90 degree
			//---------------------------------------
			if ( equal(alpha,0.0) )
			{
				type=1;
				cel_const[0]=norm_a;
			}
			//----------------------------------------
			// cos(alpha) = -1.0/3.0
			//----------------------------------------
			else if( equal(alpha, -1.0/3.0) ) 
			{
				type=2;
				cel_const[0]=norm_a*2.0/sqrt(3.0);
			}
			//----------------------------------------
			// cos(alpha) = 0.5
			//----------------------------------------
			else if( equal(alpha, 0.5) ) 
			{
				type=3;
				cel_const[0]=norm_a*sqrt(2.0);
			}
			//----------------------------------------
			// cos(alpha) = all the others
			//----------------------------------------
			else
			{
				type=7;
				cel_const[0]=norm_a;
				cel_const[3]=alpha;
			}
		}
		// Crystal classes with inequal length of lattice vectors but also with
		// A1*A2=A1*A3=A2*A3:
		// Orthogonal axes:
		else if(equal(gamma,0.0)) 
		{
			// Two axes with equal lengths means simple tetragonal: (IBRAV=5)
			// Adjustment: 'c-axis' shall be the special axis.
			if (equal(norm_a, norm_b)) 
			{
				type=5;
				cel_const[0]=norm_a;
				cel_const[2]=norm_c/norm_a;
				// No axes with equal lengths means simple orthorhombic (IBRAV=8):
				// Adjustment: Sort the axis by increasing lengths:
			}
            else if(((norm_c-norm_b)>small) && ((norm_b-norm_a)>small) ) 
			{
				type=8;
				cel_const[0]=norm_a;
				cel_const[1]=norm_b/norm_a;
				cel_const[2]=norm_c/norm_a;
			}
			// Crystal classes with A1*A3=A2*A3=/A1*A2:
		}
	}//end alpha=beta=gamma
	//-----------------------
	// TWO EQUAL ANGLES
	// alpha == beta != gamma  (gamma is special)
	//------------------------
	else if (equal(alpha-beta, 0)) 
	{
		//---------------------------------------------------------
		// alpha = beta = 90 degree
		// One axis orthogonal with respect to the other two axes:
		//---------------------------------------------------------
		if (equal(alpha, 0.0)) 
		{
			//-----------------------------------------------
			// a == b 
			// Equal length of the two nonorthogonal axes:
			//-----------------------------------------------
			if (equal(norm_a, norm_b)) 
			{
				// Cosine(alpha) equal to -1/2 means hexagonal: (IBRAV=4)
				// Adjustment: 'c-axis' shall be the special axis.
				if ( equal(gamma, -0.5))   //gamma = 120 degree
				{
					type=4;
					cel_const[0]=norm_a;
					cel_const[2]=norm_c/norm_a;
					// Other angles mean base-centered orthorhombic: (IBRAV=11)
					// Adjustment: Cosine between A1 and A2 shall be lower than zero, the
					//             'c-axis' shall be the special axis.
				}
				else if(gamma<(-1.0*small)) //gamma > 90 degree
				{
					type=11;
                    cel_const[0]=apb;
                    cel_const[1]=amb/apb;
                    cel_const[2]=norm_c/apb;
                    cel_const[5]=gamma;
				}
				// Different length of the two axes means simple monoclinic (IBRAV=12):
				// Adjustment: Cosine(gamma) should be lower than zero, special axis
				//             shall be the 'b-axis'(!!!) and |A1|<|A3|:
			}
			//----------
			// a!=b!=c
			//----------
            else if( gamma<(-1.0*small) && (norm_a-norm_b)>small) 
			{
				type=12;
				cel_const[0]=norm_b;
				cel_const[1]=norm_c/norm_b;
				cel_const[2]=norm_a/norm_b;
                cel_const[4]=gamma;
                //adjust: a->c, b->a, c->b
                ModuleBase::Vector3<double> tmp=c;
				c=a;
				a=b;
				b=tmp;
			}
		}//end gamma<small
		// Arbitrary angles between the axes:
		// |A1|=|A2|=|A3| means body-centered tetragonal (IBRAV=6):
		// Further additional criterions are: (A1+A2), (A1+A3) and (A2+A3) are
		// orthogonal to one another and (adjustment//): |A1+A3|=|A2+A3|/=|A1+A2|
		else
		{
			if( equal(norm_a, norm_b) && 
				equal(norm_b, norm_c) &&
				equal(cpa, bpc) && 
				!equal(apb, cpa) &&
				equal(norm_c*norm_c+abc,0) )
			{
				type=6;
				cel_const[0]=cpa;
				cel_const[2]=apb/cpa;
			}
			// |A1|=|A2|=/|A3| means base-centered monoclinic (IBRAV=13):
			// Adjustement: The cosine between A1 and A3 as well as the cosine
			//              between A2 and A3 should be lower than zero.
			else if( equal(norm_a,norm_b) 
					&& alpha<(-1.0*small) 
					&& beta<(-1.0*small)) 
			{
				type=13;
				cel_const[0]=apb;
				cel_const[1]=amb/apb;
				cel_const[2]=norm_c/apb;
                //cos(<a+b, c>)
                cel_const[4]=(a+b)*c/apb/norm_c;
			}
		}
	} //end alpha==beta
	//-------------------------------
	// three angles are not equal
	//-------------------------------
	else 
	{
		// Crystal classes with A1*A2=/A1*A3=/A2*A3
		// |A1|=|A2|=|A3| means body-centered orthorhombic (IBRAV=9):
		// Further additional criterions are: (A1+A2), (A1+A3) and (A2+A3) are
		// orthogonal to one another and (adjustment//): |A1+A2|>|A1+A3|>|A2+A3|
		if (equal(norm_a, norm_b) &&
				equal(norm_b, norm_c) &&
				((cpa-bpc)>small) &&
				((apb-cpa)>small) && 
				equal(norm_c*norm_c+abc, 0)) 
		{
			type=9;
			cel_const[0]=bpc;
			cel_const[1]=cpa/bpc;
			cel_const[2]=apb/bpc;
		}
		// |A1|=|A2-A3| and |A2|=|A1-A3| and |A3|=|A1-A2| means face-centered
		// orthorhombic (IBRAV=10):
		// Adjustment: |A1+A2-A3|>|A1+A3-A2|>|A2+A3-A1|
		else if(equal(amb, norm_c) &&
				equal(cma, norm_b) &&
				equal(bmc, norm_a) && 
				((apbmc-cpamb)>small) &&
				((cpamb-bpcma)>small)) 
		{
			type=10;
			cel_const[0]=bpcma;
			cel_const[1]=cpamb/bpcma;
			cel_const[2]=apbmc/bpcma;
		}
		// Now there exists only one further possibility - triclinic (IBRAV=14):
		// Adjustment: All three cosines shall be greater than zero and ordered:
		else if((gamma>beta) && (beta>alpha) && (alpha>small)) 
		{
			type=14;
			cel_const[0]=norm_a;
			cel_const[1]=norm_b/norm_a;
			cel_const[2]=norm_c/norm_a;
			cel_const[3]=alpha;
			cel_const[4]=beta;
			cel_const[5]=gamma;
		}
	}
	
	return type;
}

//---------------------------------------------------
// The lattice will be transformed to a 'standard
// cystallographic setting', the relation between
// 'origin' and 'transformed' lattice vectors will
// be givin in matrix form
// must be called before symmetry analysis
// only need to called once for each ion step
//---------------------------------------------------
void Symmetry::lattice_type(
    ModuleBase::Vector3<double> &v1,
    ModuleBase::Vector3<double> &v2,
    ModuleBase::Vector3<double> &v3,
    ModuleBase::Vector3<double> &v01,
    ModuleBase::Vector3<double> &v02,
    ModuleBase::Vector3<double> &v03,
    double *cel_const,
    double *pre_const,
    int& real_brav,
    std::string& bravname,
    const Atom* atoms,
    bool convert_atoms,
    double* newpos)const
{
    ModuleBase::TITLE("Symmetry","lattice_type");

	//----------------------------------------------
	// (1) adjustement of the basis to right hand 
	// sense by inversion of all three lattice 
	// vectors if necessary
	//----------------------------------------------
    const bool right = Symm_Other::right_hand_sense(v1, v2, v3);

	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"RIGHT HAND LATTICE",right);

	//-------------------------------------------------
	// (2) save and copy the original lattice vectors.
	//-------------------------------------------------
    v01 = v1;
    v02 = v2;
    v03 = v3;
	
	//--------------------------------------------
	// (3) calculate the 'pre_const'
	//--------------------------------------------
	ModuleBase::GlobalFunc::ZEROS(pre_const, 6);

    int pre_brav = standard_lat(v1, v2, v3, cel_const);

    for ( int i = 0; i < 6; ++i)
    {
        pre_const[i] = cel_const[i];
    }

    // find the shortest basis vectors of the lattice
    this->get_shortest_latvec(v1, v2, v3);

    Symm_Other::right_hand_sense(v1, v2, v3);

    real_brav = 15;
    double temp_const[6];

    //then we should find the best lattice vectors to make much easier the determination of the lattice symmetry
    //the method is to contrast the combination of the shortest vectors and determine their symmmetry

    ModuleBase::Vector3<double> w1, w2, w3;
    ModuleBase::Vector3<double> q1, q2, q3;
    this->get_optlat(v1, v2, v3, w1, w2, w3, real_brav, cel_const, temp_const);

    //now, the highest symmetry of the combination of the shortest vectors has been found
    //then we compare it with the original symmetry
	
    bool change_flag=false;
    for (int i = 0; i < 6; ++i) {
        if(!equal(cel_const[i], pre_const[i])) 
            {change_flag=true; break;
        }
    }

    if ( real_brav < pre_brav || change_flag )
    {
        //if the symmetry of the new vectors is higher, store the new ones
        for (int i = 0; i < 6; ++i)
        {
            cel_const[i] = temp_const[i];
        }
        q1 = w1;
        q2 = w2;
        q3 = w3;
        if(convert_atoms)
        {
            GlobalV::ofs_running <<std::endl;
            GlobalV::ofs_running <<" The lattice vectors have been changed (STRU_SIMPLE.cif)"<<std::endl;
            GlobalV::ofs_running <<std::endl;
            int at=0;
            for (int it = 0; it < this->ntype; ++it)
            {
				for (int ia = 0; ia < this->na[it]; ++ia)
				{
					ModuleBase::Mathzone::Cartesian_to_Direct(atoms[it].tau[ia].x,
							atoms[it].tau[ia].y,
							atoms[it].tau[ia].z,
							q1.x, q1.y, q1.z,
							q2.x, q2.y, q2.z,
							q3.x, q3.y, q3.z,
							newpos[3*at],newpos[3*at+1],newpos[3*at+2]);

					for(int k=0; k<3; ++k)
					{
						this->check_translation( newpos[at*3+k], -floor(newpos[at*3+k]));
						this->check_boundary( newpos[at*3+k] );
					}
					++at;
				}
			}       
        }
        // return the optimized lattice in v1, v2, v3
        v1=q1;
        v2=q2;
        v3=q3;
    }
    else
    {
        //else, store the original ones
        for (int i = 0; i < 6; ++i)
        {
            cel_const[i] = pre_const[i];
        }
        //newpos also need to be set
        if(convert_atoms)
        {
            int at=0;
            for (int it = 0; it < this->ntype; ++it)
            {
                for (int ia = 0; ia < this->na[it]; ++ia)
                {
                    ModuleBase::Mathzone::Cartesian_to_Direct(atoms[it].tau[ia].x,
                        atoms[it].tau[ia].y,
                        atoms[it].tau[ia].z,
                                    v1.x, v1.y, v1.z,
                                    v2.x, v2.y, v2.z,
                                    v3.x, v3.y, v3.z,
                                    newpos[3*at],newpos[3*at+1],newpos[3*at+2]);
                    for(int k=0; k<3; ++k)
                    {
                            this->check_translation( newpos[at*3+k], -floor(newpos[at*3+k]));
                            this->check_boundary( newpos[at*3+k] );
                    }
                    ++at;
                }
            }       
        }
    }

    /*
    bool flag3;
    if (pre_brav == temp_brav) 
	{
        flag3 = 0;
        if (!equal(temp_const[0], pre_const[0]) ||
            !equal(temp_const[1], pre_const[1]) ||
            !equal(temp_const[2], pre_const[2]) ||
            !equal(temp_const[3], pre_const[3]) ||
            !equal(temp_const[4], pre_const[4]) ||
            !equal(temp_const[5], pre_const[5])
           )
        {
            flag3 = 1;
        }
        if (flag3==0) {
            //particularly, if the symmetry of origin and new are exactly the same, we choose the original ones
            //Hey! the original vectors have been changed!!!
            v1 = s1;
            v2 = s2;
            v3 = s3;
        	change=0;
			GlobalV::ofs_running<<" The lattice vectors have been set back!"<<std::endl;
        }
    }*/
    bravname = get_brav_name(real_brav);
    return;
}
