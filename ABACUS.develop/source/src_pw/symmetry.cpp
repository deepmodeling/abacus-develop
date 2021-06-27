#include "symmetry.h"
#include "global.h"
#include "symm_other.h"

Symmetry::Symmetry()
{
    this->epsilon = 1e-6;
    this->tab = 12;
    this->available = true;
}

Symmetry::~Symmetry()
{
}


bool Symmetry::symm_flag=false;


void Symmetry::analy_sys(void)
{
    if (available == false) return;
    TITLE("Symmetry","init");
	timer::tick("Symmetry","analy_sys");

	ofs_running << "\n\n\n\n";
	ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
	ofs_running << " |                                                                    |" << endl;
	ofs_running << " | Doing symmetry analysis:                                           |" << endl;
	ofs_running << " | We calculate the norm of 3 vectors and the angles between them,    |" << endl;
	ofs_running << " | the type of Bravais lattice is given. We can judge if the unticell |" << endl;
	ofs_running << " | is a primitive cell. Finally we give the point group operation for |" << endl;
	ofs_running << " | this unitcell. We we use the point group operations to do symmetry |" << endl;
	ofs_running << " | analysis on given k-point mesh and the charge density.             |" << endl;
	ofs_running << " |                                                                    |" << endl;
	ofs_running << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
	ofs_running << "\n\n\n\n";


    this->ibrav = 0;
    // number of total atoms
    this->nat = ucell.nat;
    // number of atom species
    this->ntype = ucell.ntype;
    this->na = new int[ntype];
    this->istart = new int[ntype];
    this->ptrans = new double[(nat + 2) * 3];
    this->index = new int [nat + 2];
    ZEROS(na, ntype);
    ZEROS(istart, ntype);
    ZEROS(ptrans, (nat+2)*3);
    ZEROS(index, nat+2);

    // atom positions
    // used in checksym.
    dirpos = new double[3*nat];
	newpos = new double[3*nat];
    rotpos = new double[3*nat];
    ZEROS(dirpos, 3*nat);
	ZEROS(newpos, 3*nat);
    ZEROS(rotpos, 3*nat);

    this->a1 = ucell.a1;
    this->a2 = ucell.a2;
    this->a3 = ucell.a3;

	Matrix3 latvec1;
	latvec1.e11 = a1.x; latvec1.e12 = a1.y; latvec1.e13 = a1.z;
	latvec1.e21 = a2.x; latvec1.e22 = a2.y; latvec1.e23 = a2.z;
	latvec1.e31 = a3.x; latvec1.e32 = a3.y; latvec1.e33 = a3.z;
//  cout << "a1 = " << a1.x << " " << a1.y << " " << a1.z <<endl;
//  cout << "a1 = " << a2.x << " " << a2.y << " " << a2.z <<endl;
//  cout << "a1 = " << a3.x << " " << a3.y << " " << a3.z <<endl;

	out.printM3(ofs_running,"LATTICE VECTORS: (CARTESIAN COORDINATE: IN UNIT OF A0)",latvec1);

    int count = 0;
    istart[0] = 0;
    this->itmin_type = 0;
    this->itmin_start = 0;
    for (int it = 0; it < ntype; ++it)
    {
		Atom* atom = &ucell.atoms[it];
        this->na[it] = atom->na;
        if (it > 0) {
            istart[it] = istart[it-1] + na[it-1];
        }
        //cout << "\n istart = " << istart[it];
        if (na[it] < na[itmin_type])
        {
            this->itmin_type = it;
            this->itmin_start = istart[it];
        }

        Vector3<double> vec;
        for (int ia = 0; ia < ucell.atoms[it].na; ++ia)
        {
            dirpos[3*count + 0] = atom->taud[ia].x;
            dirpos[3*count + 1] = atom->taud[ia].y;
            dirpos[3*count + 2] = atom->taud[ia].z;
//            cout << " atom.taud = " << atom->taud[ia].x << " "<<atom->taud[ia].y<<" "<<atom->taud[ia].z<<endl;
            ++count;
        }
    }

    s1 = a1;
    s2 = a2;
    s3 = a3;

    // find the lattice type accordiing to lattice vectors.
    this->lattice_type(a1,a2,a3,ibrav,cel_const,ilattname);
  //      cout << "a1 = " << a1.x << " " << a1.y << " " << a1.z <<endl;
  //      cout << "a1 = " << a2.x << " " << a2.y << " " << a2.z <<endl;
  //      cout << "a1 = " << a3.x << " " << a3.y << " " << a3.z <<endl;


	// the atom position coordinates are changed to 
	// crystal coordinates of a1,a2,a3
	Matrix3 new_lat;
	new_lat.e11=a1.x; new_lat.e12=a1.y; new_lat.e13=a1.z;
	new_lat.e21=a2.x; new_lat.e22=a2.y; new_lat.e23=a2.z;
	new_lat.e31=a3.x; new_lat.e32=a3.y; new_lat.e33=a3.z;
	out.printM3(ofs_running,"STANDARD LATTICE VECTORS: (CARTESIAN COORDINATE: IN UNIT OF A0)",new_lat);

	int iat=0;
	for(int it=0; it<ucell.ntype; ++it)
	{
		for(int ia=0; ia<ucell.atoms[it].na; ++ia)
		{
			Mathzone::Cartesian_to_Direct(ucell.atoms[it].tau[ia].x, 
					ucell.atoms[it].tau[ia].y, 
					ucell.atoms[it].tau[ia].z,
					new_lat.e11, new_lat.e12, new_lat.e13,
					new_lat.e21, new_lat.e22, new_lat.e23,
					new_lat.e31, new_lat.e32, new_lat.e33,
					newpos[3*iat],newpos[3*iat+1],newpos[3*iat+2]);

    	//  cout << " newpos_before = " << newpos[3*iat] << " " << newpos[3*iat+1] << " " << newpos[3*iat+2] << endl;
		//	ofs_running << " newpos_before = " << newpos[3*iat] << " " << newpos[3*iat+1] << " " << newpos[3*iat+2] << endl; 
			for(int k=0; k<3; ++k)
			{
				this->check_translation( newpos[iat*3+k], -floor(newpos[iat*3+k]));
                         	this->check_boundary( this->newpos[iat*3+k] );
			}
      	// cout << " newpos = " << newpos[3*iat] << " " << newpos[3*iat+1] << " " << newpos[3*iat+2] << endl;
		// ofs_running << " newpos = " << newpos[3*iat] << " " << newpos[3*iat+1] << " " << newpos[3*iat+2] << endl; 
			++iat;
		}
	}


	Symm_Other::print1(ibrav, cel_const);

	this->change_lattice();
    //this->pricell();         // pengfei Li 2018-05-14 
         //for( iat =0 ; iat < ucell.nat ; iat++)   
//         cout << " newpos_now = " << newpos[3*iat] << " " << newpos[3*iat+1] << " " << newpos[3*iat+2] << endl;
	OUT(ofs_running,"ibrav",ibrav);
    this->setgroup(this->symop, this->nop, this->ibrav);
    //now select all symmetry operations which reproduce the lattice
    //to find those symmetry operations which reproduce the entire crystal
    this->getgroup(this->nrot, this->nrotk);
    // find the name of point group
    this->pointgroup(this->nrot, this->pgnumber, this->pgname, this->gmatrix);
    OUT(ofs_running,"POINT GROUP", this->pgname);
    //write();

    delete[] dirpos;
	delete[] newpos;
    delete[] na;
    delete[] rotpos;
    delete[] ptrans;
    delete[] index;
    delete[] istart;
    return;
}

//---------------------------------------------------
// The lattice will be transformed to a 'standard
// cystallographic setting', the relation between
// 'origin' and 'transformed' lattice vectors will
// be givin in matrix form
//---------------------------------------------------
int Symmetry::standard_lat(
    Vector3<double> &a,
    Vector3<double> &b,
    Vector3<double> &c,
    double *cel_const
)
{
    static bool first = true;
    // there are only 14 types of Bravais lattice.
    int type = 14;
	//----------------------------------------------------
    // used to calculte the volume to judge whether 
	// the lattice vectors corrispond the right-hand-sense
	//----------------------------------------------------
    double volume = 0;
    //the lattice vectors have not been changed
    change = 0;

    const double aa = a * a;
    const double bb = b * b;
    const double cc = c * c;
    const double ab = a * b; //vector: a * b * cos(alpha)
    const double bc = b * c; //vector: b * c * cos(beta)
    const double ca = c * a; //vector: c * a * cos(gamma)
    double norm_a = a.norm();
    double norm_b = b.norm();
    double norm_c = c.norm();
    double alpha = ab /( norm_a * norm_b ); // cos(alpha)
    double beta  = bc /( norm_b * norm_c ); // cos(beta)
    double gamma = ca /( norm_a * norm_c ); // cos(gamma)
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
        OUT(ofs_running,"NORM_A",norm_a);
        OUT(ofs_running,"NORM_B",norm_b);
        OUT(ofs_running,"NORM_C",norm_c);
//		OUT(ofs_running,"alpha  = ", alpha );
//		OUT(ofs_running,"beta   = " ,beta  );
//		OUT(ofs_running,"gamma  = " ,gamma );
        OUT(ofs_running,"ALPHA (DEGREE)", acos(alpha)/PI*180.0 );
        OUT(ofs_running,"BETA  (DEGREE)" ,acos(beta)/PI*180.0  );
        OUT(ofs_running,"GAMMA (DEGREE)" ,acos(gamma)/PI*180.0 );
        first = false;
    }

	
//	cout << " a=" << norm_a << endl; 
//	cout << " b=" << norm_b << endl; 
//	cout << " c=" << norm_c << endl; 
//	cout << " alpha=" << alpha << endl;
//	cout << " beta=" << beta << endl;
//	cout << " gamma=" << gamma << endl;
     

    Symm_Other::right_hand_sense(a, b, c);
	ZEROS(cel_const, 6);
	const double small = 1.0e-5;

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
		else if( equal(alpha,0.0) || equal(beta,0.0) || equal(gamma,0.0)) 
		{
			// Two axes with equal lengths means simple tetragonal: (IBRAV=5)
			// Adjustment: 'c-axis' shall be the special axis.
			if ( equal(norm_a, norm_b) || equal(norm_b, norm_c) || equal(norm_c, norm_a)) 
			{
				type=5;
				cel_const[0]=norm_a;
				cel_const[2]=norm_c/norm_a;
				// No axes with equal lengths means simple orthorhombic (IBRAV=8):
				// Adjustment: Sort the axis by increasing lengths:
			}
                        else
			//else if((abs(norm_c-norm_b)>small) && (abs(norm_b-norm_a)>small) && (abs(norm_a-norm_c)>small)) 
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
	// gamma == beta != alpha
	//------------------------
	else if (abs(gamma-beta)<small || abs(alpha-beta)<small || abs(gamma-alpha)<small ) 
	{
		//---------------------------------------------------------
		// gamma = beta = 90 degree
		// One axis orthogonal with respect to the other two axes:
		//---------------------------------------------------------
		if ( (equal(gamma, 0.0) && equal(beta, 0.0)) || (equal(gamma, 0.0) && equal(alpha, 0.0)) || (equal(alpha, 0.0) && equal(beta, 0.0))) 
		{
			//-----------------------------------------------
			// a == b 
			// Equal length of the two nonorthogonal axes:
			//-----------------------------------------------
			if ( equal(norm_a, norm_b) || equal(norm_b, norm_c) || equal(norm_c, norm_a)) 
			{
				// Cosine(alpha) equal to -1/2 means hexagonal: (IBRAV=4)
				// Adjustment: 'c-axis' shall be the special axis.
				// alpha = 60 degree or 120 degree?
				if ( equal(alpha, -0.5) || equal(alpha, 0.5) || equal(beta, -0.5) || equal(beta, 0.5) || equal(gamma, -0.5)|| equal(gamma, 0.5)) 
				{
					type=4;
					cel_const[0]=norm_a;
					cel_const[2]=norm_c/norm_a;
					// Other angles mean base-centered orthorhombic: (IBRAV=11)
					// Adjustment: Cosine between A1 and A2 shall be lower than zero, the
					//             'c-axis' shall be the special axis.
				}
				// other degree, bug alpha < 0
                                else
				//else if(alpha<(-1.0*small)) 
				{
					type=11;
					cel_const[0]=apb;
					cel_const[1]=amb/apb;
					cel_const[2]=norm_c/apb;
				}
				// Different length of the two axes means simple monoclinic (IBRAV=12):
				// Adjustment: Cosine(gamma) should be lower than zero, special axis
				//             shall be the 'b-axis'(////) and |A1|<|A3|:
			}
			//----------
			// a!=b!=c
			//----------
                        else
			//else if((alpha<(-1.0*small)) && (abs(norm_a-norm_b)>small)) 
			{
				type=12;
				cel_const[0]=norm_b;
				cel_const[1]=norm_c/norm_b;
				cel_const[2]=norm_a/norm_b;
				cel_const[4]=alpha;
				/*
				YB(1)=XB(1,3);
				YB(2)=XB(2,3);
				YB(3)=XB(3,3);
				XB(1,3)=XB(1,1);
				XB(2,3)=XB(2,1);
				XB(3,3)=XB(3,1);
				XB(1,1)=XB(1,2);
				XB(2,1)=XB(2,2);
				XB(3,1)=XB(3,2);
				XB(1,2)=YB(1);
				XB(2,2)=YB(2);
				XB(3,2)=YB(3);
				*/
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
				(abs(cpa-bpc)<small) && 
				(abs(apb-cpa)>small) &&
				(abs(norm_c*norm_c+abc)<small)) 
			{
				type=6;
				cel_const[0]=cpa;
				cel_const[2]=apb/cpa;
			}
			// |A1|=|A2|=/|A3| means base-centered monoclinic (IBRAV=13):
			// Adjustement: The cosine between A1 and A3 as well as the cosine
			//              between A2 and A3 should be lower than zero.
			/*else if((abs(norm_a-norm_b)<small) 
					&& (gamma<(-1.0*small)) 
					&& (beta<(-1.0*small))) 
			{
				type=13;
				cel_const[0]=apb;
				cel_const[1]=amb/apb;
				cel_const[2]=norm_c/apb;
				cel_const[4]=beta;
			}*/

		}
	}
	//-------------------------------
	// three angles are not equal
	//-------------------------------
	else 
	{
		// Crystal classes with A1*A2=/A1*A3=/A2*A3
		// |A1|=|A2|=|A3| means body-centered orthorhombic (IBRAV=9):
		// Further additional criterions are: (A1+A2), (A1+A3) and (A2+A3) are
		// orthogonal to one another and (adjustment//): |A1+A2|>|A1+A3|>|A2+A3|
		if ((abs(norm_a-norm_b)<small) &&
				(abs(norm_b-norm_c)<small) &&
				((cpa-bpc)>small) &&
				((apb-cpa)>small) && 
				(abs(norm_c*norm_c+abc)<small)) 
		{
			type=9;
			cel_const[0]=bpc;
			cel_const[1]=cpa/bpc;
			cel_const[2]=apb/bpc;
		}
		// |A1|=|A2-A3| and |A2|=|A1-A3| and |A3|=|A1-A2| means face-centered
		// orthorhombic (IBRAV=10):
		// Adjustment: |A1+A2-A3|>|A1+A3-A2|>|A2+A3-A1|
                else if((abs(norm_a-norm_b)<small) || (abs(norm_b-norm_c)<small) || (abs(norm_c-norm_a)<small))
                {
                        type=13;
                        cel_const[0]=apb;
                        cel_const[1]=amb/apb;
                        cel_const[2]=norm_c/apb;
                        cel_const[4]=beta;
                }

		else if((abs(amb-norm_c)<small) &&
				(abs(cma-norm_b)<small) &&
				(abs(bmc-norm_a)<small) && 
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
                else
		//else if((alpha>gamma) && (gamma>beta) && (beta>small)) 
		{
			type=14;
			cel_const[0]=norm_a;
			cel_const[1]=norm_b/norm_a;
			cel_const[2]=norm_c/norm_a;
			cel_const[3]=beta;
			cel_const[4]=gamma;
			cel_const[5]=alpha;
		}
	}
	
	return type;
}

void Symmetry::lattice_type(
    Vector3<double> &v1,
    Vector3<double> &v2,
    Vector3<double> &v3,
    int &brav,
    double *cel_const,
    string &bravname
)
{
    TITLE("Symmetry","lattice_type");
//      cout << "v1 = " << v1.x << " " << v1.y << " " << v1.z <<endl;
//      cout << "v2 = " << v2.x << " " << v2.y << " " << v2.z <<endl;
//      cout << "v3 = " << v3.x << " " << v3.y << " " << v3.z <<endl;

	//----------------------------------------------
	// (1) adjustement of the basis to right hand 
	// sense by inversion of all three lattice 
	// vectors if necessary
	//----------------------------------------------
    const bool right = Symm_Other::right_hand_sense(v1, v2, v3);
  //    cout << "v1 = " << v1.x << " " << v1.y << " " << v1.z <<endl;
  //    cout << "v2 = " << v2.x << " " << v2.y << " " << v2.z <<endl;
  //    cout << "v3 = " << v3.x << " " << v3.y << " " << v3.z <<endl;

	OUT(ofs_running,"right hand lattice",right);

	//-------------------------------------------------
	// (2) save and copy the original lattice vectors.
	//-------------------------------------------------
    s1 = v1;
    s2 = v2;
    s3 = v3;
	
	//--------------------------------------------
	// (3) calculate the 'pre_const'
	//--------------------------------------------
    double pre_const[6];
	ZEROS(pre_const, 6);
//    cout << "ATTION !!!!!!" <<endl;
//        cout << "v1 = " << v1.x << " " << v1.y << " " << v1.z <<endl;
//        cout << "v2 = " << v2.x << " " << v2.y << " " << v2.z <<endl;
//        cout << "v3 = " << v3.x << " " << v3.y << " " << v3.z <<endl;

    int pre_brav = standard_lat(v1, v2, v3, cel_const);
//    for ( int i = 0; i < 6; ++i)
//    {
//        cout << "cel = "<<cel_const[i]<<" ";
//    }
//    cout << endl;

//  cout << "pre_brav = " << pre_brav <<endl;

    for ( int i = 0; i < 6; ++i)
    {
        pre_const[i] = cel_const[i];
    }
//        cout << "v1 = " << v1.x << " " << v1.y << " " << v1.z <<endl;
//        cout << "v2 = " << v2.x << " " << v2.y << " " << v2.z <<endl;
//        cout << "v3 = " << v3.x << " " << v3.y << " " << v3.z <<endl;


    //find the shortest basis vectors of the lattice
//    shortest_vector(v1, v2, v3);
//        cout << "a1 = " << v1.x << " " << v1.y << " " << v1.z <<endl;
//        cout << "a1 = " << v2.x << " " << v2.y << " " << v2.z <<endl;
//        cout << "a1 = " << v3.x << " " << v3.y << " " << v3.z <<endl;

    Symm_Other::right_hand_sense(v1, v2, v3);
//        cout << "a1 = " << v1.x << " " << v1.y << " " << v1.z <<endl;
//        cout << "a1 = " << v2.x << " " << v2.y << " " << v2.z <<endl;
//        cout << "a1 = " << v3.x << " " << v3.y << " " << v3.z <<endl;

    int temp_brav = 15;
    double temp_const[6];

    double cos1 = 1;
    double cos2 = 1;
    double cos3 = 1;

    //then we should find the best lattice vectors to make much easier the determination of the lattice symmetry
    //the method is to contrast the combination of the shortest vectors and determine their symmmetry

    Vector3<double> r1, r2, r3;
    Vector3<double> w1, w2, w3;
    Vector3<double> q1, q2, q3;

    int nif = 0;
    for (int n33 = -2; n33 < 3; ++n33)
    {
        for (int n32 = -2; n32 < 3; ++n32)
        {
            for (int n31 = -2; n31 < 3; ++n31)
            {
                for (int n23 = -2; n23 < 3; ++n23)
                {
                    for (int n22 = -2; n22 < 3; ++n22)
                    {
                        for (int n21 = -2; n21 < 3; ++n21)
                        {
                            for (int n13 = -2; n13 < 3; ++n13)
                            {
                                for (int n12 = -2; n12 < 3; ++n12)
                                {
                                    for (int n11 = -2; n11 < 3; ++n11)
                                    {
                                        Matrix3 mat(n11, n12, n13, n21, n22, n23, n31, n32, n33);

                                        if (equal(mat.Det(),1.0))
                                        {
                                            r1.x = n11 * v1.x + n12 * v2.x + n13 * v3.x;
                                            r1.y = n11 * v1.y + n12 * v2.y + n13 * v3.y;
                                            r1.z = n11 * v1.z + n12 * v2.z + n13 * v3.z;
                                     
									        r2.x = n21 * v1.x + n22 * v2.x + n23 * v3.x;
                                            r2.y = n21 * v1.y + n22 * v2.y + n23 * v3.y;
                                            r2.z = n21 * v1.z + n22 * v2.z + n23 * v3.z;
                                     
									        r3.x = n31 * v1.x + n32 * v2.x + n33 * v3.x;
                                            r3.y = n31 * v1.y + n32 * v2.y + n33 * v3.y;
                                            r3.z = n31 * v1.z + n32 * v2.z + n33 * v3.z;
                                            //cout << "mat = " << n11 <<" " <<n12<<" "<<n13<<" "<<n21<<" "<<n22<<" "<<n23<<" "<<n31<<" "<<n32<<" "<<n33<<endl;
											
                                            brav = standard_lat(r1, r2, r3, cel_const);
//                                            if(brav == 8)
//                                            {
//                                               cout << "mat = " << n11 <<" " <<n12<<" "<<n13<<" "<<n21<<" "<<n22<<" "<<n23<<" "<<n31<<" "<<n32<<" "<<n33<<endl;
//                                            }

/*
											if(n11== 1 && n12==0 && n13==-2 && n21==2 && n22==1 && n23==-1
												&& n31==-2 && n32==-2 && n33==-1)
											{
												++nif;
												ofs_running << " " << endl;
												ofs_running << setw(8) << nif << setw(5) << n11 << setw(5) << n12
													<< setw(5) << n13 << setw(5) << n21 << setw(5) << n22
													<< setw(5) << n23 << setw(5) << n31 << setw(5) << n32
													<< setw(5) << n33 << setw(5) << ibrav << endl;
												ofs_running << " r1: " << r1.x << " " << r1.y << " " << r1.z << endl;
												ofs_running << " r2: " << r2.x << " " << r2.y << " " << r2.z << endl;
												ofs_running << " r3: " << r3.x << " " << r3.y << " " << r3.z << endl;
												ofs_running << " cel_const[3]=" << cel_const[3] << endl;
												ofs_running << " cel_const[4]=" << cel_const[4] << endl;
												ofs_running << " cel_const[5]=" << cel_const[5] << endl;
											}
											*/
//											if(brav == 14)
//											{
//												ofs_running << " ABS(CELLDM(4))=" << fabs(cel_const[3]) << endl;
//												ofs_running << " ABS(CELLDM(5))=" << fabs(cel_const[4]) << endl;
//												ofs_running << " ABS(CELLDM(6))=" << fabs(cel_const[5]) << endl;
//											}

                                            if ( brav < temp_brav || ( brav == temp_brav
                                                    && ( fabs(cel_const[3]) < (cos1-1.0e-9) )
                                                    && ( fabs(cel_const[4]) < (cos2-1.0e-9) )
                                                    && ( fabs(cel_const[5]) < (cos3-1.0e-9) )) //mohan fix bug 2012-01-15, not <=
                                               )
                                            {
												/*
												ofs_running << "\n IBRAV=" << brav << " ITYP=" << temp_brav << endl;
												ofs_running << " ABS(CELLDM(4))=" << fabs(cel_const[3]) << endl;
												ofs_running << " ABS(CELLDM(5))=" << fabs(cel_const[4]) << endl;
												ofs_running << " ABS(CELLDM(6))=" << fabs(cel_const[5]) << endl;
												ofs_running << " COS1=" << cos1 << endl;
												ofs_running << " COS2=" << cos2 << endl;
												ofs_running << " COS3=" << cos3 << endl;
												*/
												/*
												ofs_running << " r1: " << r1.x << " " << r1.y << " " << r1.z << endl;
												ofs_running << " r2: " << r2.x << " " << r2.y << " " << r2.z << endl;
												ofs_running << " r3: " << r3.x << " " << r3.y << " " << r3.z << endl;
												ofs_running << " N=" << n11 << " " << n12 << " " << n13 
												<< " " << n21 << " " << n22 << " " << n23 
												<< " " << n31 << " " << n32 << " " << n33 << endl; 
												*/
												//out.printM3("mat",mat);
                                                temp_brav = brav;
												
                                                cos1 = fabs(cel_const[3]);
                                                cos2 = fabs(cel_const[4]);
                                                cos3 = fabs(cel_const[5]);

                                                for (int i = 0; i < 6; ++i)
                                                {
                                                    temp_const[i] = cel_const[i];
                                                }
                                                w1 = r1;
                                                w2 = r2;
                                                w3 = r3;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
//        cout << "a1 = " << v1.x << " " << v1.y << " " << v1.z <<endl;
//        cout << "a1 = " << v2.x << " " << v2.y << " " << v2.z <<endl;
//        cout << "a1 = " << v3.x << " " << v3.y << " " << v3.z <<endl;

    //now, the highest symmetry of the combination of the shortest vectors has been found
    //then we compare it with the original symmetry
	
//	ofs_running << " w1" << endl;
//	ofs_running << " " << setw(15) << w1.x << setw(15) << w1.y << setw(15) << w1.z << endl;
//	ofs_running << " " << setw(15) << w2.x << setw(15) << w2.y << setw(15) << w2.z << endl;
//	ofs_running << " " << setw(15) << w3.x << setw(15) << w3.y << setw(15) << w3.z << endl;
//	ofs_running << " pre_brav=" << pre_brav << endl;
//	ofs_running << " temp_brav=" << temp_brav << endl;

    if ( temp_brav < pre_brav)
    {
        //if the symmetry of the new vectors is higher, store the new ones
        //brav = temp_brav;
        for (int i = 0; i < 6; ++i)
        {
            cel_const[i] = temp_const[i];
        }
        q1 = w1;
        q2 = w2;
        q3 = w3;
        change = 1;
        ofs_running <<endl;
        ofs_running <<" !The lattice vectors have been changed (STRU_SIMPLE.cif)"<<endl;
        ofs_running <<endl;
        int at=0;
        for(int it=0; it<ucell.ntype; ++it)
        {
                for(int ia=0; ia<ucell.atoms[it].na; ++ia)
                {
                        Mathzone::Cartesian_to_Direct(ucell.atoms[it].tau[ia].x,
                                        ucell.atoms[it].tau[ia].y,
                                        ucell.atoms[it].tau[ia].z,
                                        q1.x, q1.y, q1.z,
                                        q2.x, q2.y, q2.z,
                                        q3.x, q3.y, q3.z,
                                        newpos[3*at],newpos[3*at+1],newpos[3*at+2]);

//                        cout << " newpos_before = " << newpos[3*iat] << " " << newpos[3*iat+1] << " " << newpos[3*iat+2] << endl;
//                      ofs_running << " newpos_before = " << newpos[3*iat] << " " << newpos[3*iat+1] << " " << newpos[3*iat+2] << endl;
                        for(int k=0; k<3; ++k)
                        {
                                this->check_translation( newpos[at*3+k], -floor(newpos[at*3+k]));
                                this->check_boundary( newpos[at*3+k] );
                        }
//                        cout << " newpos = " << newpos[3*iat] << " " << newpos[3*iat+1] << " " << newpos[3*iat+2] << endl;

//                      ofs_running << " newpos = " << newpos[3*iat] << " " << newpos[3*iat+1] << " " << newpos[3*iat+2] << endl;
                        ++at;
                }
        }       
        stringstream ss;
        ss << global_out_dir << "STRU_SIMPLE.cif";

        ofstream ofs( ss.str().c_str() );
        ofs << "Lattice vector  : " << endl;
        ofs << q1.x <<"   "<<q1.y<<"  "<<q1.z<< endl;
        ofs << q2.x <<"   "<<q2.y<<"  "<<q2.z<< endl;
        ofs << q3.x <<"   "<<q3.y<<"  "<<q3.z<< endl;
        ofs << endl;
        ofs << "Direct positions : " << " " << endl;
        ofs << endl;
        at =0;
        
        for (int it=0; it<ucell.ntype; it++)
        {
            for (int ia=0; ia<ucell.atoms[it].na; ia++)
            {
                 ofs << ucell.atoms[it].label
                 << " " << newpos[3*at]
                 << " " << newpos[3*at+1]
                 << " " << newpos[3*at+2] << endl;
                 at++;
            }
        }
        ofs.close();

        
    }
    
    /*else
    {
        //else, store the original ones
        brav = pre_brav;
        for (int i = 0; i < 6; ++i)
        {
            cel_const[i] = pre_const[i];
        }
    }


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
			ofs_running<<" The lattice vectors have been set back!"<<endl;
        }
    }*/
    brav = pre_brav;
    bravname = get_brav_name(brav);

    OUT(ofs_running,"BRAVAIS TYPE",brav);
    OUT(ofs_running,"BRAVAIS LATTICE NAME",bravname);
    return;
}


void Symmetry::change_lattice(void)
{
    //if lattice vectors are changed, do the coordinates conversion
    if (test_symmetry) TITLE("Symmetry","change_lattice");

	change = 0;

	//ofs_running << "\n change = " << change;
    if (change == 1)
    {
        this->veccon(dirpos, rotpos, nat, s1, s2, s3, a1, a2, a3);
        for (int i = 0; i < nat * 3; ++i)
        {
            dirpos[i] = rotpos[i];
        }
    }
    else
    {
        for (int i = 0; i < nat * 3; ++i)
        {
            rotpos[i] = dirpos[i];
        }
    }
    return;
}

void Symmetry::pricell(void)
{
    //detect the generating cell (primitive cell) of a supercell
    if (test_symmetry) TITLE("Symmetry","pricell");

    // part 1 of pricell
    for (int it = 0; it < ntype; it++)
    {
        for (int ia = istart[it]; ia < istart[it] + na[it]; ia++)
        {
            for (int k = 0; k < 3; k++)
            {
                this->check_translation( this->dirpos[ia*3+k], -floor(dirpos[ia*3+k]));
                this->check_boundary( this->dirpos[ia*3+k] );
            }
        }
        //order original atomic positions for current species
        //cout << "\n -----------------------------------------";
        //cout << "\nCarpos it=" << it << " na=" << na[it];
        //cout << "\nBefore dirpos:" << endl;
        //print_pos(dirpos+istart*3, na[it]);
        this->atom_ordering(&dirpos[istart[it]*3], na[it], index + istart[it]);
        //cout << "\nAfter dirpos:" << endl;
        //print_pos(dirpos+istart*3, na[it]);
    }

    // part 2 of pricell
    if (na[itmin_type] == 1)
    {
        //if there is only one atom in the "itmin" species,
        //the supercell has already been a primitive cell
        ncell = 1;
        ptrans[0] = 0;
        ptrans[1] = 0;
        ptrans[2] = 0;
        p1 = a1;
        p2 = a2;
        p3 = a3;
        pbrav = ibrav;
        for (int k = 0; k < 6; ++k)
        {
            pcel_const[k] = cel_const[k];
        }
        ofs_running <<"\n The cell is a primitive cell." << endl;
        return;
    }

    int itrans = 3;
    double diff[3];
    double trans[3];

    for (int i = 1; i < 8; ++i)
    {
        ptrans[i] = 0;
    }
    ptrans[0] = 1;
    ptrans[4] = 1;
    ptrans[8] = 1;

    //All difference vectors between the first atom and all other atoms
    //of species ISMIN could be a translation vector that reproduces thelattice
    //So test all possibilities:

    //cout << "\n cartesian" << endl;
    //print_pos(dirpos, nat);

    sptmin.x = dirpos[itmin_start*3];
    sptmin.y = dirpos[itmin_start*3+1];
    sptmin.z = dirpos[itmin_start*3+2];
    for (int ia = itmin_start + 1; ia < itmin_start + na[itmin_type]; ia++)
    {
        //	cout << "\n ia=" << ia;
        //	cout << " dirpos:" << dirpos[ia*3] << " " << dirpos[ia*3+1] << " " << dirpos[ia*3+2];
        //set up the current test vector "trans"
        //and "trans" could possibly contain trivial translations:
        trans[0] = this->get_translation_vector(dirpos[ia*3+0], sptmin.x);
        trans[1] = this->get_translation_vector(dirpos[ia*3+1], sptmin.y);
        trans[2] = this->get_translation_vector(dirpos[ia*3+2], sptmin.z);

        //	cout << "\n Trans vector = " << trans[0] << " " << trans[1] << " " << trans[2];

        //translate all the atomic coordinates by "trans"
        for (int it = 0; it < ntype; it++)
        {
            for (int ia2=istart[it]; ia2<na[it]+istart[it]; ia2++)
            {
                for (int k = 0; k < 3; k++)
                {
                    this->check_translation( rotpos[ia2*3+k], trans[k] );
                    this->check_boundary( rotpos[ia2*3+k] );
                }
            }
            //order translated atomic positions for current species
            //    cout << "\n rotpos it=" << it << " na=" << na[it];
            //   cout << "\n before rotpos:" << endl;
            // print_pos(rotpos+istart[it]*3, na[it]);
            this->atom_ordering(&rotpos[istart[it]*3], na[it], index+istart[it]);
            //    cout << "\n after rotpos:" << endl;
            // print_pos(rotpos+istart[it]*3, na[it]);

            //	BLOCK_HERE("type");
        }

        bool found_tvec = true;
        //compare the two lattices 'one-by-one' whether they are identical
        for (int iat=0; iat<nat; iat++)
        {
            diff[0] = this->check_diff(dirpos[iat*3+0], rotpos[iat*3+0]);
            diff[1] = this->check_diff(dirpos[iat*3+1], rotpos[iat*3+1]);
            diff[2] = this->check_diff(dirpos[iat*3+2], rotpos[iat*3+2]);

            if (!equal(diff[0],0.0) || !equal(diff[1],0.0) || !equal(diff[2],0.0))
            {
                found_tvec = false;
                break;
            }
        }

//        cout << "\n dirpos:" << endl;
//        print_pos(dirpos, nat);
//        cout << "\n rotpos:" << endl;
//        print_pos(rotpos, nat);

        //the current test vector is a 'primitive' translation, save the detected translation vector
        if (found_tvec)
        {
            ptrans[itrans*3+0] = trans[0];
            ptrans[itrans*3+1] = trans[1];
            ptrans[itrans*3+2] = trans[2];
            ++itrans;
            //cout<<"Find "<<itrans<<" translation"<<endl;
        }
    }

	OUT(ofs_running,"Number of transition vectors",itrans);

    if (itrans == 3)
    {
        //no primitive translations are found, the supercell has already been a primitive cell
        ncell = 1;
        for (int k = 0; k < 3; k++)
        {
            ptrans[k] = 0;
        }
        p1 = a1;
        p2 = a2;
        p3 = a3;

        pbrav = ibrav;

        for (int k = 0; k < 6; ++k)
        {
            pcel_const[k] = cel_const[k];
        }
        OUT(ofs_running,"INPUT CELL IS PRIMITIVE","True");
        return;
    }

    //if some primitive translations are found, order them
    this->atom_ordering(ptrans, itrans, index);

    double first;
    int plane;

    Vector3<double> b1, b2, b3;

    //find first 'yz'-plane:
    first = ptrans[0];
    for (int i = 1; i < itrans; ++i)
    {
        if ( !equal(first,ptrans[i * 3]) )
        {
            plane = i;
            break;
        }
    }

    b1.x = ptrans[plane * 3 + 0];
    b1.y = ptrans[plane * 3 + 1];
    b1.z = ptrans[plane * 3 + 2];

    //find first 'z'-line
    first = ptrans[1];
    for (int i = 1; i < plane; ++i)
    {
        if ( !equal(first,ptrans[i * 3 + 1]) )
        {
            plane = i;
            break;
        }
    }

    b2.x = ptrans[plane * 3 + 0];
    b2.y = ptrans[plane * 3 + 1];
    b2.z = ptrans[plane * 3 + 2];

    b3.x = ptrans[0];
    b3.y = ptrans[1];
    b3.z = ptrans[2];

    p1.x = b1.x * a1.x + b1.y * a2.x + b1.z * a3.x;
    p1.y = b1.x * a1.y + b1.y * a2.y + b1.z * a3.y;
    p1.z = b1.x * a1.z + b1.y * a2.z + b1.z * a3.z;

    p2.x = b2.x * a1.x + b2.y * a2.x + b2.z * a3.x;
    p2.y = b2.x * a1.y + b2.y * a2.y + b2.z * a3.y;
    p2.z = b2.x * a1.z + b2.y * a2.z + b2.z * a3.z;

    p3.x = b3.x * a1.x + b3.y * a2.x + b3.z * a3.x;
    p3.y = b3.x * a1.y + b3.y * a2.y + b3.z * a3.y;
    p3.z = b3.x * a1.z + b3.y * a2.z + b3.z * a3.z;

    ofs_running << " a1:" << setw(20) << a1.x << setw(20) << a1.y << setw(20) << a1.z << endl;
    ofs_running << " a2:" << setw(20) << a2.x << setw(20) << a2.y << setw(20) << a2.z << endl;
    ofs_running << " a3:" << setw(20) << a3.x << setw(20) << a3.y << setw(20) << a3.z << endl;

    ofs_running << " b1:" << setw(20) << b1.x << setw(20) << b1.y << setw(20) << b1.z << endl;
    ofs_running << " b2:" << setw(20) << b2.x << setw(20) << b2.y << setw(20) << b2.z << endl;
    ofs_running << " b3:" << setw(20) << b3.x << setw(20) << b3.y << setw(20) << b3.z << endl;

    ofs_running << " p1:" << setw(20) << p1.x << setw(20) << p1.y << setw(20) << p1.z << endl;
    ofs_running << " p2:" << setw(20) << p2.x << setw(20) << p2.y << setw(20) << p2.z << endl;
    ofs_running << " p3:" << setw(20) << p3.x << setw(20) << p3.y << setw(20) << p3.z << endl;

    //analyse the data and get the symmetry infomation
//	cout<<"calculating the properties!"<<endl;

    Vector3<double> zero(0.0,0.0,0.0);
    if (p1 == zero || p2 == zero || p3 == zero)
    {
		WARNING_QUIT("Symmetry::pricell","At least one of the primitive vector is (0,0,0).");
    }

    double celvolume = 0;
    double pcelvolume = 0.0;

    this->lattice_type(p1, p2, p3, pbrav, pcel_const, plattname);

    pcelvolume = std::fabs( Symm_Other::celvol(p1,p2,p3) );
	OUT(ofs_running,"pcelvolume",pcelvolume);

    this->lattice_type(a1, a2, a3, ibrav, cel_const, ilattname);

    celvolume = std::fabs( Symm_Other::celvol(a1, a2, a3) );
	OUT(ofs_running,"celvolume",celvolume);

    ncell= celvolume / pcelvolume;

	OUT(ofs_running,"ncell",ncell);

    //output some result to the screen
	OUT(ofs_running,"BRAVAIS OF SUPERCELL",ibrav);
	OUT(ofs_running,"LATTICE NAME OF SUPERCELL",ilattname);

    ofs_running<<" LATTICE VECTORS OF SUPERCELL" << endl;
    ofs_running<<" S1:" << setw(20) << a1.x << setw(20) << a1.y << setw(20) << a1.z << endl;
    ofs_running<<" S2:" << setw(20) << a2.x << setw(20) << a2.y << setw(20) << a2.z << endl;
    ofs_running<<" S3:" << setw(20) << a3.x << setw(20) << a3.y << setw(20) << a3.z << endl;

	OUT(ofs_running,"BRAVAIS OF PRIMITIVE CELL",pbrav);
	OUT(ofs_running,"LATTICE NAME OF PRIMITIVE CELL",plattname);

    ofs_running<<" LATTICE VECTORS OF PRIMITIVE CELL:" << endl;
    ofs_running<<" P1:" << setw(20) << p1.x << setw(20) << p1.y << setw(20) << p1.z << endl;
    ofs_running<<" P2:" << setw(20) << p2.x << setw(20) << p2.y << setw(20) << p2.z << endl;
    ofs_running<<" P3:" << setw(20) << p3.x << setw(20) << p3.y << setw(20) << p3.z << endl;

	OUT(ofs_running,"PRIMITIVE CELLS",ncell);
    ofs_running<<" PRIMITIVE TRANSLATION VECTORS:"<<endl;
    for (int i = 0; i < itrans; ++i)
    {
        ofs_running << " " << setw(20) << ptrans[i*3]
		<< " " << setw(20) << ptrans[i*3+1]
		<< " " << setw(20) << ptrans[i*3+2]
		<< endl;
        //cout<<"n = "<<i*3+k<<endl;
    }
    //cout<<"pricell ends!"<<endl;
    return;
}


void Symmetry::getgroup(int &nrot, int &nrotk)
{
    TITLE("Symmetry","getgroup");

	//--------------------------------------------------------------------------------
    //GETGRP (L347 symlib.f VASP)
    //return all possible space group operators that reproduce a lattice with basis
    //out of a (maximum) pool of point group operations that is compatible with
    //the symmetry of the pure translation lattice without any basic.
	//--------------------------------------------------------------------------------

    Matrix3 zero(0,0,0,0,0,0,0,0,0);
    Matrix3 help[48];
    Vector3<double> temp[48];

    nrot = 0;
    nrotk = 0;

	//-------------------------------------------------------------------------
    //pass through the pool of (possibly allowed) symmetry operations and
    //check each operation whether it can reproduce the lattice with basis
	//-------------------------------------------------------------------------
    //cout << "nop = " <<nop <<endl;
    for (int i = 0; i < nop; ++i)
    {
    //    cout << "symop = " << symop[i].e11 <<" "<< symop[i].e12 <<" "<< symop[i].e13 <<" "<< symop[i].e21 <<" "<< symop[i].e22 <<" "<< symop[i].e23 <<" "<< symop[i].e31 <<" "<< symop[i].e32 <<" "<< symop[i].e33 << endl;
        this->checksym(this->symop[i], this->gtrans[i], this->newpos);
      //  cout << "s_flag =" <<s_flag<<endl;
        if (s_flag == 1)
        {
			//------------------------------
            // this is a symmetry operation
			// with no translation vectors
            // so ,this is pure point group 
			// operations
			//------------------------------
            if ( equal(gtrans[i].x,0.0) &&
                 equal(gtrans[i].y,0.0) &&
                 equal(gtrans[i].z,0.0))
            {
                ++nrot;
                gmatrix[nrot - 1] = symop[i];
                gtrans[nrot - 1].x = 0;
                gtrans[nrot - 1].y = 0;
                gtrans[nrot - 1].z = 0;
            }
			//------------------------------
            // this is a symmetry operation
			// with translation vectors
            // so ,this is space group 
			// operations
			//------------------------------
            else
            {
                ++nrotk;
                help[nrotk - 1] = symop[i];
                temp[nrotk - 1].x = gtrans[i].x;
                temp[nrotk - 1].y = gtrans[i].y;
                temp[nrotk - 1].z = gtrans[i].z;
            }
        }
    }

	//-----------------------------------------------------
    //If there are operations with nontrivial translations
    //then store them together in the momory
	//-----------------------------------------------------
    if (nrotk > 0)
    {
        for (int i = 0; i < nrotk; ++i)
        {
            gmatrix[nrot + i] = help[i];
            gtrans[nrot + i].x = temp[i].x;
            gtrans[nrot + i].y = temp[i].y;
            gtrans[nrot + i].z = temp[i].z;
        }
    }

	//-----------------------------------------------------
    //total number of space group operations
	//-----------------------------------------------------
    nrotk += nrot;
	OUT(ofs_running,"PURE POINT GROUP OPERATIONS",nrot);
    OUT(ofs_running,"SPACE GROUP OPERATIONS",nrotk);

	//-----------------------------------------------------
    //fill the rest of matrices and vectors with zeros
	//-----------------------------------------------------
    if (nrotk < 48)
    {
        for (int i = nrotk; i < 48; ++i)
        {
            gmatrix[i] = zero;
            gtrans[i].x = 0;
            gtrans[i].y = 0;
            gtrans[i].z = 0;
        }
    }
    return;
}

void Symmetry::checksym(Matrix3 &s, Vector3<double> &gtrans, double* pos)
{
	//----------------------------------------------
    // checks whether a point group symmetry element 
	// is a valid symmetry operation on a supercell
	//----------------------------------------------
    // the start atom index.
    bool no_diff = 0;
    Vector3<double> trans(2.0, 2.0, 2.0);
    s_flag = 0;

    for (int it = 0; it < ntype; it++)
    {
		//------------------------------------
        // impose periodic boundary condition
		// 0.5 -> -0.5
		//------------------------------------
        for (int j = istart[it]; j < istart[it] + na[it]; ++j)
        {
            this->check_boundary(pos[j*3+1]);
            this->check_boundary(pos[j*3+1]);
            this->check_boundary(pos[j*3+2]);
        }
         //for( int iat =0 ; iat < ucell.nat ; iat++)
         //cout << " newpos_now1 = " << newpos[3*iat] << " " << newpos[3*iat+1] << " " << newpos[3*iat+2] << endl;

        //order original atomic positions for current species
        this->atom_ordering(pos + istart[it] * 3, na[it], index + istart[it]);
         //for( int iat =0 ; iat < ucell.nat ; iat++)
         //cout << " newpos_now2 = " << newpos[3*iat] << " " << newpos[3*iat+1] << " " << newpos[3*iat+2] << endl;

        //Rotate atoms of current species
        for (int j = istart[it]; j < istart[it] + na[it]; ++j)
        {
            const int xx=j*3;
            const int yy=j*3+1;
            const int zz=j*3+2;
            

            rotpos[xx] = pos[xx] * s.e11
                         + pos[yy] * s.e21
                         + pos[zz] * s.e31;

            rotpos[yy] = pos[xx] * s.e12
                         + pos[yy] * s.e22
                         + pos[zz] * s.e32;

            rotpos[zz] = pos[xx] * s.e13
                         + pos[yy] * s.e23
                         + pos[zz] * s.e33;

           // cout << "pos = " << pos[xx] <<" "<<pos[yy] << " "<<pos[zz]<<endl;
           // cout << "rotpos = " << rotpos[xx] <<" "<<rotpos[yy] << " "<<rotpos[zz]<<endl;
            rotpos[xx] = fmod(rotpos[xx] + 100.5,1) - 0.5;
            rotpos[yy] = fmod(rotpos[yy] + 100.5,1) - 0.5;
            rotpos[zz] = fmod(rotpos[zz] + 100.5,1) - 0.5;
            this->check_boundary(rotpos[xx]);
            this->check_boundary(rotpos[yy]);
            this->check_boundary(rotpos[zz]);
        }
        //order rotated atomic positions for current species
        this->atom_ordering(rotpos + istart[it] * 3, na[it], index + istart[it]);
    }

	/*
	ofs_running << " ============================================= " << endl;
	ofs_running << " Matrix S " << endl;
	ofs_running << setw(5) << s.e11 << setw(5) << s.e12 << setw(5) << s.e13 << endl;
	ofs_running << setw(5) << s.e21 << setw(5) << s.e22 << setw(5) << s.e32 << endl;
	ofs_running << setw(5) << s.e23 << setw(5) << s.e23 << setw(5) << s.e33 << endl;
	ofs_running << " pos" << endl;
	print_pos(pos, nat);
	ofs_running << " rotpos" << endl;
	print_pos(rotpos, nat);
	*/

    Vector3<double> diff;

	//---------------------------------------------------------
    // itmin_start = the start atom positions of species itmin
	//---------------------------------------------------------
    sptmin.x = pos[itmin_start*3];
    sptmin.y = pos[itmin_start*3+1];
    sptmin.z = pos[itmin_start*3+2];
    for (int i = itmin_start; i < itmin_start + na[itmin_type]; ++i)
    {
        //set up the current test vector "gtrans"
        //and "gtrans" could possibly contain trivial translations:
        gtrans.x = this->get_translation_vector( pos[i*3+0], sptmin.x);
        gtrans.y = this->get_translation_vector( pos[i*3+1], sptmin.y);
        gtrans.z = this->get_translation_vector( pos[i*3+2], sptmin.z);

        //If we had already detected some translation,
        //we must only look at the vectors with coordinates smaller than those
        //of the previously detected vector (find the smallest)
        if (gtrans.x > trans.x + epsilon ||
                gtrans.y > trans.y + epsilon ||
                gtrans.z > trans.z + epsilon
           )
        {
            continue;
        }

        //translate all the atomic coordinates by "gtran"
        for (int it = 0; it < ntype; it++)
        {
            for (int ia = istart[it]; ia < na[it] + istart[it]; ia++)
            {
                this->check_boundary( rotpos[ia*3+0] );
                this->check_boundary( rotpos[ia*3+1] );
                this->check_boundary( rotpos[ia*3+2] );

                this->check_translation( rotpos[ia*3+0], gtrans.x );
                this->check_translation( rotpos[ia*3+1], gtrans.y );
                this->check_translation( rotpos[ia*3+2], gtrans.z );
            }
            //order translated atomic positions for current species
            this->atom_ordering(rotpos + istart[it] * 3, na[it], index + istart[it]);
        }

        no_diff = true;
        //compare the two lattices 'one-by-one' whether they are identical
        for (int it = 0; it < ntype; it++)
        {
            for (int ia = istart[it]; ia < na[it] + istart[it]; ia++)
            {
                //take the difference of the rotated and the original coordinates
                diff.x = this->check_diff( pos[ia*3+0], rotpos[ia*3+0]);
                diff.y = this->check_diff( pos[ia*3+1], rotpos[ia*3+1]);
                diff.z = this->check_diff( pos[ia*3+2], rotpos[ia*3+2]);

                //only if all "diff" are zero vectors, flag will remain "1"
                if (	no_diff == false||
                        !equal(diff.x,0.0)||
                        !equal(diff.y,0.0)||
                        !equal(diff.z,0.0)
                   )
                {
                    no_diff = 0;
                }
            }
        }
			

		/*
		ofs_running << " no_diff = " << no_diff << endl;
		ofs_running << " CHECK pos " << endl;
		print_pos(pos, nat);
		ofs_running << " CHECK rotpos " << endl;
		print_pos(rotpos, nat);
		*/
		//BLOCK_HERE("check symm");

				

        //the current test is successful
        if (no_diff == true)
        {
            s_flag = 1;
            //save the detected translation vector temporarily
            trans.x = gtrans.x;
            trans.y = gtrans.y;
            trans.z = gtrans.z;
        }

        //restore the original rotated coordinates by subtracting "gtrans"
        for (int it = 0; it < ntype; it++)
        {
            for (int ia = istart[it]; ia < na[it] + istart[it]; ia++)
            {
                rotpos[ia*3+0] -= gtrans.x;
                rotpos[ia*3+1] -= gtrans.y;
                rotpos[ia*3+2] -= gtrans.z;
            }
        }
    }

    if (s_flag == 1)
    {
        gtrans.x = trans.x;
        gtrans.y = trans.y;
        gtrans.z = trans.z;
    }
    return;
}


//modified by shu on 2010.01.20
void Symmetry::rho_symmetry( double *rho,
                             const int &nr1, const int &nr2, const int &nr3)
{
//  if (test_symmetry)TITLE("Symmetry","rho_symmetry");
    timer::tick("Symmetry","rho_symmetry");

    //for fft commensuration
	//nrotk : the number of space operations.
    int count_fft = 0;
    for (int i=0; i<nrotk; ++i)
    {
        if ( (gtrans[i].x * nr1 - int(gtrans[i].x * nr1) < epsilon)
           &&(gtrans[i].y * nr2 - int(gtrans[i].y * nr2) < epsilon)
           &&(gtrans[i].z * nr3 - int(gtrans[i].z * nr3) < epsilon)
           )
        {
            ++count_fft;
            this->symflag_fft[i] = true;
        }
        else
        {
            this->symflag_fft[i] = false;
        }
    }
    nrotk = count_fft;
	//cout << "\n nrotk = " << nrotk;


	// get the remaining rotation matrix.
    Matrix3 *gmatrix_fft = new Matrix3[48];

    int counter = 0;
    for (int i=0; i<48; ++i)
    {
        if (this->symflag_fft[i])
        {
            gmatrix_fft[counter] = this->gmatrix[i];
            ++counter;
        }
    }
    for (int i=0; i<48; ++i)
    {
        gmatrix[i] = gmatrix_fft[i];
    }

	// allocate flag for each FFT grid.
    bool* symflag = new bool[nr1 * nr2 * nr3];
    for (int i=0; i<nr1*nr2*nr3; i++)
    {
        symflag[i] = false;
    }

    assert(nrotk >0 );
    assert(nrotk <=48 );
    int *ri = new int[nrotk];
    int *rj = new int[nrotk];
    int *rk = new int[nrotk];

    int ci = 0;
    for (int i = 0; i< nr1; ++i)
    {
        for (int j = 0; j< nr2; ++j)
        {
            for (int k = 0; k< nr3; ++k)
            {
                if (!symflag[i * nr2 * nr3 + j * nr3 + k])
                {
                    double sum = 0;

                    for (int isym = 0; isym < nrotk; ++isym)
                    {
                        this->rotate(gmatrix[isym], gtrans[isym], i, j, k, nr1, nr2, nr3, ri[isym], rj[isym], rk[isym]);
                        const int index = ri[isym] * nr2 * nr3 + rj[isym] * nr3 + rk[isym];
                        sum += rho[ index ];
                    }
                    sum /= nrotk;

                    for (int isym = 0; isym < nrotk; ++isym)
                    {
                        const int index = ri[isym] * nr2 * nr3 + rj[isym] * nr3 + rk[isym];
                        rho[index] = sum;
                        symflag[index] = true;
                    }
                }
            }
        }
    }

    delete[] symflag;
    delete[] ri;
    delete[] rj;
    delete[] rk;
    timer::tick("Symmetry","rho_symmetry");
}

void Symmetry::force_symmetry(matrix &force , double* pos)   // pengfei 2016-12-20
{
	TITLE("Symmetry","force_symmetry");
	double *protpos;
	double *tot_force;
	int *n;
	int *start;
	double diff1,diff2,diff3;
	protpos = new double[nat*3]; ZEROS(protpos, nat*3);
	tot_force = new double[nat*3]; ZEROS(tot_force, nat*3);
	n = new int[nat]; ZEROS(n, nat);
	start = new int[ntype]; start[0] = 0;
	for(int it = 0; it < ntype; ++it)
	{
		//Atom* atom = &ucell.atoms[it];
		//na[it] = atom->na;
		if(it > 0)
		{
			start[it] = start[it-1] + ucell.atoms[it-1].na;
		}
		//cout << "na =" <<ucell.atoms[0].na<<" "<<ucell.atoms[1].na<<" "<<ucell.atoms[2].na<<endl;
		
	}
	for(int it = 0; it < ntype; it++)
	{
		for(int j = start[it]; j < start[it] + ucell.atoms[it].na; ++j)
		{
			const int xx=j*3; const int yy=j*3+1; const int zz=j*3+2;
			// cout << "xx = "<<xx<<" yy ="<<yy<<" zz = "<<zz<<endl;
			// cout << "nrotk ="<<nrotk<<endl;
			for(int k = 0 ; k < nrotk; ++k)
			{
				protpos[xx] = pos[xx] * gmatrix[k].e11 + pos[yy] * gmatrix[k].e21 + pos[zz] * gmatrix[k].e31 + gtrans[k].x;
				protpos[yy] = pos[xx] * gmatrix[k].e12 + pos[yy] * gmatrix[k].e22 + pos[zz] * gmatrix[k].e32 + gtrans[k].y;
				protpos[zz] = pos[xx] * gmatrix[k].e13 + pos[yy] * gmatrix[k].e23 + pos[zz] * gmatrix[k].e33 + gtrans[k].z;
							    			   			
				check_translation( protpos[xx], -floor(protpos[xx]));
				check_boundary( protpos[xx] );
				check_translation( protpos[yy], -floor(protpos[yy]));
				check_boundary( protpos[yy] );
				check_translation( protpos[zz], -floor(protpos[zz]));
				check_boundary( protpos[zz] );
				
				for(int l = start[it]; l < start[it] + ucell.atoms[it].na; ++l)
				{
					diff1 = check_diff( pos[l*3], protpos[xx]);
					diff2 = check_diff( pos[l*3+1], protpos[yy]);
					diff3 = check_diff( pos[l*3+2], protpos[zz]);
					if (equal(diff1,0.0) && equal(diff2,0.0) && equal(diff3,0.0))
					{
						//cout <<"nl = " << n[l]<<endl;
						tot_force[l*3] = tot_force[l*3] + force(j,0) * gmatrix[k].e11 + force(j,1) * gmatrix[k].e21 + force(j,2) * gmatrix[k].e31;
						tot_force[l*3+1] =  tot_force[l*3+1] + force(j,0) * gmatrix[k].e12 + force(j,1) * gmatrix[k].e22 + force(j,2) * gmatrix[k].e32;
						tot_force[l*3+2] =  tot_force[l*3+2] + force(j,0) * gmatrix[k].e13 + force(j,1) * gmatrix[k].e23 + force(j,2) * gmatrix[k].e33;
						n[l]++;
					}
				}
				
			}
			
		}
	}
	for(int it = 0; it < ntype; it++)
	{
		for(int j = start[it]; j < start[it] + ucell.atoms[it].na; j++)
		{
			force(j,0) = tot_force[j*3]/n[j];
			force(j,1) = tot_force[j*3+1]/n[j];
			force(j,2) = tot_force[j*3+2]/n[j];
		}
	}
	
	delete[] protpos;
	delete[] tot_force;
	delete[] n;
	delete[] start;
	
	return;
}

void Symmetry::stress_symmetry(matrix& sigma)   //zhengdy added 2017
{
	double *tot_sigma, *temp;
	tot_sigma = new double[9];
	temp = new double[9];
	ZEROS(temp, 9);
	ZEROS(tot_sigma, 9);

	temp[0]=ucell.a1.x;
	temp[1]=ucell.a1.y;
	temp[2]=ucell.a1.z;
	temp[3]=ucell.a2.x;
	temp[4]=ucell.a2.y;
	temp[5]=ucell.a2.z;
	temp[6]=ucell.a3.x;
	temp[7]=ucell.a3.y;
	temp[8]=ucell.a3.z;

	for(int i=0;i<3;i++)
	{
		for(int j= 0;j<3;j++)
		{
			for(int k=0;k<3;k++)
			{
				for(int l=0;l<3;l++)
				{
					tot_sigma[i*3 +j] += sigma(k,l) * temp[i*3+k] * temp[j*3+l];
				}
			}
		}
	}

	for(int i=0;i<3;i++)
	{
		for(int j = 0;j<3;j++)
		{
			sigma(i,j) = tot_sigma[i*3+j];
		}
	}

	ZEROS(temp, 9);
	ZEROS(tot_sigma, 9);

	for ( int k = 0 ; k < nrotk; ++k)
	{
		temp[0] = gmatrix[k].e11;
		temp[1] = gmatrix[k].e12;
		temp[2] = gmatrix[k].e13;
		temp[3] = gmatrix[k].e21;
		temp[4] = gmatrix[k].e22;
		temp[5] = gmatrix[k].e23;
		temp[6] = gmatrix[k].e31;
		temp[7] = gmatrix[k].e32;
		temp[8] = gmatrix[k].e33;

		for( int i=0; i<3; i++)
		{
			for( int j=0; j<3; j++)
			{
				for( int l=0; l<3; l++)
				{
					for( int m=0; m<3; m++)
					{
						tot_sigma[i * 3 +j] += sigma(l,m) * temp[i * 3 + l] * temp[j * 3 + m];
					}
				}
			}
		}
	}

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			sigma(i,j) = tot_sigma[i *3 + j] / nrotk;
		}
	}

	ZEROS(temp, 9);
	ZEROS(tot_sigma, 9);

	double det = ucell.a1.x*ucell.a2.y*ucell.a3.z -
		ucell.a1.x*ucell.a3.y*ucell.a2.z +
		ucell.a2.x*ucell.a3.y*ucell.a1.z -
		ucell.a2.x*ucell.a1.y*ucell.a3.z +
		ucell.a3.x*ucell.a1.y*ucell.a2.z -
		ucell.a3.x*ucell.a2.y*ucell.a1.z;

	if(det == 0)
	{
		det = 1;
	}

	temp[0] = (ucell.a2.y*ucell.a3.z - ucell.a2.z*ucell.a3.y) / det;
	temp[1] = -(ucell.a1.y*ucell.a3.z - ucell.a1.z*ucell.a3.y) / det;
	temp[2] = (ucell.a1.y*ucell.a2.z - ucell.a1.z*ucell.a2.y) / det;
	temp[3] = -(ucell.a2.x*ucell.a3.z - ucell.a2.z*ucell.a3.x) / det;
	temp[4] = (ucell.a1.x*ucell.a3.z - ucell.a1.z*ucell.a3.x) / det;
	temp[5] = -(ucell.a1.x*ucell.a2.z - ucell.a1.z*ucell.a2.x) / det;
	temp[6] = (ucell.a2.x*ucell.a3.y - ucell.a2.y*ucell.a3.x) / det;
	temp[7] = -(ucell.a1.x*ucell.a3.y - ucell.a1.y*ucell.a3.x) / det;
	temp[8] = (ucell.a1.x*ucell.a2.y - ucell.a1.y*ucell.a2.x) / det;

	for(int i=0;i<3;i++)
	{
		for(int j= 0;j<3;j++)
		{
			for(int k=0;k<3;k++)
			{
				for(int l=0;l<3;l++)
				{
					tot_sigma[i*3 +j] += sigma(k,l) * temp[i*3+k] * temp[j*3+l];
				}
			}
		}
	}

	for(int i=0;i<3;i++)
	{
		for(int j = 0;j<3;j++)
		{
			sigma(i,j) = tot_sigma[i*3+j];
		}
	}

	delete [] tot_sigma;
	delete [] temp;
	return;
}

void Symmetry::write(void)
{
    if (test_symmetry) TITLE("Symmetry","write");
    ofs_running<<endl;
    ofs_running<<"\n The point group serial number is "<<pgnumber<<"."<<endl;
    ofs_running<<"\n Its Schoenflies name is "<<pgname<<"."<<endl;
    ofs_running<<"\n There are "<<nrotk<<" rotation matrices for the pure point group symmetry:"<<endl;
    ofs_running<<"    E11 E12 E13 E21 E22 E23 E31 E32 E33"<<endl;
    for (int i = 0; i < nrot; ++i)
    {
        ofs_running << setw(2) <<i + 1<<" "
        << setw(4) << gmatrix[i].e11
        << setw(4) << gmatrix[i].e12
        << setw(4) << gmatrix[i].e13
        << setw(4) << gmatrix[i].e21
        << setw(4) << gmatrix[i].e22
        << setw(4) << gmatrix[i].e23
        << setw(4) << gmatrix[i].e31
        << setw(4) << gmatrix[i].e32
        << setw(4) << gmatrix[i].e33 << endl;

    }
    ofs_running<<endl;
    ofs_running<<"\n There are "<<nrotk<<" full space group operations:"<<endl;
    ofs_running<<"    e11 e12 e13 e21 e22 e23 e31 e32 e33"<<endl;
    for (int i = 0; i < nrotk; ++i)
    {
        ofs_running << setw(2) <<i + 1<<" "
        << setw(4) << gmatrix[i].e11
        << setw(4) << gmatrix[i].e12
        << setw(4) << gmatrix[i].e13
        << setw(4) << gmatrix[i].e21
        << setw(4) << gmatrix[i].e22
        << setw(4) << gmatrix[i].e23
        << setw(4) << gmatrix[i].e31
        << setw(4) << gmatrix[i].e32
        << setw(4) << gmatrix[i].e33
        << setw(4) << "+"
        << gtrans[i].x << " "
        << gtrans[i].y << " "
        << gtrans[i].z << endl;
    }
    ofs_running<<endl;
    return;
}

void Symmetry::print_pos(const double* pos, const int &nat)
{
    for (int i=0; i<nat; i++)
    {
        ofs_running << " pos " << i+1 << ": " << pos[i*3+0] << " " << pos[i*3+1] << " " << pos[i*3+2] << endl;
    }
    return;
}
