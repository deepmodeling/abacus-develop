//==========================================================
// AUTHOR : mohan
// DATE : 2009-03-18 ( version 1.0 )
// 2009-04-14 ( version 1.1 )
// 2009-05-19 change file name
//==========================================================
#include "Calculate_C4.h"
#include "tools.h"

int Calculate_C4::test = 0;
Calculate_C4::Calculate_C4(){}
Calculate_C4::~Calculate_C4(){}

//===================================================================
// FUNCTION : Find_Eigenvalues
// DESCRIPTION : Find 'Enumber' number of eigenvalues of 
// spherical bessel function from '0~lmax', and put the
// result eigenvalues in 'en' matrix. 
// P.S: (1) The eigenvalues are calculated '< tolerence'
// accuracy.
// P.S: (2) The rcut is important for cacluating eigenvalues
// because it makes jl( rcut * eigenvalue ) = 0.00
//====================================================================
void Calculate_C4::Find_Eigenvalues(const int &Enumber, matrix &en, 
		const int &lmax, const double &rcut, const double &tolerence)
{
	double *e = new double[Enumber];
	for(int l=0; l<lmax+1; l++)
	{
//		cout << "\n L=" << l;
		Mathzone::Spherical_Bessel_Roots(Enumber, l, tolerence, e, rcut);
		for(int i=0; i<Enumber; i++)
		{
			en(l, i) = e[i];
//			cout << "\n" << e[i];
		}
//		cout << "\n Emax = " << e[Enumber-1] << " (a.u.)";
	}
	delete[] e;
	return;
}

//===================================================================
// FUNCTION : Test_Orthogonal
// DESCRIPTION : we know the formulat
// int_{0}^{Rc} dr r^{2} jl{q1*r} * jl(q1*r) = some double value
// int_{0}^{Rc} dr r^{2} jl{q1*r} * jl(q2*r) = 0.00;
// so we choose q1 as test_ei, and calculated the 'simpson' 
// integral result to test the orthogonal relationship.
// P.S: (1) Here we need only the l as input 'L'.
// P.S: (2) choose eigenvalues q2 from matrix 'en' and 'Enumber'
// P.S: (3) We calculates the intetral points through input
// 'rcut' and 'dr'
// P.S: (4) 
//===================================================================
void Calculate_C4::Test_Orthogonal(const double &rcut, const double &dr, 
		const int &L, const matrix &en, const int &Enumber, const int &test_ei)
{
	int mesh = static_cast<int>( rcut / dr)+1;
	if(mesh % 2 == 0) mesh++;
	cout << "\n\n Mesh for test = " << mesh << " L = " << L;

	double *r = new double[mesh];
	double *rab = new double[mesh];
	double *e = new double[Enumber];
	double *jle = new double[mesh];
	double *jtest = new double[mesh];
	double *function = new double[mesh];

	for(int i=0; i<Enumber; i++)
	{
		e[i] = en(L, i);
	}

	for(int ir=0; ir< mesh; ir++)
	{
		r[ir] = ir*dr;
		rab[ir] = dr;
	}

	Mathzone::Spherical_Bessel(mesh, r, e[ test_ei ], L, jtest);

	cout << "\n r=" << r[mesh-1] << "(a.u.)  jtest=" << jtest[mesh-1];
	
	double C4=0.0;
	for(int i=0; i<Enumber; i++)
	{
		Mathzone::Spherical_Bessel(mesh, r, e[i], L, jle);
		for(int ir=0; ir<mesh; ir++)
		{
			function[ir] = jle[ir] * jtest[ir] * r[ir] * r[ir]; 
		}
		Mathzone::Simpson_Integral(mesh, function, rab, C4 );
//		if( abs(C4) < 10e-5 ) 
//		{
		//	cout << "\n" << setw(5) << i << setw(15) << "0";
//		}
//		else
//		{
			cout << "\n" << setw(5) << i << setw(15) << C4;
//		}
	}		

	delete[] function;
	delete[] jtest;
	delete[] r;
	delete[] rab;
	delete[] e;
	delete[] jle;
	return;
}

//===================================================================
// FUNCTION : key_solver
// DESCRIPTION : Calculate C4 coefficients.
// (1) Read in psi(r) (containing *r)  from file 'name'
// (2) The input L must input calfully, the same as the psi(r).
// (3) For each eigenvalue, calculate the Jl(q1*r)
// (4) use int_{0}^{Rc} psi(r)Jl(q1*r)*r dr ==> C4_firstly
// (5) ues int_{0}^{Rc} J1(q1*r)*Jl(q1*r)*r*r dr ==> norm
// (6) C4 = C4_firstly / norm
//===================================================================
void Calculate_C4::key_solver(const string &name, const int &L, const matrix &en, const int &Enumber, double *C4)
{
	ifstream ifs(name.c_str());
	if(!ifs)
	{	
		cout << "\n Can't find file : " << name;
	}	
	else
	{
		cout << "\n Open file " << name;
	}

	string title;

	int mesh;
	ifs >> title >> mesh;
	assert(mesh > 0);
	cout << "\t mesh = "<<mesh;

	double *r = new double[mesh];
	double *psi = new double[mesh];
	double *rab = new double[mesh];
	double *e = new double[Enumber];
	double *jle = new double[mesh];
	double *function = new double[mesh];
	double *function2 = new double[mesh];

	for(int i=0; i<Enumber; i++)
	{
		e[i] = en(L, i);
	}

	ifs >> title >> title >> title;

	for(int i=0; i<mesh; i++)
	{
		ifs >> r[i] >> psi[i] >> rab[i];
	}

	double norm = 0.0;
	for(int i=0; i<Enumber; i++)
	{
		Mathzone::Spherical_Bessel(mesh, r, e[i], L, jle);
		for(int ir=0; ir<mesh; ir++)
		{
			function[ir] = jle[ir] * psi[ir] * r[ir]; 
		}

		for(int ir=0; ir<mesh; ir++)
		{
			function2[ir] = jle[ir] * jle[ir] * r[ir] * r[ir]; 
		}

		Mathzone::Simpson_Integral(mesh, function, rab, C4[i] );
		Mathzone::Simpson_Integral(mesh, function2, rab, norm );
		assert(norm != 0.0);
		C4[i] /= norm;
//		cout << "\n norm = "<<norm;
	}		

	delete[] function2;
	delete[] function;
	delete[] r;
	delete[] psi;
	delete[] rab;
	delete[] e;
	delete[] jle;
	ifs.close();
	return;
}

void Calculate_C4::norm_ic( realArray &C4, const int &ntype, const int &lmax, 
		const int &nmax, const int &enumber, const double &tolerence,
		const double &rcut, const double &dr)
{
//	TITLE("Calculate_C4","norm_ic");
	assert(lmax<5);
	matrix en( lmax+1, enumber);
	Find_Eigenvalues(enumber, en, lmax, rcut, tolerence);

	// init mesh
	int mesh = static_cast<int>(rcut / dr)+1; //mohan update 2009-10-10
	if(mesh % 2 == 0) mesh++;

	if(test==1)
	{
		cout << "\n norm : mesh = " << mesh;
		cout << "\n type = " << ntype;
		cout << "\n lmax = " << lmax;
		cout << "\n enumber = " << enumber;
	}
	
	// init radial mesh
	double *r = new double[mesh];
	for(int i=0; i<mesh; i++)
	{
		r[i] = i * dr;
	}

	// init relta r, used in simpson integral
	double *rab = new double[mesh];
	for(int i=0; i<mesh; i++)
	{
		rab[i] = dr;
	}

	double *e = new double[enumber];
	double *psi = new double[mesh];
	double *function = new double[mesh];
	double *jle = new double[mesh];

	for(int it=0; it<ntype; it++)
	{
		for(int l=0;l<lmax+1; l++)
		{
			// init eigenvalue of order l
			for(int i=0; i<enumber; i++)
			{
				e[i] = en(l, i);
			}
			for(int n=0; n<nmax; n++)
			{
				ZEROS(jle, mesh);
				ZEROS(psi, mesh);
				for(int ie=0; ie<enumber; ie++)
				{
					Mathzone::Spherical_Bessel(mesh, r, e[ie], l, jle);
					for(int ir=0; ir<mesh; ir++)
					{
						psi[ir] += C4(it,l,n,ie) * jle[ir];
					}
				}

				for(int ir=0; ir<mesh; ir++)
				{
					function[ir] = psi[ir] * psi[ir] * r[ir] * r[ir];
				}
				double norm = 0.0;
				Mathzone::Simpson_Integral(mesh, function, rab, norm);

				if(test==2) cout << "\n it=" << it << " L=" << l << " N=" << n << " norm = " << norm;
				norm = sqrt(norm);
				assert(norm > 1.0e-5);
				for(int ie=0; ie<enumber; ie++)
				{
					C4(it,l,n,ie)/=norm;
				}
			}
		}
	}

	delete[] r;
	delete[] psi;
	delete[] e;
	delete[] jle;
	return;
}

void Calculate_C4::plot(ifstream &ifs)
{
	TITLE("Calculate_C4","plot");
	int nfiles;
	double dr;
	READ_VALUE( ifs, nfiles);
	READ_VALUE( ifs, dr);
	cout << "\n Number of files to plot : " << nfiles;
	cout << "\n dr = " << dr;
	assert( nfiles > 0);
	assert( dr > 0.0 );

	int l;
	double ecut;
	double rcut;
	int enumber = 0;
	double tolerence;
	double *c4;
	double *r;
	double *e;
	double *psi;
	double *jle;
	for(int i=0; i<nfiles; i++)
	{
		string filenames;
		READ_VALUE( ifs, filenames );
		ifstream ifs( filenames.c_str() );
		if(!ifs)
		{
			cout << "\n" << "Can't find file : " << filenames;
			WARNING_QUIT("Calculate_C4::plot","can't find file.");	
		}
		else
		{
			cout << "\n Open file " << filenames;
		}

		if( SCAN_BEGIN(ifs,"<HEAD>") )
		{
			string filedir;
			READ_VALUE( ifs, filedir );
			READ_VALUE( ifs, l );
			READ_VALUE( ifs, ecut );
			READ_VALUE( ifs, rcut );
			READ_VALUE( ifs, enumber );
			READ_VALUE( ifs, tolerence );
			cout << "\n l = " << l;
			cout << "\n ecut = " << ecut;
			cout << "\n enumber = " << enumber;
			SCAN_END(ifs, "</HEAD>");
		}

		// init c4
		c4 = new double[enumber];
		ZEROS(c4, enumber);
		if( SCAN_BEGIN(ifs,"<C4>") )
		{
			assert(enumber>0);
			for(int i=0; i<enumber; i++)
			{
				int id;
				ifs >> id >> c4[i];
		//		cout << "\n c4[" << i << "]=" << c4[i];
			}
			SCAN_END(ifs, "</C4>");
		}

		// init mesh
		int mesh = static_cast<int>(rcut / dr)+1;
		if(mesh % 2 == 0) mesh++;
		cout << "\n mesh = " << mesh;

		// init r
		r = new double[mesh];
		for(int i=0; i<mesh; i++)
		{
			r[i] = dr * i;
		}

		// init eigenvalue en
		assert(l<3);
		matrix en( l+1, enumber);
		Find_Eigenvalues(enumber, en, l, rcut, tolerence);
		
		// init e
		e = new double[enumber];
		for(int i=0; i<enumber; i++)
		{
			e[i] = en(l, i);
		}

		// init psi
		psi = new double[mesh];
		ZEROS(psi, mesh);

		// init jle
		jle = new double[mesh];

		// sum all jle
		for(int ie=0; ie<enumber; ie++)
		{
			ZEROS(jle, mesh);
			Mathzone::Spherical_Bessel(mesh, r, e[ie], l, jle);
			for(int ir=0; ir<mesh; ir++)
			{
				psi[ir] += c4[ie] * jle[ir];
			}
		}

		// unit
		double norm;
		double *function = new double[mesh];
		double *rab = new double[mesh];
		for(int ir=0; ir<mesh; ir++)
		{
			rab[ir] = dr;
			function[ir] = psi[ir]*psi[ir]*r[ir]*r[ir];
		}
			
		Mathzone::Simpson_Integral(mesh, function, rab, norm );
		cout << "\n norm = " << norm;
		assert(norm > 1.0e-5);// small number, 
		norm = sqrt(norm);

		stringstream ss;
		ss << filenames << ".plot.txt";

		ofstream ofs(ss.str().c_str());
		for(int ir=0; ir<mesh; ir++)
		{
			ofs << ir*dr << " " << psi[ir]/norm << endl; 
		}

		delete[] c4;
		delete[] r;
		delete[] e;
		delete[] psi;
		delete[] jle;
		delete[] function;
		delete[] rab;
	}

	exit(0);
	return;
}

void Calculate_C4::Test(ifstream &ifs, const double &ecut, const double &tolerence )
{
	TITLE("Calculate_C4","Test");
	double rcut;// rcut for input orbitals
	double dr;// dr for simpson integrals( only used in test case),
				// in output C4 case : dr is read from file as 'rab'
	int test_ei;// the eigenvalue we want to test
	int lmax; // max value of L

	READ_VALUE( ifs, rcut );
	READ_VALUE( ifs, dr );
	READ_VALUE( ifs, test_ei);
	READ_VALUE( ifs, lmax );

	cout << "\n Energy cutoff = "<< ecut;
	cout << "\n tolerence = " << tolerence;
	cout << "\n rcut = "<<rcut;
	cout << "\n dr = " <<dr;
	cout << "\n test_ei = " << test_ei;
	cout << "\n Lmax = "<<lmax;

	assert(rcut > 0.0);
	assert(dr > 0.0);
	assert(test_ei >= 0);
	assert(ecut > 0);
	assert(lmax > -1);
	assert(tolerence > 0.0);

	//================================================
	// Calculate the eigenvalues number from ecut.
	//================================================
	const int Enumber = static_cast<int>( sqrt( ecut * 2) * rcut / Mathzone::PI );
	cout << "\n Number of eigenvalues of Jl = " << Enumber;
	//================================================
	// Calculate eigenvalue of Jl
	// CALL FUNCTION : Find_Eigenvalues
	//================================================
	matrix en( lmax+1, Enumber);
	Find_Eigenvalues(Enumber, en, lmax, rcut, tolerence);

	for(int l=0; l<lmax+1; l++)
	{
		Test_Orthogonal(rcut, dr, l, en, Enumber, test_ei);
	}
	return;
}

void Calculate_C4::init(ifstream &ifs, const double &ecut, const double &rcut, const double &tolerence, const int &lmax )
{
	cout << "\n\n ==> get_C4 Program Start." << endl;

	cout << "\n rcut      = " << rcut << " (a.u.)";
	cout << "\n ecut      = " << ecut << " (Ry)";
	cout << "\n tolerence = " << tolerence;
	cout << "\n lmax      = " << lmax;

	int nw;//number of orbitals(files) to calculate C4. 
	string *name;// file name, containing psi(r)*r, rab, radial
	int *L;//angular momentum for each file( each orbital);

	READ_VALUE( ifs, nw );
	assert(nw>0);
	cout << "\n Files     = " << nw;

	name = new string[nw];
	L = new int[nw];
	//================================================
	// read in the file name and the corresponding L.
	//================================================
	for(int i=0; i<nw; i++)
	{
		READ_VALUE(ifs, name[i]);
		READ_VALUE(ifs, L[i]);
	}

	//================================================
	// Calculate the eigenvalues number from ecut.
	//================================================
	const int Enumber = static_cast<int>( sqrt( ecut * 2) * rcut / Mathzone::PI );
	cout << "\n Eigenvalue number(Jl) = " << Enumber;

	//================================================
	// Calculate eigenvalue of Jl
	// CALL FUNCTION : Find_Eigenvalues
	//================================================
	matrix en( lmax+1, Enumber);
	Find_Eigenvalues(Enumber, en, lmax, rcut, tolerence);

	for(int i=0; i<nw; i++)
	{
		double *C4 = new double[Enumber];
		//===========================================
		// CALL FUNCTION : Calculate
		// calculate C4.
		//===========================================
		this->key_solver(name[i], L[i], en, Enumber, C4);
		stringstream ss;
		ss << name[i] << ".C4";
		ofstream ofs( ss.str().c_str() );
			
		ofs << "Calculate_C4_Program";
		ofs << "\n<VERSION>";
		ofs << "\n AUTHOR : Mohan Chen";
		ofs << "\n Date : 2009-4-1";
		ofs << "\n LOCATION : LQCC, Hefei, China";
		ofs << "\n EMAIL : mohan@mail.ustc.edu.cn";
		ofs << "\n Description : Calculate the coefficients C4 of f(r) in Spherical Bessel Basis J(qr).";
		ofs << "\n Formula : C4 = integral(r)[ f(r)*jl(qr)*r^{2} ]dr ";
		ofs << "\n P.S. : We default consider f(r) read from file is in the form : ( f(r) * r ).";
		ofs << "\n</VERSION>";

		ofs << "\n\n<HEAD>";
		ofs << "\n" << name[i] << " ( f(r) read from this file)";
		ofs << "\n" << setw(8) << L[i] << " angular momentum for this orbital.";
		ofs << "\n" << setw(8) << ecut << " Energy cutoff(Ry.).";
		ofs << "\n" << setw(8) << rcut << " rcut (a.u.)";
		ofs << "\n" << setw(8) << Enumber << " eigenvalue number( sqrt(ecut*2)*rcut/PI ).";
		ofs << "\n" << setw(8) << tolerence << " tolerence to calculate eigenvalue.";
		ofs << "\n</HEAD>";

		ofs << "\n\n<C4>";
		for(int j=0; j<Enumber; j++)
		{
			ofs << "\n" << setw(5) << j << setiosflags(ios::scientific) << setprecision(8) << setw(25) << C4[j];
		}
		ofs << "\n</C4>";
		ofs.close();
		delete[] C4;
	}

	delete[] name;
	delete[] L;
}
