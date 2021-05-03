#include "global.h"
#include "atom_pseudo.h"
#include "tools.h"

Atom_pseudo::Atom_pseudo()
{
	mbl = new Vector3<int>[1];
	pseudo_fn = "not_init";
	mass = 0.0;
}

Atom_pseudo::~Atom_pseudo()
{
	delete[] mbl;
}

void Atom_pseudo::print_atom(ofstream &ofs)
{
	if(test_atom) TITLE("atom_pseudo","print_atom");

	OUT(ofs,"mass",mass);
	OUT(ofs,"pseudo_fn",pseudo_fn);

	return;
}

#ifdef __MPI
void Atom_pseudo::bcast_atom_pseudo(const int &na)
{
	TITLE("Atom_pseudo","bcast_atom_pseudo");
	Parallel_Common::bcast_double( mass );
	Parallel_Common::bcast_string( pseudo_fn );

	if(MY_RANK!=0)
	{
		mbl = new Vector3<int>[na];
	}

	for(int i=0;i<na;i++)
	{
		Parallel_Common::bcast_int( mbl[i].x );
		Parallel_Common::bcast_int( mbl[i].y );
		Parallel_Common::bcast_int( mbl[i].z );
	}
	return;
}

void Atom_pseudo::bcast_atom_pseudo2(void)
{
	TITLE("Atom_pseudo","bcast_atom_pseudo2");
// == pseudo_h ==
//int
	Parallel_Common::bcast_int( lmax );
	Parallel_Common::bcast_int( mesh );
	Parallel_Common::bcast_int( nchi );
	Parallel_Common::bcast_int( nbeta );
	Parallel_Common::bcast_int( nv );
	Parallel_Common::bcast_int( zv );

//double
	Parallel_Common::bcast_double( etotps );
	Parallel_Common::bcast_double( ecutwfc );
	Parallel_Common::bcast_double( ecutrho );

// bool
	Parallel_Common::bcast_bool( tvanp );
	Parallel_Common::bcast_bool( nlcc );
	Parallel_Common::bcast_bool( has_so );

//string
	Parallel_Common::bcast_string( psd );
	Parallel_Common::bcast_string( pp_type );
	Parallel_Common::bcast_string( dft, 4 );

	if(MY_RANK!=0)
	{
		jjj = new double [nbeta];
		els = new string[nchi];
		lchi = new int [nchi];
		oc = new double[nchi];
		jchi = new double[nchi];
		nn = new int[nchi];
	}
	
	Parallel_Common::bcast_double( jjj, nbeta );
	Parallel_Common::bcast_string( els, nchi);
	Parallel_Common::bcast_int( lchi, nchi);
	Parallel_Common::bcast_double( oc, nchi);
	Parallel_Common::bcast_double( jchi, nchi);
	Parallel_Common::bcast_int( nn, nchi);
// == end of pseudo_h

// == pseudo_atom ==
	Parallel_Common::bcast_int( msh );
	Parallel_Common::bcast_double( rcut );
	if(MY_RANK != 0)
	{
		assert(mesh!=0);
		r = new double[mesh];
		rab = new double[mesh];
		rho_atc = new double[mesh];
		rho_at = new double[mesh];
		chi.create( nchi,mesh );
	}
	
	Parallel_Common::bcast_double( r, mesh );
	Parallel_Common::bcast_double( rab, mesh );
	Parallel_Common::bcast_double( rho_atc, mesh );
	Parallel_Common::bcast_double( rho_at, mesh );
	Parallel_Common::bcast_double( chi.c, nchi * mesh );
// == end of pseudo_atom ==
	
// == pseudo_vl ==
	if(MY_RANK != 0)
	{
		vloc_at = new double[mesh];
	}
	Parallel_Common::bcast_double( vloc_at, mesh);
// == end of pseudo_vl ==

// == pseudo_nc ==
	if(MY_RANK != 0)
	{
		lll = new int[nbeta];
	}
	Parallel_Common::bcast_int( lll, nbeta );
	Parallel_Common::bcast_int( kkbeta );
	Parallel_Common::bcast_int( nh );

	int nr,nc;
	if(MY_RANK == 0)
	{
		nr = betar.nr;
		nc = betar.nc;
	}
	Parallel_Common::bcast_int( nr );
	Parallel_Common::bcast_int( nc );

	if(MY_RANK != 0)
	{
		betar.create(nr,nc);
		dion.create(nbeta, nbeta);
	}

	// below two 'bcast_double' lines of codes seem to have bugs,
	// on some computers, the code will stuck here for ever.
	// mohan note 2021-04-28
	Parallel_Common::bcast_double( dion.c , nbeta * nbeta);
	Parallel_Common::bcast_double( betar.c, nr * nc );
// == end of psesudo_nc ==

	return;
}

#endif
