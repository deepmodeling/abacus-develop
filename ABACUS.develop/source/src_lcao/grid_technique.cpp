#include "grid_technique.h"
#include "../src_pw/global.h"

#include "global_fp.h" // mohan add 2021-01-30

Grid_Technique GridT;

Grid_Technique::Grid_Technique()
{
	this->how_many_atoms = new int[1];
	this->which_atom = new int[1];
	this->which_bigcell = new int[1];
	this->which_unitcell = new int[1];
	this->bcell_start = new int[1];
	this->in_this_processor = new bool[1];
	this->trace_lo = new int[1];
	this->trace_beta = new int[1];//sun zhiyuan add

	this->total_atoms_on_grid = 0;
}


Grid_Technique::~Grid_Technique()
{
	delete[] how_many_atoms;
	delete[] which_atom;
	delete[] which_bigcell;
	delete[] which_unitcell;
	delete[] bcell_start;
	delete[] in_this_processor;
	delete[] trace_lo;
	delete[] trace_beta; //sun zhiyuan add
}


// This function is called in LOOP_ions.cpp
// after the orbital information has been read,
// this function control the routinue to generate 
// grid technique parameters.
void Grid_Technique::set_pbc_grid(
		const int &ncx_in,
		const int &ncy_in,
		const int &ncz_in,
		const int &bx_in,
		const int &by_in,
		const int &bz_in,
		const int &nbx_in,
		const int &nby_in,
		const int &nbz_in,
		const int &nbxx_in,
		const int &nbzp_start_in,
		const int &nbzp_in)
{
	TITLE("Grid_Technique","init");
	timer::tick("Grid_Technique","init",'D');

	if(OUT_LEVEL != "m") 
	{
		ofs_running << "\n SETUP EXTENDED REAL SPACE GRID FOR GRID INTEGRATION" << endl;
	}

	// (1) init_meshcell cell and big cell.
	this->set_grid_dim(
		ncx_in,ncy_in,ncz_in,
		bx_in,by_in,bz_in,
		nbx_in,nby_in,nbz_in,
		nbxx_in,nbzp_start_in,nbzp_in);

	this->init_latvec();

	this->init_big_latvec();

	this->init_meshcell_pos();

	// (2) expand the grid
	this->init_grid_expansion();

	// (3) calculate the extended grid.
	this->cal_extended_cell(this->dxe, this->dye, this->dze);

	this->init_tau_in_bigcell();

	// init meshball
	this->delete_meshball_positions(); //LiuXh add 2018-12-14

	this->init_meshball();

	this->init_atoms_on_grid();	

	this->cal_trace_lo();

	timer::tick("Grid_Technique","init",'D');
	return;
}

void Grid_Technique::init_atoms_on_grid(void)
{
	TITLE("Grid_Technique","init_atoms_on_grid");

	assert(nbxx>0);
	
	// (1) prepare data. 
	// counting the number of atoms whose orbitals have
	// values on the bigcell.
	delete[] this->how_many_atoms;
	this->how_many_atoms = new int[nbxx];
	ZEROS(how_many_atoms, nbxx);
	Memory::record("atoms_on_grid","how_many_atoms",nbxx,"int");

	// (2) start z and ended z,
	// consistent with division of FFT grid.
	// mohan add 2010-07-01
	const int zstart = nbzp_start;
	const int zend = nbzp + zstart;
	if(test_gridt)OUT(ofs_running,"zstart",zstart);
	if(test_gridt)OUT(ofs_running,"zend",zend);
	int iz_now = -1;
	
	// (3) information about gloabl grid
	// and local grid.
	// mohan add 2010-07-02
	int ix,iy,iz;	
	const int nbyz = nby * nbz;
	const int nz = nbzp;
	int f; // normal local form.
	
	// (4) Find the atoms using
	// when doing grid integration. 
	delete[] in_this_processor;
	this->in_this_processor = new bool[ucell.nat];
	for(int i=0; i<ucell.nat; i++)
	{
		in_this_processor[i] = false;
	}

	// init atoms on grid	
	assert( this->nxyze > 0);
	int* index2normal = new int[this->nxyze];
	assert( index2normal != NULL );
	Memory::record("Grid_Meshcell","index2normal",this->nxyze,"int");
	this->grid_expansion_index(1,index2normal); 

 	
	// (5) record how many atoms on
	// each local grid point (ix,iy,iz)
	int iat=0;
	int normal;
	this->total_atoms_on_grid = 0;
	for(int it=0; it<ucell.ntype; it++)
	{
		for(int ia=0; ia<ucell.atoms[it].na; ia++)
		{
			for(int im=0; im< this->meshball_ncells; im++)
			{
				// bcell[iat]: which bcell iat atom is in.
				// ball[im]: relative position of adjacent bcell.
				normal = index2normal[ this->index_atom[iat] + this->index_ball[im] ];

				if(normal >= nbxyz)
				{
					cout << " index_atom=" << index_atom[iat] << endl;
					cout << " index_ball=" << index_ball[im] << endl;
					cout << " normal=" << normal << endl;
					cout << " nbxyz=" << nbxyz << endl;
					WARNING_QUIT("Grid_Technique::init_atoms_on_grid","normal >= nbxyz");
				}

				assert(normal>=0);	

				// mohan add 2010-07-01 part1 
				iz_now = normal % nbz;
				if(iz_now<zstart)continue;
				else if(iz_now>=zend)continue;
				
				// mohan add 2010-07-01 part2
				ix = normal / nbyz;
				iy = ( normal - ix * nbyz ) / nbz;
				iz = iz_now - zstart;
				f = ix * nby * nz + iy * nz + iz;
					
				++how_many_atoms[f];
				++total_atoms_on_grid;

				this->in_this_processor[iat] = true;
			}
			++iat;
		}
	}

	if(test_gridt)OUT(ofs_running,"Total_atoms_on_grid",total_atoms_on_grid);
	
	int stop = 0;
	if(total_atoms_on_grid == 0)
	{
		ofs_running << " No atoms on this sub-FFT-mesh." << endl;
		stop = 1;
	} 	
	Parallel_Reduce::reduce_int_all( stop );
	if(stop)
	{
		WARNING("Grid_Technique::init_atoms_on_grid","No atom on this sub-FFT-mesh.");
	}

	// need how_many_atoms first.
	this->cal_grid_integration_index();
	// bcell_start is needed.
	this->init_atoms_on_grid2(index2normal);
	delete[] index2normal;	
	return;
}

void Grid_Technique::init_atoms_on_grid2(const int* index2normal)
{	
	TITLE("Grid_Techinique","init_atoms_on_grid2");

	if(total_atoms_on_grid==0) 
	{	
		WARNING("Grid_Technique::init_atoms_on_grid2","no atom on this sub FFT grid.");
		return;
	}

	int* index2ucell = new int[this->nxyze];
	assert( index2ucell != NULL );
	Memory::record("Grid_Meshcell","index2ucell",this->nxyze,"int");	
	this->grid_expansion_index(0,index2ucell);
	
	const int zstart = nbzp_start;
	const int zend = nbzp + zstart;
	int iz_now = -1;
	int ix,iy,iz;	
	const int nbyz = nby * nbz;
	const int nz = nbzp;
	int f; // normal local form.

	//--------------------------------------
	// save which atom is in the bigcell.
	//--------------------------------------
	delete[] which_atom;
	this->which_atom = new int[total_atoms_on_grid];
	assert( which_atom != 0);
	Memory::record("atoms_on_grid","which_atom",total_atoms_on_grid,"int");

	delete[] which_bigcell;
	this->which_bigcell = new int[total_atoms_on_grid];
	assert( which_bigcell != 0);
	Memory::record("atoms_on_grid","which_bigcell",total_atoms_on_grid,"int");

	delete[] which_unitcell;
	this->which_unitcell = new int[total_atoms_on_grid];
	assert( which_unitcell != 0);
	Memory::record("atoms_on_grid","which_unitcell",total_atoms_on_grid,"int");
	// for each atom, first we need to locate which cell
	// the atom is in, then we search meshball aroung this
	// grid, and record each grid's atom position.
	int count = 0;
	int iat = 0;
	ZEROS(this->how_many_atoms, nbxx);
	for(int it=0; it<ucell.ntype; it++)
	{
		for(int ia=0; ia<ucell.atoms[it].na; ia++)
		{
			// zero bigcell of meshball indicate ?
			for(int im=0; im< this->meshball_ncells; im++)
			{
				const int extgrid = this->index_atom[iat] + this->index_ball[im];
				const int normal = index2normal[ extgrid ];
			
				// mohan add 2010-07-01
				iz_now = normal % nbz;
				if(iz_now<zstart)continue;
				else if(iz_now>=zend)continue;
			
				// mohan add 2010-07-01
				ix = normal / nbyz;
				iy = ( normal - ix * nbyz ) / nbz;
				iz = iz_now - zstart;
				f = ix * nby * nz + iy * nz + iz;
				
				// it's not the normal order to calculate which_atom
				// and which_bigcell, especailly in 1D array. 
				// Each grid's adjacent atom number is different,
				// so, first we need to locate which grid, using
				// bcell_start, then we need to count which adjacent atom.
				// using how_many_atoms.
				int index = this->bcell_start[f] + this->how_many_atoms[f];
				
				// we save which_atom and which_bigcell in 1D array,
				// once you want to use this in grid integration, 
				// the only information you got is the 'normal' index, 
				// so you need to use bcell_start
				// to get the 'mesh_index', then you can you this mesh_index
				// to use which_atom or which_bigcell.
				this->which_atom[ index ] = iat;
				this->which_bigcell[ index ] = im;
				this->which_unitcell[ index ] = index2ucell[extgrid];

		//		if(im==13651)
		//		{
		//			cout << " which_unitcell=" << which_unitcell[index] << endl;
		//		}
				 
				++how_many_atoms[f];
				++count;
			}
			++iat;
		}
	}
	assert( count == total_atoms_on_grid );
	delete[] index2ucell;
	return;
}

void Grid_Technique::cal_grid_integration_index(void)
{
	// save the start 
	delete[] this->bcell_start;
	this->bcell_start = new int[nbxx];
	Memory::record("atoms_on_grid","bcell_start",nbxx,"int");
	this->bcell_start[0] = 0;
	for(int i=1; i<nbxx; i++)
	{
		this->bcell_start[i] = this->bcell_start[i-1] + this->how_many_atoms[i-1];
	}
	// calculate which grid has the largest number of atoms,
	// and how many atoms.
	this->max_atom = 0;
	for(int i=0; i<nbxx; i++)
	{
		this->max_atom = std::max( this->max_atom, this->how_many_atoms[i]);
	}

#ifdef __MPI
	int* all = new int[NPROC];
	ZEROS(all, NPROC);
	Parallel_Reduce::gather_int_all(max_atom,all);
	if(MY_RANK==0)
	{
		ofs_warning << setw(15) << "Processor" << setw(15) << "Atom" << endl;
		for(int i=0; i<NPROC; i++)
		{
			ofs_warning << setw(15) << i+1 << setw(15) << all[i] << endl; 
		}
	}
	delete[] all;
#endif

	if(test_gridt)OUT(ofs_running,"Max atom on bigcell",max_atom);
	return;
}

void Grid_Technique::cal_trace_beta(void)
{
	// save the atom information in trace_beta//
	delete[] trace_beta;
	int nkb=ORB.nkb;
	this->trace_beta = new int[nkb];
	for(int i=0; i<nkb; i++)
	{
		this->trace_beta[i] = -1;
	}
	this->lgbeta = 0;
	int iat = 0;
	int ih_all = 0;
	int ih_local = 0;

	ofs_running << "trace_beta" << endl;
	for(int it=0; it<ucell.ntype; ++it)
	{
		Atom* atom = &ucell.atoms[it];
		for(int ia=0; ia<atom->na; ++ia)
		{
			if(this->in_this_processor[iat])
			{
				for(int ih=0; ih<atom->nh; ih++)
				{
					this->trace_beta[ih_all] = ih_local;
					ofs_running << setw(5) << ih_all << setw(15) << trace_beta[ih_all] << endl;
					++ih_local;
					++ih_all;
				}
			}
			else
			{
				ih_all += atom->nh;
			}
			++iat;
		}
	}
	this->lgbeta = ih_local;
	assert(lgbeta>=0);
	OUT(ofs_running,"lgbeta",lgbeta);

	return;
}


void Grid_Technique::cal_trace_lo(void)
{	
	TITLE("Grid_Technique","cal_trace_lo");
	// save the atom information in trace_lo,
	// in fact the trace_lo dimension can be reduced
	// to ucell.nat, but I think this is another way.
	delete[] trace_lo;
	this->trace_lo = new int[NLOCAL];
	for(int i=0; i<NLOCAL; i++)
	{
		this->trace_lo[i] = -1;
	}
	Memory::record("atoms_on_grid","trace_lo",NLOCAL,"int");

	this->lnat = 0;
	this->lgd = 0;
	int iat = 0;
	int iw_all=0;
	int iw_local=0;

	for(int it=0; it<ucell.ntype; it++)
	{
		for(int ia=0; ia<ucell.atoms[it].na; ia++)
		{
			if(this->in_this_processor[iat])
			{
				++lnat;
				int nw0 = ucell.atoms[it].nw;
				if(NSPIN==4)
				{//added by zhengdy-soc, need to be double in soc
					nw0 *= 2;
					lgd += nw0;
				}
				else
					lgd += ucell.atoms[it].nw;
				
				for(int iw=0; iw<nw0; iw++)
				{
					this->trace_lo[iw_all] = iw_local;
					++iw_local; 
					++iw_all;
				}
			}
			else
			{
				iw_all += ucell.atoms[it].nw;
				if(NSPIN==4) iw_all += ucell.atoms[it].nw;
			}
			++iat;
		}
	}
	
	//------------
	// for test
	//------------
//	for(int i=0; i<NLOCAL; ++i)
//	{
//		ofs_running << " i=" << i+1 << " trace_lo=" << trace_lo[i] << endl;
//	}

	if(OUT_LEVEL != "m") 
	{
		OUT(ofs_running,"Atom number in sub-FFT-grid",lnat);
		OUT(ofs_running,"Local orbitals number in sub-FFT-grid",lgd);
	}

	assert(iw_local == lgd);
	assert(iw_all == NLOCAL);
	return;
}


