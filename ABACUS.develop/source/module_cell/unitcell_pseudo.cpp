#include "unitcell_pseudo.h"
#include "module_ORB/ORB_read.h" // to use 'ORB' -- mohan 2021-01-30
#include "src_pw/global.h"
#include <cstring>		// Peize Lin fix bug about strcmp 2016-08-02

UnitCell_pseudo::UnitCell_pseudo()
{
	if(test_pseudo_cell) TITLE("unitcell_pseudo","Constructor");
	set_atom_flag = false;
}

UnitCell_pseudo::~UnitCell_pseudo()
{
	if(set_atom_flag)
	{
		delete[] atoms;
	}
}

//==============================================================
//Calculate various lattice related quantities for given latvec
//==============================================================
void UnitCell_pseudo::setup_cell(
		const string &s_pseudopot_dir ,  
		const string &fn , ofstream &log)
{
	TITLE("UnitCell_pseudo","setup_cell");	
	// (1) init mag (global class)
	assert(ntype>0);
	delete[] mag.start_magnetization;
	mag.start_magnetization = new double[this->ntype];

	// (2) init *Atom class array.
	this->atoms = new Atom[this->ntype]; // atom species.
	this->set_atom_flag = true;


	bool ok = true;
	bool ok2 = true;

	// (3) read in atom information
	if(MY_RANK == 0)
	{
		// open "atom_unitcell" file.
		ifstream ifa(fn.c_str(), ios::in);
		if (!ifa)
		{
			ofs_warning << fn;
			ok = false;
		}

		if(ok)
		{

			ofs_running << "\n\n\n\n";
			ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
			ofs_running << " |                                                                    |" << endl;
			ofs_running << " | Reading atom information in unitcell:                              |" << endl;
			ofs_running << " | From the input file and the structure file we know the number of   |" << endl;
			ofs_running << " | different elments in this unitcell, then we list the detail        |" << endl;
			ofs_running << " | information for each element, especially the zeta and polar atomic |" << endl;
			ofs_running << " | orbital number for each element. The total atom number is counted. |" << endl;
			ofs_running << " | We calculate the nearest atom distance for each atom and show the  |" << endl;
			ofs_running << " | Cartesian and Direct coordinates for each atom. We list the file   |" << endl;
			ofs_running << " | address for atomic orbitals. The volume and the lattice vectors    |" << endl;
			ofs_running << " | in real and reciprocal space is also shown.                        |" << endl;
			ofs_running << " |                                                                    |" << endl;
			ofs_running << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
			ofs_running << "\n\n\n\n";

			ofs_running << " READING UNITCELL INFORMATION" << endl;
			//========================
			// call read_atom_species
			//========================
			this->read_atom_species(ifa);

			//==========================
			// call read_atom_positions
			//==========================
			ok2 = this->read_atom_positions(ifa);

			if(ok2)
			{
				for(int i=0;i<this->ntype;i++)
				{
					Global_File::make_dir_atom( this->atoms[i].label );
				}
			}
		}
	}
#ifdef __MPI
	Parallel_Common::bcast_bool(ok);
	Parallel_Common::bcast_bool(ok2);
	if(NSPIN==4) 
	{
		Parallel_Common::bcast_bool(DOMAG);
		Parallel_Common::bcast_bool(DOMAG_Z);
	}
#endif
	if(!ok)
	{
		WARNING_QUIT("UnitCell_pseudo::setup_cell","Can not find the file containing atom positions.!");
	}
	if(!ok2)
	{
		WARNING_QUIT("UnitCell_pseudo::setup_cell","Something wrong during read_atom_positions.");
	}

#ifdef __MPI
	this->bcast_unitcell_pseudo();

	// mohan add 2010-09-29
	ORB.bcast_files(ucell.ntype, MY_RANK);
#endif
	
	//========================================================
	// Calculate unit cell volume
	// the reason to calculate volume here is 
	// Firstly, latvec must be read in.
	//========================================================
	assert(lat0 > 0.0);
	this->omega = abs( latvec.Det() ) * this->lat0 * lat0 * lat0 ;
	if(this->omega<=0)
	{
		WARNING_QUIT("setup_cell","omega <= 0 .");
	}
	else
	{
		ofs_running << endl;
		OUT(ofs_running,"Volume (Bohr^3)", this->omega);
		OUT(ofs_running,"Volume (A^3)", this->omega * pow(BOHR_TO_A, 3));
	}
		
	//==========================================================
	// Calculate recip. lattice vectors and dot products
	// latvec have the unit of lat0, but G has the unit 2Pi/lat0
	//==========================================================
	this->GT = latvec.Inverse();
	this->G  = GT.Transpose();
	this->GGT = G * GT;
	this->invGGT = GGT.Inverse();

    //LiuXh add 20180515
    this->GT0 = latvec.Inverse();
    this->G0  = GT.Transpose();
    this->GGT0 = G * GT;
    this->invGGT0 = GGT.Inverse();

	ofs_running << endl;
	out.printM3(ofs_running,"Lattice vectors: (Cartesian coordinate: in unit of a_0)",latvec); 
	out.printM3(ofs_running,"Reciprocal vectors: (Cartesian coordinate: in unit of 2 pi/a_0)",G);
//	OUT(ofs_running,"lattice center x",latcenter.x);
//	OUT(ofs_running,"lattice center y",latcenter.y);
//	OUT(ofs_running,"lattice center z",latcenter.z);

	// read in non-local pseudopotential and ouput the projectors.

	ofs_running << "\n\n\n\n";
	ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
	ofs_running << " |                                                                    |" << endl;
	ofs_running << " | Reading pseudopotentials files:                                    |" << endl;
	ofs_running << " | The pseudopotential file is in UPF format. The 'NC' indicates that |" << endl;
	ofs_running << " | the type of pseudopotential is 'norm conserving'. Functional of    |" << endl;
	ofs_running << " | exchange and correlation is decided by 4 given parameters in UPF   |" << endl;
	ofs_running << " | file.  We also read in the 'core correction' if there exists.      |" << endl;
	ofs_running << " | Also we can read the valence electrons number and the maximal      |" << endl;
	ofs_running << " | angular momentum used in this pseudopotential. We also read in the |" << endl;
	ofs_running << " | trail wave function, trail atomic density and local-pseudopotential|" << endl;
	ofs_running << " | on logrithmic grid. The non-local pseudopotential projector is also|" << endl;
	ofs_running << " | read in if there is any.                                           |" << endl;
	ofs_running << " |                                                                    |" << endl;
	ofs_running << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
	ofs_running << "\n\n\n\n";


	this->read_cell_pseudopots(s_pseudopot_dir);
	
	if(MY_RANK == 0) 
	{
		for(int it=0; it<this->ntype; it++)
		{
			Atom* atom = &atoms[it];
			stringstream ss;
			ss << global_out_dir << atom->label 
				<< "/" << atom->label
				<< ".NONLOCAL";
			ofstream ofs(ss.str().c_str());

			ofs << "<HEADER>" << endl;
			ofs << setw(10) << atom->label << "\t" << "label" << endl;
			ofs << setw(10) << atom->pp_type << "\t" << "pseudopotential type" << endl;
			ofs << setw(10) << atom->lmax << "\t" << "lmax" << endl;
			ofs << "</HEADER>" << endl;

			ofs << "\n<DIJ>" << endl;
			ofs << setw(10) << atom->nbeta << "\t" << "nummber of projectors." << endl;
			
			for(int ib=0; ib<atom->nbeta; ib++)
			{
				for(int ib2=0; ib2<atom->nbeta; ib2++)
				{
					ofs<<setw(10) << atom->lll[ib] 
						<< " " << atom->lll[ib2]
						<< " " << atom->dion(ib,ib2)<<endl;
				}
			}
			ofs << "</DIJ>" << endl;
			
			for(int i=0; i<atom->nbeta; i++)
			{
				ofs << "<PP_BETA>" << endl;
				ofs << setw(10) << i << "\t" << "the index of projectors." <<endl;
				ofs << setw(10) << atom->lll[i] << "\t" << "the angular momentum." <<endl;

				// mohan add
				// only keep the nonzero part.
				int cut_mesh = atom->mesh; 
				for(int j=atom->mesh-1; j>=0; --j)
				{
					if( abs( atom->betar(i,j) ) > 1.0e-10 )
					{
						cut_mesh = j; 
						break;
					}
				}
				if(cut_mesh %2 == 0) ++cut_mesh;

				ofs << setw(10) << cut_mesh << "\t" << "the number of mesh points." << endl;


				for(int j=0; j<cut_mesh; ++j)
				{
					ofs << setw(15) << atom->r[j]
						<< setw(15) << atom->betar(i, j)
						<< setw(15) << atom->rab[j] << endl;
				}
				ofs << "</PP_BETA>" << endl;
			}

			ofs.close();
		}
	}

#ifdef __MPI
	this->bcast_unitcell_pseudo2();
#endif	

	for(int it=0; it<ntype; it++)
	{
		for(int j=0; j<4; j++)
		{
			if(atoms[0].dft[j]!=atoms[it].dft[j])
			{
				ofs_warning << "\n type " << atoms[0].label << " functional is " 
				<< atoms[0].dft[0] << " " << atoms[0].dft[1] << " "
				<< atoms[0].dft[2] << " " << atoms[0].dft[3];
				
				ofs_warning << "\n type " << atoms[it].label << " functional is " 
				<< atoms[it].dft[0] << " " << atoms[it].dft[1] << " "
				<< atoms[it].dft[2] << " " << atoms[it].dft[3] << endl;
				
				WARNING_QUIT("setup_cell","All DFT functional must consistent.");
			}
		}
	}

	// mohan add 2010-09-06
	// because the number of element type
	// will easily be ignored, so here
	// I warn the user again for each type.
	for(int it=0; it<ntype; it++)
	{
		xcf.which_dft(atoms[it].dft);
	}

	// setup the total number of PAOs
	this->cal_natomwfc();

	// setup NLOCAL
	this->cal_nwfc();

	// setup NBANDS
	this->cal_nelec();

	this->cal_meshx();

	// setup vdwd2 parameters
	vdwd2_para.initset(*this);		// Peize Lin add 2021.03.09

//	stringstream ss;
//	ss << global_out_dir << "unitcell_pp.log";
//	print_unitcell_pseudo( ss.str() );
	return;
}




//=========================================================
// calculate total number of electrons (nelec) and default
// number of bands (NBANDS).
//=========================================================
#include "src_pw/occupy.h"
void UnitCell_pseudo::cal_nelec(void)
{
	TITLE("UnitCell_pseudo","cal_nelec");
	//=======================================================
	// calculate the total number of electrons in the system
	// if nelec <>0; use input number (setup.f90)
	//=======================================================

	ofs_running << "\n SETUP THE ELECTRONS NUMBER" << endl;

	if (nelec == 0)
	{
		for (int it = 0; it < ntype;it++)
		{
			stringstream ss1, ss2;
			ss1 << "electron number of element " << atoms[it].label;
			const int nelec_it = atoms[it].zv * atoms[it].na;
			nelec += nelec_it;
			ss2 << "total electron number of element " << atoms[it].label; 
			
			OUT(ofs_running,ss1.str(),atoms[it].zv);
			OUT(ofs_running,ss2.str(),nelec_it);
		}
	}

	//OUT(ofs_running,"Total nelec",nelec);

	//=======================================
	// calculate number of bands (setup.f90)
	//=======================================
	double occupied_bands = static_cast<double>(nelec/DEGSPIN);	

	if( (occupied_bands - std::floor(occupied_bands)) > 0.0 )
	{
		occupied_bands = std::floor(occupied_bands) + 1.0; //mohan fix 2012-04-16
	}

	OUT(ofs_running,"occupied bands",occupied_bands);
	
	// mohan add 2010-09-04
    //cout << "nbands(ucell) = " <<NBANDS <<endl;
	if(NBANDS==occupied_bands)
	{
		if( Occupy::gauss() || Occupy::tetra() )
		{
			WARNING_QUIT("UnitCell_pseudo::cal_nelec","for smearing, num. of bands > num. of occupied bands");
		}
	}
	
	if ( CALCULATION!="scf-sto" && CALCULATION!="relax-sto" && CALCULATION!="md-sto" ) //qianrui 2021-2-20
	{
	if(NBANDS == 0)
	{
		if(NSPIN == 1)
		{
			int nbands1 = static_cast<int>(occupied_bands) + 10;
			int nbands2 = static_cast<int>(1.2 * occupied_bands);
			NBANDS = max(nbands1, nbands2);
		}
		else if (NSPIN ==2 || NSPIN == 4)
		{
			int nbands3 = nelec + 20;
			int nbands4 = 1.2 * nelec;
			NBANDS = max(nbands3, nbands4);
		}
		AUTO_SET("NBANDS",NBANDS);
	}
	//else if ( CALCULATION=="scf" || CALCULATION=="md" || CALCULATION=="relax") //pengfei 2014-10-13
	else
	{
		if(NBANDS < occupied_bands) WARNING_QUIT("unitcell","Too few bands!");
		if(NBANDS < mag.get_nelup() ) 
		{
			OUT(ofs_running,"nelup",mag.get_nelup());
			WARNING_QUIT("unitcell","Too few spin up bands!");
		}
		if(NBANDS < mag.get_neldw() )
        {
            WARNING_QUIT("unitcell","Too few spin down bands!");
        }
	}
	}

	// mohan update 2021-02-19
    // mohan add 2011-01-5
    if(BASIS_TYPE=="lcao" || BASIS_TYPE=="lcao_in_pw")
    {
        if( NBANDS > NLOCAL )
        {
            WARNING_QUIT("UnitCell_pseudo::cal_nwfc","NLOCAL < NBANDS");
        }
        else
        {
            OUT(ofs_running,"NLOCAL",NLOCAL);
            OUT(ofs_running,"NBANDS",NBANDS);
        }
    }

	OUT(ofs_running,"NBANDS",NBANDS);
	return;
}

//===========================================
// calculate the total number of local basis
// Target : nwfc, lmax,
// 			atoms[].stapos_wf
// 			NBANDS
//===========================================
void UnitCell_pseudo::cal_nwfc(void)
{
	TITLE("UnitCell_pseudo","cal_nwfc");
	assert(ntype>0);
	assert(nat>0);

	//===========================
	// (1) set iw2l, iw2n, iw2m
	//===========================
	for(int it=0; it<ntype; it++)
	{
		this->atoms[it].set_index();
	}

	//===========================
	// (2) set namax and nwmax
	//===========================
	this->namax = 0;
	this->nwmax = 0;
	for(int it=0; it<ntype; it++)
	{
		this->namax = std::max( atoms[it].na, namax );
		this->nwmax = std::max( atoms[it].nw, nwmax );
	}
	assert(namax>0);
// for tests
//		OUT(ofs_running,"max input atom number",namax);
//		OUT(ofs_running,"max wave function number",nwmax);	

	//===========================
	// (3) set nwfc and stapos_wf
	//===========================
	NLOCAL = 0;
	for(int it=0; it<ntype; it++)
	{
		atoms[it].stapos_wf = NLOCAL;
		const int nlocal_it = atoms[it].nw * atoms[it].na;
		if(NSPIN!=4) 
		{
			NLOCAL += nlocal_it;
		}
		else 
		{
			NLOCAL += nlocal_it * 2;//zhengdy-soc
		}

// for tests
//		OUT(ofs_running,ss1.str(),nlocal_it);
//		OUT(ofs_running,"start position of local orbitals",atoms[it].stapos_wf);
	}
	
	//OUT(ofs_running,"NLOCAL",NLOCAL);
	ofs_running << " " << setw(40) << "NLOCAL" << " = " << NLOCAL <<endl;
	//========================================================
	// (4) set index for iat2it, iat2ia, itia2iat, itiaiw2iwt
	//========================================================

	this->set_iat2it();
	
	delete[] iat2ia;
	this->iat2ia = new int[nat];// bug fix 2009-3-8

	// mohan add 2010-09-26
	assert(NLOCAL>0);
	delete[] iwt2iat;
	delete[] iwt2iw;
	this->iwt2iat = new int[NLOCAL];
	this->iwt2iw = new int[NLOCAL];

	this->itia2iat.create(ntype, namax);
	this->itiaiw2iwt.create(ntype, namax, nwmax*NPOL);
	int iat=0;
	int iwt=0;
	for(int it = 0;it < ntype;it++)
	{
		for(int ia=0; ia<atoms[it].na; ia++)
		{
			this->itia2iat(it, ia) = iat;
			this->iat2ia[iat] = ia;
			for(int iw=0; iw<atoms[it].nw * NPOL; iw++)
			{
				this->itiaiw2iwt(it, ia, iw) = iwt;
				this->iwt2iat[iwt] = iat;
				this->iwt2iw[iwt] = iw;
				++iwt;
			}
			++iat;
		}	
	}

	//========================
	// (5) set lmax and nmax
	//========================
	this->lmax = 0;	
	this->nmax = 0;
	for(int it=0; it<ntype; it++)
	{
		lmax = std::max( lmax, atoms[it].nwl);
		for(int l=0; l<atoms[it].nwl+1; l++)
		{
			nmax = std::max( nmax, atoms[it].l_nchi[l] );
		}

		int nchi = 0;
		for(int l=0; l<atoms[it].nwl+1; l++)
		{
			nchi += atoms[it].l_nchi[l];
		}
		this->nmax_total = std::max(nmax_total, nchi);
	}

	//=======================
	// (6) set lmax_ppwf
	//=======================
	this->lmax_ppwf = 0;
	for(int it=0; it< ntype; it++)
	{
		for(int ic=0; ic<ucell.atoms[it].nchi; ic++)
		{
			if( lmax_ppwf < ucell.atoms[it].lchi[ic] )
			{
				this->lmax_ppwf = ucell.atoms[it].lchi[ic]; 
			}
		}
	}

	/*
	for(int it=0; it< ntype; it++)
	{
		cout << " label=" << it << " nbeta=" << ucell.atoms[it].nbeta << endl;
		for(int ic=0; ic<ucell.atoms[it].nbeta; ic++)
		{
			cout << " " << ucell.atoms[it].lll[ic] << endl;
		}
	}
	*/

//	OUT(ofs_running,"lmax between L(pseudopotential)",lmax_ppwf);

	//=====================
	// Use localized basis
	//=====================
	if(BASIS_TYPE=="lcao" || BASIS_TYPE=="lcao_in_pw") //xiaohui add 2013-09-02
	{
		AUTO_SET("NBANDS",NBANDS);
	}
	else // plane wave basis
	{
		//if(winput::after_iter && winput::sph_proj)
		//{
		//	if(NBANDS < NLOCAL)
		//	{
		//		WARNING_QUIT("cal_nwfc","NBANDS must > NLOCAL !");
		//	}
		//}
	}

	return;
}

//======================
// Target : meshx
// Demand : atoms[].msh
//======================
void UnitCell_pseudo::cal_meshx()
{
	if(test_pseudo_cell) TITLE("UnitCell_pseudo","cal_meshx");
	this->meshx = 0;
	for (int it = 0;it < this->ntype;it++)
	{
		const int mesh = this->atoms[it].msh;
		if (mesh > this->meshx)
		{
			this->meshx = mesh;
		}
	}
	return;
}

//=========================
// Target : natomwfc
// Demand : atoms[].nchi
// 			atoms[].lchi
// 			atoms[].oc
// 			atoms[].na
//=========================
void UnitCell_pseudo::cal_natomwfc(void)
{
	if(test_pseudo_cell) TITLE("UnitCell_pseudo","cal_natomwfc");

	this->natomwfc = 0;
	for (int it = 0;it < ntype;it++)
	{
		//============================
		// Use pseudo-atomic orbitals
		//============================
		int tmp=0;
		for (int l = 0;l < atoms[it].nchi;l++)
		{
			if (atoms[it].oc[l] >= 0)
			{
				if(NSPIN==4)
				{
					if(atoms[it].has_so)
					{
						tmp += 2 * atoms[it].lchi[l];
						if(fabs(atoms[it].jchi[l] - atoms[it].lchi[l] - 0.5)<1e-6)
						tmp += 2 ;
					}
					else
					{
						tmp += 2 * (2 * atoms[it].lchi[l] + 1);
					}
				}
				else
					tmp += 2 * atoms[it].lchi[l] + 1;
			}
		}
		natomwfc += tmp * atoms[it].na;
	}
	OUT(ofs_running,"initial pseudo atomic orbital number",natomwfc);
	return;
}



//LiuXh add a new function here,
//20180515
void UnitCell_pseudo::setup_cell_after_vc(
        const string &s_pseudopot_dir,
        const string &fn, ofstream &log)
{
    if(MY_RANK == 0)
    {
        //ifstream ifa(fn.c_str(), ios::in);
        //this->read_atom_species_after_vc(ifa);
    }

    assert(lat0 > 0.0);
    this->omega = abs(latvec.Det()) * this->lat0 * lat0 * lat0;
    if(this->omega <= 0)
    {
        WARNING_QUIT("setup_cell_after_vc", "omega <= 0 .");
    }
    else
    {
        ofs_running << endl;
        OUT(ofs_running, "Volume (Bohr^3)", this->omega);
        OUT(ofs_running, "Volume (A^3))", this->omega * pow(BOHR_TO_A, 3));
    }

    //==========================================================
    // Calculate recip. lattice vectors and dot products
    // latvec has the unit of lat0, but G has the unit 2Pi/lat0
    //==========================================================
    this->GT = latvec.Inverse();
    this->G  = GT.Transpose();
    this->GGT = G * GT;
    this->invGGT = GGT.Inverse();

    for(int it=0; it<ucell.ntype; it++)
    {
        Atom* atom = &ucell.atoms[it];
        for(int ia =0;ia< atom->na;ia++)
        {
            atom->tau[ia] = atom->taud[ia] * ucell.latvec;
/*
#ifdef __MPI
Parallel_Common::bcast_double( atom->tau[ia].x );
Parallel_Common::bcast_double( atom->tau[ia].y );
Parallel_Common::bcast_double( atom->tau[ia].z );
Parallel_Common::bcast_double( atom->taud[ia].x );
Parallel_Common::bcast_double( atom->taud[ia].y );
Parallel_Common::bcast_double( atom->taud[ia].z );
#endif
*/
        }
    }
#ifdef __MPI
    MPI_Barrier(MPI_COMM_WORLD);
    for (int i=0;i<ucell.ntype;i++)
    {
        ucell.atoms[i].bcast_atom(); // bcast tau array
    }
#endif

    ofs_running << endl;
    out.printM3(ofs_running,"Lattice vectors: (Cartesian coordinate: in unit of a_0)",latvec);
    out.printM3(ofs_running,"Reciprocal vectors: (Cartesian coordinate: in unit of 2 pi/a_0)",G);

    return;
}
