#include "global.h"
#include "klist.h"
#include "../src_parallel/parallel_global.h"
#include "symmetry.h"

K_Vectors::K_Vectors()
{	
#ifdef _MCD_CHECK
    FILE* out;
    out=fopen("1_Memory", "w");
    if ( out == NULL )
    {
        cout << "\n Can't open file!";
        QUIT();
    }
    _MCD_RealTimeLog( out );
    _MCD_MemStatLog( out );
//	showMemStats();
#endif

    nspin = 0; // default spin.
    kc_done = false;
    kd_done = false;

    kvec_c = new Vector3<double>[1];
    kvec_d = new Vector3<double>[1];
    kvec_d_ibz = new Vector3<double>[1];

    wk = new double[1];
    wk_ibz = new double[1];
    ngk = new int[1];
    isk = new int[1];

    ibz2bz = new int[1];

    nks = 0;
    nkstot = 0;
    nkstot_ibz = 0;

    k_nkstot = 0; //LiuXh add 20180619
}

K_Vectors::~K_Vectors()
{
//	TITLE("K_Vectors","~K_Vectors");
    delete[] kvec_c;
    delete[] kvec_d;
    delete[] kvec_d_ibz;
    delete[] wk;
    delete[] wk_ibz;
    delete[] ngk;
    delete[] isk;
    delete[] ibz2bz;
#ifdef _MCD_CHECK
    showMemStats();
#endif
}

void K_Vectors::set(
    const Symmetry &symm,
    const string &k_file_name,
    const int& nspin_in,
    const Matrix3 &reciprocal_vec,
    const Matrix3 &latvec)
{
    TITLE("K_Vectors", "set");

	ofs_running << "\n\n\n\n";
	ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
	ofs_running << " |                                                                    |" << endl;
	ofs_running << " | Setup K-points                                                     |" << endl;
	ofs_running << " | We setup the k-points according to input parameters.               |" << endl;
	ofs_running << " | The reduced k-points are set according to symmetry operations.     |" << endl;
	ofs_running << " | We treat the spin as another set of k-points.                      |" << endl;
	ofs_running << " |                                                                    |" << endl;
	ofs_running << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
	ofs_running << "\n\n\n\n";

	ofs_running << "\n SETUP K-POINTS" << endl;

	// (1) set nspin, read kpoints.
	this->nspin = nspin_in;
	OUT(ofs_running,"nspin",nspin);
	if(this->nspin==4) 
	{
		this->nspin = 1;//zhengdy-soc
	}
		
	bool read_succesfully = this->read_kpoints(k_file_name);
#ifdef __MPI
	Parallel_Common::bcast_bool(read_succesfully);
#endif
	if(!read_succesfully)
	{
		WARNING_QUIT("K_Vectors::set","Something wrong while reading KPOINTS.");
	}

    // (2)
	if(Symmetry::symm_flag)
    {
        this->ibz_kpoint(symm);
        this->update_use_ibz();
        this->nks = this->nkstot = this->nkstot_ibz;
    }

    // (3)
    this->set_both_kvec(reciprocal_vec, latvec);

	int deg = 0;
	if(NSPIN == 1)
	{
		deg = 2;
	}
	else if(NSPIN == 2||NSPIN==4)
	{
		deg = 1;
	}
	else
	{
        WARNING_QUIT("K_Vectors::set", "Only available for nspin = 1 or 2 or 4");
    }
	this->normalize_wk(deg);

    // It's very important in parallel case,
    // firstly do the mpi_k() and then
    // do set_kup_and_kdw()
	Pkpoints.kinfo(nkstot);

    this->mpi_k();//2008-4-29

    this->set_kup_and_kdw();

    this->print_klists(ofs_running);

	//cout << " NUMBER OF K-POINTS   : " << nkstot << endl;

    return;
}

void K_Vectors::renew(const int &kpoint_number)
{
    delete[] kvec_c;
    delete[] kvec_d;
    delete[] wk;
    delete[] isk;
    delete[] ngk;

    kvec_c = new Vector3<double>[kpoint_number];
    kvec_d = new Vector3<double>[kpoint_number];
    wk = new double[kpoint_number];
    isk = new int[kpoint_number];
    ngk = new int[kpoint_number];

    Memory::record("K_Vectors","kvec_c",kpoint_number*3,"double");
    Memory::record("K_Vectors","kvec_d",kpoint_number*3,"double");
    Memory::record("K_Vectors","wk",kpoint_number*3,"double");
    Memory::record("K_Vectors","isk",kpoint_number*3,"int");
    Memory::record("K_Vectors","ngk",kpoint_number*3,"int");

    return;
}

bool K_Vectors::read_kpoints(const string &fn)
{
    TITLE("K_Vectors", "read_kpoints");
    if (MY_RANK != 0) return 1;

	// mohan add 2010-09-04
	if(GAMMA_ONLY_LOCAL)
	{
		ofs_warning << " Auto generating k-points file: " << fn << endl;
		ofstream ofs(fn.c_str());
		ofs << "K_POINTS" << endl;
		ofs << "0" << endl;
		ofs << "Gamma" << endl;
		ofs << "1 1 1 0 0 0" << endl;
		ofs.close();
	}

    ifstream ifk(fn.c_str());
    if (!ifk) 
	{
		ofs_warning << " Can't find File name : " << fn << endl;
		return 0;
    }

    ifk >> setiosflags(ios::uppercase);

    ifk.clear();
    ifk.seekg(0);

    string word;
    string kword;

    int ierr = 0;

    ifk.rdstate();

    while (ifk.good())
    {
        ifk >> word;
        ifk.ignore(150, '\n'); //LiuXh add 20180416, fix bug in k-point file when the first line with comments 
        if (word == "K_POINTS" || word == "KPOINTS" || word == "K" )
        {
            ierr = 1;
            break;
        }

        ifk.rdstate();
    }

    if (ierr == 0)
    {
		ofs_warning << " symbol K_POINTS not found." << endl;
		return 0;
    }

    //input k-points are in 2pi/a units
    READ_VALUE(ifk, nkstot);

    this->k_nkstot = nkstot; //LiuXh add 20180619

    //cout << " nkstot = " << nkstot << endl;
    READ_VALUE(ifk, kword);

    this->k_kword = kword; //LiuXh add 20180619

	// mohan update 2021-02-22
	int max_kpoints = 100000;
    if (nkstot > 100000)
    {
		ofs_warning << " nkstot > MAX_KPOINTS" << endl;
        return 0;
    }

    int k_type = 0;
    if (nkstot == 0) // nkstot==0, use monkhorst_pack. add by dwan
    {
        if (kword == "Gamma")
        {
            k_type = 0;
			OUT(ofs_running,"Input type of k points","Monkhorst-Pack(Gamma)");
        }
        else if (kword == "Monkhorst-Pack" || kword == "MP" || kword == "mp")
        {
            k_type = 1;
			OUT(ofs_running,"Input type of k points","Monkhorst-Pack");
        }
        else
        {
			ofs_warning << " Error: neither Gamma nor Monkhorst-Pack." << endl;
			return 0;
        }

        ifk >> nmp[0] >> nmp[1] >> nmp[2];

        ifk >> koffset[0] >> koffset[1] >> koffset[2];
        this->Monkhorst_Pack(nmp, koffset, k_type);
    }
    else if (nkstot > 0)
    {
        if (kword == "Cartesian" || kword == "C")
        {
        	this->renew(nkstot * nspin);//mohan fix bug 2009-09-01
            for (int i = 0;i < nkstot;i++)
            {
                ifk >> kvec_c[i].x >> kvec_c[i].y >> kvec_c[i].z;
				READ_VALUE(ifk, wk[i]);
            }

            this->kc_done = true;
        }
        else if (kword == "Direct" || kword == "D")
        {
        	this->renew(nkstot * nspin);//mohan fix bug 2009-09-01
            for (int i = 0;i < nkstot;i++)
            {
                ifk >> kvec_d[i].x >> kvec_d[i].y >> kvec_d[i].z;
				READ_VALUE(ifk, wk[i]);
            }
            this->kd_done = true;
        }
		else if (kword == "Line_Cartesian" )
		{
			//cout << " kword = " << kword << endl;
			if(Symmetry::symm_flag)
			{
				WARNING("K_Vectors::read_kpoints","Line mode of k-points is open, please set symmetry to 0.");
				return 0;
			}
		
		
			// how many special points.	
			int nks_special = this->nkstot;
			//cout << " nks_special = " << nks_special << endl;
			
			//------------------------------------------
			// number of points to the next k points
			//------------------------------------------
			int* nkl = new int[nks_special];
				
			//------------------------------------------
			// cartesian coordinates of special points.
			//------------------------------------------
			double *ksx = new double[nks_special];
			double *ksy = new double[nks_special];
			double *ksz = new double[nks_special];
			vector<double> kposx;
			vector<double> kposy;
			vector<double> kposz;
			ZEROS(nkl, nks_special);
			
			//recalculate nkstot.
			nkstot = 0;
			for(int iks=0; iks<nks_special; iks++)
			{
				ifk >> ksx[iks]; 
				ifk >> ksy[iks];
				ifk >> ksz[iks];
				READ_VALUE( ifk, nkl[iks] );
				//cout << " nkl[" << iks << "]=" << nkl[iks] << endl;
				assert(nkl[iks] >= 0);
				nkstot += nkl[iks];
			}
			assert( nkl[nks_special-1] == 1);

			//cout << " nkstot = " << nkstot << endl;
        	this->renew(nkstot * nspin);//mohan fix bug 2009-09-01

			int count = 0;
			for(int iks=1; iks<nks_special; iks++)
			{
				double dx = (ksx[iks] - ksx[iks-1]) / nkl[iks-1];
				double dy = (ksy[iks] - ksy[iks-1]) / nkl[iks-1];
				double dz = (ksz[iks] - ksz[iks-1]) / nkl[iks-1];
//				ofs_running << " dx=" << dx << " dy=" << dy << " dz=" << dz << endl;
				for(int is=0; is<nkl[iks-1]; is++)
				{
					kvec_c[count].x = ksx[iks-1] + is*dx;
					kvec_c[count].y = ksy[iks-1] + is*dy;
					kvec_c[count].z = ksz[iks-1] + is*dz;
					++count;
				}
			}

			// deal with the last special k point.
			kvec_c[count].x = ksx[nks_special-1];
			kvec_c[count].y = ksy[nks_special-1];
			kvec_c[count].z = ksz[nks_special-1];
			++count;
		
			//cout << " count = " << count << endl;
			assert (count == nkstot );

			for(int ik=0; ik<nkstot; ik++)
			{
				wk[ik] = 1.0;
			}
			
			ofs_warning << " Error : nkstot == -1, not implemented yet." << endl;

			delete[] nkl;
			delete[] ksx;
			delete[] ksy;
			delete[] ksz;

            this->kc_done = true;

		}

		else if (kword == "Line_Direct" || kword == "L" || kword == "Line" )
		{
			//cout << " kword = " << kword << endl;
			if(Symmetry::symm_flag)
			{
				WARNING("K_Vectors::read_kpoints","Line mode of k-points is open, please set symmetry to 0.");
				return 0;
			}
		
		
			// how many special points.	
			int nks_special = this->nkstot;
			//cout << " nks_special = " << nks_special << endl;
			
			//------------------------------------------
			// number of points to the next k points
			//------------------------------------------
			int* nkl = new int[nks_special];
				
			//------------------------------------------
			// cartesian coordinates of special points.
			//------------------------------------------
			double *ksx = new double[nks_special];
			double *ksy = new double[nks_special];
			double *ksz = new double[nks_special];
			vector<double> kposx;
			vector<double> kposy;
			vector<double> kposz;
			ZEROS(nkl, nks_special);
			
			//recalculate nkstot.
			nkstot = 0;
			for(int iks=0; iks<nks_special; iks++)
			{
				ifk >> ksx[iks]; 
				ifk >> ksy[iks];
				ifk >> ksz[iks];
				READ_VALUE( ifk, nkl[iks] );
				//cout << " nkl[" << iks << "]=" << nkl[iks] << endl;
				assert(nkl[iks] >= 0);
				nkstot += nkl[iks];
			}
			assert( nkl[nks_special-1] == 1);

			//cout << " nkstot = " << nkstot << endl;
        	this->renew(nkstot * nspin);//mohan fix bug 2009-09-01

			int count = 0;
			for(int iks=1; iks<nks_special; iks++)
			{
				double dx = (ksx[iks] - ksx[iks-1]) / nkl[iks-1];
				double dy = (ksy[iks] - ksy[iks-1]) / nkl[iks-1];
				double dz = (ksz[iks] - ksz[iks-1]) / nkl[iks-1];
//				ofs_running << " dx=" << dx << " dy=" << dy << " dz=" << dz << endl;
				for(int is=0; is<nkl[iks-1]; is++)
				{
					kvec_d[count].x = ksx[iks-1] + is*dx;
					kvec_d[count].y = ksy[iks-1] + is*dy;
					kvec_d[count].z = ksz[iks-1] + is*dz;
					++count;
				}
			}

			// deal with the last special k point.
			kvec_d[count].x = ksx[nks_special-1];
			kvec_d[count].y = ksy[nks_special-1];
			kvec_d[count].z = ksz[nks_special-1];
			++count;
		
			//cout << " count = " << count << endl;
			assert (count == nkstot );

			for(int ik=0; ik<nkstot; ik++)
			{
				wk[ik] = 1.0;
			}
			
			ofs_warning << " Error : nkstot == -1, not implemented yet." << endl;

			delete[] nkl;
			delete[] ksx;
			delete[] ksy;
			delete[] ksz;

            this->kd_done = true;

		}

        else
        {
			ofs_warning << " Error : neither Cartesian nor Direct kpoint." << endl;
			return 0;
        }
    }

    this->nks = this->nkstot;

	OUT(ofs_running,"nkstot",nkstot);
    return 1;
} // END SUBROUTINE


double K_Vectors::Monkhorst_Pack_formula( const int &k_type, const double &offset,
                                      const int& n, const int &dim)
{
    double coordinate;
    if (k_type==1) coordinate = (offset + 2.0 * (double)n - (double)dim - 1.0) / (2.0 * (double)dim);
    else coordinate = (offset + (double)n - 1.0) / (double)dim;
    return coordinate;
}

//add by dwan
void K_Vectors::Monkhorst_Pack(const int *nmp_in, const double *koffset_in, const int k_type)
{
    if (test_kpoint) TITLE("K_Vectors", "Monkhorst_Pack");
    const int mpnx = nmp_in[0];
    const int mpny = nmp_in[1];
    const int mpnz = nmp_in[2];

    this->nkstot = mpnx * mpny * mpnz;
    // only can renew after nkstot is estimated.
    this->renew(nkstot * nspin); // mohan fix bug 2009-09-01

    for (int x = 1;x <= mpnx;x++)
    {
        double v1 = Monkhorst_Pack_formula( k_type, koffset_in[0], x, mpnx);
		if( abs(v1) < 1.0e-10 ) v1 = 0.0; //mohan update 2012-06-10
        for (int y = 1;y <= mpny;y++)
        {
            double v2 = Monkhorst_Pack_formula( k_type, koffset_in[1], y, mpny);
		    if( abs(v2) < 1.0e-10 ) v2 = 0.0;
            for (int z = 1;z <= mpnz;z++)
            {
                double v3 = Monkhorst_Pack_formula( k_type, koffset_in[2], z, mpnz);
				if( abs(v3) < 1.0e-10 ) v3 = 0.0;
                // index of nks kpoint
                const int i = mpnx * mpny * (z - 1) + mpnx * (y - 1) + (x - 1);
                kvec_d[i].set(v1, v2, v3);
            }
        }
    }

    const double weight = 1.0 / static_cast<double>(nkstot);
    for (int ik=0; ik<nkstot; ik++)
    {
        wk[ik] = weight;
    }
    this->kd_done = true;

    return;
}

void K_Vectors::update_use_ibz( void )
{
    if (MY_RANK!=0) return;
    TITLE("K_Vectors","update_use_ibz");
    assert( nkstot_ibz > 0 );

	// update nkstot
    this->nkstot = this->nkstot_ibz;

	OUT(ofs_running,"nkstot now",nkstot);

    delete[] kvec_d;
    this->kvec_d = new Vector3<double>[ this->nkstot ];

    for (int i = 0; i < this->nkstot; ++i)
    {
        this->kvec_d[i] = this->kvec_d_ibz[i];

		// update weight.
        this->wk[i] = this->wk_ibz[i];
    }

    this->kd_done = true;
    this->kc_done = false;
    return;
}

void K_Vectors::ibz_kpoint(const Symmetry &symm)
{
    if (MY_RANK!=0) return;
    TITLE("K_Vectors", "ibz_kpoint");
    //===============================================
    // search in all space group operations
    // if the operations does not already included
    // inverse operation, double it.
    //===============================================
    bool include_inv = false;
    Matrix3 *kgmatrix = new Matrix3[48 * 2];
    Matrix3 inv(-1, 0, 0, 0, -1, 0, 0, 0, -1);

    int nrotkm = symm.nrotk;// change if inv not included
    for (int i = 0; i < nrotkm; ++i)
    {
        if (symm.gmatrix[i] == inv)
        {
            include_inv = true;
        }
        kgmatrix[i] = symm.gmatrix[i];
    }

    if (!include_inv)
    {
        for (int i = 0; i<symm.nrotk; ++i)
        {
            kgmatrix[i + symm.nrotk] = inv * symm.gmatrix[i];
        }
        nrotkm = 2 * symm.nrotk;
    }

    // use operation : kgmatrix to find
    // the new set kvec_d : ir_kpt
    this->nkstot_ibz = 0;

    assert(nkstot > 0 );
    delete[] kvec_d_ibz;
    delete[] wk_ibz;
    delete[] ibz2bz;
    this->kvec_d_ibz = new Vector3<double>[ this->nkstot ];
    this->wk_ibz = new double[ this->nkstot ];
    this->ibz2bz = new int[ this->nkstot ];

	// nkstot is the total input k-points number.
    const double weight = 1.0 / static_cast<double>(nkstot);

    Vector3<double> kvec_rot;
    
	
//	for(int i=0; i<nrotkm; i++)
//	{
//		out.printM3("rot matrix",kgmatrix[i]);
//	}
	
	
	// search in all k-poins.
    for (int i = 0; i < nkstot; ++i)
    {
		//cout << "\n kpoint = " << i << endl;
		//cout << "\n kvec_d = " << kvec_d[i].x << " " << kvec_d[i].y << " " << kvec_d[i].z;
        bool already_exist = false;
		int exist_number = -1;
		
        for (int j = 0; j < nrotkm; ++j)
        {
            if (!already_exist)
            {
                // rotate the kvec_d within all operations.
                // here use direct coordinates.
//                kvec_rot = kgmatrix[j] * kvec_d[i];
				// mohan modify 2010-01-30.
				// mohan modify again 2010-01-31
				// fix the bug like kvec_d * G; is wrong
				//kvec_rot = kvec_d[i] * kgmatrix[j]; //wrong for total energy, but correct for nonlocal force.
				kvec_rot = kgmatrix[j] * kvec_d[i]; //correct for total energy, but wrong for nonlocal force.
				

                kvec_rot.x = fmod(kvec_rot.x + 100, 1);
                kvec_rot.y = fmod(kvec_rot.y + 100, 1);
                kvec_rot.z = fmod(kvec_rot.z + 100, 1);

//				cout << "\n kvec_rot = " << kvec_rot.x << " " << kvec_rot.y << " " << kvec_rot.z;

                for (int k=0; k< this->nkstot_ibz; ++k)
                {
                    if (    symm.equal(kvec_rot.x, this->kvec_d_ibz[k].x) &&
                            symm.equal(kvec_rot.y, this->kvec_d_ibz[k].y) &&
                            symm.equal(kvec_rot.z, this->kvec_d_ibz[k].z))
                    {
                        already_exist = true;

						// find another ibz k point,
						// but is already in the ibz_kpoint list.
						// so the weight need to +1; 
                        this->wk_ibz[k] += weight;
						exist_number = k;
                        break;
                    }
                }
            }//end !already_exist
        }

        if (!already_exist)
        {
			//if it's a new ibz kpoint.
			//nkstot_ibz indicate the index of ibz kpoint.
            this->kvec_d_ibz[nkstot_ibz] = kvec_rot;

			//the weight should be averged k-point weight.
            this->wk_ibz[nkstot_ibz] = weight;

			//ibz2bz records the index of origin k points.
            this->ibz2bz[nkstot_ibz] = i;
            ++nkstot_ibz;
        }
		else //mohan fix bug 2010-1-30
		{
//			cout << "\n\n already exist ! ";

//			cout << "\n kvec_rot = " << kvec_rot.x << " " << kvec_rot.y << " " << kvec_rot.z;
//			cout << "\n kvec_d_ibz = " << kvec_d_ibz[exist_number].x 
//			<< " " << kvec_d_ibz[exist_number].y 
//			<< " " << kvec_d_ibz[exist_number].z;
			
			double kmol_new = kvec_d[i].norm2();
			double kmol_old = kvec_d_ibz[exist_number].norm2();

//			cout << "\n kmol_new = " << kmol_new;
//			cout << "\n kmol_old = " << kmol_old;
			

			// why we need this step?
			// because in pw_basis.cpp, while calculate ggwfc2,
			// if we want to keep the result of symmetry operation is right.
			// we need to fix the number of plane wave.
			// and the number of plane wave is depending on the |K+G|,
			// so we need to |K|max to be the same as 'no symmetry'.
			// mohan 2010-01-30
			if(kmol_new > kmol_old)
			{
				kvec_d_ibz[exist_number] = kvec_d[i];
			}
		}
//		BLOCK_HERE("check k point");
    }
	OUT(ofs_running,"nkstot_ibz",nkstot_ibz);

	ofs_running << " " << setw(8) << "IBZ" << setw(20) << "DirectX"
	<< setw(20) << "DirectY" << setw(20) << "DirectZ" 
	<< setw(20) << "Weight" << setw(10) << "ibz2bz" << endl;
    for (int ik=0; ik<nkstot_ibz; ik++)
    {
        ofs_running << " "
			<< setw(8) << ik+1
            << setw(20) << this->kvec_d_ibz[ik].x
            << setw(20) << this->kvec_d_ibz[ik].y
            << setw(20) << this->kvec_d_ibz[ik].z
            << setw(20) << this->wk_ibz[ik]
            << setw(10) << this->ibz2bz[ik] << endl;
    }
    delete[] kgmatrix;

    return;
}


void K_Vectors::set_both_kvec(const Matrix3 &G, const Matrix3 &R)
{

    if(FINAL_SCF) //LiuXh add 20180606
    {
        if(k_nkstot == 0)
        {
            kd_done = true;
            kc_done = false;
        }
        else
        {
            if(k_kword == "Cartesian" || k_kword == "C")
	    {
	        kc_done = true;
	        kd_done = false;
	    }
	    else if (k_kword == "Direct" || k_kword == "D")
	    {
	        kd_done = true;
	        kc_done = false;
	    }
            else
            {
                ofs_warning << " Error : neither Cartesian nor Direct kpoint." << endl;
            }
        }
    }

    // set cartesian k vectors.
    if (!kc_done && kd_done)
    {
        for (int i = 0;i < nkstot;i++)
        {
//wrong!!   kvec_c[i] = G * kvec_d[i];
// mohan fixed bug 2010-1-10
			if( abs(kvec_d[i].x) < 1.0e-10 ) kvec_d[i].x = 0.0;
			if( abs(kvec_d[i].y) < 1.0e-10 ) kvec_d[i].y = 0.0;
			if( abs(kvec_d[i].z) < 1.0e-10 ) kvec_d[i].z = 0.0;

			kvec_c[i] = kvec_d[i] * G;

			// mohan add2012-06-10
			if( abs(kvec_c[i].x) < 1.0e-10 ) kvec_c[i].x = 0.0;
			if( abs(kvec_c[i].y) < 1.0e-10 ) kvec_c[i].y = 0.0;
			if( abs(kvec_c[i].z) < 1.0e-10 ) kvec_c[i].z = 0.0;
        }
        kc_done = true;
    }

    // set direct k vectors
    else if (kc_done && !kd_done)
    {
        Matrix3 RT = R.Transpose();
        for (int i = 0;i < nkstot;i++)
        {
//			cout << " ik=" << i
//				<< " kvec.x=" << kvec_c[i].x
//				<< " kvec.y=" << kvec_c[i].y
//				<< " kvec.z=" << kvec_c[i].z << endl;
//wrong!            kvec_d[i] = RT * kvec_c[i];
// mohan fixed bug 2011-03-07
            kvec_d[i] = kvec_c[i] * RT;
        }
        kd_done = true;
    }

	ofs_running << "\n " << setw(8) << "KPOINTS" 
	<< setw(20) << "DIRECT_X"
	<< setw(20) << "DIRECT_Y"
	<< setw(20) << "DIRECT_Z"
	<< setw(20) << "WEIGHT" << endl;

	for(int i=0; i<nkstot; i++)
	{
        ofs_running << " "
			<< setw(8) << i+1
             << setw(20) << this->kvec_d[i].x
             << setw(20) << this->kvec_d[i].y
             << setw(20) << this->kvec_d[i].z
             << setw(20) << this->wk[i] << endl;
	}

    return;
}


void K_Vectors::normalize_wk(const int &degspin)
{
	if(MY_RANK!=0) return;
    double sum = 0.0;

    for (int ik = 0;ik < nkstot;ik++)
    {
        sum += this->wk[ik];
    }
	assert(sum>0.0);

    for (int ik = 0;ik < nkstot;ik++)
    {
        this->wk[ik] /= sum;
    }

    for (int ik = 0;ik < nkstot;ik++)
    {
        this->wk[ik] *= degspin;
    }

    return;
}

void K_Vectors::mpi_k(void)
{
#ifdef __MPI
	TITLE("K_Vectors","mpi_k");

    Parallel_Common::bcast_bool(kc_done);

    Parallel_Common::bcast_bool(kd_done);

    Parallel_Common::bcast_int(nspin);

    Parallel_Common::bcast_int(nkstot);

    Parallel_Common::bcast_int(nmp, 3);

    Parallel_Common::bcast_double(koffset, 3);

    this->nks = Pkpoints.nks_pool[MY_POOL];

	ofs_running << endl;
	OUT(ofs_running,"k-point number in this process",nks);
    int nks_minimum = this->nks;

	Parallel_Reduce::gather_min_int_all( nks_minimum );

    if (nks_minimum == 0)
    {
        WARNING_QUIT("K_Vectors::mpi_k()"," nks == 0, some processor have no k point!");
    }
    else
    {
		OUT(ofs_running,"minimum distributed K point number",nks_minimum);
    }

    int *isk_aux = new int[nkstot];

    double *wk_aux = new double[nkstot];

    double *kvec_c_aux = new double[nkstot*3];

    double *kvec_d_aux = new double[nkstot*3];

    if (MY_RANK == 0)
    {
        for (int ik = 0;ik < nkstot;ik++)
        {
            isk_aux[ik] = isk[ik];
            wk_aux[ik] = wk[ik];
            kvec_c_aux[3*ik]   = kvec_c[ik].x;
            kvec_c_aux[3*ik+1] = kvec_c[ik].y;
            kvec_c_aux[3*ik+2] = kvec_c[ik].z;
            kvec_d_aux[3*ik]   = kvec_d[ik].x;
            kvec_d_aux[3*ik+1] = kvec_d[ik].y;
            kvec_d_aux[3*ik+2] = kvec_d[ik].z;
        }
    }

    Parallel_Common::bcast_int(isk_aux, nkstot);

    Parallel_Common::bcast_double(wk_aux, nkstot);
    Parallel_Common::bcast_double(kvec_c_aux , nkstot*3);
    Parallel_Common::bcast_double(kvec_d_aux , nkstot*3);

    this->renew(this->nks * this->nspin);

    // distribute
    int k_index = 0;

    for (int i = 0;i < nks;i++)
    {
        // 3 is because each k point has three value:kx, ky, kz
        k_index = i + Pkpoints.startk_pool[MY_POOL] ;
        kvec_c[i].x = kvec_c_aux[k_index*3];
        kvec_c[i].y = kvec_c_aux[k_index*3+1];
        kvec_c[i].z = kvec_c_aux[k_index*3+2];
        kvec_d[i].x = kvec_d_aux[k_index*3];
        kvec_d[i].y = kvec_d_aux[k_index*3+1];
        kvec_d[i].z = kvec_d_aux[k_index*3+2];
        wk[i] = wk_aux[k_index];
        isk[i] = isk_aux[k_index];
    }

    delete[] isk_aux;
    delete[] wk_aux;
    delete[] kvec_c_aux;
    delete[] kvec_d_aux;
#endif
} // END SUBROUTINE



//----------------------------------------------------------
// This routine sets the k vectors for the up and down spin
//----------------------------------------------------------
// from set_kup_and_kdw.f90
void K_Vectors::set_kup_and_kdw(void)
{
    TITLE("K_Vectors", "setup_kup_and_kdw");

    //=========================================================================
    // on output: the number of points is doubled and xk and wk in the
    // first (nks/2) positions correspond to up spin
    // those in the second (nks/2) ones correspond to down spin
    //=========================================================================
    switch (nspin)
    {
    case 1:

        for (int ik = 0; ik < nks; ik++)
        {
            this->isk[ik] = 0;
        }

        break;

    case 2:

        for (int ik = 0; ik < nks; ik++)
        {
            this->kvec_c[ik+nks] = kvec_c[ik];
            this->kvec_d[ik+nks] = kvec_d[ik];
            this->wk[ik+nks]     = wk[ik];
            this->isk[ik]        = 0;
            this->isk[ik+nks]    = 1;
        }

        this->nks *= 2;
        this->nkstot *= 2;

		OUT(ofs_running,"nks(nspin=2)",nks);
		OUT(ofs_running,"nkstot(nspin=2)",nkstot);
        break;
    case 4:

        for (int ik = 0; ik < nks; ik++)
        {
            this->isk[ik] = 0;
        }

        break;
    }

    return;
} // end subroutine set_kup_and_kdw


void K_Vectors::print_klists(ofstream &ofs)
{
    TITLE("K_Vectors", "print_klists");

    if (nkstot < nks)
    {
        cout << "\n nkstot=" << nkstot;
        cout << "\n nks=" << nks;
        WARNING_QUIT("print_klists","nkstot < nks");
    }

	ofs_running << "\n " << setw(8) << "KPOINTS" 
	<< setw(20) << "CARTESIAN_X"
	<< setw(20) << "CARTESIAN_Y"
	<< setw(20) << "CARTESIAN_Z"
	<< setw(20) << "WEIGHT" << endl;
	for(int i=0; i<nks; i++)
	{
        ofs_running << " "
			<< setw(8) << i+1
             << setw(20) << this->kvec_c[i].x
             << setw(20) << this->kvec_c[i].y
             << setw(20) << this->kvec_c[i].z
             << setw(20) << this->wk[i] << endl;
	}

	ofs_running << "\n " << setw(8) << "KPOINTS" 
	<< setw(20) << "DIRECT_X"
	<< setw(20) << "DIRECT_Y"
	<< setw(20) << "DIRECT_Z"
	<< setw(20) << "WEIGHT" << endl;
	for(int i=0; i<nks; i++)
	{
        ofs_running << " "
			<< setw(8) << i+1
             << setw(20) << this->kvec_d[i].x
             << setw(20) << this->kvec_d[i].y
             << setw(20) << this->kvec_d[i].z
             << setw(20) << this->wk[i] << endl;
	}

    return;
}

//LiuXh add a new function here,
//20180515
void K_Vectors::set_after_vc(
        const Symmetry &symm,
        const string &k_file_name,
        const int& nspin_in,
        const Matrix3 &reciprocal_vec,
        const Matrix3 &latvec)
{
    TITLE("K_Vectors", "set_after_vc");

    ofs_running << "\n SETUP K-POINTS" << endl;
    this->nspin = nspin_in;
    OUT(ofs_running,"nspin",nspin);

    this->set_both_kvec_after_vc(reciprocal_vec, latvec);
    //this->set_both_kvec(reciprocal_vec, latvec);

    this->mpi_k_after_vc();

    this->set_kup_and_kdw_after_vc();

    this->print_klists(ofs_running);

    return;
}

//LiuXh add a new function here,
//20180515
void K_Vectors::mpi_k_after_vc(void)
{
#ifdef __MPI
    TITLE("K_Vectors","mpi_k_after_vc");

    Parallel_Common::bcast_bool(kc_done);
    Parallel_Common::bcast_bool(kd_done);
    Parallel_Common::bcast_int(nspin);
    Parallel_Common::bcast_int(nkstot);
    Parallel_Common::bcast_int(nmp, 3);
    Parallel_Common::bcast_double(koffset, 3);

    this->nks = Pkpoints.nks_pool[MY_POOL];
    ofs_running << endl;
    OUT(ofs_running,"k-point number in this process",nks);
    int nks_minimum = this->nks;

    Parallel_Reduce::gather_min_int_all( nks_minimum );

    if (nks_minimum == 0)
    {
        WARNING_QUIT("K_Vectors::mpi_k()"," nks == 0, some processor have no k point!");
    }
    else
    {
        OUT(ofs_running,"minimum distributed K point number",nks_minimum);
    }

    int *isk_aux = new int[nkstot];
    double *wk_aux = new double[nkstot];
    double *kvec_c_aux = new double[nkstot*3];
    double *kvec_d_aux = new double[nkstot*3];

    if (MY_RANK == 0)
    {
        for (int ik = 0;ik < nkstot;ik++)
        {
            isk_aux[ik] = isk[ik];
            wk_aux[ik] = wk[ik];
            kvec_c_aux[3*ik]   = kvec_c[ik].x;
            kvec_c_aux[3*ik+1] = kvec_c[ik].y;
            kvec_c_aux[3*ik+2] = kvec_c[ik].z;
            kvec_d_aux[3*ik]   = kvec_d[ik].x;
            kvec_d_aux[3*ik+1] = kvec_d[ik].y;
            kvec_d_aux[3*ik+2] = kvec_d[ik].z;
        }
    }

    Parallel_Common::bcast_int(isk_aux, nkstot);
    Parallel_Common::bcast_double(wk_aux, nkstot);
    Parallel_Common::bcast_double(kvec_c_aux , nkstot*3);
    Parallel_Common::bcast_double(kvec_d_aux , nkstot*3);

    int k_index = 0;
    for (int i = 0;i < nks;i++)
    {
        k_index = i + Pkpoints.startk_pool[MY_POOL] ;
        kvec_c[i].x = kvec_c_aux[k_index*3];
        kvec_c[i].y = kvec_c_aux[k_index*3+1];
        kvec_c[i].z = kvec_c_aux[k_index*3+2];
        kvec_d[i].x = kvec_d_aux[k_index*3];
        kvec_d[i].y = kvec_d_aux[k_index*3+1];
        kvec_d[i].z = kvec_d_aux[k_index*3+2];
        wk[i] = wk_aux[k_index];
        isk[i] = isk_aux[k_index];
    }

    delete[] isk_aux;
    delete[] wk_aux;
    delete[] kvec_c_aux;
    delete[] kvec_d_aux;
#endif
}

void K_Vectors::set_both_kvec_after_vc(const Matrix3 &G, const Matrix3 &R)
{
    // set cartesian k vectors.
    kd_done = true;
    kc_done = false;
    if (!kc_done && kd_done)
    {
        for (int i = 0;i < nkstot;i++)
        {
//wrong!!   kvec_c[i] = G * kvec_d[i];
// mohan fixed bug 2010-1-10
			if( abs(kvec_d[i].x) < 1.0e-10 ) kvec_d[i].x = 0.0;
			if( abs(kvec_d[i].y) < 1.0e-10 ) kvec_d[i].y = 0.0;
			if( abs(kvec_d[i].z) < 1.0e-10 ) kvec_d[i].z = 0.0;

			kvec_c[i] = kvec_d[i] * G;

			// mohan add2012-06-10
			if( abs(kvec_c[i].x) < 1.0e-10 ) kvec_c[i].x = 0.0;
			if( abs(kvec_c[i].y) < 1.0e-10 ) kvec_c[i].y = 0.0;
			if( abs(kvec_c[i].z) < 1.0e-10 ) kvec_c[i].z = 0.0;
        }
        kc_done = true;
    }

    // set direct k vectors
    else if (kc_done && !kd_done)
    {
        Matrix3 RT = R.Transpose();
        for (int i = 0;i < nkstot;i++)
        {
//			cout << " ik=" << i
//				<< " kvec.x=" << kvec_c[i].x
//				<< " kvec.y=" << kvec_c[i].y
//				<< " kvec.z=" << kvec_c[i].z << endl;
//wrong!            kvec_d[i] = RT * kvec_c[i];
// mohan fixed bug 2011-03-07
            kvec_d[i] = kvec_c[i] * RT;
        }
        kd_done = true;
    }

	ofs_running << "\n " << setw(8) << "KPOINTS" 
	<< setw(20) << "DIRECT_X"
	<< setw(20) << "DIRECT_Y"
	<< setw(20) << "DIRECT_Z"
	<< setw(20) << "WEIGHT" << endl;

	for(int i=0; i<nkstot; i++)
	{
        ofs_running << " "
			<< setw(8) << i+1
             << setw(20) << this->kvec_d[i].x
             << setw(20) << this->kvec_d[i].y
             << setw(20) << this->kvec_d[i].z
             << setw(20) << this->wk[i] << endl;
	}

    return;
}

void K_Vectors::set_kup_and_kdw_after_vc(void)
{
    TITLE("K_Vectors", "setup_kup_and_kdw_after_vc");

    //=========================================================================
    // on output: the number of points is doubled and xk and wk in the
    // first (nks/2) positions correspond to up spin
    // those in the second (nks/2) ones correspond to down spin
    //=========================================================================
    switch (nspin)
    {
    case 1:

        for (int ik = 0; ik < nks; ik++)
        {
            this->isk[ik] = 0;
        }

        break;

    case 2:

        for (int ik = 0; ik < nks; ik++)
        {
            this->kvec_c[ik+nks] = kvec_c[ik];
            this->kvec_d[ik+nks] = kvec_d[ik];
            this->wk[ik+nks]     = wk[ik];
            this->isk[ik]        = 0;
            this->isk[ik+nks]    = 1;
        }

        this->nks *= 2;
        //this->nkstot *= 2;

        OUT(ofs_running,"nks(nspin=2)",nks);
        OUT(ofs_running,"nkstot(nspin=2)",nkstot);

        break;

    case 4:

        for (int ik = 0; ik < nks; ik++)
        {
            this->isk[ik] = 0;
        }

        break;

    }

    return;
} // end subroutine set_kup_and_kdw
