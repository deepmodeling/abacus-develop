#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "klist.h"
#include "module_base/parallel_global.h"
#include "module_cell/module_symmetry/symmetry.h"
#include "module_base/parallel_reduce.h"
#include "module_base/parallel_common.h"
#include "module_base/memory.h"
#include "module_io/berryphase.h"
#include "module_base/formatter.h"
#ifdef USE_PAW
#include "module_cell/module_paw/paw_cell.h"
#endif

K_Vectors::K_Vectors()
{
#ifdef _MCD_CHECK
    FILE* out;
    out=fopen("1_Memory", "w");
    if ( out == NULL )
    {
        std::cout << "\n Can't open file!";
        ModuleBase::QUIT();
    }
    _MCD_RealTimeLog( out );
    _MCD_MemStatLog( out );
//	showMemStats();
#endif

    nspin = 0; // default spin.
    kc_done = false;
    kd_done = false;

    nks = 0;
    nkstot = 0;
    nkstot_ibz = 0;

    k_nkstot = 0; //LiuXh add 20180619
}

K_Vectors::~K_Vectors()
{
//	ModuleBase::TITLE("K_Vectors","~K_Vectors");
#ifdef _MCD_CHECK
    showMemStats();
#endif
}

void K_Vectors::set(
    const ModuleSymmetry::Symmetry &symm,
    const std::string &k_file_name,
    const int& nspin_in,
    const ModuleBase::Matrix3 &reciprocal_vec,
    const ModuleBase::Matrix3 &latvec)
{
    ModuleBase::TITLE("K_Vectors", "set");

	GlobalV::ofs_running << "\n\n\n\n";
	GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
	GlobalV::ofs_running << " |                                                                    |" << std::endl;
	GlobalV::ofs_running << " | Setup K-points                                                     |" << std::endl;
	GlobalV::ofs_running << " | We setup the k-points according to input parameters.               |" << std::endl;
	GlobalV::ofs_running << " | The reduced k-points are set according to symmetry operations.     |" << std::endl;
	GlobalV::ofs_running << " | We treat the spin as another set of k-points.                      |" << std::endl;
	GlobalV::ofs_running << " |                                                                    |" << std::endl;
	GlobalV::ofs_running << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
	GlobalV::ofs_running << "\n\n\n\n";

	GlobalV::ofs_running << "\n SETUP K-POINTS" << std::endl;

	// (1) set nspin, read kpoints.
	this->nspin = nspin_in;
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nspin",nspin);
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
		ModuleBase::WARNING_QUIT("K_Vectors::set","Something wrong while reading KPOINTS.");
	}

    // output kpoints file
    std::string skpt1="";
    std::string skpt2="";

    // (2)
    //only berry phase need all kpoints including time-reversal symmetry!
    //if symm_flag is not set, only time-reversal symmetry would be considered.
    if(!berryphase::berry_phase_flag && ModuleSymmetry::Symmetry::symm_flag != -1)
    {
        bool match = true;
        this->ibz_kpoint(symm, ModuleSymmetry::Symmetry::symm_flag, skpt1, GlobalC::ucell, match);
#ifdef __MPI
	    Parallel_Common::bcast_bool(match);
#endif
        if (!match)
        {
            std::cout<< "Optimized lattice type of reciprocal lattice cannot match the optimized real lattice. " <<std::endl;
            std::cout << "It is often because the inaccuracy of lattice parameters in STRU." << std::endl;
            if (ModuleSymmetry::Symmetry::symm_autoclose)
            {
                ModuleBase::WARNING("K_Vectors::ibz_kpoint", "Automatically set symmetry to 0 and continue ...");
                std::cout << "Automatically set symmetry to 0 and continue ..." << std::endl;
                ModuleSymmetry::Symmetry::symm_flag = 0;
                match = true;
                this->ibz_kpoint(symm, ModuleSymmetry::Symmetry::symm_flag, skpt1, GlobalC::ucell, match);
            }
            else
                ModuleBase::WARNING_QUIT("K_Vectors::ibz_kpoint", "Possible solutions: \n \
1. Refine the lattice parameters in STRU;\n \
2. Use a different`symmetry_prec`.  \n \
3. Close symemtry: set `symmetry` to 0 in INPUT. \n \
4. Set `symmetry_autoclose` to 1 in INPUT to automatically close symmetry when this error occurs.");
        }
        if (ModuleSymmetry::Symmetry::symm_flag || is_mp)
        {
            this->update_use_ibz();
            this->nks = this->nkstot = this->nkstot_ibz;
        }
    }

    // (3)
    this->set_both_kvec(reciprocal_vec, latvec, skpt2);

	if(GlobalV::MY_RANK==0)
	{
        // output kpoints file
        std::stringstream skpt;
        skpt << GlobalV::global_readin_dir << "kpoints" ;
	    std::ofstream ofkpt( skpt.str().c_str()); // clear kpoints
        ofkpt << skpt2 <<skpt1;
        ofkpt.close();
    }

	int deg = 0;
	if(GlobalV::NSPIN == 1)
	{
		deg = 2;
	}
	else if(GlobalV::NSPIN == 2||GlobalV::NSPIN==4)
	{
		deg = 1;
	}
	else
	{
        ModuleBase::WARNING_QUIT("K_Vectors::set", "Only available for nspin = 1 or 2 or 4");
    }
	this->normalize_wk(deg);

    // It's very important in parallel case,
    // firstly do the mpi_k() and then
    // do set_kup_and_kdw()
	GlobalC::Pkpoints.kinfo(nkstot);
#ifdef __MPI
    this->mpi_k();//2008-4-29
#endif

    this->set_kup_and_kdw();

    this->print_klists(GlobalV::ofs_running);

	//std::cout << " NUMBER OF K-POINTS   : " << nkstot << std::endl;

#ifdef USE_PAW
    GlobalC::paw_cell.set_isk(nks,isk.data());
#endif

    return;
}

void K_Vectors::renew(const int &kpoint_number)
{
    kvec_c.resize(kpoint_number);
    kvec_d.resize(kpoint_number);
    wk.resize(kpoint_number);
    isk.resize(kpoint_number);
    ngk.resize(kpoint_number);

    /*ModuleBase::Memory::record("KV::kvec_c",sizeof(double) * kpoint_number*3);
    ModuleBase::Memory::record("KV::kvec_d",sizeof(double) * kpoint_number*3);
    ModuleBase::Memory::record("KV::wk",sizeof(double) * kpoint_number*3);
    ModuleBase::Memory::record("KV::isk",sizeof(int) * kpoint_number*3);
    ModuleBase::Memory::record("KV::ngk",sizeof(int) * kpoint_number*3);*/

    return;
}

bool K_Vectors::read_kpoints(const std::string &fn)
{
    ModuleBase::TITLE("K_Vectors", "read_kpoints");
    if (GlobalV::MY_RANK != 0) return 1;

	// mohan add 2010-09-04
	if(GlobalV::GAMMA_ONLY_LOCAL)
	{
		GlobalV::ofs_warning << " Auto generating k-points file: " << fn << std::endl;
		std::ofstream ofs(fn.c_str());
		ofs << "K_POINTS" << std::endl;
		ofs << "0" << std::endl;
		ofs << "Gamma" << std::endl;
		ofs << "1 1 1 0 0 0" << std::endl;
		ofs.close();
	}
    else if (GlobalV::KSPACING[0] > 0.0)
    {
        if (GlobalV::KSPACING[1] <= 0 || GlobalV::KSPACING[2] <= 0){
            ModuleBase::WARNING_QUIT("K_Vectors","kspacing shold > 0");
        };
        //number of K points = max(1,int(|bi|/KSPACING+1))
        ModuleBase::Matrix3 btmp = GlobalC::ucell.G;
        double b1 = sqrt(btmp.e11 * btmp.e11 + btmp.e12 * btmp.e12 + btmp.e13 * btmp.e13);
        double b2 = sqrt(btmp.e21 * btmp.e21 + btmp.e22 * btmp.e22 + btmp.e23 * btmp.e23);
        double b3 = sqrt(btmp.e31 * btmp.e31 + btmp.e32 * btmp.e32 + btmp.e33 * btmp.e33);
        int nk1 = std::max(1,static_cast<int>(b1 * ModuleBase::TWO_PI / GlobalV::KSPACING[0] / GlobalC::ucell.lat0 + 1));
        int nk2 = std::max(1,static_cast<int>(b2 * ModuleBase::TWO_PI / GlobalV::KSPACING[1] / GlobalC::ucell.lat0 + 1));
        int nk3 = std::max(1,static_cast<int>(b3 * ModuleBase::TWO_PI / GlobalV::KSPACING[2] / GlobalC::ucell.lat0 + 1));

        GlobalV::ofs_warning << " Generate k-points file according to KSPACING: " << fn << std::endl;
		std::ofstream ofs(fn.c_str());
		ofs << "K_POINTS" << std::endl;
		ofs << "0" << std::endl;
		ofs << "Gamma" << std::endl;
		ofs << nk1 << " " << nk2 << " " << nk3 <<" 0 0 0" << std::endl;
		ofs.close();
    }

    std::ifstream ifk(fn.c_str());
    if (!ifk)
	{
		GlobalV::ofs_warning << " Can't find File name : " << fn << std::endl;
		return 0;
    }

    ifk >> std::setiosflags(std::ios::uppercase);

    ifk.clear();
    ifk.seekg(0);

    std::string word;
    std::string kword;

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
		GlobalV::ofs_warning << " symbol K_POINTS not found." << std::endl;
		return 0;
    }

    //input k-points are in 2pi/a units
    ModuleBase::GlobalFunc::READ_VALUE(ifk, nkstot);

    this->k_nkstot = nkstot; //LiuXh add 20180619

    //std::cout << " nkstot = " << nkstot << std::endl;
    ModuleBase::GlobalFunc::READ_VALUE(ifk, kword);

    this->k_kword = kword; //LiuXh add 20180619

	// mohan update 2021-02-22
	int max_kpoints = 100000;
    if (nkstot > 100000)
    {
		GlobalV::ofs_warning << " nkstot > MAX_KPOINTS" << std::endl;
        return 0;
    }

    int k_type = 0;
    if (nkstot == 0) // nkstot==0, use monkhorst_pack. add by dwan
    {
        if (kword == "Gamma")
        {
            is_mp = true;
            k_type = 0;
			ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"Input type of k points","Monkhorst-Pack(Gamma)");
        }
        else if (kword == "Monkhorst-Pack" || kword == "MP" || kword == "mp")
        {
            is_mp = true;
            k_type = 1;
			ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"Input type of k points","Monkhorst-Pack");
        }
        else
        {
			GlobalV::ofs_warning << " Error: neither Gamma nor Monkhorst-Pack." << std::endl;
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
				ModuleBase::GlobalFunc::READ_VALUE(ifk, wk[i]);
            }

            this->kc_done = true;
        }
        else if (kword == "Direct" || kword == "D")
        {
        	this->renew(nkstot * nspin);//mohan fix bug 2009-09-01
            for (int i = 0;i < nkstot;i++)
            {
                ifk >> kvec_d[i].x >> kvec_d[i].y >> kvec_d[i].z;
				ModuleBase::GlobalFunc::READ_VALUE(ifk, wk[i]);
            }
            this->kd_done = true;
        }
		else if (kword == "Line_Cartesian" )
		{
			//std::cout << " kword = " << kword << std::endl;
			if(ModuleSymmetry::Symmetry::symm_flag == 1)
			{
				ModuleBase::WARNING("K_Vectors::read_kpoints","Line mode of k-points is open, please set symmetry to 0 or -1.");
				return 0;
			}


			// how many special points.
			int nks_special = this->nkstot;
			//std::cout << " nks_special = " << nks_special << std::endl;

			//------------------------------------------
			// number of points to the next k points
			//------------------------------------------
			std::vector<int> nkl(nks_special,0);

			//------------------------------------------
			// cartesian coordinates of special points.
			//------------------------------------------
			std::vector<double> ksx(nks_special);
			std::vector<double> ksy(nks_special);
			std::vector<double> ksz(nks_special);

			//recalculate nkstot.
			nkstot = 0;
			for(int iks=0; iks<nks_special; iks++)
			{
				ifk >> ksx[iks];
				ifk >> ksy[iks];
				ifk >> ksz[iks];
				ModuleBase::GlobalFunc::READ_VALUE( ifk, nkl[iks] );
				//std::cout << " nkl[" << iks << "]=" << nkl[iks] << std::endl;
				assert(nkl[iks] >= 0);
				nkstot += nkl[iks];
			}
			assert( nkl[nks_special-1] == 1);

			//std::cout << " nkstot = " << nkstot << std::endl;
        	this->renew(nkstot * nspin);//mohan fix bug 2009-09-01

			int count = 0;
			for(int iks=1; iks<nks_special; iks++)
			{
				double dx = (ksx[iks] - ksx[iks-1]) / nkl[iks-1];
				double dy = (ksy[iks] - ksy[iks-1]) / nkl[iks-1];
				double dz = (ksz[iks] - ksz[iks-1]) / nkl[iks-1];
//				GlobalV::ofs_running << " dx=" << dx << " dy=" << dy << " dz=" << dz << std::endl;
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

			//std::cout << " count = " << count << std::endl;
			assert (count == nkstot );

			for(int ik=0; ik<nkstot; ik++)
			{
				wk[ik] = 1.0;
			}

            this->kc_done = true;

		}

		else if (kword == "Line_Direct" || kword == "L" || kword == "Line" )
		{
			//std::cout << " kword = " << kword << std::endl;
			if(ModuleSymmetry::Symmetry::symm_flag == 1)
			{
				ModuleBase::WARNING("K_Vectors::read_kpoints","Line mode of k-points is open, please set symmetry to 0 or -1.");
				return 0;
			}


			// how many special points.
			int nks_special = this->nkstot;
			//std::cout << " nks_special = " << nks_special << std::endl;

			//------------------------------------------
			// number of points to the next k points
			//------------------------------------------
			std::vector<int> nkl(nks_special,0);

			//------------------------------------------
			// cartesian coordinates of special points.
			//------------------------------------------
			std::vector<double> ksx(nks_special);
			std::vector<double> ksy(nks_special);
			std::vector<double> ksz(nks_special);

			//recalculate nkstot.
			nkstot = 0;
			for(int iks=0; iks<nks_special; iks++)
			{
				ifk >> ksx[iks];
				ifk >> ksy[iks];
				ifk >> ksz[iks];
				ModuleBase::GlobalFunc::READ_VALUE( ifk, nkl[iks] );
				//std::cout << " nkl[" << iks << "]=" << nkl[iks] << std::endl;
				assert(nkl[iks] >= 0);
				nkstot += nkl[iks];
			}
			assert( nkl[nks_special-1] == 1);

			//std::cout << " nkstot = " << nkstot << std::endl;
        	this->renew(nkstot * nspin);//mohan fix bug 2009-09-01

			int count = 0;
			for(int iks=1; iks<nks_special; iks++)
			{
				double dx = (ksx[iks] - ksx[iks-1]) / nkl[iks-1];
				double dy = (ksy[iks] - ksy[iks-1]) / nkl[iks-1];
				double dz = (ksz[iks] - ksz[iks-1]) / nkl[iks-1];
//				GlobalV::ofs_running << " dx=" << dx << " dy=" << dy << " dz=" << dz << std::endl;
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

			//std::cout << " count = " << count << std::endl;
			assert (count == nkstot );

			for(int ik=0; ik<nkstot; ik++)
			{
				wk[ik] = 1.0;
			}

            this->kd_done = true;

		}

        else
        {
			GlobalV::ofs_warning << " Error : neither Cartesian nor Direct kpoint." << std::endl;
			return 0;
        }
    }

    this->nkstot_full = this->nks = this->nkstot;

	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nkstot",nkstot);
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
    if (GlobalV::test_kpoint) ModuleBase::TITLE("K_Vectors", "Monkhorst_Pack");
    const int mpnx = nmp_in[0];
    const int mpny = nmp_in[1];
    const int mpnz = nmp_in[2];

    this->nkstot = mpnx * mpny * mpnz;
    // only can renew after nkstot is estimated.
    this->renew(nkstot * nspin); // mohan fix bug 2009-09-01

    for (int x = 1;x <= mpnx;x++)
    {
        double v1 = Monkhorst_Pack_formula( k_type, koffset_in[0], x, mpnx);
		if( std::abs(v1) < 1.0e-10 ) v1 = 0.0; //mohan update 2012-06-10
        for (int y = 1;y <= mpny;y++)
        {
            double v2 = Monkhorst_Pack_formula( k_type, koffset_in[1], y, mpny);
		    if( std::abs(v2) < 1.0e-10 ) v2 = 0.0;
            for (int z = 1;z <= mpnz;z++)
            {
                double v3 = Monkhorst_Pack_formula( k_type, koffset_in[2], z, mpnz);
				if( std::abs(v3) < 1.0e-10 ) v3 = 0.0;
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
    if (GlobalV::MY_RANK!=0) return;
    ModuleBase::TITLE("K_Vectors","update_use_ibz");
    assert( nkstot_ibz > 0 );

    // update nkstot
    this->nkstot = this->nkstot_ibz;

	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nkstot now",nkstot);

    this->kvec_d.resize(this->nkstot * nspin); //qianrui fix a bug 2021-7-13 for nspin=2 in set_kup_and_kdw()

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

void K_Vectors::ibz_kpoint(const ModuleSymmetry::Symmetry &symm, bool use_symm,std::string& skpt, const UnitCell &ucell, bool& match)
{
    if (GlobalV::MY_RANK != 0) return;
    ModuleBase::TITLE("K_Vectors", "ibz_kpoint");

    // k-lattice: "pricell" of reciprocal space
    // CAUTION: should fit into all k-input method, not only MP  !!!
    ModuleBase::Vector3<double> gb1(ucell.G.e11, ucell.G.e12, ucell.G.e13);
    ModuleBase::Vector3<double> gb2(ucell.G.e21, ucell.G.e22, ucell.G.e23);
    ModuleBase::Vector3<double> gb3(ucell.G.e31, ucell.G.e32, ucell.G.e33);
    ModuleBase::Vector3<double> gk1, gk2, gk3;
    ModuleBase::Matrix3 gk;
    if (this->is_mp)
    {
        gk1 = ModuleBase::Vector3<double>(gb1.x / nmp[0], gb1.y / nmp[0], gb1.z / nmp[0]);
        gk2 = ModuleBase::Vector3<double>(gb2.x / nmp[1], gb2.y / nmp[1], gb2.z / nmp[1]);
        gk3 = ModuleBase::Vector3<double>(gb3.x / nmp[2], gb3.y / nmp[2], gb3.z / nmp[2]);
        gk = ModuleBase::Matrix3(gk1.x, gk1.y, gk1.z, gk2.x, gk2.y, gk2.z, gk3.x, gk3.y, gk3.z);
    }


    //===============================================
    // search in all space group operations
    // if the operations does not already included
    // inverse operation, double it.
    //===============================================
    bool include_inv = false;
    std::vector<ModuleBase::Matrix3> kgmatrix(48 * 2);
    ModuleBase::Matrix3 inv(-1, 0, 0, 0, -1, 0, 0, 0, -1);
    ModuleBase::Matrix3 ind(1, 0, 0, 0, 1, 0, 0, 0, 1);

    int nrotkm;
    if(use_symm)
    {
        // bravais type of reciprocal lattice and k-lattice
        double b_const[6];
        double b0_const[6];
        double bk_const[6];
        double bk0_const[6];
        int bbrav=15;
        int bkbrav=15;
        std::string bbrav_name;
        std::string bkbrav_name;
        ModuleBase::Vector3<double>gk01 = gk1, gk02 = gk2, gk03 = gk3;

        ModuleBase::Matrix3 b_optlat = symm.optlat.Inverse().Transpose();
        //search optlat after using reciprocity relation
        ModuleBase::Vector3<double> gb01(b_optlat.e11, b_optlat.e12, b_optlat.e13);
        ModuleBase::Vector3<double> gb02(b_optlat.e21, b_optlat.e22, b_optlat.e23);
        ModuleBase::Vector3<double> gb03(b_optlat.e31, b_optlat.e32, b_optlat.e33);
        symm.lattice_type(gb1, gb2, gb3, gb01, gb02, gb03, b_const, b0_const, bbrav, bbrav_name, ucell.atoms, false, nullptr);
        GlobalV::ofs_running<<"(for reciprocal lattice: )"<<std::endl;
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"BRAVAIS TYPE", bbrav);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"BRAVAIS LATTICE NAME", bbrav_name);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "ibrav", bbrav);
        
        // the map of bravis lattice from real to reciprocal space
        // for example, 3(fcc) in real space matches 2(bcc) in reciprocal space
        std::vector<int> ibrav_a2b{ 1, 3, 2, 4, 5, 6, 7, 8, 10, 9, 11, 12, 13, 14 };
        auto ibrav_match = [&](int ibrav_b) -> bool
        {
            const int& ibrav_a = symm.real_brav;
            if (ibrav_a < 1 || ibrav_a > 14) return false;
            return (ibrav_b == ibrav_a2b[ibrav_a - 1]);
        };
        if (!ibrav_match(bbrav))
        {
            GlobalV::ofs_running << "Error: Bravais lattice type of reciprocal lattice is not compatible with that of real space lattice:" << std::endl;
            GlobalV::ofs_running << "ibrav of real space lattice: " << symm.ilattname << std::endl;
            GlobalV::ofs_running << "ibrav of reciprocal lattice: " << bbrav_name << std::endl;
            GlobalV::ofs_running << "(which should be " << ibrav_a2b[symm.real_brav - 1] << ")." << std::endl;
            match = false;
            return;
        }
        if (this->is_mp)
        {
            symm.lattice_type(gk1, gk2, gk3, gk01, gk02, gk03, bk_const, bk0_const, bkbrav, bkbrav_name, ucell.atoms, false, nullptr);
            GlobalV::ofs_running << "(for k-lattice: )" << std::endl;
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "BRAVAIS TYPE", bkbrav);
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "BRAVAIS LATTICE NAME", bkbrav_name);
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "ibrav", bkbrav);
        }
        // point-group analysis of reciprocal lattice
        ModuleBase::Matrix3 bsymop[48];
        int bnop = 0;
        // search again
        symm.lattice_type(gb1, gb2, gb3, gb1, gb2, gb3, b_const, b0_const, bbrav, bbrav_name, ucell.atoms, false, nullptr);
        ModuleBase::Matrix3 b_optlat_new(gb1.x, gb1.y, gb1.z, gb2.x, gb2.y, gb2.z, gb3.x, gb3.y, gb3.z);
        symm.setgroup(bsymop, bnop, bbrav);
        symm.gmatrix_convert(bsymop, bsymop, bnop, b_optlat_new, ucell.G);
        
        //check if all the kgmatrix are in bsymop
        auto matequal = [&symm] (ModuleBase::Matrix3 a, ModuleBase::Matrix3 b)
        {
            return (symm.equal(a.e11, b.e11) && symm.equal(a.e12, b.e12) && symm.equal(a.e13, b.e13) &&
            symm.equal(a.e21, b.e21) && symm.equal(a.e22, b.e22) && symm.equal(a.e23, b.e23) &&
            symm.equal(a.e31, b.e31) && symm.equal(a.e23, b.e23) && symm.equal(a.e33, b.e33));
        };
        for(int i=0;i<symm.nrotk;++i)
        {
            match = false;
            for(int j=0;j<bnop;++j) 
            {
                if (matequal(symm.kgmatrix[i], bsymop[j])) {match=true; break;}
            }
            if (!match) return;
        }
        nrotkm = symm.nrotk;// change if inv not included
        for (int i = 0; i < nrotkm; ++i)
        {
            if (symm.kgmatrix[i] == inv)
            {
                include_inv = true;
            }
            kgmatrix[i] = symm.kgmatrix[i];
        }

        if (!include_inv)
        {
            for (int i = 0; i<symm.nrotk; ++i)
            {
                kgmatrix[i + symm.nrotk] = inv * symm.kgmatrix[i];
            }
            nrotkm = 2 * symm.nrotk;
        }
    }
    else if(is_mp) // only include for mp grid
    {
        nrotkm = 2;
        kgmatrix[0] = ind;
        kgmatrix[1] = inv;
    }
    else
    {
        return;
    }

    // convert kgmatrix to k-lattice
    ModuleBase::Matrix3* kkmatrix = new ModuleBase::Matrix3[nrotkm];
    if (this->is_mp)symm.gmatrix_convert(kgmatrix.data(), kkmatrix, nrotkm, ucell.G, gk);
    // direct coordinates of k-points in k-lattice
    std::vector<ModuleBase::Vector3<double>> kvec_d_k(nkstot);
    if (this->is_mp) for (int i = 0;i < nkstot;++i) kvec_d_k[i] = kvec_d[i] * ucell.G * gk.Inverse();

    // use operation : kgmatrix to find
    // the new set kvec_d : ir_kpt
    this->nkstot_ibz = 0;

    assert(nkstot > 0);
    kvec_d_ibz.resize(this->nkstot);
    wk_ibz.resize(this->nkstot);
    ibz2bz.resize(this->nkstot);

	// nkstot is the total input k-points number.
    const double weight = 1.0 / static_cast<double>(nkstot);

    ModuleBase::Vector3<double> kvec_rot;
    ModuleBase::Vector3<double> kvec_rot_k;


//	for(int i=0; i<nrotkm; i++)
//	{
//		out.printM3("rot matrix",kgmatrix[i]);
//	}
    auto restrict_kpt = [&symm](ModuleBase::Vector3<double> &kvec){
        // in (-0.5, 0.5]
        kvec.x = fmod(kvec.x + 100.5-0.5*symm.epsilon, 1)-0.5+0.5*symm.epsilon;
        kvec.y = fmod(kvec.y + 100.5-0.5*symm.epsilon, 1)-0.5+0.5*symm.epsilon;
        kvec.z = fmod(kvec.z + 100.5-0.5*symm.epsilon, 1)-0.5+0.5*symm.epsilon;
        // in [0, 1)
        // kvec.x = fmod(kvec.x + 100 + symm.epsilon, 1) - symm.epsilon;
        // kvec.y = fmod(kvec.y + 100 + symm.epsilon, 1) - symm.epsilon;
        // kvec.z = fmod(kvec.z + 100 + symm.epsilon, 1) - symm.epsilon;
        if(std::abs(kvec.x) < symm.epsilon) kvec.x = 0.0;
        if(std::abs(kvec.y) < symm.epsilon) kvec.y = 0.0;
        if(std::abs(kvec.z) < symm.epsilon) kvec.z = 0.0;
        return;
    };
    // for output in kpoints file
    int ibz_index[nkstot];
	// search in all k-poins.
    for (int i = 0; i < nkstot; ++i)
    {
        //restrict to [0, 1)
        restrict_kpt(kvec_d[i]);

		//std::cout << "\n kpoint = " << i << std::endl;
		//std::cout << "\n kvec_d = " << kvec_d[i].x << " " << kvec_d[i].y << " " << kvec_d[i].z;
        bool already_exist = false;
		int exist_number = -1;
        // search over all symmetry operations
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
				kvec_rot = kvec_d[i] * kgmatrix[j]; //wrong for total energy, but correct for nonlocal force.
				//kvec_rot = kgmatrix[j] * kvec_d[i]; //correct for total energy, but wrong for nonlocal force.
                restrict_kpt(kvec_rot);
                if (this->is_mp)
                {
                    kvec_rot_k = kvec_d_k[i] * kkmatrix[j];   //k-lattice rotation
                    kvec_rot_k = kvec_rot_k * gk * ucell.G.Inverse(); //convert to recip lattice
                    restrict_kpt(kvec_rot_k);

                    assert(symm.equal(kvec_rot.x, kvec_rot_k.x));
                    assert(symm.equal(kvec_rot.y, kvec_rot_k.y));
                    assert(symm.equal(kvec_rot.z, kvec_rot_k.z));
                    // std::cout << "\n kvec_rot (in recip) = " << kvec_rot.x << " " << kvec_rot.y << " " << kvec_rot.z;
                    // std::cout << "\n kvec_rot(k to recip)= " << kvec_rot_k.x << " " << kvec_rot_k.y << " " << kvec_rot_k.z;
                    kvec_rot_k = kvec_rot_k * ucell.G * gk.Inverse(); //convert back to k-latice
                }
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
        // if really there is no equivalent k point in the list, then add it.
        if (!already_exist)
        {
			//if it's a new ibz kpoint.
			//nkstot_ibz indicate the index of ibz kpoint.
            this->kvec_d_ibz[nkstot_ibz] = kvec_rot;
            // output in kpoints file
            ibz_index[i] = nkstot_ibz;

			//the weight should be averged k-point weight.
            this->wk_ibz[nkstot_ibz] = weight;

			//ibz2bz records the index of origin k points.
            this->ibz2bz[nkstot_ibz] = i;
            ++nkstot_ibz;
        }
		else //mohan fix bug 2010-1-30
		{
//			std::cout << "\n\n already exist ! ";

//			std::cout << "\n kvec_rot = " << kvec_rot.x << " " << kvec_rot.y << " " << kvec_rot.z;
//			std::cout << "\n kvec_d_ibz = " << kvec_d_ibz[exist_number].x
//			<< " " << kvec_d_ibz[exist_number].y
//			<< " " << kvec_d_ibz[exist_number].z;

			double kmol_new = kvec_d[i].norm2();
			double kmol_old = kvec_d_ibz[exist_number].norm2();

            ibz_index[i] = exist_number;

//			std::cout << "\n kmol_new = " << kmol_new;
//			std::cout << "\n kmol_old = " << kmol_old;


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

    delete[] kkmatrix;
    // output in kpoints file
    std::stringstream ss;
    ss << " " << std::setw(40) <<"nkstot" << " = " << nkstot
        << std::setw(66) << "ibzkpt" << std::endl;
    formatter::ContextFmt fmt;
    fmt.set_context({"int_w8", "double_w12_f8", "double_w12_f8", "double_w12_f8", "int_w8", "double_w12_f8", "double_w12_f8", "double_w12_f8"});
    fmt.set_titles({"KPT", "DIRECT_X", "DIRECT_Y", "DIRECT_Z", "IBZ", "DIRECT_X", "DIRECT_Y", "DIRECT_Z"});
    fmt.set_overall_title("K-POINTS REDUCTION ACCORDING TO SYMMETRY");
    fmt.enable_title(); fmt.disable_down_frame(); fmt.disable_up_frame(); fmt.disable_mid_frame();
    fmt.right_title();

    std::vector<int> _kpt;
    std::vector<double> _direct_x;
    std::vector<double> _direct_y;
    std::vector<double> _direct_z;
    std::vector<int> _ibz;
    std::vector<double> _direct_x_ibz;
    std::vector<double> _direct_y_ibz;
    std::vector<double> _direct_z_ibz;
    // ss << " " << std::setw(8) << "KPT"
    //     << std::setw(20) << "DirectX"
	//     << std::setw(20) << "DirectY"
    //     << std::setw(20) << "DirectZ"
    //      << std::setw(8) << "IBZ"
    //     << std::setw(20) << "DirectX"
	//     << std::setw(20) << "DirectY"
    //     << std::setw(20) << "DirectZ" << std::endl;
    for (int i = 0; i < nkstot; ++i)
    {
        // ss << " "
        //     << std::setw(8) << i+1
        //     << std::setw(20) << this->kvec_d[i].x
        //     << std::setw(20) << this->kvec_d[i].y
        //     << std::setw(20) << this->kvec_d[i].z;
        // ss << std::setw(8) << ibz_index[i]+1
        //         << std::setw(20) << this->kvec_d_ibz[ibz_index[i]].x
        //         << std::setw(20) << this->kvec_d_ibz[ibz_index[i]].y
        //         << std::setw(20) << this->kvec_d_ibz[ibz_index[i]].z << std::endl;
        _kpt.push_back(i+1);
        _direct_x.push_back(this->kvec_d[i].x);
        _direct_y.push_back(this->kvec_d[i].y);
        _direct_z.push_back(this->kvec_d[i].z);
        _ibz.push_back(ibz_index[i]+1);
        _direct_x_ibz.push_back(this->kvec_d_ibz[ibz_index[i]].x);
        _direct_y_ibz.push_back(this->kvec_d_ibz[ibz_index[i]].y);
        _direct_z_ibz.push_back(this->kvec_d_ibz[ibz_index[i]].z);
    }
    fmt << _kpt << _direct_x << _direct_y << _direct_z << _ibz << _direct_x_ibz << _direct_y_ibz << _direct_z_ibz;
    ss << fmt.str() << std::endl;
    skpt = ss.str();

	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nkstot_ibz",nkstot_ibz);
    
    fmt.set_context({"int_w8", "double_w12_f8", "double_w12_f8", "double_w12_f8", "double_w6_f4", "int_w4"});
    fmt.set_titles({"IBZ", "DIRECT_X", "DIRECT_Y", "DIRECT_Z", "WEIGHT", "ibz2bz"});
    fmt.set_overall_title("K-POINTS REDUCTION ACCORDING TO SYMMETRY");
    fmt.enable_title(); fmt.disable_down_frame(); fmt.disable_up_frame(); fmt.disable_mid_frame();
    fmt.right_title();
/* 	GlobalV::ofs_running << " " << std::setw(8) << "IBZ" << std::setw(20) << "DirectX"
	<< std::setw(20) << "DirectY" << std::setw(20) << "DirectZ"
	<< std::setw(20) << "Weight" << std::setw(10) << "ibz2bz" << std::endl; */
    _ibz.clear();
    _direct_x.clear();
    _direct_y.clear();
    _direct_z.clear();
    std::vector<double> _weights;
    std::vector<int> _ibz2bz;
    for (int ik=0; ik<nkstot_ibz; ik++)
    {
/*         GlobalV::ofs_running << " "
			<< std::setw(8) << ik+1
            << std::setw(20) << this->kvec_d_ibz[ik].x
            << std::setw(20) << this->kvec_d_ibz[ik].y
            << std::setw(20) << this->kvec_d_ibz[ik].z
            << std::setw(20) << this->wk_ibz[ik]
            << std::setw(10) << this->ibz2bz[ik] << std::endl; */
        _ibz.push_back(ik+1);
        _direct_x.push_back(this->kvec_d_ibz[ik].x);
        _direct_y.push_back(this->kvec_d_ibz[ik].y);
        _direct_z.push_back(this->kvec_d_ibz[ik].z);
        _weights.push_back(this->wk_ibz[ik]);
        _ibz2bz.push_back(this->ibz2bz[ik]);
    }
    fmt << _ibz << _direct_x << _direct_y << _direct_z << _weights << _ibz2bz;
    GlobalV::ofs_running << fmt.str() << std::endl;

    return;
}


void K_Vectors::set_both_kvec(const ModuleBase::Matrix3 &G, const ModuleBase::Matrix3 &R,std::string& skpt)
{

    if(GlobalV::FINAL_SCF) //LiuXh add 20180606
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
                GlobalV::ofs_warning << " Error : neither Cartesian nor Direct kpoint." << std::endl;
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
			if( std::abs(kvec_d[i].x) < 1.0e-10 ) kvec_d[i].x = 0.0;
			if( std::abs(kvec_d[i].y) < 1.0e-10 ) kvec_d[i].y = 0.0;
			if( std::abs(kvec_d[i].z) < 1.0e-10 ) kvec_d[i].z = 0.0;

			kvec_c[i] = kvec_d[i] * G;

			// mohan add2012-06-10
			if( std::abs(kvec_c[i].x) < 1.0e-10 ) kvec_c[i].x = 0.0;
			if( std::abs(kvec_c[i].y) < 1.0e-10 ) kvec_c[i].y = 0.0;
			if( std::abs(kvec_c[i].z) < 1.0e-10 ) kvec_c[i].z = 0.0;
        }
        kc_done = true;
    }

    // set direct k vectors
    else if (kc_done && !kd_done)
    {
        ModuleBase::Matrix3 RT = R.Transpose();
        for (int i = 0;i < nkstot;i++)
        {
//			std::cout << " ik=" << i
//				<< " kvec.x=" << kvec_c[i].x
//				<< " kvec.y=" << kvec_c[i].y
//				<< " kvec.z=" << kvec_c[i].z << std::endl;
//wrong!            kvec_d[i] = RT * kvec_c[i];
// mohan fixed bug 2011-03-07
            kvec_d[i] = kvec_c[i] * RT;
        }
        kd_done = true;
    }
    formatter::ContextFmt fmt;
    fmt.set_context({"int_w8", "double_w12_f8", "double_w12_f8", "double_w12_f8", "double_w6_f4"});
    fmt.set_titles({"KPOINTS", "DIRECT_X", "DIRECT_Y", "DIRECT_Z", "WEIGHT"});
    fmt.set_overall_title("K-POINTS DIRECT COORDINATES");
    fmt.enable_title(); fmt.disable_down_frame(); fmt.disable_up_frame(); fmt.disable_mid_frame();
    fmt.right_title();
	// GlobalV::ofs_running << "\n " << std::setw(8) << "KPOINTS"
	// << std::setw(20) << "DIRECT_X"
	// << std::setw(20) << "DIRECT_Y"
	// << std::setw(20) << "DIRECT_Z"
	// << std::setw(20) << "WEIGHT" << std::endl;
    std::vector<int> _kpoints;
    std::vector<double> _direct_x;
    std::vector<double> _direct_y;
    std::vector<double> _direct_z;
    std::vector<double> _weights;
	for(int i=0; i<nkstot; i++)
	{
        //  GlobalV::ofs_running << " "
		// 	<< std::setw(8) << i+1
        //      << std::setw(20) << this->kvec_d[i].x
        //      << std::setw(20) << this->kvec_d[i].y
        //      << std::setw(20) << this->kvec_d[i].z
        //      << std::setw(20) << this->wk[i] << std::endl;
        _kpoints.push_back(i+1);
        _direct_x.push_back(this->kvec_d[i].x);
        _direct_y.push_back(this->kvec_d[i].y);
        _direct_z.push_back(this->kvec_d[i].z);
        _weights.push_back(this->wk[i]);
	}
    fmt << _kpoints << _direct_x << _direct_y << _direct_z << _weights;
    GlobalV::ofs_running << fmt.str() << std::endl;

	if(GlobalV::MY_RANK==0)
	{
        std::stringstream ss;
        ss << " " << std::setw(40) <<"nkstot now" << " = " << nkstot << std::endl;
        fmt.set_context({"int_w8", "double_w12_f8", "double_w12_f8", "double_w12_f8", "double_w6_f4"});
        fmt.set_titles({"KPOINTS", "DIRECT_X", "DIRECT_Y", "DIRECT_Z", "WEIGHT"});
        fmt.set_overall_title("K-POINTS DIRECT COORDINATES");
        fmt.enable_title(); fmt.disable_down_frame(); fmt.disable_up_frame(); fmt.disable_mid_frame();
        fmt.right_title();
        fmt << _kpoints << _direct_x << _direct_y << _direct_z << _weights;
        ss << fmt.str() << std::endl;
        // ss << " " << std::setw(8) << "KPT" << std::setw(20) << "DirectX"
	    //     << std::setw(20) << "DirectY" << std::setw(20) << "DirectZ"
	    //     << std::setw(20) << "Weight" << std::endl;
        // for (int ik=0; ik<nkstot; ik++)
        // {
        //     ss << " " << std::setw(8) << ik+1
        //         << std::setw(20) << this->kvec_d[ik].x
        //         << std::setw(20) << this->kvec_d[ik].y
        //         << std::setw(20) << this->kvec_d[ik].z
        //         << std::setw(20) << this->wk[ik] << std::endl;
        //       //  << std::setw(10) << this->ibz2bz[ik] << std::endl;
        // }
        skpt = ss.str();
	}

    return;
}


void K_Vectors::normalize_wk(const int &degspin)
{
	if(GlobalV::MY_RANK!=0) return;
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

#ifdef __MPI
void K_Vectors::mpi_k(void)
{
	ModuleBase::TITLE("K_Vectors","mpi_k");

    Parallel_Common::bcast_bool(kc_done);

    Parallel_Common::bcast_bool(kd_done);

    Parallel_Common::bcast_int(nspin);

    Parallel_Common::bcast_int(nkstot);

    Parallel_Common::bcast_int(nkstot_full);

    Parallel_Common::bcast_int(nmp, 3);

    Parallel_Common::bcast_double(koffset, 3);

    this->nks = GlobalC::Pkpoints.nks_pool[GlobalV::MY_POOL];

	GlobalV::ofs_running << std::endl;
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"k-point number in this process",nks);
    int nks_minimum = this->nks;

	Parallel_Reduce::gather_min_int_all( nks_minimum );

    if (nks_minimum == 0)
    {
        ModuleBase::WARNING_QUIT("K_Vectors::mpi_k()"," nks == 0, some processor have no k point!");
    }
    else
    {
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"minimum distributed K point number",nks_minimum);
    }

    std::vector<int> isk_aux(nkstot);
    std::vector<double> wk_aux(nkstot);
    std::vector<double> kvec_c_aux(nkstot*3);
    std::vector<double> kvec_d_aux(nkstot*3);

    if (GlobalV::MY_RANK == 0)
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

    Parallel_Common::bcast_int(isk_aux.data(), nkstot);

    Parallel_Common::bcast_double(wk_aux.data(), nkstot);
    Parallel_Common::bcast_double(kvec_c_aux.data(), nkstot*3);
    Parallel_Common::bcast_double(kvec_d_aux.data(), nkstot*3);

    this->renew(this->nks * this->nspin);

    // distribute
    int k_index = 0;

    for (int i = 0;i < nks;i++)
    {
        // 3 is because each k point has three value:kx, ky, kz
        k_index = i + GlobalC::Pkpoints.startk_pool[GlobalV::MY_POOL] ;
        kvec_c[i].x = kvec_c_aux[k_index*3];
        kvec_c[i].y = kvec_c_aux[k_index*3+1];
        kvec_c[i].z = kvec_c_aux[k_index*3+2];
        kvec_d[i].x = kvec_d_aux[k_index*3];
        kvec_d[i].y = kvec_d_aux[k_index*3+1];
        kvec_d[i].z = kvec_d_aux[k_index*3+2];
        wk[i] = wk_aux[k_index];
        isk[i] = isk_aux[k_index];
    }
} // END SUBROUTINE
#endif


//----------------------------------------------------------
// This routine sets the k vectors for the up and down spin
//----------------------------------------------------------
// from set_kup_and_kdw.f90
void K_Vectors::set_kup_and_kdw(void)
{
    ModuleBase::TITLE("K_Vectors", "setup_kup_and_kdw");

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

		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nks(nspin=2)",nks);
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nkstot(nspin=2)",nkstot);
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


void K_Vectors::print_klists(std::ofstream &ofs)
{
    ModuleBase::TITLE("K_Vectors", "print_klists");

    if (nkstot < nks)
    {
        std::cout << "\n nkstot=" << nkstot;
        std::cout << "\n nks=" << nks;
        ModuleBase::WARNING_QUIT("print_klists","nkstot < nks");
    }
    formatter::ContextFmt fmt;
    fmt.set_context({"int_w8", "double_w12_f8", "double_w12_f8", "double_w12_f8", "double_w6_f4"});
    fmt.set_titles({"KPOINTS", "CARTESIAN_X", "CARTESIAN_Y", "CARTESIAN_Z", "WEIGHT"});
    fmt.set_overall_title("K-POINTS CARTESIAN COORDINATES");
    fmt.enable_title(); fmt.disable_down_frame(); fmt.disable_up_frame(); fmt.disable_mid_frame();
    fmt.right_title();
	// GlobalV::ofs_running << "\n " << std::setw(8) << "KPOINTS"
	// << std::setw(20) << "CARTESIAN_X"
	// << std::setw(20) << "CARTESIAN_Y"
	// << std::setw(20) << "CARTESIAN_Z"
	// << std::setw(20) << "WEIGHT" << std::endl;
    std::vector<int> _kpoints;
    std::vector<double> _cartesian_x; std::vector<double> _cartesian_y; std::vector<double> _cartesian_z;
    std::vector<double> _weights;
	for(int i=0; i<nks; i++)
	{
        // GlobalV::ofs_running << " "
		// 	<< std::setw(8) << i+1
        //      << std::setw(20) << this->kvec_c[i].x
        //      << std::setw(20) << this->kvec_c[i].y
        //      << std::setw(20) << this->kvec_c[i].z
        //      << std::setw(20) << this->wk[i] << std::endl;
        _kpoints.push_back(i+1);
        _cartesian_x.push_back(this->kvec_c[i].x);
        _cartesian_y.push_back(this->kvec_c[i].y);
        _cartesian_z.push_back(this->kvec_c[i].z);
        _weights.push_back(this->wk[i]);
	}
    fmt << _kpoints << _cartesian_x << _cartesian_y << _cartesian_z << _weights;
    GlobalV::ofs_running << "\n" << fmt.str() << std::endl;

    fmt.set_context({"int_w8", "double_w12_f8", "double_w12_f8", "double_w12_f8", "double_w6_f4"});
    fmt.set_titles({"KPOINTS", "DIRECT_X", "DIRECT_Y", "DIRECT_Z", "WEIGHT"});
    fmt.set_overall_title("K-POINTS DIRECT COORDINATES");
    fmt.enable_title(); fmt.disable_down_frame(); fmt.disable_up_frame(); fmt.disable_mid_frame();
    fmt.right_title();
	// GlobalV::ofs_running << "\n " << std::setw(8) << "KPOINTS"
	// << std::setw(20) << "DIRECT_X"
	// << std::setw(20) << "DIRECT_Y"
	// << std::setw(20) << "DIRECT_Z"
	// << std::setw(20) << "WEIGHT" << std::endl;
    std::vector<double> _direct_x; std::vector<double> _direct_y; std::vector<double> _direct_z;
	for(int i=0; i<nks; i++)
	{
        // GlobalV::ofs_running << " "
		// 	<< std::setw(8) << i+1
        //      << std::setw(20) << this->kvec_d[i].x
        //      << std::setw(20) << this->kvec_d[i].y
        //      << std::setw(20) << this->kvec_d[i].z
        //      << std::setw(20) << this->wk[i] << std::endl;
        _direct_x.push_back(this->kvec_d[i].x);
        _direct_y.push_back(this->kvec_d[i].y);
        _direct_z.push_back(this->kvec_d[i].z);
	}
    fmt << _kpoints << _direct_x << _direct_y << _direct_z << _weights;
    GlobalV::ofs_running << "\n" << fmt.str() << std::endl;
    return;
}

//LiuXh add a new function here,
//20180515
void K_Vectors::set_after_vc(
        const ModuleSymmetry::Symmetry &symm,
        const std::string &k_file_name,
        const int& nspin_in,
        const ModuleBase::Matrix3 &reciprocal_vec,
        const ModuleBase::Matrix3 &latvec)
{
    ModuleBase::TITLE("K_Vectors", "set_after_vc");

    GlobalV::ofs_running << "\n SETUP K-POINTS" << std::endl;
    this->nspin = nspin_in;
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nspin",nspin);

    this->set_both_kvec_after_vc(reciprocal_vec, latvec);
    //this->set_both_kvec(reciprocal_vec, latvec);

    this->mpi_k_after_vc();

    this->set_kup_and_kdw_after_vc();

    this->print_klists(GlobalV::ofs_running);

    return;
}

//LiuXh add a new function here,
//20180515
void K_Vectors::mpi_k_after_vc(void)
{
#ifdef __MPI
    ModuleBase::TITLE("K_Vectors","mpi_k_after_vc");

    Parallel_Common::bcast_bool(kc_done);
    Parallel_Common::bcast_bool(kd_done);
    Parallel_Common::bcast_int(nspin);
    Parallel_Common::bcast_int(nkstot);
    Parallel_Common::bcast_int(nmp, 3);
    Parallel_Common::bcast_double(koffset, 3);

    this->nks = GlobalC::Pkpoints.nks_pool[GlobalV::MY_POOL];
    GlobalV::ofs_running << std::endl;
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"k-point number in this process",nks);
    int nks_minimum = this->nks;

    Parallel_Reduce::gather_min_int_all( nks_minimum );

    if (nks_minimum == 0)
    {
        ModuleBase::WARNING_QUIT("K_Vectors::mpi_k()"," nks == 0, some processor have no k point!");
    }
    else
    {
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"minimum distributed K point number",nks_minimum);
    }

    std::vector<int> isk_aux(nkstot);
    std::vector<double> wk_aux(nkstot);
    std::vector<double> kvec_c_aux(nkstot*3);
    std::vector<double> kvec_d_aux(nkstot*3);

    if (GlobalV::MY_RANK == 0)
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

    Parallel_Common::bcast_int(isk_aux.data(), nkstot);
    Parallel_Common::bcast_double(wk_aux.data(), nkstot);
    Parallel_Common::bcast_double(kvec_c_aux.data(), nkstot*3);
    Parallel_Common::bcast_double(kvec_d_aux.data(), nkstot*3);

    int k_index = 0;
    for (int i = 0;i < nks;i++)
    {
        k_index = i + GlobalC::Pkpoints.startk_pool[GlobalV::MY_POOL] ;
        kvec_c[i].x = kvec_c_aux[k_index*3];
        kvec_c[i].y = kvec_c_aux[k_index*3+1];
        kvec_c[i].z = kvec_c_aux[k_index*3+2];
        kvec_d[i].x = kvec_d_aux[k_index*3];
        kvec_d[i].y = kvec_d_aux[k_index*3+1];
        kvec_d[i].z = kvec_d_aux[k_index*3+2];
        wk[i] = wk_aux[k_index];
        isk[i] = isk_aux[k_index];
    }
#endif
}

void K_Vectors::set_both_kvec_after_vc(const ModuleBase::Matrix3 &G, const ModuleBase::Matrix3 &R)
{
    // set cartesian k vectors.
    kd_done = true;
    kc_done = false;
    if (!kc_done && kd_done)
    {
        for (int i = 0;i < nks;i++)
        {
//wrong!!   kvec_c[i] = G * kvec_d[i];
// mohan fixed bug 2010-1-10
			if( std::abs(kvec_d[i].x) < 1.0e-10 ) kvec_d[i].x = 0.0;
			if( std::abs(kvec_d[i].y) < 1.0e-10 ) kvec_d[i].y = 0.0;
			if( std::abs(kvec_d[i].z) < 1.0e-10 ) kvec_d[i].z = 0.0;

			kvec_c[i] = kvec_d[i] * G;

			// mohan add2012-06-10
			if( std::abs(kvec_c[i].x) < 1.0e-10 ) kvec_c[i].x = 0.0;
			if( std::abs(kvec_c[i].y) < 1.0e-10 ) kvec_c[i].y = 0.0;
			if( std::abs(kvec_c[i].z) < 1.0e-10 ) kvec_c[i].z = 0.0;
        }
        kc_done = true;
    }

    // set direct k vectors
    else if (kc_done && !kd_done)
    {
        ModuleBase::Matrix3 RT = R.Transpose();
        for (int i = 0;i < nks;i++)
        {
//			std::cout << " ik=" << i
//				<< " kvec.x=" << kvec_c[i].x
//				<< " kvec.y=" << kvec_c[i].y
//				<< " kvec.z=" << kvec_c[i].z << std::endl;
//wrong!            kvec_d[i] = RT * kvec_c[i];
// mohan fixed bug 2011-03-07
            kvec_d[i] = kvec_c[i] * RT;
        }
        kd_done = true;
    }
    formatter::ContextFmt fmt;
    fmt.set_context({"int_w8", "double_w12_f8", "double_w12_f8", "double_w12_f8", "double_w6_f4"});
    fmt.set_titles({"KPOINTS", "DIRECT_X", "DIRECT_Y", "DIRECT_Z", "WEIGHT"});
    fmt.set_overall_title("K-POINTS DIRECT COORDINATES");
    fmt.enable_title(); fmt.disable_down_frame(); fmt.disable_up_frame(); fmt.disable_mid_frame();
    fmt.right_title();
	// GlobalV::ofs_running << "\n " << std::setw(8) << "KPOINTS"
	// << std::setw(20) << "DIRECT_X"
	// << std::setw(20) << "DIRECT_Y"
	// << std::setw(20) << "DIRECT_Z"
	// << std::setw(20) << "WEIGHT" << std::endl;
    std::vector<int> _kpoints;
    std::vector<double> _direct_x; std::vector<double> _direct_y; std::vector<double> _direct_z;
    std::vector<double> _weights;
	for(int i=0; i<nks; i++)
	{
        // GlobalV::ofs_running << " "
		// 	<< std::setw(8) << i+1
        //      << std::setw(20) << this->kvec_d[i].x
        //      << std::setw(20) << this->kvec_d[i].y
        //      << std::setw(20) << this->kvec_d[i].z
        //      << std::setw(20) << this->wk[i] << std::endl;
        _kpoints.push_back(i+1);
        _direct_x.push_back(this->kvec_d[i].x);
        _direct_y.push_back(this->kvec_d[i].y);
        _direct_z.push_back(this->kvec_d[i].z);
        _weights.push_back(this->wk[i]);
	}
    fmt << _kpoints << _direct_x << _direct_y << _direct_z << _weights;
    GlobalV::ofs_running << fmt.str() << std::endl;
    return;
}

void K_Vectors::set_kup_and_kdw_after_vc(void)
{
    ModuleBase::TITLE("K_Vectors", "setup_kup_and_kdw_after_vc");

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

        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nks(nspin=2)",nks);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nkstot(nspin=2)",nkstot);

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
