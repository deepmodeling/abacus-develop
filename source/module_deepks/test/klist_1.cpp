#include "klist.h"

namespace Test_Deepks
{
	K_Vectors::K_Vectors()
	{	
		nspin = 0; // default spin.
		kc_done = false;
		kd_done = false;

		kvec_c = new ModuleBase::Vector3<double>[1];
		kvec_d.resize(1);

		wk = nullptr;
		isk = nullptr;

		nkstot = 0;
	}

	K_Vectors::~K_Vectors()
	{
		delete[] kvec_c;
		kvec_d.clear();
		delete[] wk;
		delete[] isk;
	}

	void K_Vectors::set(
		const std::string &k_file_name,
		const int& nspin_in,
		const ModuleBase::Matrix3 &reciprocal_vec,
		const ModuleBase::Matrix3 &latvec,
		bool &GAMMA_ONLY_LOCAL,
		std::ofstream &ofs_running,
		std::ofstream &ofs_warning)
	{

		ofs_running << "\n\n\n\n";
		ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
		ofs_running << " |                                                                    |" << std::endl;
		ofs_running << " | Setup K-points                                                     |" << std::endl;
		ofs_running << " | We setup the k-points according to input parameters.               |" << std::endl;
		ofs_running << " | The reduced k-points are set according to symmetry operations.     |" << std::endl;
		ofs_running << " | We treat the spin as another set of k-points.                      |" << std::endl;
		ofs_running << " |                                                                    |" << std::endl;
		ofs_running << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
		ofs_running << "\n\n\n\n";

		ofs_running << "\n SETUP K-POINTS" << std::endl;

		// (1) set nspin, read kpoints.
		this->nspin = nspin_in;
		ModuleBase::GlobalFunc::OUT(ofs_running,"nspin",nspin);
			
		bool read_succesfully = this->read_kpoints(
			k_file_name,
			GAMMA_ONLY_LOCAL,
			ofs_warning,
			ofs_running);
		if(!read_succesfully)
		{
			ofs_warning << "in K_Vectors::set, something wrong while reading KPOINTS." << std::endl;
			exit(1);
		}

		// (2)
		this->set_both_kvec(reciprocal_vec, latvec, ofs_running);

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
			ofs_warning << "In K_Vectors::set, Only available for nspin = 1 or 2 or 4" << std::endl;
			exit(1);
		}
		this->normalize_wk(deg);

		// It's very important in parallel case,
		// firstly do the mpi_k() and then
		// do set_kup_and_kdw()

		this->set_kup_and_kdw(ofs_running);

		this->print_klists(ofs_running);
		//std::cout << " NUMBER OF K-POINTS   : " << nkstot << std::endl;
		
		return;
	}

	void K_Vectors::renew(const int &kpoint_number)
	{
		delete[] kvec_c;
		delete[] wk;
		delete[] isk;

		kvec_c = new ModuleBase::Vector3<double>[kpoint_number];
		kvec_d.resize(kpoint_number);
		wk = new double[kpoint_number];
		isk = new int[kpoint_number];

		ModuleBase::Memory::record("KV::kvec_c",sizeof(double) * kpoint_number*3);
		ModuleBase::Memory::record("KV::kvec_d",sizeof(double) * kpoint_number*3);
		ModuleBase::Memory::record("KV::wk",sizeof(double) * kpoint_number*3);
		ModuleBase::Memory::record("KV::isk",sizeof(int) * kpoint_number*3);

		return;
	}

	bool K_Vectors::read_kpoints(const std::string &fn, bool &GAMMA_ONLY_LOCAL, std::ofstream &ofs_warning, std::ofstream &ofs_running)
	{

		std::ifstream ifk(fn.c_str());
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
			ofs_warning << " symbol K_POINTS not found." << std::endl;
			return 0;
		}

		//input k-points are in 2pi/a units
		ModuleBase::GlobalFunc::READ_VALUE(ifk, nkstot);

		//std::cout << " nkstot = " << nkstot << std::endl;
		ModuleBase::GlobalFunc::READ_VALUE(ifk, kword);

		// mohan update 2021-02-22
		int max_kpoints = 100000;
		if (nkstot > 100000)
		{
			ofs_warning << " nkstot > MAX_KPOINTS" << std::endl;
			return 0;
		}

		int k_type = 0;
		if (nkstot == 0) // nkstot==0, use monkhorst_pack. add by dwan
		{
			if (kword == "Gamma")
			{
				k_type = 0;
				ModuleBase::GlobalFunc::OUT(ofs_running,"Input type of k points","Monkhorst-Pack(Gamma)");
			}
			else if (kword == "Monkhorst-Pack" || kword == "MP" || kword == "mp")
			{
				k_type = 1;
				ModuleBase::GlobalFunc::OUT(ofs_running,"Input type of k points","Monkhorst-Pack");
			}
			else
			{
				ofs_warning << " Error: neither Gamma nor Monkhorst-Pack." << std::endl;
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
			
				// how many special points.	
				int nks_special = this->nkstot;
				//std::cout << " nks_special = " << nks_special << std::endl;
				
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
				std::vector<double> kposx;
				std::vector<double> kposy;
				std::vector<double> kposz;
				ModuleBase::GlobalFunc::ZEROS(nkl, nks_special);
				
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
				
				ofs_warning << " Error : nkstot == -1, not implemented yet." << std::endl;

				delete[] nkl;
				delete[] ksx;
				delete[] ksy;
				delete[] ksz;

				this->kc_done = true;

			}

			else if (kword == "Line_Direct" || kword == "L" || kword == "Line" )
			{
				//std::cout << " kword = " << kword << std::endl;
			
				// how many special points.	
				int nks_special = this->nkstot;
				//std::cout << " nks_special = " << nks_special << std::endl;
				
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
				std::vector<double> kposx;
				std::vector<double> kposy;
				std::vector<double> kposz;
				ModuleBase::GlobalFunc::ZEROS(nkl, nks_special);
				
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
				
				ofs_warning << " Error : nkstot == -1, not implemented yet." << std::endl;

				delete[] nkl;
				delete[] ksx;
				delete[] ksy;
				delete[] ksz;

				this->kd_done = true;

			}

			else
			{
				ofs_warning << " Error : neither Cartesian nor Direct kpoint." << std::endl;
				return 0;
			}
		}

		ModuleBase::GlobalFunc::OUT(ofs_running,"nkstot",nkstot);
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

	void K_Vectors::set_both_kvec(const ModuleBase::Matrix3 &G, const ModuleBase::Matrix3 &R, std::ofstream &ofs_running)
	{
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

		ofs_running << "\n " << std::setw(8) << "KPOINTS" 
		<< std::setw(20) << "DIRECT_X"
		<< std::setw(20) << "DIRECT_Y"
		<< std::setw(20) << "DIRECT_Z"
		<< std::setw(20) << "WEIGHT" << std::endl;

		for(int i=0; i<nkstot; i++)
		{
			ofs_running << " "
				<< std::setw(8) << i+1
				<< std::setw(20) << this->kvec_d[i].x
				<< std::setw(20) << this->kvec_d[i].y
				<< std::setw(20) << this->kvec_d[i].z
				<< std::setw(20) << this->wk[i] << std::endl;
		}

		return;
	}


	void K_Vectors::normalize_wk(const int &degspin)
	{
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

	//----------------------------------------------------------
	// This routine sets the k vectors for the up and down spin
	//----------------------------------------------------------
	// from set_kup_and_kdw.f90
	void K_Vectors::set_kup_and_kdw(std::ofstream &ofs_running)
	{
		//=========================================================================
		// on output: the number of points is doubled and xk and wk in the
		// first (nks/2) positions correspond to up spin
		// those in the second (nks/2) ones correspond to down spin
		//=========================================================================
		switch (nspin)
		{
		case 1:

			for (int ik = 0; ik < nkstot; ik++)
			{
				this->isk[ik] = 0;
			}

			break;

		case 2:

			for (int ik = 0; ik < nkstot; ik++)
			{
				this->kvec_c[ik+nkstot] = kvec_c[ik];
				this->kvec_d[ik+nkstot] = kvec_d[ik];
				this->wk[ik+nkstot]     = wk[ik];
				this->isk[ik]        = 0;
				this->isk[ik+nkstot]    = 1;
			}

			this->nkstot *= 2;

			ModuleBase::GlobalFunc::OUT(ofs_running,"nkstot(nspin=2)",nkstot);
			break;
		case 4:

			for (int ik = 0; ik < nkstot; ik++)
			{
				this->isk[ik] = 0;
			}

			break;
		}

		return;
	} // end subroutine set_kup_and_kdw


	void K_Vectors::print_klists(std::ofstream &ofs_running)
	{
		ofs_running << "\n " << std::setw(8) << "KPOINTS" 
		<< std::setw(20) << "CARTESIAN_X"
		<< std::setw(20) << "CARTESIAN_Y"
		<< std::setw(20) << "CARTESIAN_Z"
		<< std::setw(20) << "WEIGHT" << std::endl;
		for(int i=0; i<nkstot; i++)
		{
			ofs_running << " "
				<< std::setw(8) << i+1
				<< std::setw(20) << this->kvec_c[i].x
				<< std::setw(20) << this->kvec_c[i].y
				<< std::setw(20) << this->kvec_c[i].z
				<< std::setw(20) << this->wk[i] << std::endl;
		}

		ofs_running << "\n " << std::setw(8) << "KPOINTS" 
		<< std::setw(20) << "DIRECT_X"
		<< std::setw(20) << "DIRECT_Y"
		<< std::setw(20) << "DIRECT_Z"
		<< std::setw(20) << "WEIGHT" << std::endl;
		for(int i=0; i<nkstot; i++)
		{
			ofs_running << " "
				<< std::setw(8) << i+1
				<< std::setw(20) << this->kvec_d[i].x
				<< std::setw(20) << this->kvec_d[i].y
				<< std::setw(20) << this->kvec_d[i].z
				<< std::setw(20) << this->wk[i] << std::endl;
		}

		return;
	}

}
