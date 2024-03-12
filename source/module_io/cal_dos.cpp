#include "cal_dos.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/constants.h"
#include "module_base/parallel_reduce.h"

bool ModuleIO::calculate_dos
(
	const int &is,
	const std::string &fa, //file address for DOS
	const std::string &fa1, //file address for DOS_smearing
	const double &de_ev, // delta energy in ev
	const double &emax_ev,
	const double &emin_ev,// minimal energy in ev.
	const double &bcoeff,
	const int &nks,//number of k points
	const int &nkstot,
	const std::vector<double> &wk,//weight of k points
	const std::vector<int> &isk,
	const int &nbands,// number of bands
	const ModuleBase::matrix &ekb,//store energy for each k point and each band
	const ModuleBase::matrix &wg//weight of (kpoint,bands)
)
{
	ModuleBase::TITLE("ModuleIO","calculate_dos");
	std::ofstream ofs;
	std::ofstream ofs1;
	if(GlobalV::MY_RANK==0)
	{
		ofs.open(fa.c_str());//make the file clear!!
		ofs1.open(fa1.c_str());//make the file clear!!
	}
	std::vector<double> dos;
	std::vector<double> ene;
	std::vector<double> dos_smearing; //dos_smearing

#ifdef __MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	if(de_ev <= 0)
	{
		ModuleBase::WARNING("ModuleIO::calculate_dos","de <= 0 ");
		return 0; 
	}
	else if(emax_ev < emin_ev)
	{
		ModuleBase::WARNING("ModuleIO::calculate_dos","emax_ev < emin_ev");
		return 0;
	}

	// mohan fixed bug 2010-1-18
	const int npoints = static_cast<int>(std::floor ( ( emax_ev - emin_ev ) / de_ev )) ;
	dos.clear();
	ene.clear();
	if(npoints <= 0)
	{
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"npoints",npoints);
		ModuleBase::WARNING("ModuleIO::calculate_dos","npoints <= 0");
		return 0;
	}
	if(GlobalV::MY_RANK==0)
	{
		ofs << npoints << std::endl;
		ofs << nkstot << std::endl;
	}

	GlobalV::ofs_running << "\n OUTPUT DOS FILE IN: " << fa << std::endl;
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"min state energy (eV)",emin_ev);
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"max state energy (eV)",emax_ev);
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"delta energy interval (eV)",  de_ev);
	
	double *e_mod = new double[npoints];
	ModuleBase::GlobalFunc::ZEROS(e_mod, npoints);

	double sum   = 0.0;
	double e_new = emin_ev;
	double e_old = 0.0;
	while( e_new < emax_ev)
	{
//		GlobalV::ofs_running << " enew=" << e_new << std::endl;
		double count = 0.0;
		e_old = e_new ;
		e_new += de_ev;
		for(int ik=0;ik<nks;ik++)
		{
			if(is == isk[ik])
			{
				for(int ib = 0; ib < nbands; ib++)
				{
					//  compare et and e_old(e_new) in ev unit.
					if( ekb(ik, ib)*ModuleBase::Ry_to_eV >= e_old && ekb(ik, ib)*ModuleBase::Ry_to_eV < e_new)
					{
						// because count is 'double' type,so
						// we can't write count++ or ++count
						count += wk[ik]*nkstot; //mohanix bug 2012-04-23
//						GlobalV::ofs_running << " count = " << count << " wk = " << wk[ik] << " nks = " << nks << std::endl;
					}		
				}
			}
		}
#ifdef __MPI
		Parallel_Reduce::reduce_double_allpool(count);
#endif
		count = count / static_cast<double>(nkstot);
		sum += count;
		if(GlobalV::MY_RANK==0)
		{
			ofs << e_new << " " << count << std::endl;
			dos.push_back(count);
			ene.push_back(e_new);
		}

	}

	//now use Gaussian smearing to smooth the dos and write to DOS_is_smearing
	if(GlobalV::MY_RANK==0)
	{
		dos_smearing.resize(dos.size()-1);

		//double b = INPUT.dos_sigma;
		double b = sqrt(2.0)*bcoeff;
		for(int i=0;i<dos.size()-1;i++)
		{
			double Gauss=0.0;

			for(int j=0;j<dos.size()-1;j++)
			{
				double de = ene[j] - ene[i];
				double de2 = de * de;
				//----------------------------------------------------------
				// EXPLAIN : if en
				//----------------------------------------------------------
				Gauss = exp(-de2/b/b)/sqrt(3.1415926)/b;
				dos_smearing[j] += dos[i]*Gauss;
			}
		}

		//----------------------------------------------------------
		// EXPLAIN : output DOS_smearing.dat
		//----------------------------------------------------------
		double sum2=0.0;
		for(int i=0;i<dos.size()-1;i++)
		{
			sum2 += dos_smearing[i];
			ofs1 <<std::setw(20)<<ene[i]
				<<std::setw(20)<<dos_smearing[i]
				<<std::setw(20)<<sum2<<"\n";
		}
	}


	if(GlobalV::MY_RANK==0)
	{
		ofs.close();
		ofs1.close();
	}
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"number of bands",nbands);
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"sum up the states", sum);
	delete[] e_mod;

	return 1;
}
