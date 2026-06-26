#include "Out_Orbital.h"
#include "tools.h"
#include "../src_parallel/parallel_reduce.h"

Out_Orbital::Out_Orbital(){}

Out_Orbital::~Out_Orbital(){}

void Out_Orbital::write(void)
{
	TITLE(ofs_running,"Out_Orbital","write");

	ofstream ofs;

	if(MY_RANK==0)
	{
		ofs.open("ORBITAL_RESULTS.txt");
		this->version_information( ofs );
		this->INPUTs_information( ofs );
		this->spillage_information( ofs, mz.ilevel );
		this->metropolis_information( ofs);
		this->c4_information( ofs );
	}
	this->mkb_information( ofs);

	if(MY_RANK==0)
	{
		ofs.close();
	}

	return;
}

void Out_Orbital::version_information( ofstream &ofs )
{
    ofs << "\n <VERSION>";
    ofs << "\n AUTHOR : Mohan Chen";
    ofs << "\n StartDate : 2009-4-01";
	ofs << "\n LastModify: 2012-6-27";
    ofs << "\n LOCATION : LQCC, Hefei, China";
    ofs << "\n EMAIL : mohan@mail.ustc.edu.cn";
    ofs << "\n Description : Calculate the coefficients C4 of f(r) in Spherical Bessel Basis J(qr).";
    ofs << "\n Formula : C4 = integral(r)[ f(r)*jl(qr)*r^{2} ]dr ";
    ofs << "\n P.S. : We default consider f(r) read from file is in the form : ( f(r) * r ).";
    ofs << "\n</VERSION>";
}

void Out_Orbital::INPUTs_information( ofstream &ofs )
{
    ofs << "\n\n<INPUTS>";
    ofs << "\n" << setw(20) << ECUT << " Energy cutoff(Hartree.).";
    ofs << "\n" << setw(20) << RCUT << " rcut (a.u.)";
    ofs << "\n" << setw(20) << NE << " eigenvalue number( sqrt(ecut*2)*rcut/PI ).";
    ofs << "\n" << setw(20) << TOLERENCE << " tolerence to calculate eigenvalue.";
	ofs << "\n" << setw(20) << NTYPE << " Number of atom types.";
	for(int it=0; it<NTYPE; it++)
	{
		ofs << "\n" << setw(20) << LABEL[it] << " Atom Label.";
		ofs << "\n" << setw(20) << NA[it] << " Number of atoms.";
	}
	// mohan add 2010-05-02
	if(BANDS_CONTROL)
	{
		ofs << "\n" << setw(20) << BANDS_START+1 << " start band index.";
		ofs << "\n" << setw(20) << BANDS_END << " ended band index.";
	}
    ofs << "\n</INPUTS>";
}

void Out_Orbital::spillage_information( ofstream &ofs, const int &ilevel)
{
	ofs << "\n\n<Spillage>";
	ofs << "\n" << setw(20) << STRNUM << " kinds of structures.";
	
	// (1) mohan add average information 2010-04-11
	ofs << "\nAverage Spillage Value";
	for(int k=0; k<ilevel; k++)
	{
		double average = 0.0;
		for(int i=0;i<STRNUM; i++)
		{
			average+=input.SV.value_each_level(i,k);
		}
		ofs << "\n " << setw(8) << k+1 
			<< setiosflags(ios::scientific)
			<< setprecision(6) 
			<< setw(20) << average/STRNUM;
	}
	
	// (2) SV stands for 'spillage value'
	// value_old containing the newest accepted spillage value.
	for(int i=0; i<STRNUM; i++)
	{
		ofs << "\nStructureIndex " << i+1;
		for(int k=0; k<ilevel; k++)
		{
			ofs << "\n " << setw(8) << k+1 
				<< setiosflags(ios::scientific) 
				<< setprecision(6) 
				<< setw(20) << input.SV.value_each_level(i,k);
		}
	}
	ofs << "\n</Spillage>";
	ofs << resetiosflags(ios::scientific);
	return;
}

void Out_Orbital::metropolis_information( ofstream &ofs)
{
	ofs << "\n\n<Metropolis>";
	ofs << "\n" << setw(20) << mz.metro.get_spi_ini_temp() << " Start temperature (Kelvin) for spillage minimization.";
	ofs << "\n" << setw(20) << mz.metro.get_spi_rate() << " Decreasing rate of temperature.";
	ofs << "\n" << setw(20) << mz.metro.get_spi_ntemp() << " Number of different temperature (for spillage).";
	ofs << "\n" << setw(20) << mz.metro.get_spi_nsteps() << " Number of steps for each temperature (for spillage).";

	ofs << "\n" << setw(20) << mz.metro.get_kin_ini_temp() << " Start temperature (Kelvin) for kinetical energy minimization.";
	ofs << "\n" << setw(20) << mz.metro.get_kin_rate() << " Decreasing rate of temperature.";
	ofs << "\n" << setw(20) << mz.metro.get_kin_ntemp() << " Number of different temperature (for kinetical).";
	ofs << "\n" << setw(20) << mz.metro.get_kin_nsteps() << " Number of steps for each temperature (for kineitcal).";
	ofs << "\n</Metropolis>";
}

void Out_Orbital::c4_information( ofstream &ofs)
{
	ofs << "\n\n<Coefficient>";
	ofs << "\n" << setw(20) << NCHIUSED << " Total number of radial orbitals.";
	
	// use multize parameters: ntype, *lmax_type, l_nchi, ne;
	for(int it=0; it<NTYPE; it++)
	{
		for(int L=0; L<mz.lmax_type[it]+1; L++)
		{
			// last step, nmax is different if L is different
			for(int N=0; N<mz.l_nchi(it,L); N++)
			{
				ofs << "\n" << setw(20) << "Type" 
					<< setw(20) << "L"
					<< setw(20) << "Zeta-Orbital";
				ofs << "\n" << setw(20) << it+1
					<< setw(20) << L
					<< setw(20) << N+1;
				for(int ie=0; ie<mz.Level[0].ne; ie++)
				{
					if(ie%4==0) ofs << endl;
					//mohan modify 2009-09-25
//					ofs << setw(20) << input.Coef.C4_accumulate(it, L, N, ie);
					ofs << setw(25)
						<< setiosflags(ios::showpoint)
						<< setiosflags(ios::scientific)						
					    << setprecision(15) 
						<< input.Coef.C4_old(it, L, N, ie);
				}
			}
		}
	}
	ofs << "\n</Coefficient>";
	return;
}

void Out_Orbital::mkb_information( ofstream &ofs)
{
	ofs_running << "\n mkb_information." << endl;
	if(MY_RANK==0)
	{
		ofs << "\n\n<Mkb>";
		ofs << "\n" << setw(20) << mz.nlevel << " Total number of orbitals optimized levels." ;
	}

	double sum0 = 0.0;
	double *sum1 = new double[NBANDS];
	ZEROS(sum1, NBANDS);

	int bands_start;
   	int bands_end;
   	int bands_number;

   	if (BANDS_CONTROL)
   	{
       	bands_start = BANDS_START;
       	bands_end = BANDS_END;
       	bands_number = bands_end - bands_start;
   	}
   	else
   	{
       	bands_start = 0;
       	bands_end = NBANDS;
       	bands_number = bands_end - bands_start;
   	}

	ofs << "\nBands start from " << bands_start;
	ofs << "\nBands ended at " << bands_end;
	ofs << "\nOptimized bands number " << bands_number;
	ofs << "\nSpillage per band is " << 1.0/bands_number;

	//ofs << setioflags(ios::fixed);
	for(int il=0; il<mz.ilevel; il++)
	{
		if(MY_RANK==0)
		{
			ofs << "\n\nFill Left Hilbert space of each band(average) by LCAO for Level " << il+1;
			ofs << "\n" << setw(5) << "BANDS" << setw(20) << "New Fill" << setw(20) << "Total Fill" << setw(20) << "Left Spillage";
			ofs << setprecision(10);
		}
		double conb = 0.0;
		ofs_running << "\n NBANDS=" << NBANDS << " STRNUM=" << STRNUM << endl;	
	
		for(int ib=bands_start; ib<bands_end; ib++)
		{
			double mm = 0.0;
			// calculate the average fill gution from all k points.
			for(int is=0; is<STRNUM; is++)
			{
			
				for(int ik=0; ik<mz.Level[il].data[is].nks; ik++)
				{
					mm += mz.Level[il].data[is].Mkb(ik,ib) *  mz.Level[il].data[is].weight[ik] 
					/(double)bands_number
					/(double)STRNUM;
				}
#ifdef __MPI
				Parallel_Reduce::reduce_double_allpool(mm);
#endif
			}
			sum1[ib] += mm;
		
			if(MY_RANK==0)
			{
				ofs << "\n" << setw(5) << ib+1 << setw(20) << mm << setw(20) << sum1[ib] << setw(20) << 1.0/(double)bands_number - sum1[ib];
			}

			conb += mm;
		}
		sum0 += conb;	
		ofs << "\nNew   Fill Contribution = " << conb;
		ofs << "\nTotal Fill Contribution = " << sum0;

		ofs << "\nLeft spillage = " << 1.0 - sum0;
	}
	ofs << "\n<Mkb>";
	delete[] sum1;

	return;
}
