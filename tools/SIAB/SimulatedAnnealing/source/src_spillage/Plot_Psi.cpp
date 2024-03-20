#include "Plot_Psi.h"
#include "tools.h"
#include "../src_tools/Gram_Schmidt_Orth.h"

Plot_Psi::Plot_Psi()
{
	test=0;
	rcut = 0.0;
	
	dr = 0.01;

	meshU = 0; // uniform mesh.
	meshL = 0; // Log mesh.
	
	rU = new double[1];
	rL = new double[1];

	xmin = -6;
	zed = 4;//charge on the atom
	dx = 1.0e-2;
	xmax = 10; 
}

Plot_Psi::~Plot_Psi()
{
	delete[] rU;
	delete[] rL;
}

void Plot_Psi::allocate(void)
{	
	// (1) read in parameters.
    if(SCAN_BEGIN(input.ifs, "<PLOT>"))
    {
		READ_VALUE(input.ifs, this->dr);		
		assert( this->dr > 0.0);

		// mohan add 2009-11-23
		READ_VALUE(input.ifs, this->xmin);
		READ_VALUE(input.ifs, this->zed);//charge?
		READ_VALUE(input.ifs, this->dx);
		READ_VALUE(input.ifs, this->xmax);
		assert(zed>0.0);
		assert(dx>0.0);
		assert(xmax>=0.0);
        
		SCAN_END(input.ifs, "</PLOT>");
	}
	
	this->ne = NE;
	this->meshU = static_cast<int>(RCUT/dr);
	this->meshL = ( std::log( xmax * zed ) - xmin ) / dx;
	// mesh should be odd.
	if(meshU%2==0) meshU++;
	if(meshL%2==0) meshL++;
	
	/*
	cout << "\n dr = " << dr;
	cout << "\n meshU = " << meshU;
	cout << "\n meshL = " << meshL;
	
	cout << "\n xmin = " << xmin;
	cout << "\n zed = " << zed;
	cout << "\n dx = " << dx;
	cout << "\n xmax = " << xmax;
	*/
	
	// other parameters.
	this->rcut = RCUT;

	delete[] rU;
	delete[] rL;
	rU = new double[meshU];	
	rL = new double[meshL];
	
	for(int ir=0; ir<meshU; ir++)
	{
		rU[ir] = dr * ir;
	}

	for(int ir=0; ir<meshL; ir++)
	{
		rL[ir] = std::exp(xmin + dx * ir)/zed;
	}

	this->lmax = LMAXUSED;
	this->total_nchi = NCHIUSED;
	this->psiU.create( this->total_nchi, this->meshU);
	this->psiL.create( this->total_nchi, this->meshL);

	return;
}

// be called in main.cpp
void Plot_Psi::radial_wave_function(void)
{
	TITLE("Plot_Psi","radial_wave_function");
	timer::tick("Plot_Psi","radial_wf");
	
	this->allocate();

	// calculate all psi
	int ichi = 0;
	for(int it=0; it<NTYPE; it++)
	{
		if(test==2) cout << "\n mz.lmax_type[" << it << "]=" << mz.lmax_type[it];
		for(int L=0; L<mz.lmax_type[it]+1; L++)
		{
			double *eigen1 = new double[ne];
			for(int ie=0; ie<ne; ie++)
			{
				eigen1[ie] = mz.metro.Psi2.eigenvalue(L, ie);
				//cout << "\n eigen1[" << ie << "]=" << eigen1[ie]; 
			}
			
			for(int N=0; N<mz.l_nchi(it,L); N++)
			{
				if(test==2)cout << "\n T=" << it << " L=" << L << " N=" << N; 
				double *c4 = new double[ne];

				for(int ie=0; ie<ne; ie++)
				{
					c4[ie] = input.Coef.C4_old(it, L, N, ie);
				}

				// call the function to renormalize the wave functions.
				mz.metro.Psi2.norm_c4(c4, L);
				
				// calculate kinetical energy
				mz.metro.Psi2.ofsk << "\n" << "Type="<< it+1 << " L=" << L << " N=" << N+1;
				mz.metro.Psi2.ofsk << "\n" 
				<< setw(5) << "nq"
				<< setw(20) << "Coefficient" 
				<< setw(20) << "Eigen(Ry)" 
				<< setw(20) << "Kinetic_Energy";
				double kin_sum = 0.0;
				for(int ie=0; ie<ne; ie++)
				{
					//double kin = abs(c4[ie])*pow((ie+1)*PI/input.inputs.get_rcut(),2);//In Rydberg Unit.
					//double kin = abs(c4[ie])*pow(eigen1[ie]*PI/input.inputs.get_rcut(),2);//In Rydberg Unit.
					double kin = pow(c4[ie]*eigen1[ie],2)
					*mz.metro.Psi2.jjnorm(L,ie);
					kin_sum += kin;
					mz.metro.Psi2.ofsk << "\n" 
					<< setw(5) << ie
					<< setiosflags(ios::fixed)
					<< setprecision(10)
					<< setiosflags(ios::showpoint)
					//<< setiosflags(ios::scientific)
					<< setw(20) << c4[ie] 
					<< setw(20) << eigen1[ie]*eigen1[ie] 
					<< setw(20) << kin << " (Rydberg)"
					<< setprecision(6);
					
					//>>>>>>>>>>
					//The two formula are exactly same.
					//<< setw(15) << abs(c4[ie])*pow(eigen1[ie],2) 
					//<< setw(15) << abs(c4[ie])*pow((ie+1)*PI/input.inputs.get_rcut(),2);
					//<<<<<<<<<<
				}
				mz.metro.Psi2.ofsk << "\n Total kinetic energy = " << kin_sum;

				double *psi1 = new double[meshU];
				double *psi2 = new double[meshL];
				this->get_psi( meshU, rU, L, psi1, c4, eigen1, 0 ); 
				this->get_psi( meshL, rL, L, psi2, c4, eigen1, 1 ); 

				for(int ir=0; ir<meshU; ir++)
				{
					psiU(ichi, ir) = psi1[ir];
				}
				for(int ir=0; ir<meshL; ir++)
				{
					psiL(ichi, ir) = psi2[ir];
				}
				delete[] psi2;
				delete[] psi1;
				delete[] c4;
				++ichi;
			}
			delete[] eigen1;
		}
	}

	// Peize Lin add 2015-11-20
	this->orthogonalization(meshU, rU, psiU);	
	this->orthogonalization(meshL, rL, psiL);
	
	this->print_orbitalsU_used();
	this->print_orbitalsL_used();

	// mohan add 2010-04-13
	// use uniform grid to establish the energy cutoff of each orbital
	ofstream ofs("ORBITAL_ECUT.txt");
	ofs << setprecision(12);
	this->dk = 0.01;
	this->meshK = 1500;// (1500*0.01)^2 * 2 = 450 Rydberg.
	if(meshK%2==0)meshK++;
	this->psiUK.create(total_nchi,meshK);
	double *psi = new double[meshU];
	ichi = 0;
	for(int it=0; it<NTYPE; it++)
	{
		for(int L=0; L<mz.lmax_type[it]+1; L++)
		{
			for(int N=0; N<mz.l_nchi(it,L); N++)
			{
				ofs << "WaveFunctionIndex=" << ichi << endl;
				ofs << "T=" << it << " L=" << L << " N=" << N << endl;
				for(int ir=0; ir<meshU; ir++)
				{
					psi[ir] = this->psiU(ichi, ir);
				}
				this->establish_ecut(ofs, ichi, psi, L);
				++ichi;
			}
		}
	}
	assert(ichi == total_nchi);
	delete[] psi;
	ofs.close();
	
	// print the orbitals which can be viewed imediately.
	this->print_orbitals_plot();
	timer::tick("Plot_Psi","radial_wf");
	return;
}

void Plot_Psi::print_orbitals_plot(void)const
{	
	ofstream f1( "ORBITAL_PLOTU.dat");
	ofstream f2( "ORBITAL_PLOTL.dat");
	ofstream f3( "ORBITAL_PLOTUK.dat");

	f1 << setprecision(12)
	<< setiosflags(ios::fixed);
	f2 << setprecision(12)
	<< setiosflags(ios::fixed);
	f3 << setprecision(12)
	<< setiosflags(ios::fixed);
	for(int ir=0; ir<meshU; ir++)
	{
		f1 << ir*dr;
		for(int i=0; i<total_nchi; i++)
		{
			f1 << " " << psiU(i, ir);
		}
		f1 << "\n";
	}
	
	for(int ir=0; ir<meshL; ir++)
	{
		f2 << std::exp(xmin + dx * ir)/zed;
		for(int i=0; i<total_nchi; i++)
		{
			f2 << " " << psiL(i, ir);
		}
		f2 << "\n";
	}

	for(int ik=0; ik<meshK; ik++)
	{
		f3 << pow(ik*dk,2);
		for(int i=0; i<total_nchi; i++)
		{
			f3 << " " << psiUK(i, ik);
		}
		f3 << "\n";
	}
	
	f1.close();
	f2.close();
	f3.close();
	return;
}

void Plot_Psi::print_orbitalsU_used(void)const
{
	//====================================
	// (5)
	//====================================
	// output orbital for standard use form
	int ichi = 0;
	// because we can't read in r = 0 term now.
	// but liaochen suggest we should know psi(r)
	// instead of psi(r)*r at zero point,
	// I think it's correct.
	int start = 0;
	for(int it=0; it<NTYPE; it++)
	{
		// output the wave function in uniform mesh.
		stringstream ss;
		ss << "ORBITAL_" << mz.Level[0].info[it].id << "U.dat";

		ofstream f1( ss.str().c_str() );          // pengfei 2014-10-13
		f1 << "---------------------------------------------------------------------------"<<endl;
		f1 << "Element" << setw(22) << LABEL[it] <<endl;
		f1 << "Energy Cutoff(Ry)" << setw(13) << ECUT << endl;
		f1 << "Radius Cutoff(a.u.)" << setw(10) << RCUT << endl;
		f1 << "Lmax " << setw(24)<<mz.lmax_type[it]<< endl;
		for(int L=0; L<mz.lmax_type[it]+1; L++)
		{
			string orbital_type;
			switch( L )
			{
				  case 0: orbital_type = "S"; break;
				  case 1: orbital_type = "P"; break;
				  case 2: orbital_type = "D"; break;
				  case 3: orbital_type = "F"; break;
				  case 4: orbital_type = "G"; break;
			}

				f1 << "Number of "<<orbital_type << "orbital--> " << setw(7) << mz.l_nchi(it,L) << endl;
		}
		f1 << "---------------------------------------------------------------------------"<<endl;
		f1 << "SUMMARY" <<"  "<<"END"<<endl;
		f1 << endl;
		f1 << "Mesh " << setw(26) << meshU-start << endl;//mohan fix bug 2010-04-15
		f1 << "dr " << setw(29) << dr ;
		
		for(int L=0; L<mz.lmax_type[it]+1; L++)
		{
			double *eigen1 = new double[ne];
			for(int ie=0; ie<ne; ie++)
			{
				eigen1[ie] = mz.metro.Psi2.eigenvalue(L, ie);
			}
			for(int N=0; N<mz.l_nchi(it,L); N++)
			{
				f1 << "\n" << setw(20) << "Type" << setw(20) << "L" << setw(20) << "N";
				f1 << "\n" << setw(20) << it << setw(20) << L << setw(20) << N;
				f1 << setiosflags(ios::scientific) << setprecision(12);
				for(int ir=start; ir<meshU; ir++)
				{
					if(ir%4==0) f1 << endl;
					f1 << " " << this->psiU(ichi, ir);
				}
				++ichi;
			}
		}
		f1.close();
	}
	return;
}

void Plot_Psi::print_orbitalsL_used(void)const
{
	//====================================
	// (5)
	//====================================
	// output orbital for standard use form
	int ichi = 0;
	// because we can't read in r = 0 term now.
	for(int it=0; it<NTYPE; it++)
	{
		// output the wave function in uniform mesh.
		stringstream ss;
		ss << "ORBITAL_" << mz.Level[0].info[it].id << "L.dat";

		ofstream f1( ss.str().c_str() );

		f1 << "Mesh " << meshL << endl;
		f1 << "zed " << zed << endl;
		f1 << "dx " << dx << endl;
		f1 << "xmin " << xmin << endl;
		f1 << "xmax " << xmax;
		
		for(int L=0; L<mz.lmax_type[it]+1; L++)
		{
			double *eigen1 = new double[ne];
			for(int ie=0; ie<ne; ie++)
			{
				eigen1[ie] = mz.metro.Psi2.eigenvalue(L, ie);
			}
			for(int N=0; N<mz.l_nchi(it,L); N++)
			{
				f1 << "\n" << setw(20) << "Type" << setw(20) << "L" << setw(20) << "N";
				f1 << "\n" << setw(20) << it << setw(20) << L << setw(20) << N << endl;
				f1 << setiosflags(ios::scientific) << setprecision(12);
				for(int ir=0; ir<meshL; ir++)
				{
					if(ir%4==0) f1 << endl;
					f1 << " " << this->psiL(ichi, ir);
				}
				++ichi;
			}
		}
		f1.close();
	}
	return;
}

void Plot_Psi::get_psi( 
	const int &mesh,
	const double *r,
	const int &l, 
	double *psi, 
	double *c4, 
	const double *eigen1, bool logmesh)const
{
	// (1) get psi1
	double *jle = new double[mesh];
	double *g = new double[mesh]; // smooth function

	if(SMOOTH)
	{
		double sigma2 = SIGMA*SIGMA;
		for(int ir=0; ir<mesh; ir++)
		{
			g[ir] = 1.0 - std::exp(-( (r[ir]-this->rcut)*(r[ir]-this->rcut)/2.0/sigma2) );
			//if(ir==mesh-1) cout << "\n g[" << ir << "]=" << g[ir] << endl;
		}
	}

	ZEROS(psi,mesh);
	for(int ie=0; ie<ne; ie++)
	{
		ZEROS(jle, mesh);
		Mathzone::Spherical_Bessel( mesh, r, eigen1[ie], l, jle);

		// (2) normal part
		// mohan modify 2009-08-28
		// if the input file <Q|Jlq> and <Jlq|Jlq>
		// using smooth Jlq, here we output the 
		// smooth functions.
		if(SMOOTH)
		{
			for(int ir=0; ir<mesh; ir++)
			{
				jle[ir] *= g[ir];
			}
		}

		for(int ir=0; ir<mesh; ir++)
		{
			psi[ir] += c4[ie] * jle[ir];
		}
	}
	delete[] g;
	delete[] jle;

	// (2)normalization
	double norm;
	double *function = new double[mesh];
	double *rab = new double[mesh];

	// mohan fix bug 2010-04-13
	if(logmesh)
	{
		rab[0]=0.0;
		for(int ir=1; ir<mesh; ir++)
		{
			rab[ir] = r[ir] - r[ir-1];
		}
	}
	else
	{
		for(int ir=0; ir<mesh; ir++)
		{
			rab[ir] = this->dr;
		}
	}

	for(int ir=0; ir<mesh; ir++)
	{
		function[ir] = psi[ir]*psi[ir]*r[ir]*r[ir];
	}

	Mathzone::Simpson_Integral(mesh, function, rab, norm );
//	cout << "\n norm = " << norm;
	if( abs(norm) < 1.0e-5)
	{
		cout << "\n WARNING!" << " norm = " << norm; 
	}
	norm = sqrt(norm);
	for(int ir=0; ir<mesh; ir++)
	{
		psi[ir] /= norm;
	}

	for(int ie=0; ie<ne; ie++)
	{
		c4[ie] /= norm;
	}

	delete[] function;
	delete[] rab;
	return;
}

void Plot_Psi::establish_ecut(ofstream &ofs, const int &ichi, double *psi, const int &L)
{
	//(1) init the parameters
	double pi = 3.1415926535;
	double fpi = 4*pi;
	
	for(int ir=0; ir<this->meshU; ir++)
	{
		psi[ir] /= sqrt(fpi);
	}

	// (2) check the norm of psi
	double norm = 0.0;
	double* function = new double[meshU];
	double* rab = new double[meshU];
	for(int ir=0; ir<meshU; ir++)
	{
		rab[ir] = dr;
		function[ir] = std::pow(psi[ir]*this->rU[ir],2);	
	}
	Mathzone::Simpson_Integral(meshU, function, rab, norm );
	ofs << "Psi(r) norm=" << norm * fpi << endl;

	// (3) get psik
	double *psik = new double[meshK];
	double *jj = new double[meshU];
	for(int ik=0; ik<meshK; ik++)
	{
		double q = ik*this->dk;
		Mathzone::Spherical_Bessel( this->meshU, this->rU, q, L, jj);
		for(int ir=0; ir<this->meshU; ir++)
		{
			function[ir] = psi[ir]*this->rU[ir]*this->rU[ir]*jj[ir];	
		}
		Mathzone::Simpson_Integral(this->meshU, function, rab, psik[ik] );	
	}
	double prefac_k = fpi / pow(2*pi,1.5);
	for(int ik=0; ik<meshK; ik++)
	{
		psik[ik] *= prefac_k;
		this->psiUK(ichi,ik) = psik[ik];
	}

	// (4) prepare for the integral of psik.
	double normk = 0.0;
	double* rabk = new double[meshK];
	double* functionk = new double[meshK];
	for(int ik=0; ik<meshK; ik++)
	{
		rabk[ik] = dk;
	}
	for(int ik=0; ik<meshK; ik++)
	{
		functionk[ik] = std::pow(psik[ik]*ik*this->dk,2);
	}
	Mathzone::Simpson_Integral(this->meshK, functionk, rabk, normk );
	ofs << "Psi(k) norm=" << normk * fpi << endl;
	
	// (5) Find out the energy cut for each orbital.
	normk = 0.0;
	int kmesh_used = 1;

	//===================================================
	// 0.999 is the accepted value of norm, according
	// to the value we establish the energy cutoff.
	//===================================================
	while(normk*fpi< 0.99 && kmesh_used < this->meshK)
	{
		kmesh_used += 2;
		Mathzone::Simpson_Integral(kmesh_used, functionk, rabk, normk );
	}
//	ofs << "Psi(k) norm(0.99)=" << normk * fpi << endl;
	ofs << "Ecut(Ry)(norm>0.9900)=" << pow(dk*kmesh_used,2)*2 << endl;

	while(normk*fpi< 0.999 && kmesh_used < this->meshK)
	{
		kmesh_used += 2;
		Mathzone::Simpson_Integral(kmesh_used, functionk, rabk, normk );
	}
//	ofs << "Psi(k) norm(0.999)=" << normk * fpi << endl;
	ofs << "Ecut(Ry)(norm>0.9990)=" << pow(dk*kmesh_used,2)*2 << endl;

	while(normk*fpi< 0.9999 && kmesh_used < this->meshK)
	{
		kmesh_used += 2;
		Mathzone::Simpson_Integral(kmesh_used, functionk, rabk, normk );
	}
//	ofs << "Psi(k) norm(0.999)=" << normk * fpi << endl;
	ofs << "Ecut(Ry)(norm>0.9999)=" << pow(dk*kmesh_used,2)*2 << endl;

	ofs << endl;

	delete[] function;
	delete[] functionk;
	delete[] rabk;
	delete[] rab;
	delete[] psik;
	delete[] jj;
	return;
}

// Peize Lin add 2015-11-20
void Plot_Psi::orthogonalization( 
	const int &mesh,
	double *r,
	const matrix &psi
)
{
	Gram_Schmidt_Orth<double,double> gs_orth;
	gs_orth.set_coordinate(1);
	gs_orth.set_variable_num(mesh);
	gs_orth.set_rab_element(this->dr);
	gs_orth.set_r(r);

	double * tmp_index(psi.c);	

	for(int it=0; it<NTYPE; it++)
	{
		for(int L=0; L<mz.lmax_type[it]+1; L++)
		{
			gs_orth.set_func_num(mz.l_nchi(it,L));
			double **psi_TL = new double*[mz.l_nchi(it,L)];
			for(int N=0; N<mz.l_nchi(it,L); N++)
			{
				psi_TL[N] = tmp_index;
				tmp_index += mesh;
			}
			gs_orth.set_func(psi_TL);
			gs_orth.orth();
			delete[]psi_TL;
		}
	}

	return;
}
