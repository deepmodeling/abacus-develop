//==========================================================
// AUTHOR : mohan
// DATE : 2008-11-5
//==========================================================
#include "../src_pw/global.h"
#include "eximport.h"

eximport::eximport()
{
}

eximport::~eximport()
{
}

//==========================================================
// FUNCTION NAME : write_data
// DO : write the information to be compared
//==========================================================
void eximport::write_data
(
	const string &fn,
	const string &type
)
{
	TITLE("eximport","write_data");
	ofstream ofs(fn.c_str());

	this->out_input(ofs);
	this->out_wannier(ofs);
	this->out_unitcell(ofs);
	this->out_kpoints(ofs);

	if (type == "all")
	{
//==========================================================
// EXPLAIN : because wave function is memory consuming,
//           so I don't output it.
//==========================================================
//		this->out_evc(ofs);
		this->out_band(ofs);
		this->out_charge(ofs);
		this->out_energy(ofs);
	}

//==========================================================
// EXPLAIN : output wave function information
//==========================================================
	if (type == "evc")
	{
		this->out_evc(ofs);
	}

//==========================================================
// EXPLAIN : output band information 
//==========================================================
	if (type == "band")
	{
		this->out_band(ofs);
	}

//==========================================================
// EXPLAIN : output energy information 
//==========================================================
	if (type == "energy")
	{
		this->out_energy(ofs);
	}

//==========================================================
// EXPLAIN : output charge density information 
//==========================================================
	if (type == "charge")
	{
		this->out_charge(ofs);
	}

	ofs.close();
	return;
}

void eximport::read_data(const string &fn)
{
	TITLE("eximport","read_data");

	ifstream ifs(fn.c_str());
	if(!ifs)
	{
		ofs_warning << " File name : " << fn << endl;
		WARNING_QUIT("eximport::read_data","Can not find file.");
	}

	this->in_input(ifs);
	this->in_wannier(ifs);
	this->in_unitcell(ifs);
	this->in_kpoints(ifs);

	if (fn == "all")
	{
		this->in_evc(ifs);
		this->in_band(ifs);
		this->in_charge(ifs);
		this->in_energy(ifs);
	}
	return;
}

void eximport::print_data(const string &fn) const
{
	TITLE("eximport","print_data");
	ofstream ofs( fn.c_str() );

	ofs << setw(20) << "basis" << setw(20) << this->basis << endl;
	ofs << setw(20) << "latname" << setw(20) << this->latname << endl;
	ofs << setw(20) << "calculation" << setw(20) << this->calculation << endl;
	ofs << setw(20) << "ecutwfc" << setw(20) << this->ecutwfc << endl;
	ofs << setw(20) << "nband" << setw(20) << this->nband << endl;
	ofs << setw(20) << "tr2" << setw(20) << this->tr2 << endl;
	ofs << setw(20) << "nx" << setw(20) << this->nx << endl;
	ofs << setw(20) << "ny" << setw(20) << this->ny << endl;
	ofs << setw(20) << "nz" << setw(20) << this->nz << endl;
	ofs << setw(20) << "nxyz" << setw(20) << this->nxyz << endl;
	ofs << setw(20) << "startingpot" << setw(20) << this->startingpot << endl;
	ofs << setw(20) << "mixing_beta" << setw(20) << this->Mixing_beta << endl;

	ofs << setw(20) << "nks" << setw(20) << this->nks << endl;

	for (int ik = 0;ik < kv.nks;ik++)
	{ 
		ofs << setw(20) << "ngk[" << ik << "]" << setw(20) << ngk[ik] << endl; 
	}

	ofs << setw(20) << "qtot" << setw(20) << this->qtot << endl;
	ofs << setw(20) << "lat0" << setw(20) << this->lat0 << endl;
	ofs << setw(20) << "ntype" << setw(20) << this->ntype << endl;
	ofs << setw(20) << "band_energy(k-point,band)" << endl;

	for (int ik = 0;ik < this->nks;ik++)
	{
		for (int ib = 0;ib < this->nband;ib++)
		{
			ofs << setw(10) << this->band_energy[ik][ib];
		}
		ofs << endl;
	}

	ofs << setw(20) << "Omega" << setw(20) << this->omega << endl;
	ofs << setw(20) << "rho_nc" << setw(20) << this->rho_nc << endl;
	ofs << setw(20) << "iter" << setw(20) << this->iter << endl;
	ofs << setw(20) << "etot" << setw(20) << this->etot << endl;
	ofs << setw(20) << "eband" << setw(20) << this->eband << endl;
	ofs << setw(20) << "one-electron" << setw(20) << this->one_electron << endl;
	ofs << setw(20) << "hartree" << setw(20) << this->hartree << endl;
	ofs << setw(20) << "exchange-corr" << setw(20) << this->xc << endl;
	ofs << setw(20) << "ewald" << setw(20) << this->ewald << endl;
	ofs << setw(20) << "natomwfc" << setw(20) << this->natomwfc << endl;

	ofs.close();
}

void eximport::fir_wf(ComplexMatrix *psi, const int npsi, const string &fn)
{
	if(MY_RANK!=0)
	{
		cout<<"\n not the MY_RANK processor , return.";
		return;
	}

	//out put wannier function in q space, this is the file name
	//modify 2007-01-24 
	//add nwan , nwan = number of bands
	//PS: if modify here, you must modify the same name in sec_wf( 2 definition)
	ofstream ofs(fn.c_str());

    out_input(ofs);
    out_kpoints(ofs);
    out_igk(ofs);
    out_planewave(ofs);

    ofs << "\n<WAVEFUNC>";
    ofs << "\n" << npsi << " Number of wave functions." << endl;
    ofs << setprecision(6);

    for(int i=0; i<npsi; i++)
    {
        for(int ik=0; ik<kv.nks; ik++)
        {
            ofs << "\n" << ik;
            for(int ig=0; ig<kv.ngk[ik]; ig++)
            {
                if(ig%4==0) ofs << "\n";
                ofs << setw(15) << psi[ik](i, ig).real()
                    << setw(15) << psi[ik](i, ig).imag();
            }
            ofs << "\n";
        }
    }
    ofs << "\n<WAVEFUNC>";

	return;
}

void eximport::out_gspace_wan(const ComplexMatrix *psi,const int iw,const string &file_name)
{
	//cout<<"\n ==> ei.wanG_for_fit";
	ofstream out_gwan(file_name.c_str());
	int qtot = 0;
	for(int ik=0; ik<kv.nks; ik++)
	{
		for(int ig=0; ig<kv.ngk[ik]; ig++)
		{
			qtot += kv.ngk[ik];
		}
	}
	int nks = kv.nks;
	double factor = TWO_PI/ucell.lat0;
	out_gwan << qtot << endl;
	for(int ik=0;ik<nks;ik++)
	{
		//output wannier functions in G space.
		for(int ig=0;ig<kv.ngk[ik];ig++)
		{
			double g1 = pw.get_GPlusK_cartesian_projection(ik, wf.igk(ik, ig), 0);
			double g2 = pw.get_GPlusK_cartesian_projection(ik, wf.igk(ik, ig), 1);
			double g3 = pw.get_GPlusK_cartesian_projection(ik, wf.igk(ik, ig), 2);
			out_gwan 
			<< setw(15) << g1*factor 
			<< setw(15) << g2*factor
			<< setw(15) << g3*factor
			<< setw(20) << psi[ik](iw,ig).real() 
			<< setw(20) << psi[ik](iw,ig).imag() << endl;
		}
	}

	out_gwan.close();
	return;
}

bool eximport::sec_wf(ComplexMatrix *psi, const int wf_num, const string &fn)const
{
	cout<<"\n ==> ei.sec_wf()_1";
	ifstream sec(fn.c_str());

	if (sec)
	{
		cout << "Find the file :"<<endl; 
		cout << fn << endl;
	}

	if (!sec)
	{
		cout << "Can't open file : " << fn << "\treturn!" << endl;
		return 0;
	}

	for (int ik = 0;ik != kv.nks;ik++)
	{
		for (int iw = 0;iw != wf_num;iw++)
		{
			for (int ig = 0;ig != kv.ngk[ik];ig++)
			{
				// Peize Lin fix bug about rvalue 2016-08-02
				double tmp_real, tmp_imag;
				sec >> tmp_real >> tmp_imag ;
				psi[ik](iw, ig) = complex<double>(tmp_real,tmp_imag);
			}
		}
	}

	sec.close();
	return 1;
}


bool eximport::sec_wf(complex < double> ***psi, const int npsi, const string &fn)const
{
	cout<<"\n ==> ei.sec_wf()_2";
	ifstream sec(fn.c_str());

	if (sec)
	{
		cout << "\n    Find the file :"<<endl; 
		cout << fn << endl;
	}

	if (!sec)
	{
		cout << "Can't open file : " << fn << "\treturn!" << endl;
		return 0;
	}

	for (int ik = 0;ik < kv.nks;ik++)
	{
		for (int iw = 0;iw < npsi;iw++)
		{
			for (int ig = 0;ig < kv.ngk[ik];ig++)
			{
				// Peize Lin fix bug about rvalue 2016-08-02
				double tmp_real, tmp_imag;
				sec >> tmp_real >> tmp_imag ;
				psi[iw][ik][ig] = complex<double>(tmp_real,tmp_imag);				
			}
		}
	}

	//big problem!!!!!!!!
	/*
	for(int iw=0;iw<npsi;iw++)
	{
		mt.normalization(psi[iw]);
	}
	*/
	sec.close();
	return 1;
}

//******************
// wannier function
//******************
void eximport::out_wannier(ofstream &out_data)
{
	//cout<<"\n ==> out_wannier"<<endl;
	out_data << setw(20) << "WANNIER" << endl; //0
	//out_data << setw(20) << LOCAL_BASIS << endl;//1 xiaohui modify 2013-09-02
	out_data << setw(20) << BASIS_TYPE << endl; //xiaohui add 2013-09-02
	return;
}

void eximport::in_wannier(ifstream &in)
{
	cout<<"=== in_wannier ==="<<endl;
	in >> name;
	//cout<<name<<endl;		//8.0
	if (name != "WANNIER")
	{
		cout << "== wrong in in_wannier ===" << endl;
		exit(0);
	}
	in >> this->basis;	//8.1
	//cout<<"basis="<<this->basis<<endl;
	return;
}

//*************************************
// atom_unitcell
//*************************************
void eximport::out_unitcell(ofstream &out_data)
{
	//cout << "\n ==> out_unitcell" << endl;
	out_data << setw(20) << "UNITCELL" << endl;		//2.0
	out_data << setw(20) << ucell.lat0 << endl;        //2.1

	out_data << setw(20) << ucell.latvec.e11 
			 << setw(20) << ucell.latvec.e12 
			 << setw(20) << ucell.latvec.e13 << endl;	//2.2
	out_data << setw(20) << ucell.latvec.e21 
			 << setw(20) << ucell.latvec.e22 
			 << setw(20) << ucell.latvec.e23 << endl;
	out_data << setw(20) << ucell.latvec.e31 
			 << setw(20) << ucell.latvec.e32 
			 << setw(20) << ucell.latvec.e33 << endl;

	out_data << setw(20) << ucell.ntype << endl;		//2.3

	for (int i = 0;i < ucell.ntype;i++)
	{
		out_data << setw(20) << ucell.atoms[i].na ;     //2.4
	}
	out_data << endl;
	return;
}

void eximport::in_unitcell(ifstream &in)
{
	cout << "=== in_unitcell ===" << endl;
	in >> name;
	if (name != "UNITCELL")
	{
		cout << "== wrong in in_unitcell ===" << endl;
		exit(0);
	}

	int i;

	in >> this->lat0;//2.1
	latvec = new double*[3];

	for (i = 0;i < 3;i++)
	{
		latvec[i] = new double[3];
	}

	for (i = 0;i < 3;i++)
	{
		for (int j = 0;j < 3;j++)
		{
			in >> this->latvec[i][j];//2.2
		}
	}

	in >> this->ntype;//2.3

	this->na = new int[ntype];

	for (i = 0;i < ntype;i++)
	{
		in >> this->na[i];//2.4
	}
}

//***************
// out_kpoints
//***************
void eximport::out_kpoints(ofstream &out_data)
{
	//cout << "\n ==> out_k-points" << endl;
	out_data << "\n<KPOINT>";
	out_data << "\n" << kv.nkstot << " Number of total k points";      //3.1
	int sumq = 0;

	for(int ik=0; ik<kv.nks; ik++)
	{
		if(ik%10==0) out_data<<endl;
		out_data << setw(10) << kv.ngk[ik];
		sumq += kv.ngk[ik];
	}

	for (int ik = 0;ik < kv.nks;ik++)
	{
		if(ik%3==0) out_data << "\n";
		out_data << setw(10) << kv.kvec_c[ik].x
		<< setw(10) << kv.kvec_c[ik].y
		<< setw(10) << kv.kvec_c[ik].z;//3.3
	}

	out_data << "\n" << sumq << " Total Number of K+G points.";//3.4
	out_data << "\n<KPOINT>" << endl;
	return;
}

void eximport::out_planewave(ofstream &out_data)
{
	//cout << "\n ==> out_planewave" << endl;
	out_data << "\n<PLANEWAVE>";
	out_data << "\n" << ucell.lat0 << " Lattice constant";
	out_data << "\n" << pw.ngmc_g << " Number of plane waves."<<endl;
	for(int i=0; i<pw.ngmc_g; i++)
	{
		if(i%4==0) out_data<<"\n";
		out_data << setw(8) << pw.get_G_cartesian_projection(i, 0)
				 << setw(8) << pw.get_G_cartesian_projection(i, 1)
				 << setw(8) << pw.get_G_cartesian_projection(i, 2);
	}
	out_data << "\n<PLANEWAVE>";
	return;
}

void eximport::out_igk(ofstream &out_data)
{
	//cout << "\n ==> out_igk" << endl;
	out_data << "\n<KG_INDEX>";

	for(int ik=0; ik<kv.nks; ik++)
	{
		for(int ig=0; ig<kv.ngk[ik]; ig++)
		{
			if(ig%10==0) out_data<<"\n";
			out_data << setw(10) << wf.igk(ik, ig);
		}
	}
	out_data << "\n<KG_INDEX>";
	return;
}

void eximport::in_kpoints(ifstream &in)
{
	cout << "=== in_kpoints ===" << endl;
	in >> name;
	if (name != "KPOINT")
	{
		cout << "== wrong in in_input ===" << endl;
		exit(0);
	}
	int i;
	in >> this->nks;//3.1
	this->ngk = new int[nks];
	kvector = new double*[nks];

	for (i = 0;i < nks;i++)
	{
		kvector[i] = new double[3];
	}

	for (i = 0;i < nks;i++)
	{
		in >> ngk[i];//3.2
		in >> kvector[i][0] >> kvector[i][1] >> kvector[i][2];//3.3
	}

	in >> this->qtot;//3.4
}

//********************************
// input parameters
//*********************************
void eximport::out_input(ofstream &out_data)
{
	//cout << "\n ==> out_input" << endl;
	out_data << "<HEADER>";			//1.0
	out_data << "\n" << ucell.latName << " Lattice Name";//1.1

/*
	if (SCF)
	{
		out_data << setw(20) << "scf" << endl;       	//1.2
	}
	else
	{
		out_data << setw(20) << "nscf" << endl;
	}
*/

	out_data << "\n" << pw.ecutwfc << " Energy Cutoff for wave functions.";//1.3

	out_data << "\n" << ucell.nat << " Number of atoms.";

	out_data << "\n" << ucell.ntype << " Types of atoms.";

//	out_data << "\n" << NBANDS << " Number of wave functions.";//1.4
	//out_data << "\n" << LOCAL_BASIS << " 1 for Local Basis"; xiaohui modify 2013-09-02
	out_data << "\n" << BASIS_TYPE << " Basis type"; //xiaohui add 2013-09-02

	//out_data << "\n" << LINEAR_SCALING << " 1 for Linear Scaling"; xiaohui modify 2013-09-02
	out_data << "\n" << KS_SOLVER << " Diago type"; //xiaohui add 2013-09-02
	out_data << "\n" << SPARSE_MATRIX << " 1 for Sparse Matrix";

/*
	out_data << setw(20) << tr2 << endl;                 //1.5
	out_data << setw(20) << pw.nx 
			 << setw(20) << pw.ny 
			 << setw(20) << pw.nz 
			 << setw(20) << pw.nxyz << endl;//1.6
	out_data << setw(20) << pot.startingpot << endl;//1.7
	out_data << setw(20) << CHR.mixing_beta << endl;//1.8
*/
	out_data << "\n<HEADER>"<<endl;			//1.0

	return;
}

void eximport::in_input(ifstream &in)
{
	cout << " in_input ===" << endl;
	in >> name;

	if (name != "INPUT")
	{
		cout << "== wrong in in_input ===" << endl;
		exit(0);
	}

	in >> this->latname;//1.1

	in >> this->calculation;//1.2
	double ecut_in = 0.0;
	in >> ecut_in;//1.3
	if(ecut_in != pw.ecutwfc)
	{
		cout<<"Charge don't match!"<<endl;
		exit(0);
	}
	this->ecutwfc = ecut_in;

	in >> this->nband;//1.4
	in >> this->tr2;//1.5
	in >> this->nx >> this->ny >> this->nz >> this->nxyz;//1.6
	in >> this->startingpot;//1.7
	in >> this->Mixing_beta;//1.8
}

//**********
// out_band
//**********
void eximport::out_band(ofstream &out_data)
{
	//cout << "\n ==> out_band" << endl;
	out_data << setw(20) << "BAND" << endl;//6.0
	for (int ik = 0; ik < kv.nks; ik++)
	{
		for (int ib = 0; ib < NBANDS; ib++)
		{
			out_data << setw(10) << setprecision(6) << wf.ekb[ik][ib]*Ry_to_eV;//6.1
		}
		out_data << endl;
	}
	out_data << endl;
	return;
}

void eximport::in_band(ifstream &in)
{
	cout << "=== in_band ===" << endl;
	in >> name;
	if (name != "BAND")
	{
		cout << "== wrong in in_band ===" << endl;
		exit(0);
	}

	band_energy = new double*[nks];

	for (int ik = 0;ik < this->nks;ik++)
	{
		band_energy[ik] = new double[this->nband];
	}

	for (int ik = 0;ik < this->nks;ik++)
	{
		for (int ib = 0;ib < this->nband;ib++)
		{
			in >> this->band_energy[ik][ib];
		}
	}
}

//*****
// evc
//*****
void eximport::out_evc(ofstream &out_data)
{
	//cout << "=== out_evc ===" << endl;
	out_data << setw(20) << "EVC" << endl;
	out_data << setw(20) << ucell.natomwfc << endl; //4.1
	int iw;
	int ik;
	int ig;

	for (iw = 0;iw < ucell.natomwfc;iw++)
	{
		for (ik = 0;ik < kv.nks;ik++)
		{
			int npw = kv.ngk[ik];

			for (ig = 0;ig < npw;ig++)
			{
				out_data << setw(20) << wf.evc[ik](iw, ig).real() << setw(20) << wf.evc[ik](iw, ig).imag() << endl;//4.2
			}
		}
	}
}

void eximport::in_evc(ifstream &in)
{
	cout << "=== in_evc ===" << endl;
	in >> name;

	if (name != "EVC")
	{
		cout << "== wrong in in_evc ===" << endl;
		exit(0);
	}

	in >> this->natomwfc;//4.1

	this->evc = new complex <double>**[natomwfc];
	int iw = 0;
	int ik = 0;
	int ig = 0;

	for (iw = 0;iw < natomwfc;iw++)
	{
		evc[iw] = new complex <double>*[nks];

		for (ik = 0;ik < nks;ik++)
		{
			evc[iw][ik] = new complex <double>[ngk[ik]];
		}
	}

	double r;

	double i;

	for (iw = 0;iw < natomwfc;iw++)
	{
		for (ik = 0;ik < nks;ik++)
		{
			int npw = ngk[ik];

			for (ig = 0;ig < npw;ig++)
			{
				in >> r >> i;//4.2
				evc[iw][ik][ig] = complex <double>(r, i);
			}
		}
	}
}

//===========
//	energy
//===========
#include "src_pw/H_Ewald_pw.h"
#include "src_pw/H_Hartree_pw.h"
#include "src_pw/H_XC_pw.h"
void eximport::out_energy(ofstream &out_data)
{
	//cout << "\n ==> out_energy" << endl;
	out_data << setw(20) << "ENERGY" << endl;				//6.0
	out_data << setw(20) << en.etot << endl;                //6.2
	//out_data << setw(20) << elec.dE << endl;              //6.3
	out_data << setw(20) << en.eband << endl;				//6.4
	out_data << setw(20) << en.eband + en.deband << endl;   //6.5
	out_data << setw(20) << H_Hartree_pw::hartree_energy << endl;
	out_data << setw(20) << H_XC_pw::etxc - en.etxcc << endl;     //6.7
	out_data << setw(20) << H_Ewald_pw::ewald_energy << endl;                //6.8
}

void eximport::in_energy(ifstream &in)
{
	in >> name;

	if (name != "ENERGY")
	{
		cout << "== wrong in in_energy ===" << endl;
		exit(0);
	}
	in >> iter;//6.1
	in >> etot;//6.2
	in >> eband;//6.4
	in >> one_electron;//6.5
	in >> hartree;//6.6
	in >> xc;//6.7
	in >> ewald;//6.8
}

#ifdef __MPI
void eximport::in_charge_mpi(const string &dir)
{
	cout << "\n ==> eximport::in_charge()" << endl;
	double *rho_tmp = new double[pw.ncxyz]();
	assert(rho_tmp!=0);

	if(MY_RANK == 0)
	{
		ifstream in(dir.c_str());
		in >> name;

		if (name != "CHARGE")
		{
			cout << "== wrong in in_charge ===" << endl;
			exit(0);
		}
	
		in >> omega;

		in >>ncx >> ncy >> ncz;
		int ncxyz;
		in >> ncxyz;
		if(ncxyz != pw.ncxyz)
		{	
			cout<<"\n Read in ncxzy in charge file = "<<ncxyz<<endl;
			QUIT();
		}
		for (int ir = 0;ir < pw.ncxyz;ir++)
		{
			in >> rho_tmp[ir];
		}
		in.close();
	}
	MPI_Barrier(MPI_COMM_WORLD);


	//=========================
	// Read in rho , bcast it 
	//=========================
	MPI_Bcast(rho_tmp,pw.ncxyz,MPI_DOUBLE,0,MPI_COMM_WORLD);

	int ir = 0; // counters on read space mesh points
	int iz = 0; // counters on planes
	int ip = 0; // counters on processor
	//=================================================
    // Find number of planes for each cpu in this pool
    //=================================================
    int *num_z = new int[NPROC_IN_POOL]();
    for(iz=0;iz<pw.ncz;iz++)
    {
        ip = iz % NPROC_IN_POOL;
        num_z[ip]++;
    }

    //=======================================
    // Find current number of planes (nz)
    //=======================================
    int *cur_z = new int[NPROC_IN_POOL]();
    for(ip=1;ip<NPROC_IN_POOL;ip++)
    {
        cur_z[ip] = cur_z[ip-1]+num_z[ip-1];
    }

	for(ir =0;ir<pw.ncx * pw.ncy;ir++)
    {
        for(iz=0;iz<num_z[RANK_IN_POOL];iz++)
        {
            CHR.rho[0][ ir*num_z[RANK_IN_POOL]+iz ]= rho_tmp[ir*pw.ncz + cur_z[RANK_IN_POOL] + iz ];
        }
    }

	delete[] rho_tmp;
	delete[] num_z;
	delete[] cur_z;
	return;
}

void eximport::out_charge_mpi(const string &dir,double* rho_in)
{
	//cout << "\n ==> out_charge" << endl;
	//cout << dir << endl;
	ofstream out_data(dir.c_str());
	if(!out_data)
	{
		cout<<"\n Can't write charge file!"<<endl;
	}
	out_data << setw(20) << "CHARGE" << endl;	//7.0
	out_data << setw(20) << ucell.omega << endl;	//7.1
	out_data << setw(20) << pw.ncx 
			 << setw(20) << pw.ncy 
			 << setw(20) << pw.ncz << endl;		//7.2
	out_data << setw(20) << pw.ncxyz << endl;

	for (int ir = 0;ir < pw.ncxyz;ir++)
	{
		if(ir%4==0) out_data << "\n";
		out_data << setw(20) << setprecision(10) << rho_in[ir] << "\t" ;//7.4
	}
	out_data << endl;
	out_data.close();
	return;
}
#endif

//********************************
// charge
//********************************
void eximport::out_charge(ofstream &out_data)
{
	/*
	ofs_running << "\n Output charge file." << endl;
	out_data << setw(20) << "CHARGE" << endl;	//7.0
	out_data << setw(20) << pw.omega << endl;	//7.1
	out_data << setw(20) << pw.ncx 
			 << setw(20) << pw.ncy 
			 << setw(20) << pw.ncz << endl;		//7.2
	out_data << setw(20) << CHR.rho.nr			//7.3 
			 << setw(20) << CHR.rho.nc << endl;

	for (int i = 0;i < CHR.rho.nr;i++)
	{
		for (int j = 0;j < CHR.rho.nc;j++)
		{
			out_data << setw(20) << setprecision(10) << CHR.rho.c[i*CHR.rho.nr+j] << "\t" ;//7.4

			if ((i*CHR.rho.nr + j) % 4 == 3) out_data << endl;
		}
		out_data << endl;
	}
	return;
	*/
}

void eximport::in_charge(ifstream &in)
{
	in >> name;
	if (name != "CHARGE")
	{
		cout << "== wrong in in_charge ===" << endl;
		exit(0);
	}

	in >> omega;

	in >>ncx >> ncy >> ncz;
	in >> rho_nr;
	in >> rho_nc;
	rho = new double[rho_nc];

	for (int i = 0;i < rho_nr;i++)
	{
		for (int j = 0;j < rho_nc;j++)
		{
			in >> rho[j];
		}
	}
}

void eximport::nscf_chgfile(const string &chg_file)
{
	/*
	cout<<"\n ==> eximport::nscf_chgfile()";
	cout<<"       file_name : "<<chg_file<<endl;
	//int ok;
	//cin >> ok;
	ifstream in(chg_file.c_str());
	if(!in)
	{
		cerr<<"nscf_chgfile : can't find file "<<chg_file<<endl;
		exit(0);
	}
	this->in_input(in);
	this->in_wannier(in);
	this->in_unitcell(in);
	this->in_kpoints(in);
	in >> name;
	cout<<"name = "<<name<<endl;

	in >> omega;

	in >>ncx >> ncy >> ncz;
	in >> rho_nr;
	in >> rho_nc;
	cout<<"rho_nc = "<<rho_nc<<endl;
	cout<<"rho_nr = "<<rho_nr<<endl;
	rho = new double[rho_nc];
	CHR.rho.nr = rho_nr;
	CHR.rho.nc = rho_nc;

	for (int i = 0;i < rho_nr;i++)
	{
		for (int j = 0;j < rho_nc;j++)
		{
//			cout<<"j="<<j<<endl;
			in >> CHR.rho.c[j];
		}
	}
	in.close();
	*/
}

