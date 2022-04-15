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
	const std::string &fn,
	const std::string &type
)
{
	ModuleBase::TITLE("eximport","write_data");
	std::ofstream ofs(fn.c_str());

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

void eximport::read_data(const std::string &fn)
{
	ModuleBase::TITLE("eximport","read_data");

	std::ifstream ifs(fn.c_str());
	if(!ifs)
	{
		GlobalV::ofs_warning << " File name : " << fn << std::endl;
		ModuleBase::WARNING_QUIT("eximport::read_data","Can not find file.");
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

void eximport::print_data(const std::string &fn) const
{
	ModuleBase::TITLE("eximport","print_data");
	std::ofstream ofs( fn.c_str() );

	ofs << std::setw(20) << "basis" << std::setw(20) << this->basis << std::endl;
	ofs << std::setw(20) << "latname" << std::setw(20) << this->latname << std::endl;
	ofs << std::setw(20) << "calculation" << std::setw(20) << this->calculation << std::endl;
	ofs << std::setw(20) << "ecutwfc" << std::setw(20) << this->ecutwfc << std::endl;
	ofs << std::setw(20) << "nband" << std::setw(20) << this->nband << std::endl;
	ofs << std::setw(20) << "tr2" << std::setw(20) << this->tr2 << std::endl;
	ofs << std::setw(20) << "nx" << std::setw(20) << this->nx << std::endl;
	ofs << std::setw(20) << "ny" << std::setw(20) << this->ny << std::endl;
	ofs << std::setw(20) << "nz" << std::setw(20) << this->nz << std::endl;
	ofs << std::setw(20) << "nxyz" << std::setw(20) << this->nxyz << std::endl;
	ofs << std::setw(20) << "startingpot" << std::setw(20) << this->startingpot << std::endl;
	ofs << std::setw(20) << "mixing_beta" << std::setw(20) << this->Mixing_beta << std::endl;

	ofs << std::setw(20) << "nks" << std::setw(20) << this->nks << std::endl;

	for (int ik = 0;ik < GlobalC::kv.nks;ik++)
	{ 
		ofs << std::setw(20) << "ngk[" << ik << "]" << std::setw(20) << ngk[ik] << std::endl; 
	}

	ofs << std::setw(20) << "qtot" << std::setw(20) << this->qtot << std::endl;
	ofs << std::setw(20) << "lat0" << std::setw(20) << this->lat0 << std::endl;
	ofs << std::setw(20) << "ntype" << std::setw(20) << this->ntype << std::endl;
	ofs << std::setw(20) << "band_energy(k-point,band)" << std::endl;

	for (int ik = 0;ik < this->nks;ik++)
	{
		for (int ib = 0;ib < this->nband;ib++)
		{
			ofs << std::setw(10) << this->band_energy[ik][ib];
		}
		ofs << std::endl;
	}

	ofs << std::setw(20) << "Omega" << std::setw(20) << this->omega << std::endl;
	ofs << std::setw(20) << "rho_nc" << std::setw(20) << this->rho_nc << std::endl;
	ofs << std::setw(20) << "iter" << std::setw(20) << this->iter << std::endl;
	ofs << std::setw(20) << "etot" << std::setw(20) << this->etot << std::endl;
	ofs << std::setw(20) << "eband" << std::setw(20) << this->eband << std::endl;
	ofs << std::setw(20) << "one-electron" << std::setw(20) << this->one_electron << std::endl;
	ofs << std::setw(20) << "hartree" << std::setw(20) << this->hartree << std::endl;
	ofs << std::setw(20) << "exchange-corr" << std::setw(20) << this->xc << std::endl;
	ofs << std::setw(20) << "ewald" << std::setw(20) << this->ewald << std::endl;
	ofs << std::setw(20) << "natomwfc" << std::setw(20) << this->natomwfc << std::endl;

	ofs.close();
}

void eximport::fir_wf(ModuleBase::ComplexMatrix *psi, const int npsi, const std::string &fn)
{
	if(GlobalV::MY_RANK!=0)
	{
		std::cout<<"\n not the GlobalV::MY_RANK processor , return.";
		return;
	}

	//out put wannier function in q space, this is the file name
	//modify 2007-01-24 
	//add nwan , nwan = number of bands
	//PS: if modify here, you must modify the same name in sec_wf( 2 definition)
	std::ofstream ofs(fn.c_str());

    out_input(ofs);
    out_kpoints(ofs);
    out_igk(ofs);
    out_planewave(ofs);

    ofs << "\n<WAVEFUNC>";
    ofs << "\n" << npsi << " Number of wave functions." << std::endl;
    ofs << std::setprecision(6);

    for(int i=0; i<npsi; i++)
    {
        for(int ik=0; ik<GlobalC::kv.nks; ik++)
        {
            ofs << "\n" << ik;
            for(int ig=0; ig<GlobalC::kv.ngk[ik]; ig++)
            {
                if(ig%4==0) ofs << "\n";
                ofs << std::setw(15) << psi[ik](i, ig).real()
                    << std::setw(15) << psi[ik](i, ig).imag();
            }
            ofs << "\n";
        }
    }
    ofs << "\n<WAVEFUNC>";

	return;
}

void eximport::out_gspace_wan(const ModuleBase::ComplexMatrix *psi,const int iw,const std::string &file_name)
{
	//std::cout<<"\n ==> ei.wanG_for_fit";
	std::ofstream out_gwan(file_name.c_str());
	int qtot = 0;
	for(int ik=0; ik<GlobalC::kv.nks; ik++)
	{
		for(int ig=0; ig<GlobalC::kv.ngk[ik]; ig++)
		{
			qtot += GlobalC::kv.ngk[ik];
		}
	}
	int nks = GlobalC::kv.nks;
	double factor = ModuleBase::TWO_PI/GlobalC::ucell.lat0;
	out_gwan << qtot << std::endl;
	for(int ik=0;ik<nks;ik++)
	{
		//output wannier functions in G space.
		for(int ig=0;ig<GlobalC::kv.ngk[ik];ig++)
		{
			double g1 = GlobalC::pw.get_GPlusK_cartesian_projection(ik, GlobalC::wf.igk(ik, ig), 0);
			double g2 = GlobalC::pw.get_GPlusK_cartesian_projection(ik, GlobalC::wf.igk(ik, ig), 1);
			double g3 = GlobalC::pw.get_GPlusK_cartesian_projection(ik, GlobalC::wf.igk(ik, ig), 2);
			out_gwan 
			<< std::setw(15) << g1*factor 
			<< std::setw(15) << g2*factor
			<< std::setw(15) << g3*factor
			<< std::setw(20) << psi[ik](iw,ig).real() 
			<< std::setw(20) << psi[ik](iw,ig).imag() << std::endl;
		}
	}

	out_gwan.close();
	return;
}

bool eximport::sec_wf(ModuleBase::ComplexMatrix *psi, const int wf_num, const std::string &fn)const
{
	std::cout<<"\n ==> ei.sec_wf()_1";
	std::ifstream sec(fn.c_str());

	if (sec)
	{
		std::cout << "Find the file :"<<std::endl; 
		std::cout << fn << std::endl;
	}

	if (!sec)
	{
		std::cout << "Can't open file : " << fn << "\treturn!" << std::endl;
		return 0;
	}

	for (int ik = 0;ik != GlobalC::kv.nks;ik++)
	{
		for (int iw = 0;iw != wf_num;iw++)
		{
			for (int ig = 0;ig != GlobalC::kv.ngk[ik];ig++)
			{
				// Peize Lin fix bug about rvalue 2016-08-02
				double tmp_real, tmp_imag;
				sec >> tmp_real >> tmp_imag ;
				psi[ik](iw, ig) = std::complex<double>(tmp_real,tmp_imag);
			}
		}
	}

	sec.close();
	return 1;
}


bool eximport::sec_wf(std::complex < double> ***psi, const int npsi, const std::string &fn)const
{
	std::cout<<"\n ==> ei.sec_wf()_2";
	std::ifstream sec(fn.c_str());

	if (sec)
	{
		std::cout << "\n    Find the file :"<<std::endl; 
		std::cout << fn << std::endl;
	}

	if (!sec)
	{
		std::cout << "Can't open file : " << fn << "\treturn!" << std::endl;
		return 0;
	}

	for (int ik = 0;ik < GlobalC::kv.nks;ik++)
	{
		for (int iw = 0;iw < npsi;iw++)
		{
			for (int ig = 0;ig < GlobalC::kv.ngk[ik];ig++)
			{
				// Peize Lin fix bug about rvalue 2016-08-02
				double tmp_real, tmp_imag;
				sec >> tmp_real >> tmp_imag ;
				psi[iw][ik][ig] = std::complex<double>(tmp_real,tmp_imag);				
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
void eximport::out_wannier(std::ofstream &out_data)
{
	//std::cout<<"\n ==> out_wannier"<<std::endl;
	out_data << std::setw(20) << "WANNIER" << std::endl; //0
	//out_data << std::setw(20) << LOCAL_BASIS << std::endl;//1 xiaohui modify 2013-09-02
	out_data << std::setw(20) << GlobalV::BASIS_TYPE << std::endl; //xiaohui add 2013-09-02
	return;
}

void eximport::in_wannier(std::ifstream &in)
{
	std::cout<<"=== in_wannier ==="<<std::endl;
	in >> name;
	//std::cout<<name<<std::endl;		//8.0
	if (name != "WANNIER")
	{
		std::cout << "== wrong in in_wannier ===" << std::endl;
		exit(0);
	}
	in >> this->basis;	//8.1
	//std::cout<<"basis="<<this->basis<<std::endl;
	return;
}

//*************************************
// atom_unitcell
//*************************************
void eximport::out_unitcell(std::ofstream &out_data)
{
	//std::cout << "\n ==> out_unitcell" << std::endl;
	out_data << std::setw(20) << "UNITCELL" << std::endl;		//2.0
	out_data << std::setw(20) << GlobalC::ucell.lat0 << std::endl;        //2.1

	out_data << std::setw(20) << GlobalC::ucell.latvec.e11 
			 << std::setw(20) << GlobalC::ucell.latvec.e12 
			 << std::setw(20) << GlobalC::ucell.latvec.e13 << std::endl;	//2.2
	out_data << std::setw(20) << GlobalC::ucell.latvec.e21 
			 << std::setw(20) << GlobalC::ucell.latvec.e22 
			 << std::setw(20) << GlobalC::ucell.latvec.e23 << std::endl;
	out_data << std::setw(20) << GlobalC::ucell.latvec.e31 
			 << std::setw(20) << GlobalC::ucell.latvec.e32 
			 << std::setw(20) << GlobalC::ucell.latvec.e33 << std::endl;

	out_data << std::setw(20) << GlobalC::ucell.ntype << std::endl;		//2.3

	for (int i = 0;i < GlobalC::ucell.ntype;i++)
	{
		out_data << std::setw(20) << GlobalC::ucell.atoms[i].na ;     //2.4
	}
	out_data << std::endl;
	return;
}

void eximport::in_unitcell(std::ifstream &in)
{
	std::cout << "=== in_unitcell ===" << std::endl;
	in >> name;
	if (name != "UNITCELL")
	{
		std::cout << "== wrong in in_unitcell ===" << std::endl;
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
void eximport::out_kpoints(std::ofstream &out_data)
{
	//std::cout << "\n ==> out_k-points" << std::endl;
	out_data << "\n<KPOINT>";
	out_data << "\n" << GlobalC::kv.nkstot << " Number of total k points";      //3.1
	int sumq = 0;

	for(int ik=0; ik<GlobalC::kv.nks; ik++)
	{
		if(ik%10==0) out_data<<std::endl;
		out_data << std::setw(10) << GlobalC::kv.ngk[ik];
		sumq += GlobalC::kv.ngk[ik];
	}

	for (int ik = 0;ik < GlobalC::kv.nks;ik++)
	{
		if(ik%3==0) out_data << "\n";
		out_data << std::setw(10) << GlobalC::kv.kvec_c[ik].x
		<< std::setw(10) << GlobalC::kv.kvec_c[ik].y
		<< std::setw(10) << GlobalC::kv.kvec_c[ik].z;//3.3
	}

	out_data << "\n" << sumq << " Total Number of K+G points.";//3.4
	out_data << "\n<KPOINT>" << std::endl;
	return;
}

void eximport::out_planewave(std::ofstream &out_data)
{
	//std::cout << "\n ==> out_planewave" << std::endl;
	out_data << "\n<PLANEWAVE>";
	out_data << "\n" << GlobalC::ucell.lat0 << " Lattice constant";
	out_data << "\n" << GlobalC::pw.ngmc_g << " Number of plane waves."<<std::endl;
	for(int i=0; i<GlobalC::pw.ngmc_g; i++)
	{
		if(i%4==0) out_data<<"\n";
		out_data << std::setw(8) << GlobalC::pw.get_G_cartesian_projection(i, 0)
				 << std::setw(8) << GlobalC::pw.get_G_cartesian_projection(i, 1)
				 << std::setw(8) << GlobalC::pw.get_G_cartesian_projection(i, 2);
	}
	out_data << "\n<PLANEWAVE>";
	return;
}

void eximport::out_igk(std::ofstream &out_data)
{
	//std::cout << "\n ==> out_igk" << std::endl;
	out_data << "\n<KG_INDEX>";

	for(int ik=0; ik<GlobalC::kv.nks; ik++)
	{
		for(int ig=0; ig<GlobalC::kv.ngk[ik]; ig++)
		{
			if(ig%10==0) out_data<<"\n";
			out_data << std::setw(10) << GlobalC::wf.igk(ik, ig);
		}
	}
	out_data << "\n<KG_INDEX>";
	return;
}

void eximport::in_kpoints(std::ifstream &in)
{
	std::cout << "=== in_kpoints ===" << std::endl;
	in >> name;
	if (name != "KPOINT")
	{
		std::cout << "== wrong in in_input ===" << std::endl;
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
void eximport::out_input(std::ofstream &out_data)
{
	//std::cout << "\n ==> out_input" << std::endl;
	out_data << "<HEADER>";			//1.0
	out_data << "\n" << GlobalC::ucell.latName << " Lattice Name";//1.1

/*
	if (SCF)
	{
		out_data << std::setw(20) << "scf" << std::endl;       	//1.2
	}
	else
	{
		out_data << std::setw(20) << "nscf" << std::endl;
	}
*/

	out_data << "\n" << GlobalC::pw.ecutwfc << " Energy Cutoff for wave functions.";//1.3

	out_data << "\n" << GlobalC::ucell.nat << " Number of atoms.";

	out_data << "\n" << GlobalC::ucell.ntype << " Types of atoms.";

//	out_data << "\n" << GlobalV::NBANDS << " Number of wave functions.";//1.4
	//out_data << "\n" << LOCAL_BASIS << " 1 for Local Basis"; xiaohui modify 2013-09-02
	out_data << "\n" << GlobalV::BASIS_TYPE << " Basis type"; //xiaohui add 2013-09-02

	//out_data << "\n" << LINEAR_SCALING << " 1 for Linear Scaling"; xiaohui modify 2013-09-02
	out_data << "\n" << GlobalV::KS_SOLVER << " Diago type"; //xiaohui add 2013-09-02
	out_data << "\n" << GlobalV::SPARSE_MATRIX << " 1 for Sparse Matrix";

/*
	out_data << std::setw(20) << tr2 << std::endl;                 //1.5
	out_data << std::setw(20) << GlobalC::pw.nx 
			 << std::setw(20) << GlobalC::pw.ny 
			 << std::setw(20) << GlobalC::pw.nz 
			 << std::setw(20) << GlobalC::pw.nxyz << std::endl;//1.6
	out_data << std::setw(20) << GlobalC::pot.startingpot << std::endl;//1.7
	out_data << std::setw(20) << GlobalC::CHR.mixing_beta << std::endl;//1.8
*/
	out_data << "\n<HEADER>"<<std::endl;			//1.0

	return;
}

void eximport::in_input(std::ifstream &in)
{
	std::cout << " in_input ===" << std::endl;
	in >> name;

	if (name != "INPUT")
	{
		std::cout << "== wrong in in_input ===" << std::endl;
		exit(0);
	}

	in >> this->latname;//1.1

	in >> this->calculation;//1.2
	double ecut_in = 0.0;
	in >> ecut_in;//1.3
	if(ecut_in != GlobalC::pw.ecutwfc)
	{
		std::cout<<"Charge don't match!"<<std::endl;
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
void eximport::out_band(std::ofstream &out_data)
{
	//std::cout << "\n ==> out_band" << std::endl;
	out_data << std::setw(20) << "BAND" << std::endl;//6.0
	for (int ik = 0; ik < GlobalC::kv.nks; ik++)
	{
		for (int ib = 0; ib < GlobalV::NBANDS; ib++)
		{
			out_data << std::setw(10) << std::setprecision(6) << GlobalC::wf.ekb[ik][ib]*ModuleBase::Ry_to_eV;//6.1
		}
		out_data << std::endl;
	}
	out_data << std::endl;
	return;
}

void eximport::in_band(std::ifstream &in)
{
	std::cout << "=== in_band ===" << std::endl;
	in >> name;
	if (name != "BAND")
	{
		std::cout << "== wrong in in_band ===" << std::endl;
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
void eximport::out_evc(std::ofstream &out_data)
{
	//std::cout << "=== out_evc ===" << std::endl;
	out_data << std::setw(20) << "EVC" << std::endl;
	out_data << std::setw(20) << GlobalC::ucell.natomwfc << std::endl; //4.1
	int iw;
	int ik;
	int ig;

	for (iw = 0;iw < GlobalC::ucell.natomwfc;iw++)
	{
		for (ik = 0;ik < GlobalC::kv.nks;ik++)
		{
			int npw = GlobalC::kv.ngk[ik];

			for (ig = 0;ig < npw;ig++)
			{
				out_data << std::setw(20) << GlobalC::wf.evc[ik](iw, ig).real() << std::setw(20) << GlobalC::wf.evc[ik](iw, ig).imag() << std::endl;//4.2
			}
		}
	}
}

void eximport::in_evc(std::ifstream &in)
{
	std::cout << "=== in_evc ===" << std::endl;
	in >> name;

	if (name != "EVC")
	{
		std::cout << "== wrong in in_evc ===" << std::endl;
		exit(0);
	}

	in >> this->natomwfc;//4.1

	this->evc = new std::complex <double>**[natomwfc];
	int iw = 0;
	int ik = 0;
	int ig = 0;

	for (iw = 0;iw < natomwfc;iw++)
	{
		evc[iw] = new std::complex <double>*[nks];

		for (ik = 0;ik < nks;ik++)
		{
			evc[iw][ik] = new std::complex <double>[ngk[ik]];
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
				evc[iw][ik][ig] = std::complex <double>(r, i);
			}
		}
	}
}

//===========
//	energy
//===========
#include "../src_pw/H_Ewald_pw.h"
#include "../src_pw/H_Hartree_pw.h"
void eximport::out_energy(std::ofstream &out_data)
{
	//std::cout << "\n ==> out_energy" << std::endl;
	out_data << std::setw(20) << "ENERGY" << std::endl;				//6.0
	out_data << std::setw(20) << GlobalC::en.etot << std::endl;                //6.2
	//out_data << std::setw(20) << elec.dE << std::endl;              //6.3
	out_data << std::setw(20) << GlobalC::en.eband << std::endl;				//6.4
	out_data << std::setw(20) << GlobalC::en.eband + GlobalC::en.deband << std::endl;   //6.5
	out_data << std::setw(20) << H_Hartree_pw::hartree_energy << std::endl;
	out_data << std::setw(20) << GlobalC::en.etxc - GlobalC::en.etxcc << std::endl;     //6.7
	out_data << std::setw(20) << H_Ewald_pw::ewald_energy << std::endl;                //6.8
}

void eximport::in_energy(std::ifstream &in)
{
	in >> name;

	if (name != "ENERGY")
	{
		std::cout << "== wrong in in_energy ===" << std::endl;
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
void eximport::in_charge_mpi(const std::string &dir)
{
	std::cout << "\n ==> eximport::in_charge()" << std::endl;
	double *rho_tmp = new double[GlobalC::pw.ncxyz]();
	assert(rho_tmp!=0);

	if(GlobalV::MY_RANK == 0)
	{
		std::ifstream in(dir.c_str());
		in >> name;

		if (name != "CHARGE")
		{
			std::cout << "== wrong in in_charge ===" << std::endl;
			exit(0);
		}
	
		in >> omega;

		in >>ncx >> ncy >> ncz;
		int ncxyz;
		in >> ncxyz;
		if(ncxyz != GlobalC::pw.ncxyz)
		{	
			std::cout<<"\n Read in ncxzy in charge file = "<<ncxyz<<std::endl;
			ModuleBase::QUIT();
		}
		for (int ir = 0;ir < GlobalC::pw.ncxyz;ir++)
		{
			in >> rho_tmp[ir];
		}
		in.close();
	}
	MPI_Barrier(MPI_COMM_WORLD);


	//=========================
	// Read in rho , bcast it 
	//=========================
	MPI_Bcast(rho_tmp,GlobalC::pw.ncxyz,MPI_DOUBLE,0,MPI_COMM_WORLD);

	int ir = 0; // counters on read space mesh points
	int iz = 0; // counters on planes
	int ip = 0; // counters on processor
	//=================================================
    // Find number of planes for each cpu in this pool
    //=================================================
    int *num_z = new int[GlobalV::NPROC_IN_POOL]();
    for(iz=0;iz<GlobalC::pw.ncz;iz++)
    {
        ip = iz % GlobalV::NPROC_IN_POOL;
        num_z[ip]++;
    }

    //=======================================
    // Find current number of planes (nz)
    //=======================================
    int *cur_z = new int[GlobalV::NPROC_IN_POOL]();
    for(ip=1;ip<GlobalV::NPROC_IN_POOL;ip++)
    {
        cur_z[ip] = cur_z[ip-1]+num_z[ip-1];
    }

	for(ir =0;ir<GlobalC::pw.ncx * GlobalC::pw.ncy;ir++)
    {
        for(iz=0;iz<num_z[GlobalV::RANK_IN_POOL];iz++)
        {
            GlobalC::CHR.rho[0][ ir*num_z[GlobalV::RANK_IN_POOL]+iz ]= rho_tmp[ir*GlobalC::pw.ncz + cur_z[GlobalV::RANK_IN_POOL] + iz ];
        }
    }

	delete[] rho_tmp;
	delete[] num_z;
	delete[] cur_z;
	return;
}

void eximport::out_charge_mpi(const std::string &dir,double* rho_in)
{
	//std::cout << "\n ==> out_charge" << std::endl;
	//std::cout << dir << std::endl;
	std::ofstream out_data(dir.c_str());
	if(!out_data)
	{
		std::cout<<"\n Can't write charge file!"<<std::endl;
	}
	out_data << std::setw(20) << "CHARGE" << std::endl;	//7.0
	out_data << std::setw(20) << GlobalC::ucell.omega << std::endl;	//7.1
	out_data << std::setw(20) << GlobalC::pw.ncx 
			 << std::setw(20) << GlobalC::pw.ncy 
			 << std::setw(20) << GlobalC::pw.ncz << std::endl;		//7.2
	out_data << std::setw(20) << GlobalC::pw.ncxyz << std::endl;

	for (int ir = 0;ir < GlobalC::pw.ncxyz;ir++)
	{
		if(ir%4==0) out_data << "\n";
		out_data << std::setw(20) << std::setprecision(10) << rho_in[ir] << "\t" ;//7.4
	}
	out_data << std::endl;
	out_data.close();
	return;
}
#endif

//********************************
// charge
//********************************
void eximport::out_charge(std::ofstream &out_data)
{
	/*
	GlobalV::ofs_running << "\n Output charge file." << std::endl;
	out_data << std::setw(20) << "CHARGE" << std::endl;	//7.0
	out_data << std::setw(20) << GlobalC::pw.omega << std::endl;	//7.1
	out_data << std::setw(20) << GlobalC::pw.ncx 
			 << std::setw(20) << GlobalC::pw.ncy 
			 << std::setw(20) << GlobalC::pw.ncz << std::endl;		//7.2
	out_data << std::setw(20) << GlobalC::CHR.rho.nr			//7.3 
			 << std::setw(20) << GlobalC::CHR.rho.nc << std::endl;

	for (int i = 0;i < GlobalC::CHR.rho.nr;i++)
	{
		for (int j = 0;j < GlobalC::CHR.rho.nc;j++)
		{
			out_data << std::setw(20) << std::setprecision(10) << GlobalC::CHR.rho.c[i*GlobalC::CHR.rho.nr+j] << "\t" ;//7.4

			if ((i*GlobalC::CHR.rho.nr + j) % 4 == 3) out_data << std::endl;
		}
		out_data << std::endl;
	}
	return;
	*/
}

void eximport::in_charge(std::ifstream &in)
{
	in >> name;
	if (name != "CHARGE")
	{
		std::cout << "== wrong in in_charge ===" << std::endl;
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

void eximport::nscf_chgfile(const std::string &chg_file)
{
	/*
	std::cout<<"\n ==> eximport::nscf_chgfile()";
	std::cout<<"       file_name : "<<chg_file<<std::endl;
	//int ok;
	//cin >> ok;
	std::ifstream in(chg_file.c_str());
	if(!in)
	{
		std::cerr<<"nscf_chgfile : can't find file "<<chg_file<<std::endl;
		exit(0);
	}
	this->in_input(in);
	this->in_wannier(in);
	this->in_unitcell(in);
	this->in_kpoints(in);
	in >> name;
	std::cout<<"name = "<<name<<std::endl;

	in >> omega;

	in >>ncx >> ncy >> ncz;
	in >> rho_nr;
	in >> rho_nc;
	std::cout<<"rho_nc = "<<rho_nc<<std::endl;
	std::cout<<"rho_nr = "<<rho_nr<<std::endl;
	rho = new double[rho_nc];
	GlobalC::CHR.rho.nr = rho_nr;
	GlobalC::CHR.rho.nc = rho_nc;

	for (int i = 0;i < rho_nr;i++)
	{
		for (int j = 0;j < rho_nc;j++)
		{
//			std::cout<<"j="<<j<<std::endl;
			in >> GlobalC::CHR.rho.c[j];
		}
	}
	in.close();
	*/
}

