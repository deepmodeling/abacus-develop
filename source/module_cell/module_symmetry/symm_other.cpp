#include "symm_other.h"

namespace ModuleSymmetry
{
void Symm_Other::print1(const int &ibrav, const double *cel_const, std::ofstream &ofs_running)
{
	ModuleBase::TITLE("Symm_Other","print1");

	ModuleBase::GlobalFunc::OUT(ofs_running,"IBRAV",ibrav);
	if(ibrav==1)
	{
		ModuleBase::GlobalFunc::OUT(ofs_running,"BRAVAIS","SIMPLE CUBIC");
    	ModuleBase::GlobalFunc::OUT(ofs_running,"LATTICE CONSTANT A",cel_const[0]);
	}
	else if(ibrav==2)
	{
		ModuleBase::GlobalFunc::OUT(ofs_running,"BRAVAIS","BODY CENTERED CUBIC");
    	ModuleBase::GlobalFunc::OUT(ofs_running,"LATTICE CONSTANT A",cel_const[0]);
	}
	else if(ibrav==3)
	{
		ModuleBase::GlobalFunc::OUT(ofs_running,"BRAVAIS","FACE CENTERED CUBIC");
    	ModuleBase::GlobalFunc::OUT(ofs_running,"LATTICE CONSTANT A",cel_const[0]);
	}
	else if(ibrav==4)
	{
		ModuleBase::GlobalFunc::OUT(ofs_running,"BRAVAIS","HEXAGONAL CELL");
    	ModuleBase::GlobalFunc::OUT(ofs_running,"LATTICE CONSTANT A",cel_const[0]);
		ModuleBase::GlobalFunc::OUT(ofs_running,"C/A RATIO",cel_const[2]);
	}
	else if(ibrav==5)
	{
		ModuleBase::GlobalFunc::OUT(ofs_running,"BRAVAIS","SIMPLE TETROGONAL CELL");
    	ModuleBase::GlobalFunc::OUT(ofs_running,"LATTICE CONSTANT A",cel_const[0]);
		ModuleBase::GlobalFunc::OUT(ofs_running,"C/A RATIO",cel_const[2]);
	}
	else if(ibrav==6)
	{
		ModuleBase::GlobalFunc::OUT(ofs_running,"BRAVAIS","BODY CENTERED TETROGONAL CELL");
    	ModuleBase::GlobalFunc::OUT(ofs_running,"LATTICE CONSTANT A",cel_const[0]);
		ModuleBase::GlobalFunc::OUT(ofs_running,"C/A RATIO",cel_const[2]);
	}
	else if(ibrav==7)
	{
		ModuleBase::GlobalFunc::OUT(ofs_running,"BRAVAIS","TRIGONAL (RHOMBOEDRIC) CELL");
    	ModuleBase::GlobalFunc::OUT(ofs_running,"LATTICE CONSTANT A",cel_const[0]);
		ModuleBase::GlobalFunc::OUT(ofs_running,"COS(ALPHA)",cel_const[3]);
	}
	else if(ibrav==8)
	{
		ModuleBase::GlobalFunc::OUT(ofs_running,"BRAVAIS","SIMPLE ORTHORHOMBIC CELL");
    	ModuleBase::GlobalFunc::OUT(ofs_running,"LATTICE CONSTANT A",cel_const[0]);
		ModuleBase::GlobalFunc::OUT(ofs_running,"B/A RATIO",cel_const[1]);
		ModuleBase::GlobalFunc::OUT(ofs_running,"C/A RATIO",cel_const[2]);
	}
	else if(ibrav==9)
	{
		ModuleBase::GlobalFunc::OUT(ofs_running,"BRAVAIS","BODY CENTERED ORTHORHOMBIC CELL");
    	ModuleBase::GlobalFunc::OUT(ofs_running,"LATTICE CONSTANT A",cel_const[0]);
		ModuleBase::GlobalFunc::OUT(ofs_running,"B/A RATIO",cel_const[1]);
		ModuleBase::GlobalFunc::OUT(ofs_running,"C/A RATIO",cel_const[2]);
	}
	else if(ibrav==10)
	{
		ModuleBase::GlobalFunc::OUT(ofs_running,"BRAVAIS","FACE CENTERED ORTHORHOMBIC CELL");
    	ModuleBase::GlobalFunc::OUT(ofs_running,"LATTICE CONSTANT A",cel_const[0]);
		ModuleBase::GlobalFunc::OUT(ofs_running,"B/A RATIO",cel_const[1]);
		ModuleBase::GlobalFunc::OUT(ofs_running,"C/A RATIO",cel_const[2]);
	}
	else if(ibrav==11)
	{
		ModuleBase::GlobalFunc::OUT(ofs_running,"BRAVAIS","BASE CENTERED ORTHORHOMBIC CELL");
    	ModuleBase::GlobalFunc::OUT(ofs_running,"LATTICE CONSTANT A",cel_const[0]);
		ModuleBase::GlobalFunc::OUT(ofs_running,"B/A RATIO",cel_const[1]);
		ModuleBase::GlobalFunc::OUT(ofs_running,"C/A RATIO",cel_const[2]);
	}
	else if(ibrav==12)
	{
		ModuleBase::GlobalFunc::OUT(ofs_running,"BRAVAIS","SIMPLE MONOLINIC CELL");
    	ModuleBase::GlobalFunc::OUT(ofs_running,"LATTICE CONSTANT A",cel_const[0]);
		ModuleBase::GlobalFunc::OUT(ofs_running,"B/A RATIO",cel_const[1]);
		ModuleBase::GlobalFunc::OUT(ofs_running,"C/A RATIO",cel_const[2]);
		ModuleBase::GlobalFunc::OUT(ofs_running,"COS(BETA)",cel_const[4]);
	}
	else if(ibrav==13)
	{
		ModuleBase::GlobalFunc::OUT(ofs_running,"BRAVAIS","BASE CENTERED MONOLINIC CELL");
    	ModuleBase::GlobalFunc::OUT(ofs_running,"LATTICE CONSTANT A",cel_const[0]);
		ModuleBase::GlobalFunc::OUT(ofs_running,"B/A RATIO",cel_const[1]);
		ModuleBase::GlobalFunc::OUT(ofs_running,"C/A RATIO",cel_const[2]);
		ModuleBase::GlobalFunc::OUT(ofs_running,"COS(BETA)",cel_const[4]);
	}
	else if(ibrav==14)
	{
		ModuleBase::GlobalFunc::OUT(ofs_running,"BRAVAIS","TRICLINIC CELL");
    	ModuleBase::GlobalFunc::OUT(ofs_running,"LATTICE CONSTANT A",cel_const[0]);
		ModuleBase::GlobalFunc::OUT(ofs_running,"B/A RATIO",cel_const[1]);
		ModuleBase::GlobalFunc::OUT(ofs_running,"C/A RATIO",cel_const[2]);
		ModuleBase::GlobalFunc::OUT(ofs_running,"COS(ALPHA)",cel_const[3]);	
		ModuleBase::GlobalFunc::OUT(ofs_running,"COS(BETA)",cel_const[4]);	
		ModuleBase::GlobalFunc::OUT(ofs_running,"COS(GAMMA)",cel_const[5]);	
	}
	else
	{
		ModuleBase::WARNING_QUIT("Symm_Other::print1","ibrav is wrong.");
	}
	return;
}

bool Symm_Other::right_hand_sense(ModuleBase::Vector3<double> &v1,ModuleBase::Vector3<double> &v2,ModuleBase::Vector3<double> &v3)
{
	double volume = Symm_Other::celvol(v1,v2,v3);
	//OUT(ofs_running,"volume = ",volume);
	if(volume < 0)
	{
		v1.reverse();
		v2.reverse();
		v3.reverse();
		return false;
	}
	return true;
}

//calculate the volume of the cell spanned by the vectors
double Symm_Other::celvol(const ModuleBase::Vector3<double> &a, 
const ModuleBase::Vector3<double> &b, const ModuleBase::Vector3<double> &c)
{
	return a.x * ( b.y * c.z - b.z * c.y ) + a.y * ( b.z * c.x - b.x * c.z ) 
	+ a.z * ( b.x * c.y - b.y * c.x );
}
}

