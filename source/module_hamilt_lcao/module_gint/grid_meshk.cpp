#include "grid_meshk.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

Grid_MeshK::Grid_MeshK()
{
	ucell_index2x = nullptr;
	ucell_index2y = nullptr;
	ucell_index2z = nullptr;
}

Grid_MeshK::~Grid_MeshK()
{
	delete[] ucell_index2x;
	delete[] ucell_index2y;
	delete[] ucell_index2z;
}

int Grid_MeshK::cal_Rindex(const int &u1, const int &u2, const int &u3)const
{
	const int x1 = u1 - this->minu1;
	const int x2 = u2 - this->minu2;
	const int x3 = u3 - this->minu3;
	
	if(x1<0 || x2<0 || x3<0)
	{
		std::cout << " u1=" << u1 << " minu1=" << minu1 << std::endl;
		std::cout << " u2=" << u2 << " minu2=" << minu2 << std::endl;
		std::cout << " u3=" << u3 << " minu3=" << minu3 << std::endl;
		ModuleBase::WARNING_QUIT("Grid_MeshK::cal_Rindex","x1<0 || x2<0 || x3<0 !");
	}

	assert(x1>=0);
	assert(x2>=0);
	assert(x3>=0);
	return (x3 + x2 * this->nu3 + x1 * this->nu2 * this->nu3);
}

void Grid_MeshK::cal_extended_cell(const int &dxe, const int &dye, const int &dze)
{
	ModuleBase::TITLE("Grid_MeshK","cal_extended_cell");

	//--------------------------------------
	// max and min unitcell in expaned grid.
	//--------------------------------------
	this->maxu1 = dxe / GlobalC::bigpw->nbx + 1;
	this->maxu2 = dye / GlobalC::bigpw->nby + 1;
	this->maxu3 = dze / GlobalC::bigpw->nbz + 1;

	this->minu1 = (-dxe+1) / GlobalC::bigpw->nbx - 1; 
	this->minu2 = (-dye+1) / GlobalC::bigpw->nby - 1; 
	this->minu3 = (-dze+1) / GlobalC::bigpw->nbz - 1; 

	if(GlobalV::test_gridt)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"MaxUnitcell",maxu1,maxu2,maxu3);
	if(GlobalV::test_gridt)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"MinUnitcell",minu1,minu2,minu3);

	//--------------------------------------
	// number of unitcell in each direction.
	//--------------------------------------
	this->nu1 = maxu1 - minu1 + 1;
	this->nu2 = maxu2 - minu2 + 1;
	this->nu3 = maxu3 - minu3 + 1;
	this->nutot = nu1 * nu2 * nu3;

	if(GlobalV::test_gridt)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"UnitCellNumber",nu1,nu2,nu3);
	//xiaohui add 'GlobalV::OUT_LEVEL' line, 2015-09-16
	if(GlobalV::OUT_LEVEL != "m") ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"UnitCellTotal",nutot);

//	std::cout << " nu1 = " << nu1 << " nu2 = " << nu2 << " nu3 = " << nu3 << std::endl;
//	std::cout << " nutot = " << nutot << std::endl;

	delete[] ucell_index2x;
	delete[] ucell_index2y;
	delete[] ucell_index2z;
	this->ucell_index2x = new int[nutot];
	this->ucell_index2y = new int[nutot];
	this->ucell_index2z = new int[nutot];
	ModuleBase::GlobalFunc::ZEROS(ucell_index2x, nutot);
	ModuleBase::GlobalFunc::ZEROS(ucell_index2y, nutot);
	ModuleBase::GlobalFunc::ZEROS(ucell_index2z, nutot);

	this->nutot = nu1 * nu2 * nu3;

	for(int i=minu1; i<=maxu1; i++)
	{
		for(int j=minu2; j<=maxu2; j++)
		{
			for(int k=minu3; k<=maxu3; k++)
			{
				const int cell = cal_Rindex(i,j,k);	
				assert(cell<nutot);

				this->ucell_index2x[cell] = i-minu1;
				this->ucell_index2y[cell] = j-minu2;
				this->ucell_index2z[cell] = k-minu3;

			}
		}
	}

	return;
}
