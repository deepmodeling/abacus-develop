#include "grid_meshk.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_parameter/parameter.h"

Grid_MeshK::Grid_MeshK()
{
}

Grid_MeshK::~Grid_MeshK()
{
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

void Grid_MeshK::init_ucell_para(void)
{
	this->max_ucell_para=std::vector<int>(3,0);
    this->max_ucell_para[0]=this->maxu1;
    this->max_ucell_para[1]=this->maxu2;
    this->max_ucell_para[2]=this->maxu3;
	
	this->min_ucell_para=std::vector<int>(3,0);
    this->min_ucell_para[0]=this->minu1;
    this->min_ucell_para[1]=this->minu2;
    this->min_ucell_para[2]=this->minu3;
    
	this->num_ucell_para=std::vector<int>(4,0);
	this->num_ucell_para[0]=this->nu1;
    this->num_ucell_para[1]=this->nu2;
    this->num_ucell_para[2]=this->nu3;
    this->num_ucell_para[3]=this->nutot;
}


void Grid_MeshK::cal_extended_cell(const int &dxe, const int &dye, const int &dze,const int& nbx, const int& nby, const int& nbz)
{
	ModuleBase::TITLE("Grid_MeshK","cal_extended_cell");

	//--------------------------------------
	// max and min unitcell in expaned grid.
	//--------------------------------------
	this->maxu1 = dxe / nbx + 1;
	this->maxu2 = dye / nby + 1;
	this->maxu3 = dze / nbz + 1;

	this->minu1 = (-dxe+1) / nbx - 1; 
	this->minu2 = (-dye+1) / nby - 1; 
	this->minu3 = (-dze+1) / nbz - 1; 

	if(GlobalV::test_gridt)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"MaxUnitcell",maxu1,maxu2,maxu3);
	if(GlobalV::test_gridt)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"MinUnitcell",minu1,minu2,minu3);

	//--------------------------------------
	// number of unitcell in each direction.
	//--------------------------------------
	this->nu1 = maxu1 - minu1 + 1;
	this->nu2 = maxu2 - minu2 + 1;
	this->nu3 = maxu3 - minu3 + 1;
	this->nutot = nu1 * nu2 * nu3;

	init_ucell_para();
	if(GlobalV::test_gridt)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"UnitCellNumber",nu1,nu2,nu3);
	//xiaohui add 'PARAM.inp.out_level' line, 2015-09-16
	if(PARAM.inp.out_level != "m") ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"UnitCellTotal",nutot);


    this->ucell_index2x = std::vector<int>(nutot, 0);
    this->ucell_index2y = std::vector<int>(nutot, 0);
    this->ucell_index2z = std::vector<int>(nutot, 0);

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