#include "symmetry.h"
using namespace ModuleSymmetry;

void Symmetry::getgroup(int& nrot, int& nrotk, std::ofstream& ofs_running, 
		const int& nop, const ModuleBase::Matrix3* symop, ModuleBase::Matrix3* gmatrix, 
		ModuleBase::Vector3<double>* gtrans, double* pos, double* rotpos, 
		int* index, const int ntype, const int itmin_type, 
		const int itmin_start, int* istart, int* na)const
{
    ModuleBase::TITLE("Symmetry", "getgroup");

	//--------------------------------------------------------------------------------
    //return all possible space group operators that reproduce a lattice with basis
    //out of a (maximum) pool of point group operations that is compatible with
    //the symmetry of the pure translation lattice without any basic.
	//--------------------------------------------------------------------------------

    ModuleBase::Matrix3 zero(0,0,0,0,0,0,0,0,0);
    ModuleBase::Matrix3 help[48];
    ModuleBase::Vector3<double> temp[48];

    nrot = 0;
    nrotk = 0;

	//-------------------------------------------------------------------------
    //pass through the pool of (possibly allowed) symmetry operations and
    //check each operation whether it can reproduce the lattice with basis
	//-------------------------------------------------------------------------
    //std::cout << "nop = " <<nop <<std::endl;
    for (int i = 0; i < nop; ++i)
    {
        bool s_flag = this->checksym(symop[i], gtrans[i], pos, rotpos, index, ntype, itmin_type, itmin_start, istart, na);
        if (s_flag == 1)
        {
			//------------------------------
            // this is a symmetry operation
			// with no translation vectors
            // so ,this is pure point group 
			// operations
			//------------------------------
            if ( equal(gtrans[i].x,0.0) &&
                 equal(gtrans[i].y,0.0) &&
                 equal(gtrans[i].z,0.0))
            {
                ++nrot;
                gmatrix[nrot - 1] = symop[i];
                gtrans[nrot - 1].x = 0;
                gtrans[nrot - 1].y = 0;
                gtrans[nrot - 1].z = 0;
            }
			//------------------------------
            // this is a symmetry operation
			// with translation vectors
            // so ,this is space group 
			// operations
			//------------------------------
            else
            {
                ++nrotk;
                help[nrotk - 1] = symop[i];
                temp[nrotk - 1].x = gtrans[i].x;
                temp[nrotk - 1].y = gtrans[i].y;
                temp[nrotk - 1].z = gtrans[i].z;
            }
        }
    }

	//-----------------------------------------------------
    //If there are operations with nontrivial translations
    //then store them together in the momory
	//-----------------------------------------------------
    if (nrotk > 0)
    {
        for (int i = 0; i < nrotk; ++i)
        {
            gmatrix[nrot + i] = help[i];
            gtrans[nrot + i].x = temp[i].x;
            gtrans[nrot + i].y = temp[i].y;
            gtrans[nrot + i].z = temp[i].z;
        }
    }

	//-----------------------------------------------------
    //total number of space group operations
	//-----------------------------------------------------
    nrotk += nrot;

    if(test_brav)
    {
	    ModuleBase::GlobalFunc::OUT(ofs_running,"PURE POINT GROUP OPERATIONS",nrot);
        ModuleBase::GlobalFunc::OUT(ofs_running,"SPACE GROUP OPERATIONS",nrotk);
    }

	//-----------------------------------------------------
    //fill the rest of matrices and vectors with zeros
	//-----------------------------------------------------
    if (nrotk < 48)
    {
        for (int i = nrotk; i < 48; ++i)
        {
            gmatrix[i] = zero;
            gtrans[i].x = 0;
            gtrans[i].y = 0;
            gtrans[i].z = 0;
        }
    }
    return;
}
