#include "symmetry.h"
#include "source_base/mymath.h"
#include "source_base/formatter.h"

namespace ModuleSymmetry
{
void Symmetry_Basic::matrigen(ModuleBase::Matrix3 *symgen, const int ngen, ModuleBase::Matrix3* symop, int &nop) const
{
    int m1 = 0;
    int m2 = 0;
    int n = 0;

    // allocate memory for the symmetry operations
    ModuleBase::Matrix3  iden(1,0,0,0,1,0,0,0,1);
    ModuleBase::Matrix3   sig(1,0,0,0,1,0,0,0,1);
    ModuleBase::Matrix3 temp1(1,0,0,0,1,0,0,0,1);
    ModuleBase::Matrix3 temp2(1,0,0,0,1,0,0,0,1);

    bool flag = false; // mark whether the symmetry operation is a new one
    int order = 0;
    int now = 0;

    symop[0] = iden; //identity (the trivial element)
    nop = 1; // counter of the symmetry operations

    // take all generators
    for (int i = 0; i < ngen; ++i)
    {
        sig = symgen[i];
        flag = true; // assume it is a new symmetry operation
        // search if the symmetry operation already exists among the found symmetry operations
        // if so, skip it
        for (int j = 0; j < nop; ++j)
        {
            if (symop[j] == sig)
            {
                flag = 0; // not a new symmetry operation
                break;
            }
        }
        if (flag == 0) // if old, return
        {
            continue;
        }
        // otherwise

        // determine the order of the operation: by which power will the operation return
        // to the identity operation.
        temp1 = sig;
        for (int j = 1; j < 100; ++j)
        {
            order = j;
            if (temp1 == iden)
            {
                break;
            }
            temp1 = sig * temp1;
        }
        now = nop;
        for (int j = 0; j < nop; ++j)
        {
            temp1 = symop[j];
            for (int k = 1; k < order; ++k)
            {
                temp1 = sig * temp1;

                for (int l = 0; l < nop; ++l)
                {
                    temp2 = symop[l] * temp1;
                    flag = 1;
                    for (int m = 0; m < now; ++m)
                    {
                        if (symop[m] == temp2)
                        {
                            flag = 0;
                            break;
                        }
                    }
                    if (flag == 0)
                    {
                        continue;    //the newly-found element has already existed.
                    }

                    ++now; // the number of elements we found
                    if (now > 48) // number of symm_op cannot be more than 48 (of O_h point group)
                    {
                        std::cout << "\n a: now= "<<now<<std::endl;
                        std::cout << "\n There are too many symmetrical matrices!"<<std::endl;
                        return;
                    }
                    symop[now - 1] = temp2;
                }
            }
            if(j == 0)
            {
                n = now;
            }
        }

        m1 = nop;
        m2 = now;
        for (int j = 1; j < 50; ++j)
        {
            for (int k = nop; k < n; ++k)
            {
                for (int m = m1; m < m2; ++m)
                {
                    temp1 = symop[k] * symop[m];
                    flag = 1;
                    for (int l = 0; l < now; ++l)
                    {
                        if (symop[l] == temp1)
                        {
                            flag = 0;
                            break;
                        }
                    }
                    if (flag == 0)
                    {
                        continue;    //the new-found element has already existed
                    }

                    ++now;
                    if (now > 48)
                    {
                        std::cout << "\n b: now= "<<now<<std::endl;
                        std::cout << "\n There are too many symmetrical matrices!"<<std::endl;
                        return;
                    }
                    symop[now - 1] = temp1;
                }
            }
            if (now == m2)
            {
                break;    //if no more new element could be found, stop the loop
            }
            m1 = m2;
            m2 = now;
        }
        nop = now;
    }
}

//--------------------------------------------------------------
// set up all possible space group operators
// (integer rotation matrices and nontrivial translations
// given in crystal coordinates)
// of a lattice with some arbitrary basis (atomic arrangement).
//--------------------------------------------------------------
void Symmetry_Basic::setgroup(ModuleBase::Matrix3* symop, int &nop, const int &ibrav,
                               const int* cal_symm_repr) const
{
    if(cal_symm_repr != nullptr && cal_symm_repr[0] > 1) {
        ModuleBase::TITLE("Symmetry_Basic", "setgroup");
    }
    ModuleBase::Matrix3 symgen[3]; // the number of generators is up to 3

    ModuleBase::Matrix3    inv(-1, 0, 0, 0,-1, 0, 0, 0,-1); // (x, y, z) -> (-x, -y, -z)
    ModuleBase::Matrix3    r3d( 0, 1, 0, 0, 0, 1, 1, 0, 0); // (x, y, z) -> (y, z, x)
    ModuleBase::Matrix3    r6z( 1, 1, 0,-1, 0, 0, 0, 0, 1); // (x, y, z) -> (x+y, -x, z)
    ModuleBase::Matrix3  r2hex( 1, 0, 0,-1,-1, 0, 0, 0,-1); // (x, y, z) -> (x, -x-y, -z)
    ModuleBase::Matrix3  r2tri(-1, 0, 0, 0, 0,-1, 0,-1, 0); // (x, y, z) -> (-x, -z, -y)
    ModuleBase::Matrix3   r4zp( 0, 1, 0,-1, 0, 0, 0, 0, 1); // (x, y, z) -> (y, -x, z)
    ModuleBase::Matrix3   r2yp(-1, 0, 0, 0, 1, 0, 0, 0,-1); // (x, y, z) -> (-x, y, -z)
    ModuleBase::Matrix3  r4zbc( 0, 0,-1, 1, 1, 1, 0,-1, 0); // (x, y, z) -> (-z, x+y+z, -y)
    ModuleBase::Matrix3  r4zfc( 1, 0,-1, 1, 0, 0, 1,-1, 0); // (x, y, z) -> (x-z, x, x-y)
    ModuleBase::Matrix3   r2zp(-1, 0, 0, 0,-1, 0, 0, 0, 1); // (x, y, z) -> (-x, -y, z)
    ModuleBase::Matrix3  r2ybc( 0, 0, 1,-1,-1,-1, 1, 0, 0); // (x, y, z) -> (z, -x-y-z, x)
    ModuleBase::Matrix3  r2zbc( 0, 1, 0, 1, 0, 0,-1,-1,-1); // (x, y, z) -> (y, x, -x-y-z)
    ModuleBase::Matrix3 r2ybas( 0,-1, 0,-1, 0, 0, 0, 0,-1); // (x, y, z) -> (-y, -x, -z)
    ModuleBase::Matrix3  r2yfc( 0,-1, 1, 0,-1, 0, 1,-1, 0); // (x, y, z) -> (-y+z, -y, x-y)
    ModuleBase::Matrix3  r2zfc( 0, 1,-1, 1, 0,-1, 0, 0,-1); // (x, y, z) -> (y-z, x-z, -z)

    //the pure translation lattice (bravais lattice) has some maximum symmetry
    //set first up the point group operations for this symmetry.
    symgen[0] = inv;
    // generate the point group operations for the bravais lattice
    // rewrite with switch-case to get better performance and readability
    switch (ibrav) {
        case 1:
            symgen[1] = r3d;
            symgen[2] = r4zp;
            this->matrigen(symgen, 3, symop, nop);
        break;
        case 2:
            symgen[1] = r3d;
            symgen[2] = r4zbc;
            this->matrigen(symgen, 3, symop, nop);
        break;
        case 3:
            symgen[1] = r3d;
            symgen[2] = r4zfc;
            this->matrigen(symgen, 3, symop, nop);
        break;
        case 4:
            symgen[1] = r6z;
            symgen[2] = r2hex;
            this->matrigen(symgen, 3, symop, nop);
        break;
        case 5:
            symgen[1] = r4zp;
            symgen[2] = r2yp;
            this->matrigen(symgen, 3, symop, nop);
        break;
        case 6:
            symgen[1] = r4zbc;
            symgen[2] = r2ybc;
            this->matrigen(symgen, 3, symop, nop);
        break;
        case 7:
            symgen[1] = r2tri;
            symgen[2] = r3d;
            this->matrigen(symgen, 3, symop, nop);
        break;
        case 8:
            symgen[1] = r2zp;
            symgen[2] = r2yp;
            this->matrigen(symgen, 3, symop, nop);
        break;
        case 9:
            symgen[1] = r2zbc;
            symgen[2] = r2ybc;
            this->matrigen(symgen, 3, symop, nop);
        break;
        case 10:
            symgen[1] = r2zfc;
            symgen[2] = r2yfc;
            this->matrigen(symgen, 3, symop, nop);
        break;
        case 11:
            symgen[1] = r2zp;
            symgen[2] = r2ybas;
            this->matrigen(symgen, 3, symop, nop);
        break;
        case 12:
            symgen[1] = r2yp;
            this->matrigen(symgen, 2, symop, nop);
        break;
        case 13:
            symgen[1] = r2ybas;
            this->matrigen(symgen, 2, symop, nop);
        break;
        case 14:
            this->matrigen(symgen, 1, symop, nop);
        break;
        default:
            ModuleBase::WARNING_QUIT("Symmetry_Basic::setgroup",
                "ibrav = " + std::to_string(ibrav) + " is not supported.");
        break;
    }

    // print
    if (test_brav)
    {
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Number of rotation matrices", nop);
    }

    // print the symmetry operations
    if (cal_symm_repr != nullptr && cal_symm_repr[0] > 0)
    {
        GlobalV::ofs_running << std::endl
                             << " ======================================================================\n"
                             << " MATRIX REPRESENTATION OF SYMMETRY OPERATION\n"
                             << " ======================================================================\n"
                             << " There are " << nop << " symmetry operation representation matrices.\n"
                             << " For each matrix, the elements are arranged like: \n"
                             << " [[e11, e12, e13], [e21, e22, e23], [e31, e32, e33]].reshape(3, 3)\n"
                             << std::endl;

        // control the digits
        const int precision = cal_symm_repr[1];
        const int width = precision + 4;
        std::string fmtstr = " %" + std::to_string(width) + "." + std::to_string(precision) + "f";
        fmtstr += fmtstr + fmtstr + "\n";

        // print the symmetry operations
        std::string mat;
        for (int i = 0; i < nop; ++i)
        {
            mat = " " + FmtCore::format("No. %3d", i + 1) + "\n"
                + FmtCore::format(fmtstr.c_str(), symop[i].e11, symop[i].e12, symop[i].e13)
                + FmtCore::format(fmtstr.c_str(), symop[i].e21, symop[i].e22, symop[i].e23)
                + FmtCore::format(fmtstr.c_str(), symop[i].e31, symop[i].e32, symop[i].e33);
            GlobalV::ofs_running << mat << std::endl;
        }
        GlobalV::ofs_running << " ======================================================================\n";
    }

    return;
}
}
