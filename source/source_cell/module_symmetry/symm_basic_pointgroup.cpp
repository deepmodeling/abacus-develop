#include "symmetry.h"
#include "source_base/mymath.h"

namespace ModuleSymmetry
{
int Symmetry_Basic::subgroup(const int& nrot, const int& ninv,
    const int& nc2, const int& nc3, const int& nc4, const int& nc6,
    const int& ns1, const int& ns3, const int& ns4, const int& ns6)const
{
    if (nrot > 24)
    {
        if (ninv)
        {
            // if (nc2 >= 7 && nc3 >= 2 && nc6 >= 2 && ns1 >= 7 && ns3 >= 2 && ns6 >= 2) { return 27; } //D_6h
            if (nc2 >= 3 && nc3 >= 8 && ns1 >= 3 && ns6 >= 8) { return 29; } //T_h
        }
        else
        {
            if (nc2 >= 9 && nc3 >= 8 && nc4 >= 6) { return 30; } //O
            if (nc2 >= 3 && nc3 >= 8 && ns1 >= 6 && ns4 >= 6) { return 31; } //T_d
        }
    }
    if (nrot > 16)//not else if: nrot>24 can also fall in this part and below
    {
        if (ninv && nc2 >= 5 && nc4 >= 2 && ns1 >= 5 && ns4 >= 2) { return 20; } //D_4h
    }
    if (nrot > 12)
    {
        if (ninv)
        {
            if (nc2 >= 1 && nc3 >= 2 && nc6 >= 2 && ns1 >= 1 && ns3 >= 2 && ns6 >= 2) { return 23; } //C_6h
            if (nc2 >= 3 && nc3 >= 2 && ns1 >= 3 && ns6 >= 2) { return 13; } //D_3d
        }
        else
        {
            if (nc2 >= 3 && nc3 >= 8) { return 28; } //T
            if (nc2 >= 3 && nc3 >= 2 && ns1 >= 4 && ns3 >= 2) { return 26; } //D_3h
            if (nc2 >= 1 && nc3 >= 2 && nc6 >= 2 && ns1 >= 6) { return 25; } //C_6v
            if (nc2 >= 7 && nc3 >= 2 && nc6 >= 2) { return 24; } //D_6
        }
    }
    if (nrot > 8)
    {
        if (ninv)
        {
            if (nc2 >= 1 && nc4 >= 2 && ns1 >= 1 && ns4 >= 2) { return 16; } //C_4h
            if (nc2 >= 3 && ns1 >= 3) { return 8; } //D_2h
        }
        else
        {
            if (nc2 >= 3 && ns1 >= 2 && ns4 >= 2) { return 19; } //D_2d
            if (nc2 >= 1 && nc4 >= 2 && ns1 >= 4) { return 18; } //C_4v
            if (nc2 >= 5 && nc4 >= 2) { return 17; } //D_4
        }
    }
    if (nrot > 6)
    {
        if (nc3 >= 2 && ns1 >= 1 && ns3 >= 2) { return 22; } //C_3h
        if (nc2 >= 1 && nc3 >= 2 && nc6 >= 2) { return 21; } //C_6
        if (nc3 >= 2 && ns1 >= 3) { return 12; } //C_3v
        if (nc2 >= 3 && nc3 >= 2) { return 11; } //D_3
        if (ninv && nc3 >= 2 && ns3 >= 2) { return 10; }//S_6
    }
    if (nrot > 4)
    {
        if (nc2 >= 1 && ns4 >= 2) { return 15; } //S_4
        if (nc2 >= 1 && nc4 >= 2) { return 14; } //C_4
        if (nc2 >= 1 && ns1 >= 2) { return 7; } //C_2v
        if (nc2 >= 3) { return 6; } //D_2
        if (ninv && nc2 >= 1 && ns1 >= 1) { return 5; } //C_2h
    }
    if (nrot > 3)
    {
        if (nc3 >= 2) { return 9; } //C_3
    }
    if (nrot > 2)
    {
        if (ns1 >= 1) { return 4; } //C_1h
        if (nc2 >= 1) { return 3; } //C_2
        if (ninv) { return 2; } //S_2
    }
    return 1; //C_1
}


bool Symmetry_Basic::pointgroup(const int& nrot, int& pgnumber,
  std::string& pgname, const ModuleBase::Matrix3* gmatrix, std::ofstream& ofs_running,
  const int* cal_symm_repr)const
{
    //-------------------------------------------------------------------------
    //return the name of the point group
    //the "name" (Schoenflies mark) of the group defined by following key:
    //       1 --> C_1       9 --> C_3      17 --> D_4      25 --> C_6v     *
    //       2 --> S_2      10 --> S_6      18 --> C_4v     26 --> D_3h     *
    //       3 --> C_2      11 --> D_3      19 --> D_2d     27 --> D_6h     *
    //       4 --> C_1h     12 --> C_3v     20 --> D_4h     28 --> T        *
    //       5 --> C_2h     13 --> D_3d     21 --> C_6      29 --> T_h      *
    //       6 --> D_2      14 --> C_4      22 --> C_3h     30 --> O        *
    //       7 --> C_2v     15 --> S_4      23 --> C_6h     31 --> T_d      *
    //       8 --> D_2h     16 --> C_4h     24 --> D_6      32 --> O_h      *
    //-------------------------------------------------------------------------

    //there are four trivial cases which could be easily determined
    //because the number of their elements are exclusive
    if (cal_symm_repr != nullptr && cal_symm_repr[0] > 1) {
        ModuleBase::TITLE("Symmetry_Basic", "pointgroup");
    }

    std::vector<std::string> pgdict = { "none", "C_1", "S_2", "C_2", "C_1h", "C_2h",
    "D_2", "C_2v", "D_2h", "C_3", "S_6", "D_3", "C_3v", "D_3d", "C_4", "S_4", "C_4h",
    "D_4", "C_4v", "D_2d", "D_4h", "C_6", "C_3h", "C_6h", "D_6", "C_6v", "D_3h", "D_6h",
    "T", "T_h", "O", "T_d", "O_h" };

    if(nrot == 1)
    {
        pgnumber = 1;
        pgname="C_1";
        return true;
    }
    if(nrot == 3)
    {
        pgnumber = 9;
        pgname="C_3";
        return true;
    }
    if(nrot == 16)
    {
        pgnumber = 20;
        pgname="D_4h";
        return true;
    }
    if(nrot == 48)
    {
        pgnumber = 32;
        pgname="O_h";
        return true;
    }

    //-------------------------------------------------------------------------------
    //all other groups need further investigations and detailed analysis
    //first determine the type of elements and count them
    //Possible elements are E, I, C_2, C_3, C_4, C_6 and S_1, S_3, S_4, S_6 (S_1 = m)
    //The type of a symmetry operation can be identified simply by
    //calculating the trace and the determinant of the rotation matrix. The
    //combination of these two quantities is specific for specific elements:
    //-------------------------------------------------------------------------------

    // Element:         E    I  C_2  C_3  C_4  C_6  S_1  S_6  S_4  S_3
    // Trace:          +3   -3   -1    0   +1   +2   +1    0   -1   -2
    // Determinant:    +1   -1   +1   +1   +1   +1   -1   -1   -1   -1

    int trace = 0;
    int det = 0;
    int ninv = 0;

    int nc2 = 0;
    int nc3 = 0;
    int nc4 = 0;
    int nc6 = 0;
    int ns1 = 0;
    int ns3 = 0; //mohan add 2012-01-15
    int ns4 = 0;
    int ns6 = 0; //mohan add 2012-01-15

    for(int i = 0; i < nrot; ++i)
    {
        //calculate the trace of a matrix
        trace = int(gmatrix[i].e11+gmatrix[i].e22+gmatrix[i].e33);
        //calculate the determinant of a matrix
        det = int(gmatrix[i].Det());

        if(trace == 3)
        {
            continue;    //found unity operator (trivial)
        }
        //found inversion
        if(trace == -3)
        {
            ninv = 1;
            continue;
        }

        if(trace == -1 && det == 1)       { ++nc2; }
        else if(trace == 0 && det == 1)   { ++nc3; }
        else if(trace == 1 && det == 1)   { ++nc4; }
        else if(trace == 2 && det == 1)   { ++nc6; }
        else if(trace == 1 && det == -1)  { ++ns1; }
        else if(trace == 0 && det == -1)  { ++ns6; } //mohan add 2012-01-15
        else if(trace == -1 && det == -1) { ++ns4; }
        else if(trace == -2 && det == -1) { ++ns3; } //mohan add 2012-01-15
    }

    if(test_brav)
    {
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "C2", nc2);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "C3", nc3);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "C4", nc4);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "C6", nc6);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "S1", ns1);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "S3", ns3);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "S4", ns4);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "S6", ns6);
    }

    if(nrot == 2)
    {
        if(ninv == 1)
        {
            pgnumber = 2;
            pgname="S_2";
            return true;
        }
        if(nc2 == 1)
        {
            pgnumber = 3;
            pgname="C_2";
            return true;
        }
        if(ns1 == 1)
        {
            pgnumber = 4;
            pgname="C_1h";
            return true;
        }
    }
    if(nrot == 4)
    {
        if(ninv == 1)
        {
            pgnumber = 5;
            pgname="C_2h";
            return true;
        }
        if(nc2 == 3)
        {
            pgnumber = 6;
            pgname="D_2";
            return true;
        }
        if(ns1 == 2)
        {
            pgnumber = 7;
            pgname="C_2v";
            return true;
        }
        if(nc4 == 2)
        {
            pgnumber = 14;
            pgname="C_4";
            return true;
        }
        if(ns4 == 2)
        {
            pgnumber = 15;
            pgname="S_4";
            return true;
        }
    }
    if(nrot == 6)
    {
        if(ninv == 1)
        {
            pgnumber = 10;
            pgname="S_6";
            return true;
        }
        if(nc2 == 3)
        {
            pgnumber = 11;
            pgname="D_3";
            return true;
        }
        if(ns1 == 3)
        {
            pgnumber = 12;
            pgname="C_3v";
            return true;
        }
        if(nc2 == 1)
        {
            pgnumber = 21;
            pgname="C_6";
            return true;
        }
        if(ns1 == 1)
        {
            pgnumber = 22;
            pgname="C_3h";
            return true;
        }
    }
    if(nrot == 8)
    {
        if(ns1 == 3)
        {
            pgnumber = 8;
            pgname="D_2h";
            return true;
        }
        if(ns1 == 1)
        {
            pgnumber = 16;
            pgname="C_4h";
            return true;
        }
        if(ns1 == 0)
        {
            pgnumber = 17;
            pgname="D_4";
            return true;
        }
        if(ns1 == 4)
        {
            pgnumber = 18;
            pgname="C_4v";
            return true;
        }
        if(ns1 == 2)
        {
            pgnumber = 19;
            pgname="D_2d";
            return true;
        }
    }
    if(nrot == 12)
    {
        if(ns1 == 3)
        {
            pgnumber = 13;
            pgname="D_3d";
            return true;
        }
        if(ns1 == 1)
        {
            pgnumber = 23;
            pgname="C_6h";
            return true;
        }
        if(nc2 == 7)
        {
            pgnumber = 24;
            pgname="D_6";
            return true;
        }
        if(ns1 == 6)
        {
            pgnumber = 25;
            pgname="C_6v";
            return true;
        }
        if(ns1 == 4)
        {
            pgnumber = 26;
            pgname="D_3h";
            return true;
        }
        if(nc3 == 8)
        {
            pgnumber = 28;
            pgname="T";
            return true;
        }
    }
    if(nrot == 24)
    {
        if(nc6 == 2)
        {
            pgnumber = 27;
            pgname="D_6h";
            return true;
        }
        if(ninv == 1)
        {
            pgnumber = 29;
            pgname="T_h";
            return true;
        }
        if(nc4 == 6)
        {
            pgnumber = 30;
            pgname="O";
            return true;
        }
        if(ns4 == 6)
        {
            pgnumber = 31;
            pgname="T_d";
            return true;
        }
    }
    GlobalV::ofs_running << "\n WARNING: Symmetry operations cannot completely constitute a point group.\n\
    It'll be better to try another `symmetry_prec`.\n  Now search the subgroups ..." << std::endl;
    pgnumber = this->subgroup(nrot, ninv, nc2, nc3, nc4, nc6, ns1, ns3, ns4, ns6);
    pgname = pgdict[pgnumber];
    return false;
}
}
