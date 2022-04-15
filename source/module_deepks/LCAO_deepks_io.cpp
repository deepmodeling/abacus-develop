//wenfei 2022-1-11
//This file contains subroutines that contains interface with libnpy
//since many arrays must be saved in numpy format
//It also contains subroutines for printing density matrices
//which is used in unit tests

//There are 2 subroutines for printing density matrices:
//1. print_dm : for gamma only
//2. print_dm_k : for multi-k

//There are 4 subroutines in this file that prints to npy file:
//1. save_npy_d : descriptor ->dm_eig.npy
//2. save_npy_gvx : gvx ->grad_vx.npy
//3. save_npy_e : energy
//4. save_npy_f : force
#ifdef __DEEPKS

#include "LCAO_deepks.h"
#include "npy.hpp"

void LCAO_Deepks::print_dm(const ModuleBase::matrix &dm)
{
    ofstream ofs("dm");
    ofs << std::setprecision(15);
    for (int mu=0;mu<GlobalV::NLOCAL;mu++)
    {
        for (int nu=0;nu<GlobalV::NLOCAL;nu++)
        {
            ofs << dm(mu,nu) << " ";
        }
        ofs << std::endl;
    }
}

void LCAO_Deepks::print_dm_k(const int nks, const std::vector<ModuleBase::ComplexMatrix>& dm)
{
    stringstream ss;
    for(int ik=0;ik<nks;ik++)
    {
        ss.str("");
        ss<<"dm_"<<ik;
        ofstream ofs(ss.str().c_str());
        ofs << std::setprecision(15);
        
        for (int mu=0;mu<GlobalV::NLOCAL;mu++)
        {
            for (int nu=0;nu<GlobalV::NLOCAL;nu++)
            {
                ofs << dm[ik](mu,nu) << " ";
            }
            ofs << std::endl;
        }
    }
}

//saves descriptor into dm_eig.npy
void LCAO_Deepks::save_npy_d(const int nat)
{
    ModuleBase::TITLE("LCAO_Deepks", "save_npy_d");
    if(GlobalV::MY_RANK!=0) return;
    //save descriptor in .npy format
    vector<double> npy_des;
    for (int inl = 0;inl < inlmax;++inl)
    {
        int nm = 2*inl_l[inl] + 1;
        for(int im=0;im<nm;im++)
        {
            npy_des.push_back(this->d_tensor[inl].index({im}).item().toDouble());
        }
    }
    const long unsigned dshape[] = {(long unsigned) nat, (long unsigned) this->des_per_atom };
    if (GlobalV::MY_RANK == 0)
    {
        npy::SaveArrayAsNumpy("dm_eig.npy", false, 2, dshape, npy_des);
    }
    return;
}

//saves gvx into grad_vx.npy
void LCAO_Deepks::save_npy_gvx(const int nat)
{
    ModuleBase::TITLE("LCAO_Deepks", "save_npy_gvx");
    if(GlobalV::MY_RANK!=0) return;
    //save grad_vx.npy (when  force label is in use)
    //unit: /Bohr
    const long unsigned gshape[] = {(long unsigned) nat, 3, nat, this->des_per_atom};
    vector<double> npy_gvx;
    for (int ibt = 0;ibt < nat;++ibt)
    {
        for (int i = 0;i < 3;i++)
        {
            for (int iat = 0;iat < nat;++iat)
            {
                for(int p=0;p<this->des_per_atom;++p)
                {
                    npy_gvx.push_back(this->gvx_tensor.index({ ibt, i, iat, p }).item().toDouble());
                }
            }
        }
    }
    npy::SaveArrayAsNumpy("grad_vx.npy", false, 4, gshape, npy_gvx);
    return;
}

//saves energy in numpy format
void LCAO_Deepks::save_npy_e(const double &e, const std::string &e_file)
{
    ModuleBase::TITLE("LCAO_Deepks", "save_npy_e");
    if(GlobalV::MY_RANK!=0) return;
    //save e_base
    const long unsigned eshape[] = { 1 };
    vector<double> npy_e;
    npy_e.push_back(e);
    npy::SaveArrayAsNumpy(e_file, false, 1, eshape, npy_e);
    return;
}

//saves force in numpy format
void LCAO_Deepks::save_npy_f(const ModuleBase::matrix &f, const std::string &f_file, const int nat)
{
    ModuleBase::TITLE("LCAO_Deepks", "save_npy_f");
    if(GlobalV::MY_RANK!=0) return;
    //save f_base
    //caution: unit: Rydberg/Bohr
    const long unsigned fshape[] = {(long unsigned) nat, 3 };
    vector<double> npy_f;
    for (int iat = 0;iat < nat;++iat)
    {
        for (int i = 0;i < 3;i++)
        {
            npy_f.push_back(f(iat, i));
        }
    }
    npy::SaveArrayAsNumpy(f_file, false, 2, fshape, npy_f);
    return;
}

void LCAO_Deepks::save_npy_o(const double &bandgap, const std::string &o_file)
{
    ModuleBase::TITLE("LCAO_Deepks", "save_npy_o");
    if(GlobalV::MY_RANK!=0) return;
    //save o_base
    const long unsigned oshape[] = { 1 };
    vector<double> npy_o;
    npy_o.push_back(bandgap);
    npy::SaveArrayAsNumpy(o_file, false, 1, oshape, npy_o);
    return;
}

void LCAO_Deepks::save_npy_orbital_precalc(const int nat)
{
    ModuleBase::TITLE("LCAO_Deepks", "save_npy_orbital_precalc");
    if(GlobalV::MY_RANK!=0) return;
    //save orbital_precalc.npy (when bandgap label is in use)
    //unit: a.u.
    const long unsigned gshape[] = {(long unsigned) 1, nat, this->des_per_atom};
    vector<double> npy_orbital_precalc;
    for (int hl = 0;hl < 1; ++hl)
    {
        
        for (int iat = 0;iat < nat;++iat)
        {
            for(int p=0; p<this->des_per_atom; ++p)
            {
                npy_orbital_precalc.push_back(this->orbital_precalc_tensor.index({ hl, iat, p }).item().toDouble());
            }
        }
        
    }
    npy::SaveArrayAsNumpy("orbital_precalc.npy", false, 3, gshape, npy_orbital_precalc);
    return;
}

#endif