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
//5. save_npy_s : stress
//6. save_npy_o : bandgap
//7. save_npy_orbital_precalc : orbital_precalc

#ifdef __DEEPKS

#include "LCAO_deepks.h"
#include "npy.hpp"

void LCAO_Deepks::print_dm(const std::vector<double> &dm)
{
    std::ofstream ofs("dm");
    ofs << std::setprecision(15);
    for (int mu=0;mu<GlobalV::NLOCAL;mu++)
    {
        for (int nu=0;nu<GlobalV::NLOCAL;nu++)
        {
            //ofs << dm(mu,nu) << " ";
            ofs << dm[mu * this->pv->nrow + nu] << " ";
        }
        ofs << std::endl;
    }
}

void LCAO_Deepks::print_dm_k(const int nks, const std::vector<std::vector<std::complex<double>>>& dm)
{
    std::stringstream ss;
    for(int ik=0;ik<nks;ik++)
    {
        ss.str("");
        ss<<"dm_"<<ik;
        std::ofstream ofs(ss.str().c_str());
        ofs << std::setprecision(15);

        for (int mu=0;mu<GlobalV::NLOCAL;mu++)
        {
            for (int nu=0;nu<GlobalV::NLOCAL;nu++)
            {
                //ofs << dm[ik](mu,nu) << " ";
                ofs << dm[ik][mu * this->pv->nrow + nu] << " ";
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
    if(!if_equiv)
    {
        vector<double> npy_des;
        for (int inl = 0;inl < inlmax;++inl)
        {
            int nm = 2*inl_l[inl] + 1;
            for(int im=0;im<nm;im++)
            {
                npy_des.push_back(this->d_tensor[inl].index({im}).item().toDouble());
            }
        }
        const long unsigned dshape[] = {static_cast<unsigned long>(nat), static_cast<unsigned long>(this->des_per_atom)};
        if (GlobalV::MY_RANK == 0)
        {
            npy::SaveArrayAsNumpy("dm_eig.npy", false, 2, dshape, npy_des);
        }
    }
    else
    {
        // a rather unnecessary way of writing this, but I'll do it for now
        std::vector<double> npy_des;
        for(int iat = 0; iat < nat; iat ++)
        {
            for(int i = 0; i < this->des_per_atom; i++)
            {
                npy_des.push_back(this->d_tensor[iat].index({i}).item().toDouble());
            }
        }
        const long unsigned dshape[] = {static_cast<unsigned long>(nat), static_cast<unsigned long>(this->des_per_atom)};
        if (GlobalV::MY_RANK == 0)
        {
            npy::SaveArrayAsNumpy("dm_eig.npy", false, 2, dshape, npy_des);
        }        
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
    const long unsigned gshape[]
        = {static_cast<unsigned long>(nat), 3UL, static_cast<unsigned long>(nat), static_cast<unsigned long>(this->des_per_atom)};
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

//saves gvx into grad_vepsl.npy
void LCAO_Deepks::save_npy_gvepsl(const int nat)
{
    ModuleBase::TITLE("LCAO_Deepks", "save_npy_gvepsl");
    if(GlobalV::MY_RANK!=0) return;
    //save grad_vepsl.npy (when  stress label is in use)
    //unit: none
    const long unsigned gshape[] = {6UL, static_cast<unsigned long>(nat), static_cast<unsigned long>(this->des_per_atom)};
    vector<double> npy_gvepsl;
    for (int i = 0;i < 6;i++)
    {
        for (int ibt = 0;ibt < nat;++ibt)
        {

            for(int p=0;p<this->des_per_atom;++p)
            {
                npy_gvepsl.push_back(this->gvepsl_tensor.index({ i, ibt, p }).item().toDouble());
            }

        }
    }
    npy::SaveArrayAsNumpy("grad_vepsl.npy", false, 3, gshape, npy_gvepsl);
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
    const long unsigned fshape[] = {static_cast<unsigned long>(nat), 3 };
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

//saves stress in numpy format
void LCAO_Deepks::save_npy_s(const ModuleBase::matrix &s, const std::string &s_file, const double omega)
{
    ModuleBase::TITLE("LCAO_Deepks", "save_npy_s");
    if(GlobalV::MY_RANK!=0) return;
    //save f_base
    //caution: unit: Rydberg/Bohr
    const long unsigned sshape[] = { 6 };
    vector<double> npy_s;

    for (int ipol = 0;ipol < 3;++ipol)
    {
        for (int jpol = ipol;jpol < 3;jpol++)
        {
            npy_s.push_back(s(ipol, jpol)*omega);
        }
    }
    npy::SaveArrayAsNumpy(s_file, false, 1, sshape, npy_s);
    return;
}

void LCAO_Deepks::save_npy_o(const ModuleBase::matrix &bandgap, const std::string &o_file, const int nks)
{
    ModuleBase::TITLE("LCAO_Deepks", "save_npy_o");
    if(GlobalV::MY_RANK!=0) return;
    //save o_base
    const long unsigned oshape[] = {static_cast<unsigned long>(nks), 1 };
    vector<double> npy_o;
    for (int iks = 0; iks < nks; ++iks)
    {
        for (int hl = 0;hl < 1;hl++)
        {
            npy_o.push_back(bandgap(iks,hl));
        }
    }

    npy::SaveArrayAsNumpy(o_file, false, 2, oshape, npy_o);
    return;
}

void LCAO_Deepks::save_npy_orbital_precalc(const int nat, const int nks)
{
    ModuleBase::TITLE("LCAO_Deepks", "save_npy_orbital_precalc");
    if(GlobalV::MY_RANK!=0) return;
    //save orbital_precalc.npy (when bandgap label is in use)
    //unit: a.u.
    const long unsigned gshape[] = {static_cast<unsigned long>(nks),
                                    1,
                                    static_cast<unsigned long>(nat),
                                    static_cast<unsigned long>(this->des_per_atom)};
    vector<double> npy_orbital_precalc;
    for (int iks = 0; iks < nks; ++iks)
    {
        for (int hl = 0; hl < 1; ++hl)
        {
            for (int iat = 0;iat < nat;++iat)
            {
                for(int p=0; p<this->des_per_atom; ++p)
                {
                    npy_orbital_precalc.push_back(this->orbital_precalc_tensor.index({iks, hl, iat, p }).item().toDouble());
                }
            }
        }
    }
    npy::SaveArrayAsNumpy("orbital_precalc.npy", false, 4, gshape, npy_orbital_precalc);
    return;
}

#endif