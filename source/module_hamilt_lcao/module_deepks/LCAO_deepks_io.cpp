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

void LCAO_Deepks::print_dm(const ModuleBase::matrix &dm)
{
    std::ofstream ofs("dm");
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
    const long unsigned dshape[] = {static_cast<unsigned long>(nat), static_cast<unsigned long>(this->des_per_atom)};
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

//just for gamma only
void LCAO_Deepks::save_npy_h(const ModuleBase::matrix &H,const std::string &h_file,const int nlocal)
{
    ModuleBase::TITLE("LCAO_Deepks", "save_npy_h");
    int nks=1;
    const long unsigned hshape[] = {static_cast<unsigned long>(nks),static_cast<unsigned long>(nlocal), static_cast<unsigned long>(nlocal) };
    vector<double> npy_h;
    for(int k=0; k<nks; k++)
    {
        for (int i=0; i<nlocal; i++)
        {
            for (int j=0; j<nlocal; j++)
            {
                npy_h.push_back(H(i,j));
            }
        }         
    }
    npy::SaveArrayAsNumpy(h_file, false, 3, hshape, npy_h);
    return;    
}

void LCAO_Deepks::save_npy_v_delta_precalc(const int nat, const int nks,const int nlocal)
{
    ModuleBase::TITLE("LCAO_Deepks", "save_npy_v_delta_precalc");
    if(GlobalV::MY_RANK!=0) return;
    // timeval t_start;
    // gettimeofday(&t_start,NULL);
    //save v_delta_precalc.npy (when v_delta label is in use)
    //unit: a.u.
    const long unsigned gshape[] = {static_cast<unsigned long>(nks),
                                    static_cast<unsigned long>(nlocal),
                                    static_cast<unsigned long>(nlocal),
                                    static_cast<unsigned long>(nat),
                                    static_cast<unsigned long>(this->des_per_atom)};
    vector<double> npy_v_delta_precalc;
    for (int iks = 0; iks < nks; ++iks)
    {
        for (int mu = 0; mu < nlocal; ++mu)
        {
            for (int nu = 0; nu < nlocal; ++nu)
            {
                for (int iat = 0;iat < nat;++iat)
                {
                    for(int p=0; p<this->des_per_atom; ++p)
                    {
                        npy_v_delta_precalc.push_back(this->v_delta_precalc_tensor.index({iks, mu, nu, iat, p }).item().toDouble());
                    }
                }                
            }
        }
    }
    npy::SaveArrayAsNumpy("v_delta_precalc.npy", false, 5, gshape, npy_v_delta_precalc);
    // timeval t_end;
    // gettimeofday(&t_end,NULL);
    // std::cout<<"save v_delta_precalc time:\t"<<(double)(t_end.tv_sec-t_start.tv_sec) + (double)(t_end.tv_usec-t_start.tv_usec)/1000000.0<<std::endl;
    return;
}

void LCAO_Deepks::save_npy_psialpha(const int nat, const int nks,const int nlocal)
{
    ModuleBase::TITLE("LCAO_Deepks", "save_npy_psialpha");
    if(GlobalV::MY_RANK!=0) return;
    // timeval t_start;
    // gettimeofday(&t_start,NULL);
    //save psialpha.npy (when v_delta label == 2)
    //unit: a.u.
    int nlmax = this->inlmax/nat;
    int mmax = 2*this->lmaxd+1;
    const long unsigned gshape[] = {static_cast<unsigned long>(nat),
                                    static_cast<unsigned long>(nlmax),
                                    static_cast<unsigned long>(nks),
                                    static_cast<unsigned long>(nlocal),
                                    static_cast<unsigned long>(mmax)};
    vector<double> npy_psialpha;
    for(int iat=0; iat< nat ; iat++) 
    {
        for(int nl = 0; nl < nlmax; nl++)
        {
            for (int iks = 0; iks < nks ; iks++)
            {
                for(int mu = 0; mu < nlocal ; mu++)
                {
                    for(int m=0; m< mmax; m++)
                    {
                        npy_psialpha.push_back(this->psialpha_tensor.index({ iat,nl, iks, mu, m }).item().toDouble());
                    }
                }                
            }
        }
    }
    npy::SaveArrayAsNumpy("psialpha.npy", false, 5, gshape, npy_psialpha);
    // timeval t_end;
    // gettimeofday(&t_end,NULL);
    // std::cout<<"save psialpha time:\t"<<(double)(t_end.tv_sec-t_start.tv_sec) + (double)(t_end.tv_usec-t_start.tv_usec)/1000000.0<<std::endl;
    return;
}

void LCAO_Deepks::save_npy_gevdm(const int nat)
{
    ModuleBase::TITLE("LCAO_Deepks", "save_npy_gevdm");
    if(GlobalV::MY_RANK!=0) return;
    // timeval t_start;
    // gettimeofday(&t_start,NULL);
    //save grad_evdm.npy (when v_delta label == 2)
    //unit: a.u.
    int nlmax = this->inlmax/nat;
    int mmax = 2*this->lmaxd+1;
    const long unsigned gshape[] = {static_cast<unsigned long>(nat),
                                    static_cast<unsigned long>(nlmax),
                                    static_cast<unsigned long>(mmax),
                                    static_cast<unsigned long>(mmax),
                                    static_cast<unsigned long>(mmax)};
    vector<double> npy_gevdm;
    for(int iat=0; iat< nat ; iat++)
    {
        for(int nl = 0; nl < nlmax; nl++)
        {
            for(int v=0; v< mmax; v++)
            {
                for(int m=0; m< mmax; m++)
                {
                    for(int n=0; n< mmax; n++)
                    {
                        npy_gevdm.push_back(this->gevdm_tensor.index({ iat, nl, v, m, n}).item().toDouble());
                    }         
                }
            }                
        }
    }
    npy::SaveArrayAsNumpy("grad_evdm.npy", false, 5, gshape, npy_gevdm);
    // timeval t_end;
    // gettimeofday(&t_end,NULL);
    // std::cout<<"save gevdm time:\t"<<(double)(t_end.tv_sec-t_start.tv_sec) + (double)(t_end.tv_usec-t_start.tv_usec)/1000000.0<<std::endl;
    return;
}

#endif