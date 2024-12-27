//wenfei 2022-1-11
//This file contains subroutines that contains interface with libnpy
#include "module_parameter/parameter.h"
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
//8. save_npy_h : Hamiltonian
//9. save_npy_v_delta_precalc : v_delta_precalc
//10. save_npy_phialpha : phialpha
//11. save_npy_gevdm : grav_evdm , can use phialpha and gevdm to calculate v_delta_precalc

#ifdef __DEEPKS

#include <mpi.h>
#include "LCAO_deepks_io.h"
#include "npy.hpp"

template <typename TK>
void LCAO_deepks_io::print_dm(const int nks, 
                                const int nlocal,
                                const int nrow,
                                const std::vector<std::vector<TK>>& dm)
{
    std::stringstream ss;
    for(int ik=0;ik<nks;ik++)
    {
        ss.str("");
        ss<<"dm_"<<ik;
        std::ofstream ofs(ss.str().c_str());
        ofs << std::setprecision(15);

        for (int mu=0;mu<nlocal;mu++)
        {
            for (int nu=0;nu<nlocal;nu++)
            {
                ofs << dm[ik][mu * nrow + nu] << " ";
            }
            ofs << std::endl;
        }
    }
}


void LCAO_deepks_io::load_npy_gedm(const int nat, 
                                   const int des_per_atom,
                                   double** gedm,
                                   double& e_delta,
                                   const int rank)
{
    ModuleBase::TITLE("LCAO_deepks_io", "load_npy_gedm");

    if(rank==0)
    {
        //load gedm.npy
        std::vector<double> npy_gedm;
        std::vector<unsigned long> dshape = {static_cast<unsigned long>(nat), 
                                             static_cast<unsigned long>(des_per_atom)};

        std::string gedm_file = "gedm.npy";

        npy::LoadArrayFromNumpy(gedm_file, dshape, npy_gedm);

        for (int iat = 0; iat < nat; iat++)
        {
            for(int ides = 0; ides < des_per_atom; ides++)
            {
                gedm[iat][ides] = npy_gedm[iat*des_per_atom + ides] * 2.0; //Ha to Ry
            }
        }

        //load ec.npy
        std::vector<double> npy_ec;
        std::vector<unsigned long> eshape = { 1ul };
        std::string ec_file = "ec.npy";
        npy::LoadArrayFromNumpy(ec_file, eshape, npy_ec);
        e_delta = npy_ec[0] * 2.0; //Ha to Ry
    }

#ifdef __MPI
    for(int iat = 0; iat < nat; iat++)
    {
        MPI_Bcast(gedm[iat], des_per_atom, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    MPI_Bcast(&e_delta, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

}


//saves descriptor into dm_eig.npy
void LCAO_deepks_io::save_npy_d(const int nat, 
                                const int des_per_atom,
                                const int inlmax,
                                const int *inl_l,
                                const bool deepks_equiv,
                                const std::vector<torch::Tensor> &d_tensor,
                                const std::string& out_dir,
                                const int rank)
{
    ModuleBase::TITLE("LCAO_deepks_io", "save_npy_d");

	if(rank!=0) 
	{ 
		return;
	}

    //save descriptor in .npy format
    // deepks_equiv was PARAM.inp.deepks_equiv
    if(!deepks_equiv)
    {
        std::vector<double> npy_des;
        for (int inl = 0;inl < inlmax;++inl)
        {
            int nm = 2*inl_l[inl] + 1;
            for(int im=0;im<nm;im++)
            {
                npy_des.push_back(d_tensor[inl].index({im}).item().toDouble());
            }
        }
        const long unsigned dshape[] = {static_cast<unsigned long>(nat), static_cast<unsigned long>(des_per_atom)};
        if (rank == 0)
        {
            std::string file_dm_eig = out_dir + "deepks_dm_eig.npy";
            //std::string file_dm_eig = "dm_eig.npy";
            npy::SaveArrayAsNumpy(file_dm_eig, false, 2, dshape, npy_des);
        }
    }
    else
    {
        // a rather unnecessary way of writing this, but I'll do it for now
        std::vector<double> npy_des;
        for(int iat = 0; iat < nat; iat ++)
        {
            for(int i = 0; i < des_per_atom; i++)
            {
                npy_des.push_back(d_tensor[iat].index({i}).item().toDouble());
            }
        }
        const long unsigned dshape[] = {static_cast<unsigned long>(nat), static_cast<unsigned long>(des_per_atom)};
        if (rank == 0)
        {
            std::string file_dm_eig = out_dir + "deepks_dm_eig.npy";
            //std::string file_dm_eig = "dm_eig.npy";
            npy::SaveArrayAsNumpy(file_dm_eig, false, 2, dshape, npy_des);
        }        
    }
    return;
}


//saves gvx into grad_vx.npy
void LCAO_deepks_io::save_npy_gvx(const int nat, 
                                  const int des_per_atom,
                                  const torch::Tensor &gvx_tensor,
                                  const std::string &out_dir,
                                  const int rank)
{
    ModuleBase::TITLE("LCAO_deepks_io", "save_npy_gvx");

	if(rank!=0)
	{
		return;
	}

    assert(nat>0);

    //save grad_vx.npy (when  force label is in use)
    //unit: /Bohr
    const long unsigned gshape[]
        = {static_cast<unsigned long>(nat), 
           3UL, 
           static_cast<unsigned long>(nat), 
           static_cast<unsigned long>(des_per_atom)};

    std::vector<double> npy_gvx;
    for (int ibt = 0;ibt < nat;++ibt)
    {
        for (int i = 0;i < 3;i++)
        {
            for (int iat = 0;iat < nat;++iat)
            {
                for(int p=0;p<des_per_atom;++p)
                {
                    npy_gvx.push_back(gvx_tensor.index({ ibt, i, iat, p }).item().toDouble());
                }
            }
        }
    }
    
    std::string file_gradvx = out_dir + "deepks_gradvx.npy";
    npy::SaveArrayAsNumpy(file_gradvx, false, 4, gshape, npy_gvx);
    return;
}

//saves gvx into grad_vepsl.npy
void LCAO_deepks_io::save_npy_gvepsl(const int nat,
                                     const int des_per_atom,
                                     const torch::Tensor &gvepsl_tensor,
                                     const std::string& out_dir,
                                     const int rank)
{
    ModuleBase::TITLE("LCAO_deepks_io", "save_npy_gvepsl");

	if(rank!=0) 
	{ 
		return;
	}

    //save grad_vepsl.npy (when  stress label is in use)
    //unit: none
    const long unsigned gshape[] = {6UL, 
                                    static_cast<unsigned long>(nat), 
                                    static_cast<unsigned long>(des_per_atom)};
 
    std::vector<double> npy_gvepsl;

    for (int i = 0;i < 6;i++)
    {
        for (int ibt = 0;ibt < nat;++ibt)
        {

            for(int p=0;p<des_per_atom;++p)
            {
                npy_gvepsl.push_back(gvepsl_tensor.index({ i, ibt, p }).item().toDouble());
            }

        }
    }

    // change the name from grad_vepsl.npy to deepks_gvepsl.npy
    const std::string file = out_dir + "deepks_gvepsl.npy";
    npy::SaveArrayAsNumpy(file, false, 3, gshape, npy_gvepsl);
    return;
}


//saves energy in numpy format
void LCAO_deepks_io::save_npy_e(const double &e, 
                                const std::string &e_file,
                                const int rank)
{
    ModuleBase::TITLE("LCAO_deepks_io", "save_npy_e");
	if(rank!=0) 
	{ 
		return;
	}

	//save e_base
    const long unsigned eshape[] = { 1 };
    std::vector<double> npy_e;
    npy_e.push_back(e);
    npy::SaveArrayAsNumpy(e_file, false, 1, eshape, npy_e);
    return;
}


//saves force in numpy format
void LCAO_deepks_io::save_npy_f(const ModuleBase::matrix &f, 
                                const std::string &f_file, 
                                const int nat,
                                const int rank)
{
    ModuleBase::TITLE("LCAO_deepks_io", "save_npy_f");

	if(rank!=0) 
	{ 
		return;
	}

    assert(nat>0);

	//save f_base
    //caution: unit: Rydberg/Bohr
    const long unsigned fshape[] = {static_cast<unsigned long>(nat), 3};
    std::vector<double> npy_f;
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
void LCAO_deepks_io::save_npy_s(const ModuleBase::matrix &stress, 
                                const std::string &s_file, 
                                const double &omega,
                                const int rank)
{
    ModuleBase::TITLE("LCAO_deepks_io", "save_npy_s");
	if(rank!=0) 
	{
		return;
	}

    const long unsigned sshape[] = { 6 };
    std::vector<double> npy_s;

    for (int ipol = 0;ipol < 3;++ipol)
    {
        for (int jpol = ipol;jpol < 3;jpol++)
        {
            npy_s.push_back(stress(ipol, jpol)*omega);
        }
    }
    npy::SaveArrayAsNumpy(s_file, false, 1, sshape, npy_s);
    return;
}


void LCAO_deepks_io::save_npy_o(const ModuleBase::matrix &bandgap, 
                                const std::string &o_file, 
                                const int nks,
                                const int rank)
{
	ModuleBase::TITLE("LCAO_deepks_io", "save_npy_o");
	if(rank!=0) 
	{ 
		return;
	}

    //save o_base
    const long unsigned oshape[] = {static_cast<unsigned long>(nks), 1 };

    std::vector<double> npy_o;
    for (int iks = 0; iks < nks; ++iks)
    {
        for (int hl = 0;hl < 1; ++hl)
        {
            npy_o.push_back(bandgap(iks,hl));
        }
    }

    npy::SaveArrayAsNumpy(o_file, false, 2, oshape, npy_o);
    return;
}


void LCAO_deepks_io::save_npy_orbital_precalc(const int nat, 
                                              const int nks, 
                                              const int des_per_atom,
                                              const torch::Tensor& orbital_precalc_tensor,
                                              const std::string& out_dir,
                                              const int rank)
{
    ModuleBase::TITLE("LCAO_deepks_io", "save_npy_orbital_precalc");
	if(rank!=0) 
	{ 
		return;
	}

    //save orbital_precalc.npy (when bandgap label is in use)
    //unit: a.u.
    const long unsigned gshape[] = {static_cast<unsigned long>(nks),
                                    1,
                                    static_cast<unsigned long>(nat),
                                    static_cast<unsigned long>(des_per_atom)};

    std::vector<double> npy_orbital_precalc;
    for (int iks = 0; iks < nks; ++iks)
    {
        for (int hl = 0; hl < 1; ++hl)
        {
            for (int iat = 0;iat < nat;++iat)
            {
                for(int p=0; p<des_per_atom; ++p)
                {
                    npy_orbital_precalc.push_back(orbital_precalc_tensor.index({iks, hl, iat, p }).item().toDouble());
                }
            }
        }
    }

    const std::string file_orbpre = out_dir + "deepks_orbpre.npy";
    npy::SaveArrayAsNumpy(file_orbpre, false, 4, gshape, npy_orbital_precalc);
    return;
}

template <typename TK, typename TH>
void LCAO_deepks_io::save_npy_h(const std::vector<TH> &hamilt,
                                const std::string &h_file,
                                const int nlocal,
                                const int nks,
                                const int rank)
{
    ModuleBase::TITLE("LCAO_deepks_io", "save_npy_h");
	if(rank!=0)
	{
		return;
	}

    const long unsigned hshape[] = {static_cast<unsigned long>(nks),
                                    static_cast<unsigned long>(nlocal), 
                                    static_cast<unsigned long>(nlocal) };

    std::vector<TK> npy_h;
    for(int k=0; k<nks; k++)
    {
        for (int i=0; i<nlocal; i++)
        {
            for (int j=0; j<nlocal; j++)
            {
                npy_h.push_back(hamilt[k](i,j));
            }
        }         
    }

    npy::SaveArrayAsNumpy(h_file, false, 3, hshape, npy_h);
    return;    
}

template <typename TK>
void LCAO_deepks_io::save_npy_v_delta_precalc(const int nat, 
                                              const int nks,
                                              const int nlocal, 
                                              const int des_per_atom,
                                              const torch::Tensor& v_delta_precalc_tensor,
                                              const std::string& out_dir,
                                              const int rank)
{
    ModuleBase::TITLE("LCAO_deepks_io", "save_npy_v_delta_precalc");
	if(rank!=0)
	{
		return;
	}

	// timeval t_start;
    // gettimeofday(&t_start,NULL);
    //save v_delta_precalc.npy (when v_delta label is in use)
    //unit: a.u.
    const long unsigned gshape[] = {static_cast<unsigned long>(nks),
                                    static_cast<unsigned long>(nlocal),
                                    static_cast<unsigned long>(nlocal),
                                    static_cast<unsigned long>(nat),
                                    static_cast<unsigned long>(des_per_atom)};

    std::vector<TK> npy_v_delta_precalc;    
    for (int iks = 0; iks < nks; ++iks)
    {
        for (int mu = 0; mu < nlocal; ++mu)
        {
            for (int nu = 0; nu < nlocal; ++nu)
            {
                for (int iat = 0;iat < nat;++iat)
                {
                    for(int p=0; p<des_per_atom; ++p)
                    {
                        if constexpr (std::is_same<TK, double>::value)
                        {
                            npy_v_delta_precalc.push_back(v_delta_precalc_tensor.index({iks, mu, nu, iat, p }).item().toDouble());
                        }
                        else
                        {
                            std::complex<double> value(torch::real(v_delta_precalc_tensor.index({iks, mu, nu, iat, p})).item<double>(),
                                                        torch::imag(v_delta_precalc_tensor.index({iks, mu, nu, iat, p})).item<double>());
                            npy_v_delta_precalc.push_back(value);
                        }
                    }
                }                
            }
        }
    }
    const std::string file_vdpre = out_dir + "deepks_vdpre.npy";
    npy::SaveArrayAsNumpy(file_vdpre, false, 5, gshape, npy_v_delta_precalc);
    return;
}

template <typename TK>
void LCAO_deepks_io::save_npy_phialpha(const int nat, 
                                       const int nks,
                                       const int nlocal,
                                       const int inlmax,
                                       const int lmaxd,
                                       const torch::Tensor &phialpha_tensor,
                                       const std::string& out_dir,
                                       const int rank)
{
    ModuleBase::TITLE("LCAO_deepks_io", "save_npy_phialpha");
	if(rank!=0)
	{ 
		return;
	}

    //save phialpha.npy (when v_delta label == 2)
    //unit: a.u.
    const int nlmax = inlmax/nat;
    const int mmax = 2*lmaxd+1;
    const long unsigned gshape[] = {static_cast<unsigned long>(nat),
                                    static_cast<unsigned long>(nlmax),
                                    static_cast<unsigned long>(nks),
                                    static_cast<unsigned long>(nlocal),
                                    static_cast<unsigned long>(mmax)};
    std::vector<TK> npy_phialpha;
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
                        if constexpr (std::is_same<TK, double>::value)
                        {
                            npy_phialpha.push_back(phialpha_tensor.index({ iat,nl, iks, mu, m }).item().toDouble());
                        }
                        else
                        {
                            std::complex<double> value(torch::real(phialpha_tensor.index({ iat, nl, iks, mu, m })).item<double>(), 
                                                       torch::imag(phialpha_tensor.index({ iat, nl, iks, mu, m })).item<double>());
                            npy_phialpha.push_back(value);
                        }
                    }
                }                
            }
        }
    }
    const std::string file_phialpha = out_dir + "deepks_phialpha.npy";
    npy::SaveArrayAsNumpy(file_phialpha, false, 5, gshape, npy_phialpha);
    return;
}


void LCAO_deepks_io::save_npy_gevdm(const int nat,
                                    const int inlmax,
                                    const int lmaxd,
                                    const torch::Tensor& gevdm_tensor,
                                    const std::string& out_dir,
                                    const int rank)
{
    ModuleBase::TITLE("LCAO_deepks_io", "save_npy_gevdm");
	if(rank!=0) 
	{ 
		return;
	}

    assert(nat>0);

	//save grad_evdm.npy (when v_delta label == 2)
    //unit: a.u.
    const int nlmax = inlmax/nat;
    const int mmax = 2*lmaxd+1;
    const long unsigned gshape[] = {static_cast<unsigned long>(nat),
                                    static_cast<unsigned long>(nlmax),
                                    static_cast<unsigned long>(mmax),
                                    static_cast<unsigned long>(mmax),
                                    static_cast<unsigned long>(mmax)};
    std::vector<double> npy_gevdm;
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
                        npy_gevdm.push_back(gevdm_tensor.index({ iat, nl, v, m, n}).item().toDouble());
                    }         
                }
            }                
        }
    }
    const std::string file_gevdm = out_dir + "deepks_gevdm.npy";
    npy::SaveArrayAsNumpy(file_gevdm, false, 5, gshape, npy_gevdm);
    return;
}


template void LCAO_deepks_io::print_dm<double>(const int nks, 
                                               const int nlocal,
                                               const int nrow,
                                               const std::vector<std::vector<double>>& dm);

template void LCAO_deepks_io::print_dm<std::complex<double>>(const int nks, 
                                                             const int nlocal,
                                                             const int nrow,
                                                             const std::vector<std::vector<std::complex<double>>>& dm);

template void LCAO_deepks_io::save_npy_h<double>(const std::vector<ModuleBase::matrix> &hamilt,
                                                 const std::string &h_file,
                                                 const int nlocal,
                                                 const int nks,
                                                 const int rank);

template void LCAO_deepks_io::save_npy_h<std::complex<double>>(const std::vector<ModuleBase::ComplexMatrix> &hamilt,
                                                               const std::string &h_file,
                                                               const int nlocal,
                                                               const int nks,
                                                               const int rank);

template void LCAO_deepks_io::save_npy_v_delta_precalc<double>(const int nat, 
                                                               const int nks,
                                                               const int nlocal, 
                                                               const int des_per_atom,
                                                               const torch::Tensor& v_delta_precalc_tensor,
                                                               const std::string& out_dir,
                                                               const int rank);

template void LCAO_deepks_io::save_npy_v_delta_precalc<std::complex<double>>(const int nat, 
                                                                             const int nks,
                                                                             const int nlocal, 
                                                                             const int des_per_atom,
                                                                             const torch::Tensor& v_delta_precalc_tensor,
                                                                             const std::string& out_dir,
                                                                             const int rank);

template void LCAO_deepks_io::save_npy_phialpha<double>(const int nat,
                                                        const int nks,
                                                        const int nlocal,
                                                        const int inlmax,
                                                        const int lmaxd,
                                                        const torch::Tensor &phialpha_tensor,
                                                        const std::string& out_dir,
                                                        const int rank);

template void LCAO_deepks_io::save_npy_phialpha<std::complex<double>>(const int nat,
                                                                      const int nks,
                                                                      const int nlocal,
                                                                      const int inlmax,
                                                                      const int lmaxd,
                                                                      const torch::Tensor &phialpha_tensor,
                                                                      const std::string& out_dir,
                                                                      const int rank);

#endif
