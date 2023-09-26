//wenfei 2022-1-5
//This file contains constructor and destructor of the class LCAO_deepks, 
//as well as subroutines for initializing and releasing relevant data structures 

//Other than the constructor and the destructor, it contains 3 types of subroutines:
//1. subroutines that are related to calculating descriptors:
//  - init : allocates some arrays
//  - init_index : records the index (inl)
//  - allocate_nlm : allocates data structures (nlm_save) which is used to store <chi|alpha>
//2. subroutines that are related to calculating force label:
//  - init_gdmx : allocates gdmx; it is a private subroutine
//  - del_gdmx : releases gdmx
//3. subroutines that are related to calculating force label:
//  - init_gdmepsl : allocates gdm_epsl; it is a private subroutine
//  - del_gdmepsl : releases gdm_epsl
//4. subroutines that are related to V_delta:
//  - allocate_V_delta : allocates H_V_delta; if calculating force, it also calls
//      init_gdmx, as well as allocating F_delta
//  - allocate_V_deltaR : allcoates H_V_deltaR, for multi-k calculations

#ifdef __DEEPKS

#include "LCAO_deepks.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

namespace GlobalC
{
    LCAO_Deepks ld;
}

//Constructor of the class
LCAO_Deepks::LCAO_Deepks()
{
    alpha_index = new ModuleBase::IntArray[1];
    inl_index = new ModuleBase::IntArray[1];
    inl_l = nullptr;
    H_V_deltaR = nullptr;
    gedm = nullptr;
}

//Desctructor of the class
LCAO_Deepks::~LCAO_Deepks()
{
    delete[] alpha_index;
    delete[] inl_index;
    delete[] inl_l;
    delete[] H_V_deltaR;

    //=======1. to use deepks, pdm is required==========
    //delete pdm**
    for (int inl = 0;inl < this->inlmax;inl++)
    {
        delete[] pdm[inl];
    }
    delete[] pdm;
    //=======2. "deepks_scf" part==========
    //if (GlobalV::deepks_scf)
    if (gedm)
    {
        //delete gedm**
        for (int inl = 0;inl < this->inlmax;inl++)
        {
            delete[] gedm[inl];
        }
        delete[] gedm;
    }

    del_gdmx();
}

void LCAO_Deepks::init(
    const LCAO_Orbitals& orb,
    const int nat,
    const int ntype,
    const Parallel_Orbitals& pv_in,
    std::vector<int> na)
{
    ModuleBase::TITLE("LCAO_Deepks", "init");

    GlobalV::ofs_running << " Initialize the descriptor index for DeePKS (lcao line)" << std::endl;

    const int lm = orb.get_lmax_d();
    const int nm = orb.get_nchimax_d();
    const int tot_inl_per_atom = orb.Alpha[0].getTotal_nchi();

    assert(lm >= 0);
    assert(nm >= 0);
    assert(tot_inl_per_atom >= 0);
    
    const int tot_inl = tot_inl_per_atom * nat;

    this->lmaxd = lm;
    this->nmaxd = nm;
    this->inlmax = tot_inl;
    GlobalV::ofs_running << " lmax of descriptor = " << this->lmaxd << std::endl;
    GlobalV::ofs_running << " nmax of descriptor= " << nmaxd << std::endl;
	GlobalV::ofs_running << " total basis (all atoms) for descriptor= " << std::endl;
    
    //init pdm**
    const int pdm_size = (this->lmaxd * 2 + 1) * (this->lmaxd * 2 + 1);
    this->pdm = new double* [this->inlmax];
    for (int inl = 0;inl < this->inlmax;inl++)
    {
        this->pdm[inl] = new double[pdm_size];
        ModuleBase::GlobalFunc::ZEROS(this->pdm[inl], pdm_size);
    }

    // cal n(descriptor) per atom , related to Lmax, nchi(L) and m. (not total_nchi!)
	this->des_per_atom=0; // mohan add 2021-04-21
    for (int l = 0; l <= this->lmaxd; l++)
    {
        this->des_per_atom += orb.Alpha[0].getNchi(l) * (2 * l + 1);
    }

    this->n_descriptor = nat * this->des_per_atom;

    this->init_index(ntype, nat, na, tot_inl, orb);
    this->allocate_nlm(nat);

    this->pv = &pv_in;

    return;
}

void LCAO_Deepks::init_index(const int ntype, const int nat, std::vector<int> na, const int Total_nchi, const LCAO_Orbitals &orb)
{
    delete[] this->alpha_index;
    this->alpha_index = new ModuleBase::IntArray[ntype];
    delete[] this->inl_index;
    this->inl_index = new ModuleBase::IntArray[ntype];
    delete[] this->inl_l;
    this->inl_l = new int[this->inlmax];
    ModuleBase::GlobalFunc::ZEROS(this->inl_l, this->inlmax);

    int inl = 0;
    int alpha = 0;
    for (int it = 0; it < ntype; it++)
    {
        this->alpha_index[it].create(
            na[it],
            this->lmaxd + 1, // l starts from 0
            this->nmaxd,
            2 * this->lmaxd + 1); // m ==> 2*l+1

        this->inl_index[it].create(
            na[it],
            this->lmaxd + 1,
            this->nmaxd); 

        GlobalV::ofs_running << " Type " << it + 1
                    << " number_of_atoms " << na[it] << std::endl;

        for (int ia = 0; ia < na[it]; ia++)
        {
            //alpha
            for (int l = 0; l < this->lmaxd + 1; l++)
            {
                for (int n = 0; n < orb.Alpha[0].getNchi(l); n++)
                {
                    for (int m = 0; m < 2 * l + 1; m++)
                    {
                        this->alpha_index[it](ia, l, n, m) = alpha;
                        alpha++;
                    }
                    this->inl_index[it](ia, l, n) = inl;
                    this->inl_l[inl] = l;
                    inl++;
                }
            }
        }//end ia
    }//end it
    assert(this->n_descriptor == alpha);
    assert(Total_nchi == inl);
    GlobalV::ofs_running << " descriptors_per_atom " << this->des_per_atom << std::endl;
    GlobalV::ofs_running << " total_descriptors " << this->n_descriptor << std::endl;
	return;
}

void LCAO_Deepks::allocate_nlm(const int nat)
{
    if(GlobalV::GAMMA_ONLY_LOCAL)
    {
        this->nlm_save.resize(nat);
    }
    else
    {
        this->nlm_save_k.resize(nat);
    }
}

void LCAO_Deepks::init_gdmx(const int nat)
{
    this->gdmx = new double** [nat];
    this->gdmy = new double** [nat];
    this->gdmz = new double** [nat];
    for (int iat = 0;iat < nat;iat++)
    {
        this->gdmx[iat] = new double* [inlmax];
        this->gdmy[iat] = new double* [inlmax];
        this->gdmz[iat] = new double* [inlmax];
        for (int inl = 0;inl < inlmax;inl++)
        {
            this->gdmx[iat][inl] = new double [(2 * lmaxd + 1) * (2 * lmaxd + 1)];
            this->gdmy[iat][inl] = new double [(2 * lmaxd + 1) * (2 * lmaxd + 1)];
            this->gdmz[iat][inl] = new double[(2 * lmaxd + 1) * (2 * lmaxd + 1)];
            ModuleBase::GlobalFunc::ZEROS(gdmx[iat][inl], (2 * lmaxd + 1) * (2 * lmaxd + 1));
            ModuleBase::GlobalFunc::ZEROS(gdmy[iat][inl], (2 * lmaxd + 1) * (2 * lmaxd + 1));
            ModuleBase::GlobalFunc::ZEROS(gdmz[iat][inl], (2 * lmaxd + 1) * (2 * lmaxd + 1));
        }
    }
    this->nat_gdm = nat;
    return;
}

//void LCAO_Deepks::del_gdmx(const int nat)
void LCAO_Deepks::del_gdmx()
{
    for (int iat = 0;iat < nat_gdm;iat++)
    {
        for (int inl = 0;inl < inlmax;inl++)
        {
            delete[] this->gdmx[iat][inl];
            delete[] this->gdmy[iat][inl];
            delete[] this->gdmz[iat][inl];
        }
        delete[] this->gdmx[iat];
        delete[] this->gdmy[iat];
        delete[] this->gdmz[iat];
    }
    delete[] this->gdmx;
    delete[] this->gdmy;
    delete[] this->gdmz;
    return;
}

void LCAO_Deepks::init_gdmepsl()
{
    this->gdm_epsl = new double** [6];
    
    for (int ipol = 0;ipol < 6;ipol++)
    {
        this->gdm_epsl[ipol] = new double* [inlmax];
        for (int inl = 0;inl < inlmax;inl++)
        {
            this->gdm_epsl[ipol][inl] = new double [(2 * lmaxd + 1) * (2 * lmaxd + 1)];
            ModuleBase::GlobalFunc::ZEROS(gdm_epsl[ipol][inl], (2 * lmaxd + 1) * (2 * lmaxd + 1));
        }
    }
    return;
}

void LCAO_Deepks::del_gdmepsl()
{
    for (int ipol = 0;ipol < 6;ipol++)
    {
        for (int inl = 0;inl < inlmax;inl++)
        {
            delete[] this->gdm_epsl[ipol][inl];
        }
        delete[] this->gdm_epsl[ipol];
    }
    delete[] this->gdm_epsl;
    return;
}


void LCAO_Deepks::allocate_V_delta(const int nat, const int nks)
{
    ModuleBase::TITLE("LCAO_Deepks", "allocate_V_delta");
    nks_V_delta = nks;

    //initialize the H matrix H_V_delta
    if(GlobalV::GAMMA_ONLY_LOCAL)
    {
        this->H_V_delta.resize(pv->nloc);
        ModuleBase::GlobalFunc::ZEROS(this->H_V_delta.data(), pv->nloc);
    }
    else
    {
        H_V_delta_k.resize(nks);
        for(int ik=0;ik<nks;ik++)
        {
            this->H_V_delta_k[ik].resize(pv->nloc);
            ModuleBase::GlobalFunc::ZEROS(this->H_V_delta_k[ik].data(), pv->nloc);
        }
    }

    //init gedm**
    const int pdm_size = (this->lmaxd * 2 + 1) * (this->lmaxd * 2 + 1);
    this->gedm = new double* [this->inlmax];
    for (int inl = 0;inl < this->inlmax;inl++)
    {
        this->gedm[inl] = new double[pdm_size];
        ModuleBase::GlobalFunc::ZEROS(this->gedm[inl], pdm_size);
    }
    if (GlobalV::CAL_FORCE)
    {
        //init F_delta
        F_delta.create(nat, 3);
        if(GlobalV::deepks_out_labels) 
        { 
            this->init_gdmx(nat);
            this->init_gdmepsl();
        }
        //gdmx is used only in calculating gvx
    }

    if (GlobalV::deepks_bandgap)
    {
        //init o_delta
        o_delta.create(nks, 1);
        
    }

    return;
}

void LCAO_Deepks::allocate_V_deltaR(const int nnr)
{
    ModuleBase::TITLE("LCAO_Deepks", "allocate_V_deltaR");
    GlobalV::ofs_running << nnr << std::endl;
    delete[] H_V_deltaR;
    H_V_deltaR = new double[nnr];
    ModuleBase::GlobalFunc::ZEROS(H_V_deltaR, nnr);
}

void LCAO_Deepks::init_orbital_pdm_shell(const int nks)
{
    
    this->orbital_pdm_shell = new double*** [nks];

    for (int iks=0; iks<nks; iks++)
    {
        this->orbital_pdm_shell[iks] = new double** [1];
        for (int hl=0; hl < 1; hl++)
        {
            this->orbital_pdm_shell[iks][hl] = new double* [this->inlmax];

            for(int inl = 0; inl < this->inlmax; inl++)
            {
                this->orbital_pdm_shell[iks][hl][inl] = new double [(2 * this->lmaxd + 1) * (2 * this->lmaxd + 1)];
                ModuleBase::GlobalFunc::ZEROS(orbital_pdm_shell[iks][hl][inl], (2 * this->lmaxd + 1) * (2 * this->lmaxd + 1));
            }
        }
    }

    return;
}


void LCAO_Deepks::del_orbital_pdm_shell(const int nks)
{
    for (int iks=0; iks<nks; iks++)
    {
        for (int hl=0; hl<1; hl++)
        {
            for (int inl = 0;inl < this->inlmax; inl++)
            {
                delete[] this->orbital_pdm_shell[iks][hl][inl];
            }
            delete[] this->orbital_pdm_shell[iks][hl];
        }
        delete[] this->orbital_pdm_shell[iks];
    }
     delete[] this->orbital_pdm_shell;    

    return;
}

#endif
