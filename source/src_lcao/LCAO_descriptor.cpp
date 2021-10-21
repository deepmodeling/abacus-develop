//caoyu add 2021-03-2
#ifdef __DEEPKS

#include "LCAO_descriptor.h"
#include "LCAO_matrix.h"
#include "../module_base/lapack_connector.h"
#include "../module_base/intarray.h"
#include "../module_base/complexmatrix.h"
#include "global_fp.h"
#include "../src_pw/global.h"
#include "../src_io/winput.h"

#include <torch/script.h>
#include <torch/csrc/autograd/autograd.h>
#include <npy.hpp>
#include <torch/csrc/api/include/torch/linalg.h>

namespace GlobalC
{
    LCAO_Descriptor ld;
}
LCAO_Descriptor::LCAO_Descriptor()
{
    alpha_index = new ModuleBase::IntArray[1];
    inl_index = new ModuleBase::IntArray[1];
    inl_l = new int[1];
    d = new double[1];
    H_V_delta = new double[1];
    dm_double = new double[1];
    DH_V_delta_x = new double[1];
    DH_V_delta_y = new double[1];
    DH_V_delta_z = new double[1];
}


LCAO_Descriptor::~LCAO_Descriptor()
{
    delete[] alpha_index;
    delete[] inl_index;
    delete[] inl_l;
    delete[] d;
    delete[] H_V_delta;
    delete[] dm_double;
    delete[] DH_V_delta_x;
    delete[] DH_V_delta_y;
    delete[] DH_V_delta_z;

    //=======1. "out_descriptor" part==========
    //delete S_mu_alpha**
    for (int inl = 0;inl < this->inlmax;inl++)
    {
        delete[] S_mu_alpha[inl];
    }
    delete[] S_mu_alpha;
    //delete pdm**
    for (int inl = 0;inl < this->inlmax;inl++)
    {
        delete[] pdm[inl];
    }
    delete[] pdm;
    //=======2. "deepks_scf" part==========
    if (GlobalV::deepks_scf)
    {
        //delete gedm**
        for (int inl = 0;inl < this->inlmax;inl++)
        {
            delete[] gedm[inl];
        }
        delete[] gedm;
        if (GlobalV::FORCE)
        {
            //delete DS_mu_alpha**
            for (int inl = 0;inl < this->inlmax;inl++)
            {
                delete[] DS_mu_alpha_x[inl];
                delete[] DS_mu_alpha_y[inl];
                delete[] DS_mu_alpha_z[inl];
            }
            delete[] DS_mu_alpha_x;
            delete[] DS_mu_alpha_y;
            delete[] DS_mu_alpha_z;
        }
    }
}


void LCAO_Descriptor::init(
	const int lm, // max L for descriptor 
	const int nm, // max n for descriptor
	const int tot_inl) // total number of atomic orbital basis for all of the atoms 
{
    ModuleBase::TITLE("LCAO_Descriptor", "init");

    GlobalV::ofs_running << " Initialize the descriptor index for DeePKS (lcao line)" << std::endl;

    assert(lm >= 0);
    assert(nm >= 0);
    assert(tot_inl >= 0);
    
    this->lmaxd = lm;
    this->nmaxd = nm;
    this->inlmax = tot_inl;
    GlobalV::ofs_running << " lmax of descriptor = " << this->lmaxd << std::endl;
    GlobalV::ofs_running << " nmax of descriptor= " << nmaxd << std::endl;
	GlobalV::ofs_running << " total basis (all atoms) for descriptor= " << std::endl;
    
    //initialize the density matrix (dm)
    delete[] this->dm_double;
    this->dm_double = new double[GlobalV::NLOCAL * GlobalV::NLOCAL];
    ModuleBase::GlobalFunc::ZEROS(this->dm_double, GlobalV::NLOCAL * GlobalV::NLOCAL);
    
    //init S_mu_alpha**
    this->S_mu_alpha = new double* [this->inlmax];    //inl
    for (int inl = 0;inl < this->inlmax;inl++)
    {
        this->S_mu_alpha[inl] = new double[GlobalV::NLOCAL * (2 * this->lmaxd + 1)];     //GlobalV::NLOCAL*nm
        ModuleBase::GlobalFunc::ZEROS(S_mu_alpha[inl], GlobalV::NLOCAL * (2 * this->lmaxd+ 1));
    }

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
        this->des_per_atom += GlobalC::ORB.Alpha[0].getNchi(l) * (2 * l + 1);
    }

    this->n_descriptor = GlobalC::ucell.nat * this->des_per_atom;

    this->init_index();
    
    return;
}


void LCAO_Descriptor::init_index(void)
{
    delete[] this->alpha_index;
    this->alpha_index = new ModuleBase::IntArray[GlobalC::ucell.ntype];
    delete[] this->inl_index;
    this->inl_index = new ModuleBase::IntArray[GlobalC::ucell.ntype];
    delete[] this->inl_l;
    this->inl_l = new int[this->inlmax];
    ModuleBase::GlobalFunc::ZEROS(this->inl_l, this->inlmax);

    int inl = 0;
    int alpha = 0;
    for (int it = 0; it < GlobalC::ucell.ntype; it++)
    {
        this->alpha_index[it].create(
            GlobalC::ucell.atoms[it].na,
            this->lmaxd + 1, // l starts from 0
            this->nmaxd,
            2 * this->lmaxd + 1); // m ==> 2*l+1

        this->inl_index[it].create(
            GlobalC::ucell.atoms[it].na,
            this->lmaxd + 1,
            this->nmaxd); 

        GlobalV::ofs_running << " Type " << it + 1
                    << " number_of_atoms " << GlobalC::ucell.atoms[it].na << std::endl;

        for (int ia = 0; ia < GlobalC::ucell.atoms[it].na; ia++)
        {
            //alpha
            for (int l = 0; l < this->lmaxd + 1; l++)
            {
                for (int n = 0; n < GlobalC::ORB.Alpha[0].getNchi(l); n++)
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
    assert(GlobalC::ucell.nat * GlobalC::ORB.Alpha[0].getTotal_nchi() == inl);
    GlobalV::ofs_running << " descriptors_per_atom " << this->des_per_atom << std::endl;
    GlobalV::ofs_running << " total_descriptors " << this->n_descriptor << std::endl;
	return;
}


void LCAO_Descriptor::build_S_descriptor(const bool& calc_deri)
{
    ModuleBase::TITLE("LCAO_Descriptor", "build_S_descriptor");
    //array to store data
    double olm[3] = {0.0, 0.0, 0.0};

    //\sum{T} e**{ikT} <\phi_{ia}|d\phi_{k\beta}(T)>	//???
    ModuleBase::Vector3<double> tau1, tau2, dtau;
    ModuleBase::Vector3<double> dtau1, dtau2, tau0;
    for (int T1 = 0; T1 < GlobalC::ucell.ntype; ++T1)
    {
        Atom *atom1 = &GlobalC::ucell.atoms[T1];
        for (int I1 = 0; I1 < atom1->na; ++I1)
        {
            tau1 = atom1->tau[I1];
            //GlobalC::GridD.Find_atom(tau1);
            GlobalC::GridD.Find_atom(GlobalC::ucell, tau1, T1, I1);

            for (int ad = 0; ad < GlobalC::GridD.getAdjacentNum() + 1; ++ad)
            {
                const int T2 = GlobalC::GridD.getType(ad);
                const int I2 = GlobalC::GridD.getNatom(ad);
                //Atom *atom2 = &GlobalC::ucell.atoms[T2];
                tau2 = GlobalC::GridD.getAdjacentTau(ad);
                dtau = tau2 - tau1;
                double distance = dtau.norm() * GlobalC::ucell.lat0;
                double rcut = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Alpha[0].getRcut(); //Rcut is subject to GlobalC::ORB.Phi to keep dimension of S_mu_alpha same as Sloc
                if (distance < rcut)
                {
                    int iw1_all = GlobalC::ucell.itiaiw2iwt(T1, I1, 0); //iw1_all = combined index (it, ia, iw)

                    for (int jj = 0; jj < atom1->nw * GlobalV::NPOL; ++jj)
                    {
                        const int jj0 = jj / GlobalV::NPOL;
                        const int L1 = atom1->iw2l[jj0];
                        const int N1 = atom1->iw2n[jj0];
                        const int m1 = atom1->iw2m[jj0];

                        for (int L2 = 0; L2 <= GlobalC::ORB.Alpha[0].getLmax(); ++L2)
                        {
                            for (int N2 = 0; N2 < GlobalC::ORB.Alpha[0].getNchi(L2); ++N2)
                            {
                                for (int m2 = 0; m2 < 2 * L2 + 1; ++m2)
                                {
                                    olm[0] = olm[1] = olm[2] = 0.0;
                                    if (!calc_deri)
                                    {
                                        GlobalC::UOT.snap_psipsi(GlobalC::ORB, olm, 0, 'D', tau1,
                                                T1, L1, m1, N1, GlobalC::GridD.getAdjacentTau(ad),
                                                T2, L2, m2, N2, GlobalV::NSPIN);
                                        if (GlobalV::GAMMA_ONLY_LOCAL)
                                        {
                                            this->set_S_mu_alpha(iw1_all, inl_index[T2](I2,L2,N2), m2, olm[0]);
                                        }
                                    }
                                    else
                                    {
                                        GlobalC::UOT.snap_psipsi(GlobalC::ORB, olm, 1, 'D', tau1,
                                            T1, L1, m1, N1, GlobalC::GridD.getAdjacentTau(ad),
                                            T2, L2, m2, N2, GlobalV::NSPIN);
                                        if (GlobalV::GAMMA_ONLY_LOCAL)
                                        {
                                            this->set_DS_mu_alpha(iw1_all, inl_index[T2](I2,L2,N2), m2, olm[0], olm[1], olm[2]);
                                        }
                                    }

                                } //m2
                            }     //N2
                        }         //nw2(L2)
                        ++iw1_all;
                    } // nw1
                }     // distance
            }         // ad
        } // I1
    }     // T1
    if (!GlobalV::GAMMA_ONLY_LOCAL)
    {
        ModuleBase::WARNING_QUIT("LCAO_Descriptor::build_S_descriptor", "muti-kpoint method for descriptor is not implemented yet! ");
    }
    return;
}


void LCAO_Descriptor::set_S_mu_alpha(
	const int &iw1_all, 
	const int &inl, 
	const int &im, 
	const double &v)
{
    //const int ir = GlobalC::ParaO.trace_loc_row[iw1_all];
    //const int ic = GlobalC::ParaO.trace_loc_col[iw2_all];
    //no parellel yet
    const int ir = iw1_all;
    const int ic = im;
    //const int index = ir * GlobalC::ParaO.ncol + ic;
    int index;
    if (GlobalV::KS_SOLVER == "genelpa" || GlobalV::KS_SOLVER == "scalapack_gvx") // save the matrix as column major format
    {
        index = ic * GlobalV::NLOCAL + ir;
    }
    else
    {
        index = ir * (2*inl_l[inl]+1)  + ic; //row: lcao orbitals; col: descriptor basis
    }
    this->S_mu_alpha[inl][index] += v;
    return;
}



// compute the full projected density matrix for each atom
// save the matrix for each atom in order to minimize the usage of memory
// --mohan 2021-08-04
void LCAO_Descriptor::cal_dm_as_descriptor(const ModuleBase::matrix &dm)
{
	ModuleBase::TITLE("LCAO_Descriptor", "cal_proj_dm");

	for(int it=0; it<GlobalC::ucell.ntype; ++it)
	{
		for(int ia=0; ia<GlobalC::ucell.atoms[it].na; ++ia)
		{
			// compute S^T * dm * S to obtain descriptor
			// of each atom 
			// and then diagonalize it
		}
	}

	return;
}


void LCAO_Descriptor::cal_projected_DM(const ModuleBase::matrix &dm)
{
    ModuleBase::TITLE("LCAO_Descriptor", "cal_projected_DM");
    //step 1: get dm: the coefficient of wfc, not charge density
    //now,  dm is an input arg of this func, but needed converting to double*
    this->getdm_double(dm);

    //step 2: get S_alpha_mu and S_nu_beta
    double **ss = this->S_mu_alpha;

    //step 3 : multiply: cal ST*DM*S
    
    //init tmp_pdm*
    const int tmp_pdm_size = GlobalV::NLOCAL * (lmaxd*2+1);
    double* tmp_pdm = new double[tmp_pdm_size];
    ModuleBase::GlobalFunc::ZEROS(tmp_pdm, tmp_pdm_size);

    for (int inl = 0;inl < inlmax;inl++)
    {   
        int nm = 2 * inl_l[inl] + 1;   //1,3,5,...
        const char t = 'T';  //transpose
        const char nt = 'N'; //non transpose
        const double alpha = 1;
        const double beta = 0;
        double *a = this->dm_double;
        double *b = ss[inl];
        double *c = tmp_pdm;
        dgemm_(&nt, &nt, &GlobalV::NLOCAL, &nm, &GlobalV::NLOCAL, &alpha, a, &GlobalV::NLOCAL, b, &GlobalV::NLOCAL, &beta, c, &GlobalV::NLOCAL); //DM*S
        a = ss[inl];
        b = c;
        c = this->pdm[inl];
        dgemm_(&t, &nt, &nm, &nm, &GlobalV::NLOCAL, &alpha, a, &GlobalV::NLOCAL, b, &GlobalV::NLOCAL, &beta, c, &nm); //ST*DM*S
    }
    
    delete[] tmp_pdm;
    return;
}


void LCAO_Descriptor::cal_descriptor(void)
{
    ModuleBase::TITLE("LCAO_Descriptor", "cal_descriptor");
    delete[] d;
    d = new double[this->n_descriptor];
    const int lmax = GlobalC::ORB.get_lmax_d();
    int id = 0;
    for (int it = 0; it < GlobalC::ucell.ntype; it++)
    {
        for (int ia = 0; ia < GlobalC::ucell.atoms[it].na; ia++)
        {
            for (int l = 0; l <= lmax; l++)
            {
                int nmax = GlobalC::ORB.Alpha[0].getNchi(l);
                for (int n = 0; n < nmax; n++)
                {
                    const int dim = 2 * l + 1;
                    const int inl = inl_index[it](ia, l, n);
                    // descriptor for atom (it, ia)
                    ModuleBase::ComplexMatrix des(dim, dim);
                    for (int m = 0; m < dim; m++)
                    {
                        for (int m2 = 0; m2 < dim; m2++)
                        {
                            int index = m * dim + m2;
                            complex<double> tmp(this->pdm[inl][index], 0);
                            des(m, m2) += tmp;
                        }
                    }
                    /*
                    stringstream ss;
                    ss << winput::spillage_outdir << "/"<< "PDM.dat";
                    if (GlobalV::MY_RANK == 0)
                    {
                        ofstream ofs;
                        ofs.open(ss.str().c_str());
                        this->print_projected_DM(ofs, des, it, ia, l, n);
                        ofs.close();
                    }
                    */
                    if (l == 0)
                    {
                        this->d[id] = des(0, 0).real();
                        ++id;
                    }
                    else
                    {
                        // diagonalizae
                        // assume des matrix is Hermitian
                        char jobz = 'N'; // eigenvalues only
                        char uplo = 'U'; // upper matrix is stored
                        int ndim = des.nr;
                        double *tmpd = new double[ndim]();
                        const int lwork = 2 * ndim;
                        complex<double> *work = new complex<double>[lwork]();
                        double *rwork = new double[3 * ndim - 2]();
                        int infor = 0;
                        // diag by calling zheev
                        LapackConnector::zheev(jobz, uplo, ndim, des, ndim, tmpd, work, lwork, rwork, &infor);
                        // put the eigenvalues into d (descriptor)
                        for (int idim = 0; idim < ndim; ++idim)
                        {
                            this->d[id] = tmpd[idim];
                            ++id;
                        }
                        delete[] tmpd;
                        delete[] rwork;
                        delete[] work;
                    }
                } //n
            }     //l
        }         //ia
    }             //it
    this->print_descriptor();
    return;
}



void LCAO_Descriptor::print_projected_DM(
	ofstream& ofs, 
	ModuleBase::ComplexMatrix& des, 
	const int& it, // index for atom type
	const int& ia, // index for atoms
	const int& l, // index for angular momentum quantum number L
	const int& n) // index for principal quantum number n
{
    ofs << "L=" << l << "   N=" << n << std::endl;
    for (int i = 0; i < 2 * l + 1; i++)
    {
        for (int j = 0; j < 2 * l + 1; j++)
        {
            ofs << des(i, j).real() << " ";
        }
        ofs << std::endl;
    }
    return;
}


void LCAO_Descriptor::print_descriptor(void)
{
    ModuleBase::TITLE("LCAO_Descriptor", "print_descriptor");
    ofstream ofs;
    stringstream ss;
    // the parameter 'winput::spillage_outdir' is read from INPUTw.
    ss << winput::spillage_outdir << "/"
       << "descriptor.dat";
    if (GlobalV::MY_RANK == 0)
    {
        ofs.open(ss.str().c_str());
    }
    for (int it = 0; it < GlobalC::ucell.ntype; it++)
    {
        for (int ia = 0; ia < GlobalC::ucell.atoms[it].na; ia++)
        {
            ofs << GlobalC::ucell.atoms[it].label << " atom_index " << ia + 1 << " n_descriptor " << this->des_per_atom << std::endl;
            int id0 = this->alpha_index[it](ia, 0, 0, 0);
            for (int id = id0; id < id0 + this->des_per_atom; ++id)
            {
                if ((id - id0) > 0 && (id - id0) % 8 == 0)
                    ofs << std::endl;
                ofs << d[id] << " ";
            }
            ofs << std::endl << std::endl;
        }
    }
    GlobalV::ofs_running << " Descriptors have been printed to " << ss.str() << std::endl;

    return;
}


void LCAO_Descriptor::set_DS_mu_alpha(
	const int& iw1_all, 
	const int& inl, 
	const int& im,
    const double& vx, 
	const double& vy, 
	const double& vz)
{
    //const int ir = GlobalC::ParaO.trace_loc_row[iw1_all];
    //const int ic = GlobalC::ParaO.trace_loc_col[iw2_all];
    //no parellel yet
    const int ir = iw1_all;
    const int ic = im;
    //const int index = ir * GlobalC::ParaO.ncol + ic;
    int index;
    if (GlobalV::KS_SOLVER == "genelpa" || GlobalV::KS_SOLVER == "scalapack_gvx") // save the matrix as column major format
    {
        index = ic * GlobalV::NLOCAL + ir;
    }
    else
    {
        index = ir * (2*inl_l[inl]+1)  + ic; //row: lcao orbitals; col: descriptor basis
    }
    this->DS_mu_alpha_x[inl][index] += vx;
    this->DS_mu_alpha_y[inl][index] += vy;
    this->DS_mu_alpha_z[inl][index] += vz;
    return;
}


void LCAO_Descriptor::getdm_double(const ModuleBase::matrix &dm)
{
    for (int i = 0; i < dm.nr; i++)
    {
        for (int j = 0; j < dm.nc; j++)
        {
            this->dm_double[i * GlobalV::NLOCAL + j] = dm(i, j); //only consider default GlobalV::NSPIN = 1
        }
    }
	return;
}


void LCAO_Descriptor::cal_gdmx(const ModuleBase::matrix &dm)
{
    ModuleBase::TITLE("LCAO_Descriptor", "cal_gdmx");
    //get DS_alpha_mu and S_nu_beta
    double** ss = this->S_mu_alpha;
    double** dsx = this->DS_mu_alpha_x;
    double** dsy = this->DS_mu_alpha_y;
    double** dsz = this->DS_mu_alpha_z;
    for (int inl = 0;inl < inlmax;inl++)
    {
        //dE/dD will be multiplied in cal_f_delta, here only calculate dD/dx_I
        int nm = 2 * inl_l[inl] + 1;   //1,3,5,...
        for (int mu =0; mu < GlobalV::NLOCAL;++mu) 
        {
            const int iat = GlobalC::ucell.iwt2iat[mu];//the atom whose force being calculated
            for (int nu= 0;nu < GlobalV::NLOCAL; ++nu)
            {
                //const int mu = GlobalC::ParaO.trace_loc_row[j];
                //const int nu = GlobalC::ParaO.trace_loc_col[i];
                if (mu >= 0 && nu >= 0)
                {
                    for (int m1 = 0;m1 < nm;++m1)
                    {
                        for (int m2 = 0;m2 < nm;++m2)
                        {
                            for (int is = 0;is < GlobalV::NSPIN;++is)
                            {
								//  save the matrix as column major format
                                if (GlobalV::KS_SOLVER == "genelpa" || GlobalV::KS_SOLVER == "scalapack_gvx")
                                {
                                    gdmx[iat][inl][m1*nm + m2] += 
									2 * dsx[inl][m1*GlobalV::NLOCAL + mu] * dm(mu, nu) * ss[inl][m2*GlobalV::NLOCAL + nu];
                                    gdmy[iat][inl][m1*nm + m2] += 
									2 * dsy[inl][m1*GlobalV::NLOCAL + mu] * dm(mu, nu) * ss[inl][m2*GlobalV::NLOCAL + nu];
                                    gdmz[iat][inl][m1*nm + m2] += 
									2 * dsz[inl][m1*GlobalV::NLOCAL + mu] * dm(mu, nu) * ss[inl][m2*GlobalV::NLOCAL + nu];
                                }
                                else
                                {
                                    gdmx[iat][inl][m1*nm + m2] += 
									2 * dsx[inl][mu*nm + m1] * dm(mu, nu) * ss[inl][nu*nm + m2];
                                    gdmy[iat][inl][m1*nm + m2] += 
									2 * dsy[inl][mu*nm + m1] * dm(mu, nu) * ss[inl][nu*nm + m2];
                                    gdmz[iat][inl][m1*nm + m2] += 
									2 * dsz[inl][mu*nm + m1] * dm(mu, nu) * ss[inl][nu*nm + m2];
                                }
                            }
                        }//end m2
                    } //end m1
                }//end if
            }//end nu
        }//end mu
    }//end inl
    //print GDMX
    return;
}


void LCAO_Descriptor::init_gdmx(void)
{
    this->gdmx = new double** [GlobalC::ucell.nat];
    this->gdmy = new double** [GlobalC::ucell.nat];
    this->gdmz = new double** [GlobalC::ucell.nat];
    for (int iat = 0;iat < GlobalC::ucell.nat;iat++)
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
    return;
}


void LCAO_Descriptor::del_gdmx(void)
{
    for (int iat = 0;iat < GlobalC::ucell.nat;iat++)
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


void LCAO_Descriptor::deepks_pre_scf(const string& model_file)
{
    ModuleBase::TITLE("LCAO_Descriptor", "deepks_pre_scf");

	// load the DeePKS model from deep neural network
    this->load_model(model_file);
    
    //initialize the H matrix H_V_delta
    delete[] this->H_V_delta;
    this->H_V_delta = new double[GlobalV::NLOCAL * GlobalV::NLOCAL];
    ModuleBase::GlobalFunc::ZEROS(this->H_V_delta, GlobalV::NLOCAL * GlobalV::NLOCAL);

    //init gedm**
    const int pdm_size = (this->lmaxd * 2 + 1) * (this->lmaxd * 2 + 1);
    this->gedm = new double* [this->inlmax];
    for (int inl = 0;inl < this->inlmax;inl++)
    {
        this->gedm[inl] = new double[pdm_size];
        ModuleBase::GlobalFunc::ZEROS(this->gedm[inl], pdm_size);
    }
    if (GlobalV::FORCE)
    {
        //init F_delta
        F_delta.create(GlobalC::ucell.nat, 3);
        //init DS_mu_alpha**
        this->DS_mu_alpha_x = new double* [this->inlmax];
        this->DS_mu_alpha_y = new double* [this->inlmax];
        this->DS_mu_alpha_z = new double* [this->inlmax];
        for (int inl = 0;inl < this->inlmax;inl++)
        {
            this->DS_mu_alpha_x[inl] = new double[GlobalV::NLOCAL * (2 * this->lmaxd + 1)];
            this->DS_mu_alpha_y[inl] = new double[GlobalV::NLOCAL * (2 * this->lmaxd + 1)];
            this->DS_mu_alpha_z[inl] = new double[GlobalV::NLOCAL * (2 * this->lmaxd + 1)];
            ModuleBase::GlobalFunc::ZEROS(DS_mu_alpha_x[inl], GlobalV::NLOCAL * (2 * this->lmaxd + 1));
            ModuleBase::GlobalFunc::ZEROS(DS_mu_alpha_y[inl], GlobalV::NLOCAL * (2 * this->lmaxd + 1));
            ModuleBase::GlobalFunc::ZEROS(DS_mu_alpha_z[inl], GlobalV::NLOCAL * (2 * this->lmaxd + 1));
        }
        //init DH_V_delta*
        delete[] DH_V_delta_x;
        delete[] DH_V_delta_y;
        delete[] DH_V_delta_z;
        this->DH_V_delta_x = new double[GlobalV::NLOCAL * GlobalV::NLOCAL];
        this->DH_V_delta_y = new double [GlobalV::NLOCAL * GlobalV::NLOCAL];
        this->DH_V_delta_z = new double[GlobalV::NLOCAL * GlobalV::NLOCAL];
        ModuleBase::GlobalFunc::ZEROS(DH_V_delta_x, GlobalV::NLOCAL * GlobalV::NLOCAL);
        ModuleBase::GlobalFunc::ZEROS(DH_V_delta_y, GlobalV::NLOCAL * GlobalV::NLOCAL);
        ModuleBase::GlobalFunc::ZEROS(DH_V_delta_z, GlobalV::NLOCAL * GlobalV::NLOCAL);
    }
    return;
}


void LCAO_Descriptor::cal_v_delta(const ModuleBase::matrix& dm)
{
    ModuleBase::TITLE("LCAO_Descriptor", "cal_v_delta");
    //1.  (dE/dD)<alpha_m'|psi_nv> (descriptor changes in every scf iter)
    this->cal_gedm(dm);
    
    //2. multiply overlap matrice and sum
    double* tmp_v1 = new double[(2 * lmaxd + 1) * GlobalV::NLOCAL];
    double* tmp_v2 = new double[GlobalV::NLOCAL *GlobalV::NLOCAL];

    ModuleBase::GlobalFunc::ZEROS(this->H_V_delta, GlobalV::NLOCAL * GlobalV::NLOCAL); //init before calculate
    
    for (int inl = 0;inl < inlmax;inl++)
    {
        ModuleBase::GlobalFunc::ZEROS(tmp_v1, (2 * lmaxd + 1) * GlobalV::NLOCAL);
        ModuleBase::GlobalFunc::ZEROS(tmp_v2, GlobalV::NLOCAL * GlobalV::NLOCAL);
        int nm = 2 * inl_l[inl] + 1;   //1,3,5,...
        const char t = 'T';  //transpose
        const char nt = 'N'; //non transpose
        const double alpha = 1;
        const double beta = 0;
        double* a = this->gedm[inl];//[nm][nm]
        double* b = S_mu_alpha[inl];//[GlobalV::NLOCAL][nm]--trans->[nm][GlobalV::NLOCAL]
        double* c = tmp_v1;
        
        //2.1  (dE/dD)*<alpha_m'|psi_nv>
        dgemm_(&nt, &t, &nm, &GlobalV::NLOCAL, &nm, &alpha, a, &nm, b, &GlobalV::NLOCAL, &beta, c, &nm);

        //2.2  <psi_mu|alpha_m>*(dE/dD)*<alpha_m'|psi_nv>
        a = b; //[GlobalV::NLOCAL][nm]
        b = c;//[nm][GlobalV::NLOCAL]
        c = tmp_v2;//[GlobalV::NLOCAL][GlobalV::NLOCAL]
        dgemm_(&nt, &nt, &GlobalV::NLOCAL, &GlobalV::NLOCAL, &nm, &alpha, a, &GlobalV::NLOCAL, b, &nm, &beta, c, &GlobalV::NLOCAL);

        //3. sum of Inl
        for (int i = 0;i < GlobalV::NLOCAL * GlobalV::NLOCAL;++i)
        {
            this->H_V_delta[i] += c[i];
        }
    }
    delete[] tmp_v1;
    delete[] tmp_v2;

    GlobalV::ofs_running << " Finish calculating H_V_delta" << std::endl;
    return;
}

//for GAMMA_ONLY, search adjacent atoms from I0
void LCAO_Descriptor::build_v_delta_alpha(const bool& calc_deri)
{
    ModuleBase::TITLE("LCAO_Descriptor", "build_v_delta_alpha");
    ModuleBase::GlobalFunc::ZEROS(this->H_V_delta, GlobalV::NLOCAL * GlobalV::NLOCAL); //init before calculate

    std::vector<double> Rcut;
	for(int it1=0; it1<GlobalC::ucell.ntype; ++it1)
        Rcut.push_back(GlobalC::ORB.Phi[it1].getRcut() + GlobalC::ORB.Alpha[0].getRcut());
    
    for (int T0 = 0;T0 < GlobalC::ucell.ntype;++T0)
    {
        for (int I0 = 0;I0 < GlobalC::ucell.atoms[T0].na;++I0)
        {
            const ModuleBase::Vector3<double> tau0 = GlobalC::ucell.atoms[T0].tau[I0];
            //Rcut in this function may need to be changed ?! (I think the range of adjacent atoms here should be rcut(phi)+rcut(alpha))
            GlobalC::GridD.Find_atom(GlobalC::ucell, tau0, T0, I0);

            //adj atom pairs
            for (int ad1 = 0;ad1 < GlobalC::GridD.getAdjacentNum() + 1;++ad1)
            {
                const int T1 = GlobalC::GridD.getType(ad1);
                const int I1 = GlobalC::GridD.getNatom(ad1);
				//const int iat1 = GlobalC::ucell.itia2iat(T1, I1);
                const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
                const ModuleBase::Vector3<double> tau1 = GlobalC::GridD.getAdjacentTau(ad1);
                const Atom* atom1 = &GlobalC::ucell.atoms[T1];
                const int nw1_tot = atom1->nw * GlobalV::NPOL;
                for (int ad2 = 0;ad2 < GlobalC::GridD.getAdjacentNum() + 1;++ad2)
                {
                    const int T2 = GlobalC::GridD.getType(ad2);
					const int I2 = GlobalC::GridD.getNatom(ad2);
					const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);
					const ModuleBase::Vector3<double> tau2 = GlobalC::GridD.getAdjacentTau(ad2);
					const Atom* atom2 = &GlobalC::ucell.atoms[T2];
                    const int nw2_tot = atom2->nw * GlobalV::NPOL;

                    ModuleBase::Vector3<double> dtau10 = tau1 - tau0;
                    ModuleBase::Vector3<double> dtau20 = tau2 - tau0;
                    double distance10 = dtau10.norm() * GlobalC::ucell.lat0;
                    double distance20 = dtau20.norm() * GlobalC::ucell.lat0;
                    
                    if (distance10 < Rcut[T1] && distance20 < Rcut[T2])
                    {
                        for (int iw1=0; iw1<nw1_tot; ++iw1)
						{
							const int iw1_all = start1 + iw1;
							const int iw1_local = GlobalC::ParaO.trace_loc_row[iw1_all];
							if(iw1_local < 0)continue;
							const int iw1_0 = iw1/GlobalV::NPOL;

							for (int iw2=0; iw2<nw2_tot; ++iw2)
							{
								const int iw2_all = start2 + iw2;
								const int iw2_local = GlobalC::ParaO.trace_loc_col[iw2_all];
								if(iw2_local < 0)continue;
								const int iw2_0 = iw2/GlobalV::NPOL;

								double nlm[3];
								nlm[0] = nlm[1] = nlm[2] = 0.0;

								if(!calc_deri)
								{
									GlobalC::UOT.snap_psialpha(
											nlm, 0, tau1, T1,
											atom1->iw2l[ iw1_0 ], // L1
											atom1->iw2m[ iw1_0 ], // m1
											atom1->iw2n[ iw1_0 ], // N1
											tau2, T2,
											atom2->iw2l[ iw2_0 ], // L2
											atom2->iw2m[ iw2_0 ], // m2
											atom2->iw2n[ iw2_0 ], // n2
											GlobalC::ucell.atoms[T0].tau[I0], T0, I0, 
                                            this->inl_index,
                                            this->gedm);
                                    //GlobalC::LM.set_HSgamma(iw1_all, iw2_all, nlm[0], 'L');
                                    int index = iw2_all * GlobalV::NLOCAL + iw1_all;     //for genelpa
                                    this->H_V_delta[index] += nlm[0];
                                }
								else  // calculate force
								{
									GlobalC::UOT.snap_psialpha(
											nlm, 1, tau1, T1,
											atom1->iw2l[ iw1_0 ], // L1
											atom1->iw2m[ iw1_0 ], // m1
											atom1->iw2n[ iw1_0 ], // N1
											tau2, T2,
											atom2->iw2l[ iw2_0 ], // L2
											atom2->iw2m[ iw2_0 ], // m2
											atom2->iw2n[ iw2_0 ], // n2
											GlobalC::ucell.atoms[T0].tau[I0], T0, I0, 
                                            this->inl_index,
                                            this->gedm);
                                    //for Pulay Force
                                    //GlobalC::LM.set_force(iw1_all, iw2_all, nlm[0], nlm[1], nlm[2], 'N');
                                    const int ir = GlobalC::ParaO.trace_loc_row[ iw1_all ];
                                    const int ic = GlobalC::ParaO.trace_loc_col[ iw2_all ];
                                    const long index = ir * GlobalC::ParaO.ncol + ic;
                                    this->DH_V_delta_x[index] += nlm[0];
                                    this->DH_V_delta_y[index] += nlm[1];
                                    this->DH_V_delta_z[index] += nlm[2];
                                }
							}// end iw2
						}// end iw1
					} // end distance
                }//end ad2
            }//end ad1
        }//end I0
    }//end T0
    return;
}

//for multi-k, search adjacent atoms from mu
void LCAO_Descriptor::build_v_delta_mu(const bool& calc_deri)
{
    ModuleBase::TITLE("LCAO_Descriptor", "build_v_delta_mu");
    ModuleBase::GlobalFunc::ZEROS(this->H_V_delta, GlobalV::NLOCAL * GlobalV::NLOCAL); //init before calculate
    //timer::tick ("LCAO_gen_fixedH","build_Nonlocal_mu");

    // < phi1 | beta > < beta | phi2 >
	// phi1 is within the unitcell.
	// while beta is in the supercell.
	// while phi2 is in the supercell.

	int nnr = 0;
	ModuleBase::Vector3<double> tau1, tau2, dtau;
	ModuleBase::Vector3<double> dtau1, dtau2, tau0;
	double distance = 0.0;
	double distance1, distance2;
	double rcut = 0.0;
	double rcut1, rcut2;
		
//	Record_adj RA;
//	RA.for_2d();

	// psi1
    for (int T1 = 0; T1 < GlobalC::ucell.ntype; ++T1)
    {
		const Atom* atom1 = &GlobalC::ucell.atoms[T1];
        for (int I1 =0; I1< atom1->na; ++I1)
        {
            //GlobalC::GridD.Find_atom( atom1->tau[I1] );
            GlobalC::GridD.Find_atom(GlobalC::ucell, atom1->tau[I1] ,T1, I1);
			//const int iat1 = GlobalC::ucell.itia2iat(T1, I1);
			const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
            tau1 = atom1->tau[I1];

			// psi2
            for (int ad2=0; ad2<GlobalC::GridD.getAdjacentNum()+1; ++ad2)
			{
				const int T2 = GlobalC::GridD.getType(ad2);
				const Atom* atom2 = &GlobalC::ucell.atoms[T2];
                
				const int I2 = GlobalC::GridD.getNatom(ad2);
				//const int iat2 = GlobalC::ucell.itia2iat(T2, I2);
                const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);
                tau2 = GlobalC::GridD.getAdjacentTau(ad2);

				bool is_adj = false;
					
				dtau = tau2 - tau1;
				distance = dtau.norm() * GlobalC::ucell.lat0;
				// this rcut is in order to make nnr consistent 
				// with other matrix.
				rcut = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut();
				if(distance < rcut) is_adj = true;
				else if(distance >= rcut)
				{
                    for (int ad0 = 0; ad0 < GlobalC::GridD.getAdjacentNum()+1; ++ad0)
                    {
						const int T0 = GlobalC::GridD.getType(ad0);
						//const int I0 = GlobalC::GridD.getNatom(ad0);
						//const int T0 = RA.info[iat1][ad0][3];
						//const int I0 = RA.info[iat1][ad0][4];
                        //const int iat0 = GlobalC::ucell.itia2iat(T0, I0);
                        //const int start0 = GlobalC::ucell.itiaiw2iwt(T0, I0, 0);

                        tau0 = GlobalC::GridD.getAdjacentTau(ad0);
                        dtau1 = tau0 - tau1;
                        dtau2 = tau0 - tau2;

                        double distance1 = dtau1.norm() * GlobalC::ucell.lat0;
                        double distance2 = dtau2.norm() * GlobalC::ucell.lat0;

                        rcut1 = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Alpha[0].getRcut();
                        rcut2 = GlobalC::ORB.Phi[T2].getRcut() + GlobalC::ORB.Alpha[0].getRcut();

                        if( distance1 < rcut1 && distance2 < rcut2 )
                        {
                            is_adj = true;
                            break;
                        }
                    }
				}


				if(is_adj)
				{
					// < psi1 | all projectors | psi2 >
					// ----------------------------- enter the nnr increaing zone -------------------------
					for (int j=0; j<atom1->nw*GlobalV::NPOL; j++)
					{
						const int j0 = j/GlobalV::NPOL;//added by zhengdy-soc
						const int iw1_all = start1 + j;
						const int mu = GlobalC::ParaO.trace_loc_row[iw1_all];
						if(mu < 0)continue; 

						// fix a serious bug: atom2[T2] -> atom2
						// mohan 2010-12-20
						for (int k=0; k<atom2->nw*GlobalV::NPOL; k++)
						{
							const int k0 = k/GlobalV::NPOL;
							const int iw2_all = start2 + k;
							const int nu = GlobalC::ParaO.trace_loc_col[iw2_all];						
							if(nu < 0)continue;


							//(3) run over all projectors in nonlocal pseudopotential.
							for (int ad0=0; ad0 < GlobalC::GridD.getAdjacentNum()+1 ; ++ad0)
							{
								const int T0 = GlobalC::GridD.getType(ad0);
                                const int I0 = GlobalC::GridD.getNatom(ad0);
								tau0 = GlobalC::GridD.getAdjacentTau(ad0);

								dtau1 = tau0 - tau1;
								dtau2 = tau0 - tau2;
								distance1 = dtau1.norm() * GlobalC::ucell.lat0;
								distance2 = dtau2.norm() * GlobalC::ucell.lat0;

								// seems a bug here!! mohan 2011-06-17
								rcut1 = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Alpha[0].getRcut();
								rcut2 = GlobalC::ORB.Phi[T2].getRcut() + GlobalC::ORB.Alpha[0].getRcut();

								if(distance1 < rcut1 && distance2 < rcut2)
								{
									//const Atom* atom0 = &GlobalC::ucell.atoms[T0];
									double nlm[3]={0,0,0};
									if(!calc_deri)
									{
										GlobalC::UOT.snap_psialpha(
												nlm, 0, tau1, T1,
												atom1->iw2l[ j0 ], // L1
												atom1->iw2m[ j0 ], // m1
												atom1->iw2n[ j0 ], // N1
												tau2, T2,
												atom2->iw2l[ k0 ], // L2
												atom2->iw2m[ k0 ], // m2
												atom2->iw2n[ k0 ], // n2
												tau0, T0, I0,
                                                this->inl_index,
                                                this->gedm);


										if(GlobalV::GAMMA_ONLY_LOCAL)
										{
											// mohan add 2010-12-20
											if( nlm[0]!=0.0 )
											{
                                                //GlobalC::LM.set_HSgamma(iw1_all,iw2_all,nlm[0],'N');//N stands for nonlocal.
                                                int index = iw2_all * GlobalV::NLOCAL + iw1_all;     //for genelpa
                                                this->H_V_delta[index] += nlm[0];
                                            }
										}
										else
                                        {
                                            //for multi-k, not prepared yet 
                                        }
									}// calc_deri
									else // calculate the derivative
									{
										if(GlobalV::GAMMA_ONLY_LOCAL)
										{
                                            GlobalC::UOT.snap_psialpha(
                                                    nlm, 1, tau1, T1,
                                                    atom1->iw2l[ j0 ], // L1
                                                    atom1->iw2m[ j0 ], // m1
                                                    atom1->iw2n[ j0 ], // N1
                                                    tau2, T2,
                                                    atom2->iw2l[ k0 ], // L2
                                                    atom2->iw2m[ k0 ], // m2
                                                    atom2->iw2n[ k0 ], // n2
                                                    tau0, T0, I0,
                                                    this->inl_index,
                                                    this->gedm);

											// sum all projectors for one atom.
											//GlobalC::LM.set_force (iw1_all, iw2_all,	nlm[0], nlm[1], nlm[2], 'N');
										}
										else
										{
											// mohan change the order on 2011-06-17
											// origin: < psi1 | beta > < beta | dpsi2/dtau >
											//now: < psi1/dtau | beta > < beta | psi2 >
											GlobalC::UOT.snap_psialpha(
													nlm, 1, tau2, T2,
													atom2->iw2l[ k0 ], // L2
													atom2->iw2m[ k0 ], // m2
													atom2->iw2n[ k0 ], // n2
													tau1, T1,
													atom1->iw2l[ j0 ], // L1
													atom1->iw2m[ j0 ], // m1
													atom1->iw2n[ j0 ], // N1
                                                    tau0, T0, I0,
                                                    this->inl_index,
                                                    this->gedm);

											//GlobalC::LM.DHloc_fixedR_x[nnr] += nlm[0];
											//GlobalC::LM.DHloc_fixedR_y[nnr] += nlm[1];
											//GlobalC::LM.DHloc_fixedR_z[nnr] += nlm[2];
										}
									}//!calc_deri
								}// distance
							} // ad0
							++nnr;
						}// k
					} // j 
				}// end is_adj
			} // ad2
		} // I1
	} // T1

    //timer::tick ("LCAO_gen_fixedH","build_Nonlocal_mu");
	return;
}


void LCAO_Descriptor::add_v_delta(void)
{
    ModuleBase::TITLE("LCAO_DESCRIPTOR", "add_v_delta");

    if (GlobalV::GAMMA_ONLY_LOCAL)
    {
        for (int iw1 = 0;iw1 < GlobalV::NLOCAL;++iw1)
        {
            for (int iw2 = 0;iw2 < GlobalV::NLOCAL;++iw2)
            {
				if (!GlobalC::ParaO.in_this_processor(iw1,iw2))
				{
					continue;
				}
                GlobalC::LM.set_HSgamma(iw1, iw2, this->H_V_delta[iw1 * GlobalV::NLOCAL + iw2], 'L');
            }
        }
    }
    else
    {
		ModuleBase::WARNING_QUIT("add_v_delta","not implemented yet.");
        //call set_HSk, complex Matrix
    }
	return;
}


void LCAO_Descriptor::cal_f_delta(const ModuleBase::matrix &dm)
{
    ModuleBase::TITLE("LCAO_Descriptor", "cal_f_delta");
    this->F_delta.zero_out();
    //1. cal gedm
    this->cal_gedm(dm);
    
    //2. cal gdmx
    this->init_gdmx();
    this->cal_gdmx(dm);
    
    //3.multiply and sum for each atom
    //3.1 Pulay term 
    // \sum_{Inl}\sum_{mm'} <gedm, gdmx>_{mm'}
    //notice: sum of multiplied corresponding element(mm') , not matrix multiplication !
    int iat = 0;    //check if the index same as GlobalC::ucell.iw2iat or not !!
    for (int it = 0;it < GlobalC::ucell.ntype;++it)
    {
        for (int ia = 0;ia < GlobalC::ucell.atoms[it].na;++ia)
        {
            for (int inl = 0;inl < this->inlmax;++inl)
            {
                int nm = 2 * inl_l[inl] + 1;
                for (int m1 = 0;m1 < nm;++m1)
                {
                    for (int m2 = 0; m2 < nm;++m2)
                    {
                        this->F_delta(iat, 0) += this->gedm[inl][m1 * nm + m2] * gdmx[iat][inl][m1 * nm + m2];
                        this->F_delta(iat, 1) += this->gedm[inl][m1 * nm + m2] * gdmy[iat][inl][m1 * nm + m2];
                        this->F_delta(iat, 2) += this->gedm[inl][m1 * nm + m2] * gdmz[iat][inl][m1 * nm + m2];
                    }
                }
            }//end inl
            ++iat;
        }
    }
    this->print_F_delta("F_delta_pulay_old.dat");
    this->F_delta.zero_out();
    iat = 0;
    for (int it = 0;it < GlobalC::ucell.ntype;++it)
    {
        for (int ia = 0;ia < GlobalC::ucell.atoms[it].na;++ia)
        {
            //3.2 HF term
            double** ss = this->S_mu_alpha;
            double** dsx = this->DS_mu_alpha_x;
            double** dsy = this->DS_mu_alpha_y;
            double** dsz = this->DS_mu_alpha_z;
            for (int mu = 0;mu < GlobalV::NLOCAL;++mu)
            {
                for (int nu = 0;nu < GlobalV::NLOCAL;++nu)
                {
                    for (int l = 0;l <= GlobalC::ORB.Alpha[0].getLmax();++l)
                    {
                        for (int n = 0;n < GlobalC::ORB.Alpha[0].getNchi(l);++n)
                        {
                            for (int m1 = 0;m1 < 2 * l + 1;++m1)
                            {
                                for (int m2 = 0;m2 < 2 * l + 1;++m2)
                                {
                                    if (GlobalV::KS_SOLVER == "genelpa" || GlobalV::KS_SOLVER == "scalapack_gvx")
                                    {
                                        this->F_delta(iat, 0) -= 2*dm(mu, nu) * dsx[inl_index[it](ia, l, n)][m1 * GlobalV::NLOCAL + mu]
                                            * this->gedm[inl_index[it](ia, l, n)][m1 * (2 * l + 1) + m2] * ss[inl_index[it](ia, l, n)][m2 * GlobalV::NLOCAL + nu];
                                        this->F_delta(iat, 1) -= 2*dm(mu, nu) * dsy[inl_index[it](ia, l, n)][m1 * GlobalV::NLOCAL + mu]
                                            * this->gedm[inl_index[it](ia, l, n)][m1 * (2 * l + 1) + m2] * ss[inl_index[it](ia, l, n)][m2 * GlobalV::NLOCAL + nu];
                                        this->F_delta(iat, 2) -= 2*dm(mu, nu) * dsz[inl_index[it](ia, l, n)][m1 * GlobalV::NLOCAL + mu]
                                            * this->gedm[inl_index[it](ia, l, n)][m1 * (2 * l + 1) + m2] * ss[inl_index[it](ia, l, n)][m2 * GlobalV::NLOCAL + nu];
                                    }
                                    else
                                    {
                                        this->F_delta(iat, 0) -= 2*dm(mu, nu) * dsx[inl_index[it](ia, l, n)][mu* (2*l+1) + m1]
                                            * this->gedm[inl_index[it](ia, l, n)][m1 * (2 * l + 1) + m2] * ss[inl_index[it](ia, l, n)][nu* (2*l+1) + m2];
                                        this->F_delta(iat, 1) -= 2*dm(mu, nu) * dsy[inl_index[it](ia, l, n)][mu* (2*l+1) + m1]
                                            * this->gedm[inl_index[it](ia, l, n)][m1 * (2 * l + 1) + m2] * ss[inl_index[it](ia, l, n)][nu* (2*l+1) + m2];
                                        this->F_delta(iat, 2) -= 2*dm(mu, nu) * dsz[inl_index[it](ia, l, n)][mu* (2*l+1) + m1]
                                            * this->gedm[inl_index[it](ia, l, n)][m1 * (2 * l + 1) + m2] * ss[inl_index[it](ia, l, n)][nu* (2*l+1) + m2];
                                    }
                                }//end m2
                            }//end m1
                        }//end n
                    }//end l
                }//end nu
            }//end mu
            ++iat;
        }//end ia
    }//end it
    this->print_F_delta("F_delta_hf_old.dat");
    //3.3 Overlap term
    //somthing in NN, which not included in Hamiltonian
    /*
    for (int mu = 0;mu < GlobalV::NLOCAL;++mu)
    {
        const int iat = GlobalC::ucell.iwt2iat[mu];
        for (int nu = 0;nu < GlobalV::NLOCAL;++nu)
        {
            this->F_delta(iat, 0) += 2*(this->E_delta - this->e_delta_band)* dm(mu, nu) * GlobalC::LM.DSloc_x[mu * GlobalV::NLOCAL + nu];
            this->F_delta(iat, 1) += 2*(this->E_delta - this->e_delta_band) * dm(mu, nu) * GlobalC::LM.DSloc_y[mu * GlobalV::NLOCAL + nu];
            this->F_delta(iat, 2) += 2*(this->E_delta - this->e_delta_band) * dm(mu, nu) * GlobalC::LM.DSloc_z[mu * GlobalV::NLOCAL + nu];
        }
    }*/
    this->del_gdmx();
    return;
}

void LCAO_Descriptor::cal_f_delta_hf(const ModuleBase::matrix& dm)
{
    ModuleBase::TITLE("LCAO_Descriptor", "cal_f_delta_hf");
    this->F_delta.zero_out();
    for (int iat = 0; iat < GlobalC::ucell.nat; ++iat)
    {
        const int it = GlobalC::ucell.iat2it[iat];
        const int ia = GlobalC::ucell.iat2ia[iat];
        const ModuleBase::Vector3<double> tau0 = GlobalC::ucell.atoms[it].tau[ia];
		GlobalC::GridD.Find_atom(GlobalC::ucell, GlobalC::ucell.atoms[it].tau[ia] ,it, ia);
		const double Rcut_Alpha = GlobalC::ORB.Alpha[0].getRcut();

        //FOLLOWING ARE CONTRIBUTIONS FROM
        //VNL DUE TO PROJECTOR'S DISPLACEMENT
        for (int ad1 =0 ; ad1 < GlobalC::GridD.getAdjacentNum()+1; ad1++)
        {
            const int T1 = GlobalC::GridD.getType (ad1);
            const Atom* atom1 = &GlobalC::ucell.atoms[T1];
            const int I1 = GlobalC::GridD.getNatom (ad1);
            const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
			const ModuleBase::Vector3<double> tau1 = GlobalC::GridD.getAdjacentTau (ad1);
			const double Rcut_AO1 = GlobalC::ORB.Phi[T1].getRcut();

            for (int ad2 =0 ; ad2 < GlobalC::GridD.getAdjacentNum()+1; ad2++)
            {
                const int T2 = GlobalC::GridD.getType (ad2);
                const Atom* atom2 = &GlobalC::ucell.atoms[T2];
                const int I2 = GlobalC::GridD.getNatom (ad2);
                const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);
                const ModuleBase::Vector3<double> tau2 = GlobalC::GridD.getAdjacentTau (ad2);
                const double Rcut_AO2 = GlobalC::ORB.Phi[T2].getRcut();

                const double dist1 = (tau1-tau0).norm() * GlobalC::ucell.lat0;
                const double dist2 = (tau2-tau0).norm() * GlobalC::ucell.lat0;

                if (dist1 > Rcut_Alpha + Rcut_AO1
                        || dist2 > Rcut_Alpha + Rcut_AO2)
                {
                    continue;
                }

                for (int jj = 0; jj < GlobalC::ucell.atoms[T1].nw; jj++)
                {
                    const int iw1_all = start1 + jj;
                    const int mu = GlobalC::ParaO.trace_loc_row[iw1_all];
                    if(mu<0) continue;
                    for (int kk = 0; kk < GlobalC::ucell.atoms[T2].nw; kk++)
                    {
                        const int iw2_all = start2 + kk;
                        const int nu = GlobalC::ParaO.trace_loc_col[iw2_all];
                        if(nu<0) continue;
                    
                        double nlm[3] = {0,0,0};
                                
                        GlobalC::UOT.snap_psialpha(
                            nlm, 1,
                            tau1, T1,
                            atom1->iw2l[jj], // L2
                            atom1->iw2m[jj], // m2
                            atom1->iw2n[jj], // N2
                            tau2, T2,
                            atom2->iw2l[kk], // L1
                            atom2->iw2m[kk], // m1
                            atom2->iw2n[kk], // n1
                            tau0, it, ia, 
                            this->inl_index,
                            this->gedm); // mohan  add 2021-05-07

                        double nlm1[3] = {0,0,0};

                        const int index = mu * GlobalC::ParaO.ncol + nu;

                        // HF term is minus, only one projector for each atom force.

                        double sum_dm = 0.0;
                        //remaining: sum for is
                        this->F_delta(iat, 0) -= 2 * dm(mu, nu) * nlm[0];
                        this->F_delta(iat, 1) -= 2 * dm(mu, nu) * nlm[1];
                        this->F_delta(iat, 2) -= 2 * dm(mu, nu) * nlm[2];
                        //this->F_delta(iat, 0) -= 4 * dm(mu, nu) * nlm[0];   //2 for v_delta(not calculated togethor), 2 for e_delta
                        //this->F_delta(iat, 1) -= 4 * dm(mu, nu) * nlm[1];
                        //this->F_delta(iat, 2) -= 4 * dm(mu, nu) * nlm[2];
                    }//!kk
                }//!ad2
            }//!jj
        }//!ad1
    }//!iat
    return;
}
void LCAO_Descriptor::cal_f_delta_pulay(const ModuleBase::matrix& dm)
{
    ModuleBase::TITLE("LCAO_Descriptor", "cal_f_delta_pulay");
    //this->F_delta.zero_out();
    this->build_v_delta_alpha(1);
    //this->build_v_delta_mu(1);    //, if multi-k
    for (int i = 0;i < GlobalV::NLOCAL;++i)   //col, diff
    {
        const int iat = GlobalC::ucell.iwt2iat[i];//the atom whose force being calculated
        for (int j = 0;j < GlobalV::NLOCAL;++j)   //row
        {
            const int mu = GlobalC::ParaO.trace_loc_row[j];
            const int nu = GlobalC::ParaO.trace_loc_col[i];
            if (mu >= 0 && nu >= 0)
            {
                const int index = mu * GlobalC::ParaO.ncol + nu;
                this->F_delta(iat, 0) += 2 * dm(mu, nu) * this->DH_V_delta_x[index];
                this->F_delta(iat, 1) += 2 * dm(mu, nu) * this->DH_V_delta_y[index];
                this->F_delta(iat, 2) += 2 * dm(mu, nu) * this->DH_V_delta_z[index];
            }
        }
    }
    return;
}
void LCAO_Descriptor::cal_descriptor_tensor(void)
{
    ModuleBase::TITLE("LCAO_Descriptor", "cal_descriptor_tensor");
    //init pdm_tensor and d_tensor
    torch::Tensor tmp;

    //if pdm_tensor and d_tensor is not empty, clear it !!
    if (!this->d_tensor.empty())
    {
        this->d_tensor.erase(this->d_tensor.begin(), this->d_tensor.end());
    }
    if (!this->pdm_tensor.empty())
    {
        this->pdm_tensor.erase(this->pdm_tensor.begin(), this->pdm_tensor.end());
    }

    for (int inl = 0;inl < this->inlmax;++inl)
    {
        int nm = 2 * inl_l[inl] + 1;
        tmp = torch::ones({ nm, nm }, torch::TensorOptions().dtype(torch::kFloat64));
        for (int m1 = 0;m1 < nm;++m1)
        {
            for (int m2 = 0;m2 < nm;++m2)
            {
                tmp.index_put_({m1, m2}, this->pdm[inl][m1 * nm + m2]);
            }
        }
        //torch::Tensor tmp = torch::from_blob(this->pdm[inl], { nm, nm }, torch::requires_grad());
        tmp.requires_grad_(true);
        this->pdm_tensor.push_back(tmp);
        this->d_tensor.push_back(torch::ones({ nm }, torch::requires_grad(true)));
    }

    //cal d_tensor
    for (int inl = 0;inl < inlmax;++inl)
    {
        torch::Tensor vd;
        std::tuple<torch::Tensor, torch::Tensor> d_v(this->d_tensor[inl], vd);
        //d_v = torch::symeig(pdm_tensor[inl], /*eigenvalues=*/true, /*upper=*/true);
        d_v = torch::linalg::eigh(pdm_tensor[inl], /*uplo*/"U");
        d_tensor[inl] = std::get<0>(d_v);
    }
    return;
}

void LCAO_Descriptor::load_model(const string& model_file)
{
    ModuleBase::TITLE("LCAO_Descriptor", "load_model");

    try
	{
        this->module = torch::jit::load(model_file);
    }
    catch (const c10::Error& e)

	{
        std::cerr << "error loading the model" << std::endl;
        return;
    }
	return;
}


void LCAO_Descriptor::cal_gedm(const ModuleBase::matrix &dm)
{
    //using this->pdm_tensor
    ModuleBase::TITLE("LCAO_Descriptor", "cal_gedm");
    //-----prepare for autograd---------
    this->cal_projected_DM(dm);
    this->cal_descriptor();
    this->cal_descriptor_tensor();  //use torch::linalg::eigh
    //-----prepared-----------------------
    //forward
    std::vector<torch::jit::IValue> inputs;
    //input_dim:(natom, des_per_atom)
    inputs.push_back(torch::cat(this->d_tensor, /*dim=*/0).reshape({ GlobalC::ucell.nat, this->des_per_atom }));
    std::vector<torch::Tensor> ec;
    ec.push_back(module.forward(inputs).toTensor());    //Hartree
    this->E_delta = ec[0].item().toDouble() * 2;//Ry; *2 is for Hartree to Ry
    
    //cal gedm
    std::vector<torch::Tensor> gedm_shell;
    gedm_shell.push_back(torch::ones_like(ec[0]));
    this->gedm_tensor = torch::autograd::grad(ec, this->pdm_tensor, gedm_shell, /*retain_grad=*/true, /*create_graph=*/false, /*allow_unused=*/true);

    //gedm_tensor(Hartree) to gedm(Ry)
    for (int inl = 0;inl < inlmax;++inl)
    {
        int nm = 2 * inl_l[inl] + 1;
        for (int m1 = 0;m1 < nm;++m1)
        {
            for (int m2 = 0;m2 < nm;++m2)
            {
                int index = m1 * nm + m2;
                //*2 is for Hartree to Ry
                this->gedm[inl][index] = this->gedm_tensor[inl].index({ m1,m2 }).item().toDouble() * 2;
            }
        }
    }
    return;
}

void LCAO_Descriptor::cal_gvdm()
{
    ModuleBase::TITLE("LCAO_Descriptor", "cal_gvdm");
    //cal gevdm(d(EigenValue(D))/dD)
    int nlmax = inlmax/GlobalC::ucell.nat;
    for (int nl=0;nl<nlmax;++nl)
    {
        std::vector<torch::Tensor> avmmv;
        for (int iat = 0;iat<GlobalC::ucell.nat;++iat)
        {
            int inl = iat*nlmax+nl;
            int nm = 2*this->inl_l[inl]+1;
            //repeat each block for nm times in an additional dimension
            torch::Tensor tmp_x = this->pdm_tensor[inl].reshape({nm, nm}).unsqueeze(0).repeat({nm, 1, 1});
            //torch::Tensor tmp_y = std::get<0>(torch::symeig(tmp_x, true));
            torch::Tensor tmp_y = std::get<0>(torch::linalg::eigh(tmp_x, "U"));
            torch::Tensor tmp_yshell = torch::eye(nm, torch::TensorOptions().dtype(torch::kFloat64));
            std::vector<torch::Tensor> tmp_rpt;     //repeated-pdm-tensor (x)
            std::vector<torch::Tensor> tmp_rdt; //repeated-d-tensor (y)
            std::vector<torch::Tensor> tmp_gst; //gvx-shell
            tmp_rpt.push_back(tmp_x);
            tmp_rdt.push_back(tmp_y);
            tmp_gst.push_back(tmp_yshell);
            std::vector<torch::Tensor> tmp_res;
            tmp_res = torch::autograd::grad(tmp_rdt, tmp_rpt, tmp_gst, false, false, /*allow_unused*/true); //nm(v)**nm*nm
            avmmv.push_back(tmp_res[0]);
        }
        torch::Tensor avmm = torch::stack(avmmv, 0); //nat*nv**nm*nm
        this->gevdm_vector.push_back(avmm);
    }
    assert(this->gevdm_vector.size() == nlmax);
    return;
}

void LCAO_Descriptor::print_H_V_delta(void)
{
    ModuleBase::TITLE("LCAO_Descriptor", "print_H_V_delta");

    ofstream ofs;
    stringstream ss;

    // the parameter 'winput::spillage_outdir' is read from INPUTw.
    ss << winput::spillage_outdir << "/"
       << "H_V_delta.dat";

    if (GlobalV::MY_RANK == 0)
    {
        ofs.open(ss.str().c_str());
    }

    ofs << "E_delta(Ry) from deepks model: " << this->E_delta << std::endl;
    ofs << "E_delta(eV) from deepks model: " << this->E_delta * ModuleBase::Ry_to_eV << std::endl;
    ofs << "H_delta(Ry)(gamma only)) from deepks model: " << std::endl;

    for (int i = 0;i < GlobalV::NLOCAL;++i)
    {
        for (int j = 0;j < GlobalV::NLOCAL;++j)
        {
            ofs<< std::setw(12)<< this->H_V_delta[i * GlobalV::NLOCAL + j] << " ";
        }
        ofs << std::endl;
    }
    ofs << "H_delta(eV)(gamma only)) from deepks model: " << std::endl;

    for (int i = 0;i < GlobalV::NLOCAL;++i)
    {
        for (int j = 0;j < GlobalV::NLOCAL;++j)
        {
            ofs<< std::setw(12)<< this->H_V_delta[i * GlobalV::NLOCAL + j] *ModuleBase::Ry_to_eV<< " ";
        }
        ofs << std::endl;
    }

    GlobalV::ofs_running << " H_delta has been printed to " << ss.str() << std::endl;
    return;
}


void LCAO_Descriptor::print_F_delta(const string& fname)
{
    ModuleBase::TITLE("LCAO_Descriptor", "print_F_delta");

    ofstream ofs;
    stringstream ss;
    // the parameter 'winput::spillage_outdir' is read from INPUTw.
    ss << winput::spillage_outdir << "/"<< fname ;

    if (GlobalV::MY_RANK == 0)
    {
        ofs.open(ss.str().c_str());
    }

    ofs << "F_delta(Hatree/Bohr) from deepks model: " << std::endl;
    ofs << std::setw(12) << "type" << std::setw(12) << "atom" << std::setw(15) << "dF_x" << std::setw(15) << "dF_y" << std::setw(15) << "dF_z" << std::endl;

    for (int it = 0;it < GlobalC::ucell.ntype;++it)
    {
        for (int ia = 0;ia < GlobalC::ucell.atoms[it].na;++ia)
        {
            int iat = GlobalC::ucell.itia2iat(it, ia);
            ofs << std::setw(12) << GlobalC::ucell.atoms[it].label << std::setw(12) << ia
                << std::setw(15) << this->F_delta(iat, 0) / 2 << std::setw(15) << this->F_delta(iat, 1) / 2
                << std::setw(15) << this->F_delta(iat, 2) / 2 << std::endl;
        }
    }

    ofs << "F_delta(eV/Angstrom) from deepks model: " << std::endl;
    ofs << std::setw(12) << "type" << std::setw(12) << "atom" << std::setw(15) << "dF_x" << std::setw(15) << "dF_y" << std::setw(15) << "dF_z" << std::endl;

    for (int it = 0;it < GlobalC::ucell.ntype;++it)
    {
        for (int ia = 0;ia < GlobalC::ucell.atoms[it].na;++ia)
        {
            int iat = GlobalC::ucell.itia2iat(it, ia);
            ofs << std::setw(12) << GlobalC::ucell.atoms[it].label << std::setw(12)
                << ia << std::setw(15) << this->F_delta(iat, 0) * ModuleBase::Ry_to_eV/ModuleBase::BOHR_TO_A
                << std::setw(15) << this->F_delta(iat, 1) * ModuleBase::Ry_to_eV/ModuleBase::BOHR_TO_A
                << std::setw(15) << this->F_delta(iat, 2) * ModuleBase::Ry_to_eV/ModuleBase::BOHR_TO_A << std::endl;
        }
    }

    GlobalV::ofs_running << " F_delta has been printed to " << ss.str() << std::endl;
    ofs.close();

    /*
    //============for test: double check of 2 methods=============== 
    //1. same as old method
    ModuleBase::matrix F_delta_old;
    F_delta_old.create(GlobalC::ucell.nat, 3);
    F_delta_old.zero_out();
    for (int inl = 0;inl < this->inlmax;++inl)
    {
        int nm = 2 * inl_l[inl] + 1;
        for (int m1 = 0;m1 < nm;++m1)
        {
            for (int m2 = 0;m2 < nm;++m2)
            {
                for (int mu = 0;mu < GlobalV::NLOCAL;++mu)
                {
                    int iat = GlobalC::ucell.iwt2iat[mu];
                    for (int nu = 0;nu < GlobalV::NLOCAL;++nu)
                    {
                        F_delta_old(iat, 0) += 2*LOC.wfc_dm_2d.dm_gamma[0](mu, nu)*gedm[inl][m1 * nm + m2] * this->DS_mu_alpha_x[inl][m1 * GlobalV::NLOCAL + mu] * this->S_mu_alpha[inl][m2 * GlobalV::NLOCAL + nu];
                        F_delta_old(iat, 1) += 2*LOC.wfc_dm_2d.dm_gamma[0](mu, nu)*gedm[inl][m1 * nm + m2] * this->DS_mu_alpha_y[inl][m1 * GlobalV::NLOCAL + mu] * this->S_mu_alpha[inl][m2 * GlobalV::NLOCAL + nu];
                        F_delta_old(iat,2) += 2*LOC.wfc_dm_2d.dm_gamma[0](mu, nu)*gedm[inl][m1 * nm + m2] * this->DS_mu_alpha_z[inl][m1 * GlobalV::NLOCAL + mu] * this->S_mu_alpha[inl][m2 * GlobalV::NLOCAL + nu];
                    }
                }
            }
        }
    }
    //print F_new
    stringstream ss1;
    ss1 << winput::spillage_outdir << "/"
       << "F_delta_old.dat";
    if (GlobalV::MY_RANK == 0)
    {
        ofs.open(ss1.str().c_str());
    }
    for (int iat = 0;iat < GlobalC::ucell.nat;++iat)
    {
        for (int i = 0;i < 3;++i)
        {
            ofs<< std::setw(8)<< F_delta_old(iat,i)/2<< " ";
        }
        ofs << std::endl;
    }
    ofs.close();

    //2. same as new method
    ModuleBase::matrix F_delta_new;
    F_delta_new.create(GlobalC::ucell.nat, 3);
    F_delta_new.zero_out();
    for (int inl = 0;inl < this->inlmax;++inl)
    {
        int nm = 2 * inl_l[inl] + 1;
        for (int m1 = 0;m1 < nm;++m1)
        {
            for (int m2 = 0;m2 < nm;++m2)
            {
                for (int mu = 0;mu < GlobalV::NLOCAL;++mu)
                {
                    for (int nu = 0;nu < GlobalV::NLOCAL;++nu)
                    {
                        int iat = GlobalC::ucell.iwt2iat[nu];
                        F_delta_new(iat, 0) += 2 * LOC.wfc_dm_2d.dm_gamma[0](mu, nu) * gedm[inl][m1 * nm + m2] * this->DS_mu_alpha_x[inl][m1 * GlobalV::NLOCAL + mu] * this->S_mu_alpha[inl][m2 * GlobalV::NLOCAL + nu];
                        F_delta_new(iat, 1) += 2*LOC.wfc_dm_2d.dm_gamma[0](mu, nu)*gedm[inl][m1 * nm + m2] * this->DS_mu_alpha_y[inl][m1 * GlobalV::NLOCAL + mu] * this->S_mu_alpha[inl][m2 * GlobalV::NLOCAL + nu];
                        F_delta_new(iat,2) += 2*LOC.wfc_dm_2d.dm_gamma[0](mu, nu)*gedm[inl][m1 * nm + m2] * this->DS_mu_alpha_z[inl][m1 * GlobalV::NLOCAL + mu] * this->S_mu_alpha[inl][m2 * GlobalV::NLOCAL + nu];
                    }
                }
            }
        }
    }
    //print F_new
    stringstream ss2;
    ss2 << winput::spillage_outdir << "/"
       << "F_delta_new.dat";
    if (GlobalV::MY_RANK == 0)
    {
        ofs.open(ss2.str().c_str());
    }
    for (int iat = 0;iat < GlobalC::ucell.nat;++iat)
    {
        for (int i = 0;i < 3;++i)
        {
            ofs<< std::setw(8)<< F_delta_new(iat,i)/2<< " ";
        }
        ofs << std::endl;
    }
    ofs.close();
*/
    return;
}


void LCAO_Descriptor::save_npy_d(void)
{
    ModuleBase::TITLE("LCAO_Descriptor", "save_npy_d");
    //save descriptor in .npy format
    vector<double> npy_des;
    for (int i = 0;i < this->n_descriptor;++i)
    {
        npy_des.push_back(this->d[i]);
    }
    const long unsigned dshape[] = {(long unsigned) GlobalC::ucell.nat, (long unsigned) this->des_per_atom };
    npy::SaveArrayAsNumpy("dm_eig.npy", false, 2, dshape, npy_des);
    return;
}


void LCAO_Descriptor::save_npy_e(const double &e, const std::string &e_file)
{
    ModuleBase::TITLE("LCAO_Descriptor", "save_npy_e");
    //save e_base
    const long unsigned eshape[] = { 1 };
    vector<double> npy_e;
    npy_e.push_back(e);
    npy::SaveArrayAsNumpy(e_file, false, 1, eshape, npy_e);
    return;
}

void LCAO_Descriptor::save_npy_f(const ModuleBase::matrix &f, const std::string &f_file)
{
    ModuleBase::TITLE("LCAO_Descriptor", "save_npy_f");
    //save f_base
    //caution: unit: Rydberg/Bohr
    const long unsigned fshape[] = {(long unsigned) GlobalC::ucell.nat, 3 };
    vector<double> npy_f;
    for (int iat = 0;iat < GlobalC::ucell.nat;++iat)
    {
        for (int i = 0;i < 3;i++)
        {
            npy_f.push_back(f(iat, i));
        }
    }
    npy::SaveArrayAsNumpy(f_file, false, 2, fshape, npy_f);
    return;
}

void LCAO_Descriptor::save_npy_gvx(void)
{
    ModuleBase::TITLE("LCAO_Descriptor", "save_npy_gvx");
    //save grad_vx.npy (when  force label is in use)
    //unit: /Bohr
    const long unsigned gshape[] = {(long unsigned) GlobalC::ucell.nat, 3, GlobalC::ucell.nat, this->des_per_atom};
    vector<double> npy_gvx;
    for (int ibt = 0;ibt < GlobalC::ucell.nat;++ibt)
    {
        for (int i = 0;i < 3;i++)
        {
            for (int iat = 0;iat < GlobalC::ucell.nat;++iat)
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

void LCAO_Descriptor::cal_e_delta_band(const std::vector<ModuleBase::matrix> &dm)
{
    ModuleBase::TITLE("LCAO_Descriptor", "cal_e_delta_band");
    this->e_delta_band = 0;
    for (int i = 0; i < GlobalV::NLOCAL; ++i)
    {
        for (int j = 0; j < GlobalV::NLOCAL; ++j)
        {
            const int mu = GlobalC::ParaO.trace_loc_row[j];
            const int nu = GlobalC::ParaO.trace_loc_col[i];
            if (mu >= 0 && nu >= 0)
            {
                const int index = mu * GlobalC::ParaO.ncol + nu;
                for (int is = 0; is < GlobalV::NSPIN; ++is)
                {
                    this->e_delta_band += dm[is](nu, mu) * this->H_V_delta[i * GlobalV::NLOCAL + j];
                }
            }
        }
    }
    return;
}

void LCAO_Descriptor::cal_gvx(const ModuleBase::matrix &dm)
{
    ModuleBase::TITLE("LCAO_Descriptor","cal_gvx");
    //preconditions
    this->cal_gvdm();

    this->build_S_descriptor(1);
    this->init_gdmx();
    this->cal_gdmx(dm); //checked

    //make gdmx as tensor
    int nlmax = this->inlmax/GlobalC::ucell.nat;
    for (int nl=0;nl<nlmax;++nl)
    {
        std::vector<torch::Tensor> bmmv;
        for (int ibt=0;ibt<GlobalC::ucell.nat;++ibt)
        {
            std::vector<torch::Tensor> xmmv;
            for (int i=0;i<3;++i)
            {
                std::vector<torch::Tensor> ammv;
                for (int iat=0; iat<GlobalC::ucell.nat; ++iat)
                {
                    int inl = iat*nlmax + nl;
                    int nm = 2*this->inl_l[inl]+1;
                    std::vector<double> mmv;
                    for (int m1=0;m1<nm;++m1)
                    {
                        for(int m2=0;m2<nm;++m2)
                        {
                            if(i==0) mmv.push_back(this->gdmx[ibt][inl][m1*nm+m2]);
                            if(i==1) mmv.push_back(this->gdmy[ibt][inl][m1*nm+m2]);
                            if(i==2) mmv.push_back(this->gdmz[ibt][inl][m1*nm+m2]);
                        }
                    }//nm^2
                    torch::Tensor mm = torch::tensor(mmv, torch::TensorOptions().dtype(torch::kFloat64) ).reshape({nm, nm});    //nm*nm
                    ammv.push_back(mm);
                }
                torch::Tensor amm = torch::stack(ammv, 0);  //nat*nm*nm
                xmmv.push_back(amm);
            }
            torch::Tensor bmm = torch::stack(xmmv, 0);  //3*nat*nm*nm
            bmmv.push_back(bmm); 
        }
        this->gdmr_vector.push_back(torch::stack(bmmv, 0)); //nbt*3*nat*nm*nm
    }
    assert(this->gdmr_vector.size()==nlmax);

    std::cout<<"gdmr-ok"<<std::endl;
    std::cout << nlmax <<" " << this->gdmr_vector.size()<<" "<<this->gevdm_vector.size()<<std::endl;
    //einsum for each inl: 
    std::vector<torch::Tensor> gvx_vector;
    for (int nl = 0;nl<nlmax;++nl)
    {
        gvx_vector.push_back(at::einsum("bxamn, avmn->bxav", {this->gdmr_vector[nl], this->gevdm_vector[nl]}));
    }//
    
    // cat nv-> \sum_nl(nv) = \sum_nl(nm_nl)=des_per_atom
    this->gvx_tensor = torch::cat(gvx_vector, -1);

    assert(this->gvx_tensor.size(0) == GlobalC::ucell.nat);
    assert(this->gvx_tensor.size(1) == 3);
    assert(this->gvx_tensor.size(2) == GlobalC::ucell.nat);
    assert(this->gvx_tensor.size(3) == this->des_per_atom);

    return;
}
#endif