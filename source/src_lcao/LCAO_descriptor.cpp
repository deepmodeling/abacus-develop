//caoyu add 2021-03-29

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


LCAO_Descriptor::LCAO_Descriptor(int lm, int inlm):lmaxd(lm), inlmax(inlm)
{
    alpha_index = new IntArray[1];
    inl_index = new IntArray[1];
    inl_l = new int[1];
    d = new double[1];
    H_V_delta = new double[1];
    F_delta.create(ucell.nat, 3);

    //init S_mu_alpha**
    this->S_mu_alpha = new double* [this->inlmax];    //inl
    for (int inl = 0;inl < this->inlmax;inl++)
    {
        this->S_mu_alpha[inl] = new double[NLOCAL * (2 * this->lmaxd + 1)];     //NLOCAL*nm
        ZEROS(S_mu_alpha[inl], NLOCAL * (2 * this->lmaxd+ 1));
    }
    
    //init DS_mu_alpha**
    this->DS_mu_alpha_x = new double* [this->inlmax];
    this->DS_mu_alpha_y = new double* [this->inlmax];
    this->DS_mu_alpha_z = new double* [this->inlmax];
    for (int inl = 0;inl < this->inlmax;inl++)
    {
        this->DS_mu_alpha_x[inl] = new double[NLOCAL * (2 * this->lmaxd + 1)];
        this->DS_mu_alpha_y[inl] = new double[NLOCAL * (2 * this->lmaxd + 1)];
        this->DS_mu_alpha_z[inl] = new double[NLOCAL * (2 * this->lmaxd + 1)];
        ZEROS(DS_mu_alpha_x[inl], NLOCAL * (2 * this->lmaxd + 1));
        ZEROS(DS_mu_alpha_y[inl], NLOCAL * (2 * this->lmaxd + 1));
        ZEROS(DS_mu_alpha_z[inl], NLOCAL * (2 * this->lmaxd + 1));
    }
    //init pdm**
    const int PDM_size = (lmaxd * 2 + 1) * (lmaxd * 2 + 1);
    this->pdm = new double* [inlmax];
    for (int inl = 0;inl < inlmax;inl++)
    {
        this->pdm[inl] = new double[PDM_size];
        ZEROS(this->pdm[inl], PDM_size);
    }
    //init gedm**
    this->gedm = new double* [inlmax];
    for (int inl = 0;inl < inlmax;inl++)
    {
        this->gedm[inl] = new double[PDM_size];
        ZEROS(this->gedm[inl], PDM_size);
    }
}
LCAO_Descriptor::~LCAO_Descriptor()
{
    delete[] alpha_index;
    delete[] inl_index;
    delete[] inl_l;
    delete[] d;
    delete[] H_V_delta;

    //delete S_mu_alpha**
    for (int inl = 0;inl < this->inlmax;inl++)
    {
        delete[] S_mu_alpha[inl];
    }
    delete[] S_mu_alpha;

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

    //delete pdm**
    for (int inl = 0;inl < this->inlmax;inl++)
    {
        delete[] pdm[inl];
    }
    delete[] pdm;
    //delete gedm**
    for (int inl = 0;inl < this->inlmax;inl++)
    {
        delete[] gedm[inl];
    }
    delete[] gedm;
}

void LCAO_Descriptor::build_S_descriptor(const bool &calc_deri)
{
    TITLE("LCAO_Descriptor", "build_S_descriptor");
    // =======init==============
    // cal n(descriptor) per atom , related to Lmax, nchi(L) and m. (not total_nchi!)
	this->des_per_atom=0; // mohan add 2021-04-21
    for (int l = 0; l <= this->lmaxd; l++)
    {
        this->des_per_atom += ORB.Alpha[0].getNchi(l) * (2 * l + 1);
    }
    this->n_descriptor = ucell.nat * this->des_per_atom;

    this->init_index();
    
    // =======init==============

    //array to store data
    double olm[3] = {0.0, 0.0, 0.0};

    //\sum{T} e**{ikT} <\phi_{ia}|d\phi_{k\beta}(T)>	//???
    Vector3<double> tau1, tau2, dtau;
    Vector3<double> dtau1, dtau2, tau0;
    for (int T1 = 0; T1 < ucell.ntype; ++T1)
    {
        Atom *atom1 = &ucell.atoms[T1];
        for (int I1 = 0; I1 < atom1->na; ++I1)
        {
            tau1 = atom1->tau[I1];
            //GridD.Find_atom(tau1);
            GridD.Find_atom(ucell, tau1, T1, I1);

            for (int ad = 0; ad < GridD.getAdjacentNum() + 1; ++ad)
            {
                const int T2 = GridD.getType(ad);
                const int I2 = GridD.getNatom(ad);
                //Atom *atom2 = &ucell.atoms[T2];
                tau2 = GridD.getAdjacentTau(ad);
                dtau = tau2 - tau1;
                double distance = dtau.norm() * ucell.lat0;
                double rcut = ORB.Phi[T1].getRcut() + ORB.Alpha[0].getRcut(); //Rcut is subject to ORB.Phi to keep dimension of S_mu_alpha same as Sloc
                if (distance < rcut)
                {
                    int iw1_all = ucell.itiaiw2iwt(T1, I1, 0); //iw1_all = combined index (it, ia, iw)

                    for (int jj = 0; jj < atom1->nw * NPOL; ++jj)
                    {
                        const int jj0 = jj / NPOL;
                        const int L1 = atom1->iw2l[jj0];
                        const int N1 = atom1->iw2n[jj0];
                        const int m1 = atom1->iw2m[jj0];

                        for (int L2 = 0; L2 <= ORB.Alpha[0].getLmax(); ++L2)
                        {
                            for (int N2 = 0; N2 < ORB.Alpha[0].getNchi(L2); ++N2)
                            {
                                for (int m2 = 0; m2 < 2 * L2 + 1; ++m2)
                                {
                                    olm[0] = olm[1] = olm[2] = 0.0;
                                    if (!calc_deri)
                                    {
                                        UOT.snap_psipsi(olm, 0, 'D', tau1,
                                                T1, L1, m1, N1, GridD.getAdjacentTau(ad),
                                                T2, L2, m2, N2, NSPIN);
                                        if (GAMMA_ONLY_LOCAL)
                                        {
                                            this->set_S_mu_alpha(iw1_all, inl_index[T2](I2,L2,N2), m2, olm[0]);
                                        }
                                    }
                                    else
                                    {
                                        UOT.snap_psipsi(olm, 1, 'D', tau1,
                                            T1, L1, m1, N1, GridD.getAdjacentTau(ad),
                                            T2, L2, m2, N2, NSPIN);
                                        if (GAMMA_ONLY_LOCAL)
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
    if (!GAMMA_ONLY_LOCAL)
    {
        WARNING_QUIT("LCAO_Descriptor::build_S_descriptor", "muti-kpoint method for descriptor is not implemented yet! ");
    }
    return;
}

void LCAO_Descriptor::set_S_mu_alpha(const int &iw1_all, const int &inl, const int &im, const double &v)
{
    //const int ir = ParaO.trace_loc_row[iw1_all];
    //const int ic = ParaO.trace_loc_col[iw2_all];
    //no parellel yet
    const int ir = iw1_all;
    const int ic = im;
    //const int index = ir * ParaO.ncol + ic;
    int index;
    if (KS_SOLVER == "genelpa" || KS_SOLVER == "scalapack_gvx") // save the matrix as column major format
    {
        index = ic * NLOCAL + ir;
    }
    else
    {
        index = ir * (2*inl_l[inl]+1)  + ic; //row: lcao orbitals; col: descriptor basis
    }
    this->S_mu_alpha[inl][index] += v;
    return;
}

void LCAO_Descriptor::cal_projected_DM()
{
    TITLE("LCAO_Descriptor", "cal_projected_DM");
    //step 1: get dm: the coefficient of wfc, not charge density
    double *dm = new double[NLOCAL * NLOCAL];
    ZEROS(dm, NLOCAL * NLOCAL);
    this->getdm(dm);

    //step 2: get S_alpha_mu and S_nu_beta
    double **ss = this->S_mu_alpha;

    //step 3 : multiply: cal ST*DM*S
    
    //init tmp_pdm*
    const int tmp_PDM_size = NLOCAL * (lmaxd*2+1);
    double* tmp_pdm = new double[tmp_PDM_size];
    ZEROS(tmp_pdm, tmp_PDM_size);

    for (int inl = 0;inl < inlmax;inl++)
    {   
        int nm = 2 * inl_l[inl] + 1;   //1,3,5,...
        const char t = 'T';  //transpose
        const char nt = 'N'; //non transpose
        const double alpha = 1;
        const double beta = 0;
        double *a = dm;
        double *b = ss[inl];
        double *c = tmp_pdm;
        dgemm_(&nt, &nt, &NLOCAL, &nm, &NLOCAL, &alpha, a, &NLOCAL, b, &NLOCAL, &beta, c, &NLOCAL); //DM*S
        a = ss[inl];
        b = c;
        c = this->pdm[inl];
        dgemm_(&t, &nt, &nm, &nm, &NLOCAL, &alpha, a, &NLOCAL, b, &NLOCAL, &beta, c, &nm); //ST*DM*S
    }
    
    delete[] tmp_pdm;
    delete[] dm;
    return;
}

void LCAO_Descriptor::cal_descriptor()
{
    TITLE("LCAO_Descriptor", "cal_descriptor");
    delete[] d;
    d = new double[this->n_descriptor];
    //==========print preparation=============
    ofs_running << " print out each DM_inl" << endl;
    ofstream ofs;
    stringstream ss;
    ss << winput::spillage_outdir << "/"
       << "projected_DM.dat";
    if (MY_RANK == 0)
    {
        ofs.open(ss.str().c_str());
    }
    //==========print preparation=============
    const int lmax = ORB.get_lmax_d();
    int id = 0;
    for (int it = 0; it < ucell.ntype; it++)
    {
        for (int ia = 0; ia < ucell.atoms[it].na; ia++)
        {
            ofs << ucell.atoms[it].label << " atom_index " << ia + 1 << " n_descriptor " << this->des_per_atom << endl;
            for (int l = 0; l <= lmax; l++)
            {
                int nmax = ORB.Alpha[0].getNchi(l);
                for (int n = 0; n < nmax; n++)
                {
                    const int dim = 2 * l + 1;
                    const int inl = inl_index[it](ia, l, n);
                    // descriptor for atom (it, ia)
                    ComplexMatrix des(dim, dim);
                    for (int m = 0; m < dim; m++)
                    {
                        for (int m2 = 0; m2 < dim; m2++)
                        {
                            int index = m * dim + m2;
                            complex<double> tmp(this->pdm[inl][index], 0);
                            des(m, m2) += tmp;
                        }
                    }

                    this->print_projected_DM(ofs, des, it, ia, l, n);

                    //ofs_running << "dimension of des is " << 2 * l + 1 << endl;
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
    this->cal_descriptor_tensor();  //use torch::symeig
    this->print_descriptor();
    return;
}

void LCAO_Descriptor::init_index(void)
{
    TITLE("LCAO_Descriptor", "init_index");
    ofs_running << " Initialize the descriptor index for deepks (lcao line)" << endl;
    const int lmaxd = ORB.get_lmax_d();
    const int nmaxd = ORB.get_nchimax_d();
    assert(lmaxd >= 0);
    assert(nmaxd >= 0);
    ofs_running << " lmax of descriptor = " << lmaxd << endl;
    ofs_running << " nmax of descriptor= " << nmaxd << endl;

    delete[] this->alpha_index;
    this->alpha_index = new IntArray[ucell.ntype];
    delete[] this->inl_index;
    this->inl_index = new IntArray[ucell.ntype];
    delete[] this->inl_l;
    this->inl_l = new int[inlmax];
    ZEROS(inl_l, inlmax);

    int inl = 0;
    int alpha = 0;
    for (int it = 0; it < ucell.ntype; it++)
    {
        this->alpha_index[it].create(
            ucell.atoms[it].na,
            lmaxd + 1, // l starts from 0
            nmaxd,
            2 * lmaxd + 1); // m ==> 2*l+1

        this->inl_index[it].create(
            ucell.atoms[it].na,
            lmaxd + 1,
            nmaxd); 

        ofs_running << "Type " << it + 1
                    << " number_of_atoms " << ucell.atoms[it].na << endl;

        for (int ia = 0; ia < ucell.atoms[it].na; ia++)
        {
            //alpha
            for (int l = 0; l < lmaxd + 1; l++)
            {
                for (int n = 0; n < ORB.Alpha[0].getNchi(l); n++)
                {
                    for (int m = 0; m < 2 * l + 1; m++)
                    {
                        this->alpha_index[it](ia, l, n, m) = alpha;
                        alpha++;
                    }
                    this->inl_index[it](ia, l, n) = inl;
                    inl_l[inl] = l;
                    inl++;
                }
            }
        }//end ia
    }//end it
    assert(this->n_descriptor == alpha);
    assert(ucell.nat * ORB.Alpha[0].getTotal_nchi() == inl);
    ofs_running << "descriptors_per_atom " << this->des_per_atom << endl;
    ofs_running << "total_descriptors " << this->n_descriptor << endl;

    return;
}


void LCAO_Descriptor::print_projected_DM(ofstream& ofs, ComplexMatrix& des, const int& it, const int& ia, const int& l, const int& n)
{
    ofs << "L=" << l << "   N=" << n << endl;
    for (int i = 0; i < 2 * l + 1; i++)
    {
        for (int j = 0; j < 2 * l + 1; j++)
        {
            ofs << des(i, j).real() << " ";
        }
        ofs << endl;
    }
    return;
}
void LCAO_Descriptor::print_descriptor()
{
    TITLE("LCAO_Descriptor", "print_descriptor");
    ofstream ofs;
    stringstream ss;
    // the parameter 'winput::spillage_outdir' is read from INPUTw.
    ss << winput::spillage_outdir << "/"
       << "descriptor.dat";
    if (MY_RANK == 0)
    {
        ofs.open(ss.str().c_str());
    }
    for (int it = 0; it < ucell.ntype; it++)
    {
        for (int ia = 0; ia < ucell.atoms[it].na; ia++)
        {
            ofs << ucell.atoms[it].label << " atom_index " << ia + 1 << " n_descriptor " << this->des_per_atom << endl;
            int id0 = this->alpha_index[it](ia, 0, 0, 0);
            for (int id = id0; id < id0 + this->des_per_atom; ++id)
            {
                if ((id - id0) > 0 && (id - id0) % 8 == 0)
                    ofs << endl;
                ofs << d[id] << " ";
            }
            ofs << endl << endl;
        }
    }
    ofs_running << "descriptors are printed" << endl;

    //print descriptor in .npy format
    vector<double>npy_des;
    for (int i = 0;i < n_descriptor;++i)
    {
        npy_des.push_back(this->d[i]);
    }
    const long unsigned shape [] = {ucell.nat, des_per_atom};
    npy::SaveArrayAsNumpy("dm_eig.npy", false, 2, shape, npy_des);

    return;
}

void LCAO_Descriptor::set_DS_mu_alpha(const int& iw1_all, const int& inl, const int& im,
    const double& vx, const double& vy, const double& vz)
{
    //const int ir = ParaO.trace_loc_row[iw1_all];
    //const int ic = ParaO.trace_loc_col[iw2_all];
    //no parellel yet
    const int ir = iw1_all;
    const int ic = im;
    //const int index = ir * ParaO.ncol + ic;
    int index;
    if (KS_SOLVER == "genelpa" || KS_SOLVER == "scalapack_gvx") // save the matrix as column major format
    {
        index = ic * NLOCAL + ir;
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

void LCAO_Descriptor::getdm(double* dm)
{
    for (int i = 0; i < LOC.wfc_dm_2d.dm_gamma[0].nr; i++)
    {
        for (int j = 0; j < LOC.wfc_dm_2d.dm_gamma[0].nc; j++)
        {
            dm[i * NLOCAL + j] = LOC.wfc_dm_2d.dm_gamma[0](i, j); //only consider default NSPIN = 1
        }
    }
}

void LCAO_Descriptor::cal_gdmx(matrix &dm)
{
    TITLE("LCAO_Descriptor", "cal_gdmx");
    //get DS_alpha_mu and S_nu_beta
    double** ss = this->S_mu_alpha;
    double** dsx = this->DS_mu_alpha_x;
    double** dsy = this->DS_mu_alpha_y;
    double** dsz = this->DS_mu_alpha_z;
    for (int inl = 0;inl < inlmax;inl++)
    {
        //dE/dD will be multiplied in cal_f_delta, here only calculate dD/dx_I
        int nm = 2 * inl_l[inl] + 1;   //1,3,5,...
        for (int i =0; i < NLOCAL;++i) 
        {
            const int iat = ucell.iwt2iat[i];//the atom whose force being calculated
            for (int j= 0;j < NLOCAL; ++j)
            {
                const int mu = ParaO.trace_loc_row[j];
                const int nu = ParaO.trace_loc_col[i];
                if (mu >= 0 && nu >= 0)
                {
                    for (int m1 = 0;m1 < nm;++m1)
                    {
                        for (int m2 = 0;m2 < nm;++m2)
                        {
                            for (int is = 0;is < NSPIN;++is)
                            {
                                if (KS_SOLVER == "genelpa" || KS_SOLVER == "scalapack_gvx") // save the matrix as column major format
                                {
                                    gdmx[iat][inl][m1 * nm + m2] += 4 * dsx[inl][m1 * NLOCAL + mu] * dm(mu, nu) * ss[inl][m2 * NLOCAL + nu];
                                    gdmy[iat][inl][m1 * nm + m2] += 4 * dsy[inl][m1 * NLOCAL + mu] *  dm(mu, nu)  * ss[inl][m2 * NLOCAL + nu];
                                    gdmz[iat][inl][m1 * nm + m2] += 4 * dsz[inl][m1 * NLOCAL + mu] * dm(mu, nu)  * ss[inl][m2 * NLOCAL + nu];
                                }
                                else
                                {
                                    gdmx[iat][inl][m1 * nm + m2] += 4 * dsx[inl][mu * nm + m1] * dm(mu, nu) * ss[inl][nu * nm + m2];
                                    gdmy[iat][inl][m1 * nm + m2] += 4 * dsy[inl][mu * nm + m1] * dm(mu, nu) * ss[inl][nu * nm + m2];
                                    gdmz[iat][inl][m1 * nm + m2] += 4 * dsz[inl][mu * nm + m1] *  dm(mu, nu) * ss[inl][nu * nm + m2];
                                }
                            }
                        }//end m2
                    } //end m1
                }//end if
            }//end nu
        }//end mu
    }//end inl
    return;
}

void LCAO_Descriptor::init_gdmx()
{
    this->gdmx = new double** [ucell.nat];
    this->gdmy = new double** [ucell.nat];
    this->gdmz = new double** [ucell.nat];
    for (int iat = 0;iat < ucell.nat;iat++)
    {
        this->gdmx[iat] = new double* [inlmax];
        this->gdmy[iat] = new double* [inlmax];
        this->gdmz[iat] = new double* [inlmax];
        for (int inl = 0;inl < inlmax;inl++)
        {
            this->gdmx[iat][inl] = new double [(2 * lmaxd + 1) * (2 * lmaxd + 1)];
            this->gdmy[iat][inl] = new double [(2 * lmaxd + 1) * (2 * lmaxd + 1)];
            this->gdmz[iat][inl] = new double[(2 * lmaxd + 1) * (2 * lmaxd + 1)];
            ZEROS(gdmx[iat][inl], (2 * lmaxd + 1) * (2 * lmaxd + 1));
            ZEROS(gdmy[iat][inl], (2 * lmaxd + 1) * (2 * lmaxd + 1));
            ZEROS(gdmz[iat][inl], (2 * lmaxd + 1) * (2 * lmaxd + 1));
        }
    }
    return;
}

void LCAO_Descriptor::del_gdmx()
{
    for (int iat = 0;iat < ucell.nat;iat++)
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

void LCAO_Descriptor::cal_v_delta(const string& model_name)
{
    TITLE("LCAO_Descriptor", "cal_v_delta");
    //1.  (dE/dD)<alpha_m'|psi_nv>
    this->load_model(model_name);
    this->cal_gedm();

    //2. multiply and sum
    double* tmp_v1 = new double[(2 * lmaxd + 1) * NLOCAL];
    ZEROS(tmp_v1, (2 * lmaxd + 1) * NLOCAL);
    double* tmp_v2 = new double[NLOCAL *NLOCAL];
    ZEROS(tmp_v2, NLOCAL * NLOCAL);
    //init H_V_delta
    delete[] this->H_V_delta;
    this->H_V_delta = new double[NLOCAL * NLOCAL];
    
    for (int inl = 0;inl < inlmax;inl++)
    {
        int nm = 2 * inl_l[inl] + 1;   //1,3,5,...
        const char t = 'T';  //transpose
        const char nt = 'N'; //non transpose
        const double alpha = 1;
        const double beta = 0;
        double* a = this->gedm[inl];//[nm][nm]
        double* b = S_mu_alpha[inl];//[NLOCAL][nm]--trans->[nm][NLOCAL]
        double* c = tmp_v1;
        
        dgemm_(&nt, &t, &nm, &NLOCAL, &nm, &alpha, a, &nm, b, &NLOCAL, &beta, c, &nm);

        //2. <psi_mu|alpha_m>*(dE/dD)*<alpha_m'|psi_nv>
        a = b; //[NLOCAL][nm]
        b = c;//[nm][NLOCAL]
        c = tmp_v2;//[NLOCAL][NLOCAL]
        dgemm_(&nt, &nt, &NLOCAL, &NLOCAL, &nm, &alpha, a, &NLOCAL, b, &nm, &beta, c, &NLOCAL);

        //3. sum of Inl
        for (int i = 0;i < NLOCAL * NLOCAL;++i)
        {
            this->H_V_delta[i] += c[i];
        }
    }
    delete[] tmp_v1;
    delete[] tmp_v2;
    ofs_running << "finish calculating H_V_delta" << endl;
    return;
}

void LCAO_Descriptor::cal_f_delta(matrix& dm)
{
    TITLE("LCAO_Descriptor", "cal_f_delta");
    int iat = 0;    //check if the index same as ucell.iw2iat or not !!
    for (int it = 0;it < ucell.ntype;++it)
    {
        for (int ia = 0;ia < ucell.atoms[it].na;++ia)
        {
            for (int inl = 0;inl < inlmax;++inl)
            {
                int nm = 2 * inl_l[inl] + 1;
                
                //1. cal gdmx
                this->init_gdmx();
                this->cal_gdmx(dm);

                //2.multiply and sum for each atom
                // \sum_{Inl}\sum_{mm'} <gedm, gdmx>_{mm'}
                //notice: sum of multiplied corresponding element(mm') , not matrix multiplication !
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
        }//end ia
    }//end it
    this->del_gdmx();
    return;
}

void LCAO_Descriptor::cal_descriptor_tensor()
{
    TITLE("LCAO_Descriptor", "cal_descriptor_tensor");
    //init pdm_tensor and d_tensor
    torch::Tensor tmp;
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
        d_v = torch::symeig(pdm_tensor[inl], /*eigenvalues=*/true, /*upper=*/true);
        d_tensor[inl] = std::get<0>(d_v);
    }
    return;
}

void LCAO_Descriptor::load_model(const string& model_name)
{
    try {
        module = torch::jit::load(model_name);
    }
    catch (const c10::Error& e) {
        std::cerr << "error loading the model" << std::endl;
        return;
    }
}
void LCAO_Descriptor::cal_gedm()
{
    TITLE("LCAO_Descriptor", "cal_gedm");
    //forward
    std::vector<torch::jit::IValue> inputs;
    //input_dim:(natom, des_per_atom)
    inputs.push_back(torch::cat(d_tensor, /*dim=*/0).reshape({ ucell.nat, des_per_atom }));
    std::vector<torch::Tensor> ec;
    ec.push_back(module.forward(inputs).toTensor());
    this->E_delta = ec[0].item().toDouble();
    
    //cal gedm
    std::vector<torch::Tensor> grad_shell;
    grad_shell.push_back(torch::ones_like(ec[0]));
    this->gedm_tensor = torch::autograd::grad(ec, this->pdm_tensor, grad_shell);

    //gedm_tensor to gedm
    for (int inl = 0;inl < inlmax;++inl)
    {
        int nm = 2 * inl_l[inl] + 1;
        for (int m1 = 0;m1 < nm;++m1)
        {
            for (int m2 = 0;m2 < nm;++m2)
            {
                int index = m1 * nm + m2;
                this->gedm[inl][index] = this->gedm_tensor[inl].index({ m1,m2 }).item().toDouble();
            }
        }
    }
    return;
}

void LCAO_Descriptor::print_H_V_delta()
{
    TITLE("LCAO_Descriptor", "print_H_V_delta");
    ofstream ofs;
    stringstream ss;
    // the parameter 'winput::spillage_outdir' is read from INPUTw.
    ss << winput::spillage_outdir << "/"
       << "H_V_delta.dat";
    if (MY_RANK == 0)
    {
        ofs.open(ss.str().c_str());
    }
    ofs << "E_delta(Hartree) from deepks model: " << this->E_delta << endl;
    ofs << "E_delta(eV) from deepks model: " << this->E_delta * Hartree_to_eV << endl;
    ofs << "H_delta(Hartree)(gamma only)) from deepks model: " << endl;
    for (int i = 0;i < NLOCAL;++i)
    {
        for (int j = 0;j < NLOCAL;++j)
        {
            ofs<< setw(12)<< this->H_V_delta[i * NLOCAL + j] << " ";
        }
        ofs << endl;
    }
    ofs << "H_delta(eV)(gamma only)) from deepks model: " << endl;
    for (int i = 0;i < NLOCAL;++i)
    {
        for (int j = 0;j < NLOCAL;++j)
        {
            ofs<< setw(12)<< this->H_V_delta[i * NLOCAL + j] *Hartree_to_eV<< " ";
        }
        ofs << endl;
    }
    ofs_running << "H_delta is printed" << endl;
    return;
}

void LCAO_Descriptor::print_F_delta()
{
    TITLE("LCAO_Descriptor", "print_F_delta");
    ofstream ofs;
    stringstream ss;
    // the parameter 'winput::spillage_outdir' is read from INPUTw.
    ss << winput::spillage_outdir << "/"
       << "F_delta.dat";
    if (MY_RANK == 0)
    {
        ofs.open(ss.str().c_str());
    }
    ofs << "F_delta(Hatree/Angstrom) from deepks model: " << endl;
    ofs << setw(12) << "type" << setw(12) << "atom" << setw(15) << "dF_x" << setw(15) << "dF_y" << setw(15) << "dF_z" << endl;
    for (int it = 0;it < ucell.ntype;++it)
    {
        for (int ia = 0;ia < ucell.atoms[it].na;++ia)
        {
            int iat = ucell.itia2iat(it, ia);
            ofs << setw(12) << ucell.atoms[it].label << setw(12) << ia
                << setw(15) << this->F_delta(iat, 0) << setw(15) << this->F_delta(iat, 1)
                << setw(15) << this->F_delta(iat, 2) << endl;
        }
    }
    ofs << "F_delta(eV/Angstrom) from deepks model: cd" << endl;
    ofs << setw(12) << "type" << setw(12) << "atom" << setw(15) << "dF_x" << setw(15) << "dF_y" << setw(15) << "dF_z" << endl;
    for (int it = 0;it < ucell.ntype;++it)
    {
        for (int ia = 0;ia < ucell.atoms[it].na;++ia)
        {
            int iat = ucell.itia2iat(it, ia);
            ofs << setw(12) << ucell.atoms[it].label << setw(12)
                << ia << setw(15) << this->F_delta(iat, 0) * Hartree_to_eV
                << setw(15) << this->F_delta(iat, 1) * Hartree_to_eV
                << setw(15) << this->F_delta(iat, 2) * Hartree_to_eV << endl;
        }
    }
    ofs_running << "F_delta is printed" << endl;
    return;
}