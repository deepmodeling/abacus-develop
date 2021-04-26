//caoyu add 2021-03-29
#include "LCAO_descriptor.h"
#include "LCAO_matrix.h"
#include "../src_global/lapack_connector.h"
#include "../src_global/intarray.h"
#include "../src_global/complexmatrix.h"
#include "global_fp.h"
#include "../src_pw/global.h"
#include "../src_io/winput.h"

LCAO_Descriptor::LCAO_Descriptor()
{
    S_mu_alpha = new double[1];
    PDM = new double[1];
    mu_index = new IntArray[1];
    d = new double[1];
}
LCAO_Descriptor::~LCAO_Descriptor()
{
    delete[] S_mu_alpha;
    delete[] PDM;
    delete[] mu_index;
    delete[] d;
}

void LCAO_Descriptor::build_S_descriptor(const bool &calc_deri)
{
    TITLE("LCAO_Descriptor", "build_S_descriptor");

    // =======init==============
    // cal n(descriptor) per atom , related to Lmax, nchi(L) and m. (not total_nchi!)
	this->des_per_atom=0; // mohan add 2021-04-21
    for (int l = 0; l <= ORB.get_lmax_d(); l++)
    {
        this->des_per_atom += ORB.Alpha[0].getNchi(l) * (2 * l + 1);
    }
    this->n_descriptor = ucell.nat * this->des_per_atom;
    const long Ssize = this->n_descriptor * NLOCAL;
    delete[] S_mu_alpha;
    S_mu_alpha = new double[Ssize];
    ZEROS(S_mu_alpha, Ssize);

    this->init_mu_index();
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
            GridD.Find_atom(tau1, T1, I1);
            int *T2arr = new int[GridD.getAdjacentNum() + 1];
            int *I2arr = new int[GridD.getAdjacentNum() + 1];
            for (int ad = 0; ad < GridD.getAdjacentNum() + 1; ++ad)
            {
                T2arr[ad] = GridD.getType(ad);
                I2arr[ad] = GridD.getNatom(ad);
            }
            for (int ad = 0; ad < GridD.getAdjacentNum() + 1; ++ad)
            {
                //const int T2 = GridD.getType(ad);
                //const int I2 = GridD.getNatom(ad);
                Atom *atom2 = &ucell.atoms[T2arr[ad]];
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

                        //init iw2_all
                        int iw2_all = 0;
                        int iatom = 0;
                        for (int it = 0; it < T2arr[ad]; it++)
                        {
                            for (int ia = 0; ia < ucell.atoms[it].na; ia++)
                            {
                                iatom++; // cal how many atoms before ad in ucell
                            }
                        }
                        iatom += I2arr[ad];
                        iw2_all = iatom * this->des_per_atom;

                        for (int L2 = 0; L2 <= ORB.Alpha[0].getLmax(); ++L2)
                        {
                            for (int N2 = 0; N2 < ORB.Alpha[0].getNchi(L2); ++N2)
                            {
                                for (int m2 = 0; m2 < 2 * L2 + 1; ++m2)
                                {
                                    olm[0] = olm[1] = olm[2] = 0.0;
                                    if (!calc_deri)
                                    {
                                        UOT.snap_psialpha(olm, 0, tau1,
                                                          T1, L1, m1, N1, GridD.getAdjacentTau(ad),
                                                          T2arr[ad], L2, m2, N2);
                                        if (GAMMA_ONLY_LOCAL)
                                        {
                                            this->set_S_mu_alpha(iw1_all, iw2_all, olm[0]);
                                        }
                                    }
                                    ++iw2_all;
                                } //m2
                            }     //N2
                        }         //nw2(L2)
                        ++iw1_all;
                    } // nw1
                }     // distance
            }         // ad
            delete[] T2arr;
            delete[] I2arr;
        } // I1
    }     // T1
    if (!GAMMA_ONLY_LOCAL)
    {
        WARNING_QUIT("LCAO_Descriptor::build_S_descriptor", "muti-kpoint method for descriptor is not implemented yet! ");
    }
    return;
}

void LCAO_Descriptor::set_S_mu_alpha(const int &iw1_all, const int &iw2_all, const double &v)
{
    //const int ir = ParaO.trace_loc_row[iw1_all];
    //const int ic = ParaO.trace_loc_col[iw2_all];
    //no parellel yet
    const int ir = iw1_all;
    const int ic = iw2_all;
    //const int index = ir * ParaO.ncol + ic;
    long index;
    if (KS_SOLVER == "genelpa" || KS_SOLVER == "scalapack_gvx") // save the matrix as column major format
    {
        index = ic * NLOCAL + ir;
    }
    else
    {
        index = ir * this->n_descriptor + ic; //row: lcao orbitals; col: descriptor basis
    }
    this->S_mu_alpha[index] += v;
    return;
}

void LCAO_Descriptor::cal_projected_DM()
{
    //step 1: get dm: the coefficient of wfc, not charge density
    double *dm = new double[NLOCAL * NLOCAL];
    ZEROS(dm, NLOCAL * NLOCAL);
    for (int i = 0; i < LOC.wfc_dm_2d.dm_gamma[0].nr; i++)
    {
        for (int j = 0; j < LOC.wfc_dm_2d.dm_gamma[0].nc; j++)
        {
            dm[i * NLOCAL + j] = LOC.wfc_dm_2d.dm_gamma[0](i, j); //only consider default NSPIN = 1
        }
    }
/*
    //===============test==============
    cout << "test: out wfc_dm_2d.dm_gamma[0](i, j)" << endl;
    for (int ir = 0; ir < NLOCAL; ir++)
    {
        for (int ic = 0; ic < NLOCAL; ic++)
        {
            cout << dm[ir * NLOCAL + ic] << " ";
        }
        cout << endl;
    }
    //===============\test==============
*/
    //step 2: get SS_alpha_mu and SS_nu_beta
    double *ss = this->S_mu_alpha; //SS_nu_beta
/*
    //===============test==============
    cout << "test: out S_nu_beta" << endl;
    for (int ir = 0; ir < NLOCAL; ir++)
    {
        for (int ic = 0; ic < this->n_descriptor; ic++)
        {
            cout << ss[ir * this->n_descriptor + ic] << " ";
        }
        cout << endl;
    }
    //===============\test==============
*/
    //step 3 : multiply
    //cal ssT*DM*ss

    const long tmp_PDM_size = NLOCAL * this->n_descriptor;
    double *tmp_PDM = new double[tmp_PDM_size];
    ZEROS(tmp_PDM, tmp_PDM_size);
    const long PDM_size = this->n_descriptor * this->n_descriptor;
    delete[] this->PDM;
    this->PDM = new double[PDM_size];
    ZEROS(this->PDM, PDM_size);

    const char t = 'T';  //transpose
    const char nt = 'N'; //non transpose
    const double alpha = 1;
    const double beta = 0;
    double *a = dm;
    double *b = ss;
    double *c = tmp_PDM;
    dgemm_(&nt, &nt, &NLOCAL, &n_descriptor, &NLOCAL, &alpha, a, &NLOCAL, b, &NLOCAL, &beta, c, &NLOCAL); //S_nu_nu*SS_nu_beta
    a = ss;
    b = c;
    c = this->PDM;
    dgemm_(&t, &nt, &n_descriptor, &n_descriptor, &NLOCAL, &alpha, a, &NLOCAL, b, &NLOCAL, &beta, c, &n_descriptor); //SS_alpha_mu*S_mu_mu*DM*S_nu_nu*SS_nu_beta
/*
    //===============test==============
    cout << "test: out PDM" << endl;
    for (int ir = 0; ir < n_descriptor; ir++)
    {
        for (int ic = 0; ic < n_descriptor; ic++)
        {
            cout << this->PDM[ir * n_descriptor + ic] << " ";
        }
        cout << endl;
    }
    //===============\test==============
*/
    delete[] tmp_PDM;
    delete[] dm;
    return;
}

void LCAO_Descriptor::cal_descriptor()
{
    delete[] d;
    d = new double[this->n_descriptor];
    //==========print preparation=============
    ofs_running << " print out each DM_Inl" << endl;
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
                    // descriptor for atom (it, ia)
                    ComplexMatrix des(dim, dim);
                    for (int m = 0; m < 2 * l + 1; m++)
                    {
                        const int ii = mu_index[it](ia, l, n, m);
                        for (int m2 = 0; m2 < 2 * l + 1; m2++)
                        {
                            const int jj = mu_index[it](ia, l, n, m2);

                            long index = ii * this->n_descriptor + jj;
                            assert(index >= 0);
                            assert(index < this->n_descriptor * this->n_descriptor);
                            complex<double> tmp(this->PDM[index], 0);
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
    this->print_descriptor();
    return;
}

void LCAO_Descriptor::init_mu_index(void)
{
    ofs_running << " Initialize the mu index for deepks (lcao line)" << endl;
    const int lmax = ORB.get_lmax_d();
    const int nmax = ORB.get_nchimax_d();
    assert(lmax >= 0);
    assert(nmax >= 0);
    ofs_running << " lmax = " << lmax << endl;
    ofs_running << " nmax = " << nmax << endl;

    delete[] this->mu_index;
    this->mu_index = new IntArray[ucell.ntype];

    int mu = 0;
    for (int it = 0; it < ucell.ntype; it++)
    {
        this->mu_index[it].create(
            ucell.atoms[it].na,
            lmax + 1, // l starts from 0
            nmax,
            2 * lmax + 1); // m ==> 2*l+1

        ofs_running << "Type " << it + 1
                    << " number_of_atoms " << ucell.atoms[it].na << endl;

        for (int ia = 0; ia < ucell.atoms[it].na; ia++)
        {
            for (int l = 0; l < lmax + 1; l++)
            {
                for (int n = 0; n < ORB.Alpha[0].getNchi(l); n++)
                {
                    for (int m = 0; m < 2 * l + 1; m++)
                    {
                        this->mu_index[it](ia, l, n, m) = mu;
                        mu++;
                    }
                }
            }
        }
    }
    assert(this->n_descriptor == mu);
    ofs_running << "descriptors_per_atom " << this->des_per_atom << endl;
    ofs_running << "total_descriptors " << this->n_descriptor << endl;

    return;
}

void LCAO_Descriptor::print_projected_DM(ofstream &ofs, ComplexMatrix &des, const int &it, const int &ia, const int &l, const int &n)
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
            int id0 = this->mu_index[it](ia, 0, 0, 0);
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
    return;
}
