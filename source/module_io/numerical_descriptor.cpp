#include "numerical_descriptor.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_cell/module_symmetry/symmetry.h"
#include "winput.h"
#include "module_base/math_ylmreal.h"
#include "module_base/lapack_connector.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"

Numerical_Descriptor::Numerical_Descriptor() 
{
	this->init_label = false;
	this->lmax = -1;
	this->nmax = -1;
}

Numerical_Descriptor::~Numerical_Descriptor() 
{
	if(init_label==true)
	{
		delete[] mu_index;
	}
	return;
}


void Numerical_Descriptor::output_descriptor(const psi::Psi<std::complex<double>> &psi, const int &lmax_in, const double &rcut_in, const double &tol_in, const int nks_in)
{
	ModuleBase::TITLE("Numerical_Descriptor","output_descriptor");
	ModuleBase::GlobalFunc::NEW_PART("DeepKS descriptor: D_{Inl}");

	//-----------------------------------
	// 1. Initialize parameters
	//-----------------------------------

	//GlobalV::ofs_running << "D_{Inl}_m_m'=sum_{i}<J_Inl_m|Psi_i><Psi_i|J_Inl_m'>" << std::endl;
	GlobalV::ofs_running << "input lmax = " << lmax_in << std::endl;
	GlobalV::ofs_running << "input rcut = " << rcut_in << std::endl;
	GlobalV::ofs_running << "input tolerence = " << tol_in << std::endl;
	this->lmax = lmax_in;
	assert(lmax>=0);

    const int nks = nks_in;
    int ne = 0; 
	
	// Peize Lin change 2022.12.15
    // 0 stands for : 'Faln' is not used.
    this->bessel_basis.init(
		0,
		std::stod(INPUT.bessel_descriptor_ecut),
		GlobalC::ucell.ntype,
		this->lmax,
		INPUT.bessel_descriptor_smooth,
		INPUT.bessel_descriptor_sigma,
		rcut_in,
		tol_in,
        GlobalC::ucell
        );
	this->nmax = Numerical_Descriptor::bessel_basis.get_ecut_number();
    this->init_mu_index();
    this->init_label = true;

	assert(nmax>0);

	// Currently we are not considering doing DeePKS in PW basis
	// hence this subroutine is used only for generating projectors and save to jle.orb
	// As a result, I will return here and the rest of the code is saved for future use
	return;

/*
	//-----------------------------------
	// 2. Open the file
	//-----------------------------------
    std::ofstream ofs;
    std::stringstream ss;
    // the parameter 'winput::spillage_outdir' is read from INPUTw.
    ss << winput::spillage_outdir << "/" << "descriptor.dat";
    if (GlobalV::MY_RANK==0)
    {
        ofs.open(ss.str().c_str());
    }


	//-------------------------------------
	// 3. Initialize overlap_Q1 and Q2 
	//-------------------------------------
	// OVERLAP : < J_mu | Psi >
    ModuleBase::realArray overlap_Q1(nks, GlobalV::NBANDS, this->nlocal );
    ModuleBase::realArray overlap_Q2(nks, GlobalV::NBANDS, this->nlocal );

    ModuleBase::GlobalFunc::ZEROS(overlap_Q1.ptr, overlap_Q1.getSize() );
    ModuleBase::GlobalFunc::ZEROS(overlap_Q2.ptr, overlap_Q2.getSize() );

	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"number of k points",overlap_Q1.getBound1());
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"number of bands",overlap_Q1.getBound2());
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"number of local orbitals",overlap_Q1.getBound3());


	//-------------------------------------
	// 4. Compute overlap_Q1 and Q2 
	//-------------------------------------
    // nks now is the reduced k-points.
    for (int ik=0; ik<nks; ik++)
    {
        const int npw= p_kv->ngk[ik];
		GlobalV::ofs_running << " --------------------------------------------------------" << std::endl;
		GlobalV::ofs_running << " Print the overlap matrixs Q and S for this kpoint";
        GlobalV::ofs_running << "\n " << std::setw(8) << "ik" << std::setw(8) << "npw";
        GlobalV::ofs_running << "\n " << std::setw(8) << ik+1 << std::setw(8) << npw << std::endl;
		GlobalV::ofs_running << " --------------------------------------------------------" << std::endl;
        // search for all k-points.
		psi.fix_k(ik);
        this->jlq3d_overlap(overlap_Q1, overlap_Q2, ik, ik, npw, psi);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"jlq3d_overlap");
	}

#ifdef __MPI
    Parallel_Reduce::reduce_double_pool( overlap_Q1.ptr, overlap_Q1.getSize() );
    Parallel_Reduce::reduce_double_pool( overlap_Q2.ptr, overlap_Q2.getSize() );
#endif

	// do not need to output <J|psi> here
    //this->output_overlap_Q( ofs, overlap_Q1, overlap_Q2 );


	
	//-------------------------------------
	// 5. Generate descriptors for each atom 
	//-------------------------------------
	
	for (int it=0; it<GlobalC::ucell.ntype; it++)
	{
		GlobalV::ofs_running << GlobalC::ucell.atoms[it].label << " label" << std::endl;
		for (int ia=0; ia<GlobalC::ucell.atoms[it].na; ia++)
		{
			//--------------------------------------------------
			// compute the number of descriptors for each atom
			// this is a fixed number for all of the atoms
			//--------------------------------------------------
			int l_mu = mu_index[it](ia, 0, 0, 0); //min mu_index for each atom
			int r_mu = mu_index[it](ia, lmax, nmax-1, 2*lmax);//max mu_index for each atom
			const int nd =r_mu - l_mu + 1;
			
			GlobalV::ofs_running << " atom_index: " << ia+1 << " dimension of descritpor: " << nd << std::endl;

			double* d = new double[nd]; //descriptor for each atom
				
			// (it, ia) we know the index 'I' for descriptor
			// each atom has channel up to 'lmax',
			// for each 'lmax' we have 'n' up to 'ecut_number'
			this->generate_descriptor(overlap_Q1, overlap_Q2, it ,ia, d, nd);

			ofs << GlobalC::ucell.atoms[it].label << " atom_index " << ia+1 << " n_descriptor " << nd << std::endl;
			for(int id=0; id<nd; ++id)
			{
				if(id>0 && id%8==0) ofs << std::endl;
			//	if(abs(d[id]>1.0e-9)) ofs << d[id] << " ";
			//	else ofs << "0 ";
				ofs << d[id] << " ";
			}
			ofs << std::endl;

			delete[] d;
		}
	}



    if (GlobalV::MY_RANK==0) ofs.close();
    return;
*/
}

/*
void Numerical_Descriptor::generate_descriptor(ModuleBase::realArray &overlap_Q1, ModuleBase::realArray &overlap_Q2,
const int &it, const int &ia, double *d, const int &nd)
{
    int nbands = overlap_Q1.getBound2();
    // nwfc = nd * ne
    int nwfc = overlap_Q1.getBound3();
    int start0 = mu_index[it](ia, 0, 0, 0); //min mu_index for each atom


    // for 1 k-point only now
    int ik=0;


    GlobalV::ofs_running << " print out each descriptor" << std::endl;
    int id=0;
    for (int l=0; l<lmax+1; l++)
    {
        for (int n=0; n<nmax; n++)
        {
                const int dim = 2*l+1;
                // descriptor for atom (it, ia)
                ModuleBase::ComplexMatrix des(dim, dim);
                for (int m=0; m<2*l+1; m++)
                {
                        const int ii=mu_index[it](ia,l,n,m);
                        for (int m2=0; m2<2*l+1; m2++)
                        {
                                const int jj=mu_index[it](ia,l,n,m2);
                                for(int ib=0; ib<nbands; ib++) // sum for nbands
                                {
                                        std::complex<double> c1(overlap_Q1(ik, ib, ii),  overlap_Q2(ik, ib, ii));
                                        std::complex<double> c2(overlap_Q1(ik, ib, jj), -overlap_Q2(ik, ib, jj));
                                        des(m,m2) += c1*c2;
                                }
            //					GlobalV::ofs_running << std::setw(15) << des(m,m2);
                        }
            //			GlobalV::ofs_running << std::endl;
                }
                GlobalV::ofs_running << "dimension of des is " << 2*l+1 << std::endl;

                if(l==0)
                {
                    d[id]=des(0,0).real();
                    ++id;
                }
                else
                {
                    // diagonalizae
                    // assume des matrix is Hermitian
                    char jobz = 'N'; // eigenvalues only
                    char uplo = 'U'; // upper matrix is stored
                    int ndim = des.nr;
                    double* tmpd = new double[ndim]();
                    const int lwork = 2*ndim;
                    std::complex<double>* work = new std::complex<double>[lwork]();
                    double* rwork = new double[3 * ndim - 2]();
                    int infor = 0;
                    // diag by calling zheev
                    LapackConnector::zheev(jobz, uplo, ndim, des, ndim, tmpd, work, lwork, rwork, &infor);
                    // put the eigenvalues into d (descriptor)
                    for(int idim=0; idim<ndim; ++idim)
                    {
                        d[id] = tmpd[idim];
                        ++id;
                    }
                    delete[] tmpd;
                    delete[] rwork;
                    delete[] work;
                }
        }
    }

    return;
}



void Numerical_Descriptor::jlq3d_overlap(
    ModuleBase::realArray &overlap_Q1,
    ModuleBase::realArray &overlap_Q2,
    const int &ik_ibz,
    const int &ik,
    const int &np,
    const psi::Psi<std::complex<double>> &psi)
{
    ModuleBase::TITLE("Numerical_Descriptor","jlq3d_overlap");
    ModuleBase::timer::tick("Numerical_Descriptor","jlq3d_overlap");

    GlobalV::ofs_running << " OUTPUT THE OVERLAP BETWEEN SPHERICAL BESSEL FUNCTIONS AND BLOCH WAVE FUNCTIONS" <<
std::endl; GlobalV::ofs_running << " Q = < J_it_ia_il_in_im | Psi_n, k > " << std::endl;

    const double normalization = (4 * ModuleBase::PI) / sqrt(GlobalC::ucell.omega);// Peize Lin add normalization
2015-12-29

    const int total_lm = ( this->lmax + 1) * ( this->lmax + 1);
    ModuleBase::matrix ylm(total_lm, np);

    ModuleBase::Vector3<double> *gk = new ModuleBase::Vector3 <double> [np];
    for (int ig=0; ig<np; ig++)
    {
        gk[ig] = wfcpw->getgpluskcar(ik,ig);
    }

    ModuleBase::YlmReal::Ylm_Real(total_lm, np, gk, ylm);

    GlobalV::ofs_running << "\n " << std::setw(5) << "ik"
    << std::setw(8) << "Type1"
    << std::setw(8) << "Atom1"
    << std::setw(8) << "L"
    << std::endl;

    double *flq = new double[np];
    std::complex<double> overlapQ = ModuleBase::ZERO;
    for (int T1 = 0; T1 < GlobalC::ucell.ntype; T1++)
    {
        for (int I1 = 0; I1 < GlobalC::ucell.atoms[T1].na; I1++)
        {
            std::complex<double> *sk = sf.get_sk(ik, T1, I1,wfcpw);
            for (int L=0; L< lmax+1; L++)
            {
                GlobalV::ofs_running << " " << std::setw(5) << ik+1
                            << std::setw(8) << GlobalC::ucell.atoms[T1].label
                            << std::setw(8) << I1+1
                            << std::setw(8) << L
                            << std::endl;
                //OUT("l",l);
                std::complex<double> lphase = normalization * pow(ModuleBase::IMAG_UNIT, L);			// Peize Lin add
normalization 2015-12-29 for (int ie=0; ie < nmax; ie++)
                {
                    for (int ig=0; ig<np; ig++)
                    {
                        flq[ig] = Numerical_Descriptor::bessel_basis.Polynomial_Interpolation2
                                  (L, ie, gk[ig].norm() * GlobalC::ucell.tpiba );
                    }

                    for (int m=0; m<2*L+1; m++)
                    {
                        const int lm = L*L+m;
                        for (int ib=0; ib<GlobalV::NBANDS; ib++)
                        {
                            std::complex<double> overlap_tmp = ModuleBase::ZERO;
                            for (int ig=0; ig<np; ig++)
                            {
                                const std::complex<double> local_tmp = lphase * sk[ig] * ylm(lm, ig) * flq[ig];
                                overlap_tmp += conj( local_tmp ) * psi(ib, ig); // psi is bloch orbitals
                            }
                            overlap_Q1(ik_ibz, ib, mu_index[T1](I1, L, ie, m)) = overlap_tmp.real();
                            overlap_Q2(ik_ibz, ib, mu_index[T1](I1, L, ie, m)) = overlap_tmp.imag();
                        }
                    }
                }//end ie
            }//end l
            delete[] sk;
        }
    }

    delete[] flq;
    delete[] gk;
    ModuleBase::timer::tick("Numerical_Descriptor","jlq3d_overlap");
    return;
}
*/

void Numerical_Descriptor::init_mu_index(void)
{
	GlobalV::ofs_running << " Initialize the mu index for deepks" << std::endl;
	GlobalV::ofs_running << " lmax = " << this->lmax << std::endl;
	GlobalV::ofs_running << " nmax = " << this->nmax << std::endl;
    Numerical_Descriptor::mu_index = new ModuleBase::IntArray[GlobalC::ucell.ntype];

	assert(lmax>=0);
	assert(nmax>0);
	
	int mu=0;
	for (int it=0; it<GlobalC::ucell.ntype; ++it)
	{
		this->mu_index[it].create(
			GlobalC::ucell.atoms[it].na,
			lmax+1, // l starts from 0
			nmax,
			2*lmax+1); // m ==> 2*l+1

		GlobalV::ofs_running << "Type " << it+1 
		<< " number_of_atoms " << GlobalC::ucell.atoms[it].na
		<< " number_of_L " << lmax+1
		<< " number_of_n " << nmax
		<< " number_of_m " << 2*lmax+1 << std::endl;

        for (int ia=0; ia<GlobalC::ucell.atoms[it].na; ia++)
		{
				for (int l=0; l<lmax+1; l++)
				{
						for (int n=0; n<nmax; n++)
						{
								for (int m=0; m<2*l+1; m++)
								{
										this->mu_index[it](ia,l,n,m) = mu;
										mu++;
								}
						}
				}
		}

	}

	this->nlocal = mu;
	GlobalV::ofs_running << " total number of atomic orbitals " << nlocal << std::endl;

	return;
}
