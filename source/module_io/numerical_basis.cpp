#include "numerical_basis.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_cell/module_symmetry/symmetry.h"
#include "winput.h"
#include "module_base/math_ylmreal.h"
#include "module_base/parallel_reduce.h"
#include <cstring>
#include <functional>
#include <algorithm>
#include "module_base/timer.h"

Numerical_Basis::Numerical_Basis() {}
Numerical_Basis::~Numerical_Basis() {}

//============================================================
// MEMBER FUNCTION :
// NAME : init
// DESCRIPTION : Two main functions:
// (1) start_from_file = true;
// Firstly, use check(1) to call bessel_basis.init
// to generate TableOne.
// Secondly readin C4 from file.
// Thirdly generate 3D atomic wfc in G space, put the
// results in psi.
//
// (2) If output overlap Q, start_from_file = false;
// Firstly, use check(0) to call bessel_basis,init
// to generate TableOne
// Secondly output overlap, use psi(evc) and jlq3d.
//============================================================
void Numerical_Basis::start_from_file_k(const int& ik, ModuleBase::ComplexMatrix& psi, const Structure_Factor& sf, const ModulePW::PW_Basis_K* wfcpw)
{
    ModuleBase::TITLE("Numerical_Basis","start_from_file_k");

    if (!this->init_label)
    {
        // true stands for : start_from_file
        this->bessel_basis.init(
			true,
			std::stod(INPUT.bessel_nao_ecut),
			GlobalC::ucell.ntype,
			GlobalC::ucell.lmax,
			INPUT.bessel_nao_smooth,
			INPUT.bessel_nao_sigma,
			INPUT.bessel_nao_rcut,
			INPUT.bessel_nao_tolerence,
            GlobalC::ucell
            );
        this->mu_index = this->init_mu_index();
        this->init_label = true;
    }
    this->numerical_atomic_wfc(ik, wfcpw, psi, sf);
}

// The function is called in run_fp.cpp.
void Numerical_Basis::output_overlap(const psi::Psi<std::complex<double>>& psi, const Structure_Factor& sf, const K_Vectors& kv, const ModulePW::PW_Basis_K* wfcpw)
{
    ModuleBase::TITLE("Numerical_Basis","output_overlap");
    ModuleBase::GlobalFunc::NEW_PART("Overlap Data For Spillage Minimization");

	//---------------------------------------------------------
	// if the numerical_basis hasn't been initialized yet,
    // then we initial here.
	//---------------------------------------------------------
    if (!this->init_label)
    {
        // false stands for : 'Faln' is not used.
        this->bessel_basis.init(
			false,
			std::stod(INPUT.bessel_nao_ecut),
			GlobalC::ucell.ntype,
			GlobalC::ucell.lmax,
			INPUT.bessel_nao_smooth,
			INPUT.bessel_nao_sigma,
			INPUT.bessel_nao_rcut,
			INPUT.bessel_nao_tolerence,
            GlobalC::ucell
            );
        this->mu_index = this->init_mu_index();
        this->init_label = true;
    }
    ModuleBase::GlobalFunc::MAKE_DIR(winput::spillage_outdir);
	for(int derivative_order=0; derivative_order<=1; ++derivative_order)            // Peize Lin add 2020.04.23
	{
        std::ofstream ofs;
        std::stringstream ss;
        // the parameter 'winput::spillage_outdir' is read from INPUTw.
        ss << winput::spillage_outdir << "/" <<  "orb_matrix." << derivative_order << ".dat";
        if (GlobalV::MY_RANK==0)
        {
            ofs.open(ss.str().c_str());
        }

        // OVERLAP : < J_mu | Psi >
        std::vector<ModuleBase::ComplexArray> overlap_Q(kv.nks);

        // OVERLAP : < J_mu | J_nu >
        std::vector<ModuleBase::ComplexArray> overlap_Sq(kv.nks);

        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "number of k points", kv.nks);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"number of bands", GlobalV::NBANDS);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"number of local orbitals", GlobalV::NLOCAL);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"number of eigenvalues of Jl(x)", this->bessel_basis.get_ecut_number());

        // nks now is the reduced k-points.
        for (int ik = 0; ik < kv.nks; ik++)
        {
            const int npw = kv.ngk[ik];
            GlobalV::ofs_running << " --------------------------------------------------------" << std::endl;
            GlobalV::ofs_running << " Print the overlap matrixs Q and S for this kpoint" << std::endl;
            GlobalV::ofs_running << std::setw(8) << "ik" << std::setw(8) << "npw" << std::endl;
            GlobalV::ofs_running << std::setw(8) << ik+1 << std::setw(8) << npw << std::endl;
            GlobalV::ofs_running << " --------------------------------------------------------" << std::endl;

            // search for all k-points.
            psi.fix_k(ik);
            overlap_Q[ik]
                = this->cal_overlap_Q(ik, npw, wfcpw, psi, static_cast<double>(derivative_order), sf);
            ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"cal_overlap_Q");

            // (2) generate Sq matrix if necessary.
            if (winput::out_spillage == 2)
            {
                overlap_Sq[ik] = this->cal_overlap_Sq( ik, npw, static_cast<double>(derivative_order), sf, wfcpw);
                ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"cal_overlap_Sq");
            }
        }

        const ModuleBase::matrix overlap_V = this->cal_overlap_V(wfcpw,
                                                                 psi,
                                                                 static_cast<double>(derivative_order),
                                                                 kv); // Peize Lin add 2020.04.23

#ifdef __MPI
        for (int ik = 0; ik < kv.nks; ik++)
        {
            Parallel_Reduce::reduce_pool(overlap_Q[ik].ptr, overlap_Q[ik].getSize());
            Parallel_Reduce::reduce_pool(overlap_Sq[ik].ptr, overlap_Sq[ik].getSize());
        }
        Parallel_Reduce::reduce_pool(overlap_V.c, overlap_V.nr * overlap_V.nc);		// Peize Lin add 2020.04.23
    #endif

        this->output_info(ofs, bessel_basis, kv);

        this->output_k(ofs, kv);

        this->output_overlap_Q(ofs, overlap_Q, kv);

        if (winput::out_spillage == 2)
        {
            this->output_overlap_Sq(ss.str(), ofs, overlap_Sq, kv);
        }

        this->output_overlap_V(ofs, overlap_V); // Peize Lin add 2020.04.23

        if (GlobalV::MY_RANK==0) ofs.close();
    }
    return;
}

ModuleBase::ComplexArray Numerical_Basis::cal_overlap_Q(const int& ik,
                                                        const int& np,
                                                        const ModulePW::PW_Basis_K* wfcpw,
                                                        const psi::Psi<std::complex<double>>& psi,
                                                        const double derivative_order,
                                                        const Structure_Factor& sf) const
{
    ModuleBase::TITLE("Numerical_Basis","cal_overlap_Q");
    ModuleBase::timer::tick("Numerical_Basis","cal_overlap_Q");

	GlobalV::ofs_running << " OUTPUT THE OVERLAP BETWEEN SPHERICAL BESSEL FUNCTIONS AND BLOCH WAVE FUNCTIONS" << std::endl;
	GlobalV::ofs_running << " Q = < J_mu, q | Psi_n, k > " << std::endl;

    ModuleBase::ComplexArray overlap_Q(GlobalV::NBANDS, GlobalV::NLOCAL, this->bessel_basis.get_ecut_number() );
    overlap_Q.zero_out();

	const double normalization = (4 * ModuleBase::PI) / sqrt(GlobalC::ucell.omega);			// Peize Lin add normalization 2015-12-29

    std::vector<ModuleBase::Vector3<double>> gk(np);
    for (int ig = 0; ig < np; ig++)
        gk[ig] = wfcpw->getgpluskcar(ik, ig);

    const std::vector<double> gpow = Numerical_Basis::cal_gpow(gk, derivative_order);

	const ModuleBase::realArray flq = this->cal_flq(ik, gk);

    const ModuleBase::matrix ylm = Numerical_Basis::cal_ylm(gk);

    GlobalV::ofs_running << "\n " << std::setw(5)
        << "ik" << std::setw(8)
        << "Type1" << std::setw(8)
        << "Atom1" << std::setw(8)
        << "L" << std::endl;

    for (int T = 0; T < GlobalC::ucell.ntype; T++)
    {
        //OUT("T",T);
        for (int I = 0; I < GlobalC::ucell.atoms[T].na; I++)
        {
            //OUT("I",I);
            std::complex<double>* sk = sf.get_sk(ik, T, I, wfcpw);
            for (int L=0; L< GlobalC::ucell.atoms[T].nwl+1; L++)
            {
                GlobalV::ofs_running << " " << std::setw(5) << ik+1
                            << std::setw(8) << GlobalC::ucell.atoms[T].label
                            << std::setw(8) << I+1
							<< std::setw(8) << L
							<< std::endl;
                //OUT("l",l);
                std::complex<double> lphase = normalization * pow(ModuleBase::IMAG_UNIT, L);			// Peize Lin add normalization 2015-12-29
                for (int ie=0; ie < this->bessel_basis.get_ecut_number(); ie++)
                {
                    const int N = 0;
                    assert( GlobalC::ucell.nmax == 1);
                    for (int m=0; m<2*L+1; m++)
                    {
                        const int lm = L*L+m;
                        for (int ib=0; ib<GlobalV::NBANDS; ib++)
                        {
                            std::complex<double> overlap_tmp = ModuleBase::ZERO;
                            for (int ig=0; ig<np; ig++)
                            {
//                              const std::complex<double> local_tmp = lphase * sk[ig] * ylm(lm, ig) * flq[ig];
                                const std::complex<double> local_tmp = lphase * sk[ig] * ylm(lm, ig) * flq(L,ie,ig) * gpow[ig];		// Peize Lin add for dpsi 2020.04.23
                                overlap_tmp += conj( local_tmp ) * psi(ib, ig); // psi is bloch orbitals
                            }
                            overlap_Q(ib, this->mu_index[T](I, L, N, m), ie) = overlap_tmp;
                        }
                    }
                }//end ie
            }//end l
            delete[] sk;    sk=nullptr;
        }
    }

    ModuleBase::timer::tick("Numerical_Basis","cal_overlap_Q");
    return overlap_Q;
}

ModuleBase::ComplexArray Numerical_Basis::cal_overlap_Sq(const int& ik,
                                                         const int& np,
                                                         const double derivative_order,
                                                         const Structure_Factor& sf,
                                                         const ModulePW::PW_Basis_K* wfcpw) const
{
    ModuleBase::TITLE("Numerical_Basis","cal_overlap_Sq");
    ModuleBase::timer::tick("Numerical_Basis","cal_overlap_Sq");

	GlobalV::ofs_running << " OUTPUT THE OVERLAP BETWEEN SPHERICAL BESSEL FUNCTIONS"  << std::endl;
	GlobalV::ofs_running << " S = < J_mu,q1 | J_nu,q2 >" << std::endl;

    const int enumber = this->bessel_basis.get_ecut_number();
    ModuleBase::ComplexArray overlap_Sq( GlobalV::NLOCAL, GlobalV::NLOCAL, enumber, enumber );
    overlap_Sq.zero_out();

	const double normalization = (4 * ModuleBase::PI) * (4 * ModuleBase::PI) / GlobalC::ucell.omega;			// Peize Lin add normalization 2015-12-29

    std::vector<ModuleBase::Vector3<double>> gk(np);
    for (int ig = 0; ig < np; ig++)
        gk[ig] = wfcpw->getgpluskcar(ik, ig);

    const std::vector<double> gpow = Numerical_Basis::cal_gpow(gk, derivative_order);

	const ModuleBase::realArray flq = this->cal_flq(ik, gk);

    const ModuleBase::matrix ylm = Numerical_Basis::cal_ylm(gk);

    GlobalV::ofs_running << "\n " << std::setw(5)
        << "ik" << std::setw(8)
        << "Type1" << std::setw(8)
        << "Atom1" << std::setw(8)
        << "L1" << std::setw(8)
        << "Type2" << std::setw(8)
        << "Atom2" << std::setw(8)
        << "L2" << std::endl;

    for (int T1 = 0; T1 < GlobalC::ucell.ntype; T1++) // 1.1
    {
        for (int I1 = 0; I1 < GlobalC::ucell.atoms[T1].na; I1++) // 1.2
        {
            std::complex<double>* sk1 = sf.get_sk(ik, T1, I1, wfcpw);
            for (int T2=0; T2<GlobalC::ucell.ntype; T2++) // 2.1
            {
                for (int I2=0; I2<GlobalC::ucell.atoms[T2].na; I2++) // 2.2
                {
                    std::complex<double>* sk2 = sf.get_sk(ik, T2, I2, wfcpw);
                    for (int l1 = 0; l1 < GlobalC::ucell.atoms[T1].nwl+1; l1++) // 1.3
                    {
                        const std::complex<double> lphase1 = normalization * pow(ModuleBase::IMAG_UNIT, l1);			// Peize Lin add normalization 2015-12-29
                        for (int l2 = 0; l2 < GlobalC::ucell.atoms[T2].nwl+1; l2++) // 2.3
                        {
                            GlobalV::ofs_running << " " << std::setw(5)
                                << ik+1 << std::setw(8)
                                << GlobalC::ucell.atoms[T1].label << std::setw(8)
                                << I1+1 << std::setw(8)
                                << l1 << std::setw(8)
                                << GlobalC::ucell.atoms[T2].label << std::setw(8)
                                << I2+1 << std::setw(8)
                                << l2 << std::setw(8) << std::endl;

                            const std::complex<double> lphase2 = pow(ModuleBase::IMAG_UNIT, l2);
                            for (int ic1=0; ic1 < GlobalC::ucell.nmax; ic1++) // 1.5
                            {
                                for (int ic2=0; ic2 < GlobalC::ucell.nmax; ic2++) // 2.5
                                {
                                    for (int m1=0; m1<2*l1+1; m1++) // 1.6
                                    {
                                        const int lm1 = l1*l1+m1;
                                        const int iwt1 = this->mu_index[T1](I1,l1,ic1,m1);

                                        std::vector<std::complex<double>> about_ig1(np, std::complex<double>(0.0,0.0));
                                        for (int ig=0; ig<np; ig++)
                                            about_ig1[ig] = conj( lphase1 * sk1[ig] * ylm(lm1, ig) ) * gpow[ig];		// Peize Lin add for dpsi 2020.04.23

                                        for (int m2=0; m2<2*l2+1; m2++) // 2.6
                                        {
                                            const int lm2 = l2*l2+m2;
                                            const int iwt2 = this->mu_index[T2](I2,l2,ic2,m2);

                                            std::vector<std::complex<double>> about_ig2(np, std::complex<double>(0.0,0.0));
                                            for (int ig=0; ig<np; ++ig)
                                                about_ig2[ig] = lphase2 * sk2[ig] * ylm(lm2, ig) * about_ig1[ig];

                                            /* same as:
                                            for (int ig=0; ig<np; ig++)
                                                for (int ie1=0; ie1 < enumber; ie1++)
                                                    for (int ie2=0; ie2 < enumber; ie2++)
                                                        overlap_Sq( iwt1, iwt2, ie1, ie2) +=
                                                            about_ig2[ig] * flq(l1,ie1,ig) * flq(l2,ie2,ig);
                                            */

                                            ModuleBase::ComplexMatrix about_ig3_1(enumber,np);
                                            std::copy( &flq(l1,0,0), &flq(l1,0,0)+enumber*np, about_ig3_1.c );

                                            ModuleBase::ComplexMatrix about_ig3_2(enumber,np);
                                            for(int ie2=0; ie2<enumber; ++ie2)
                                                std::transform( &flq(l2,ie2,0), &flq(l2,ie2,0)+np, about_ig2.data(), about_ig3_2.c+ie2*np, std::multiplies<std::complex<double>>() );

                                            BlasConnector::gemm('N', 'T', enumber, enumber, np,
                                                1.0, about_ig3_1.c, np, about_ig3_2.c, np,
                                                1.0, &overlap_Sq(iwt1,iwt2,0,0), enumber);
                                        }
                                    }
                                }
                            }
                        }
                    }
                    delete[] sk2;   sk2=nullptr;
                }
            }
            delete[] sk1;   sk1=nullptr;
        }
    }

    ModuleBase::timer::tick("Numerical_Basis","cal_overlap_Sq");
    return overlap_Sq;
}

// Peize Lin add for dpsi 2020.04.23
ModuleBase::matrix Numerical_Basis::cal_overlap_V(const ModulePW::PW_Basis_K* wfcpw,
                                                  const psi::Psi<std::complex<double>>& psi,
                                                  const double derivative_order,
                                                  const K_Vectors& kv)
{
    ModuleBase::matrix overlap_V(kv.nks, GlobalV::NBANDS);
    for (int ik = 0; ik < kv.nks; ++ik)
    {
        std::vector<ModuleBase::Vector3<double>> gk(kv.ngk[ik]);
        for (int ig=0; ig<gk.size(); ig++)
            gk[ig] = wfcpw->getgpluskcar(ik,ig);

        const std::vector<double> gpow = Numerical_Basis::cal_gpow(gk, derivative_order);

		for(int ib=0; ib<GlobalV::NBANDS; ++ib)
            for (int ig = 0; ig < kv.ngk[ik]; ++ig)
                overlap_V(ik,ib)+= norm(psi(ik,ib,ig)) * gpow[ig];
	}
	return overlap_V;
}

ModuleBase::realArray Numerical_Basis::cal_flq(const int ik, const std::vector<ModuleBase::Vector3<double>> &gk) const
{
	const int np = gk.size();
	const int enumber = this->bessel_basis.get_ecut_number();

    // get flq(G) = \int f(r)jl(G*r) from interpolation table.
    ModuleBase::realArray flq(GlobalC::ucell.lmax+1, enumber, np);
    for (int il=0; il<GlobalC::ucell.lmax+1; il++)
        for (int ie=0; ie<enumber; ie++)
            for (int ig=0; ig<np; ig++)
                flq(il,ie,ig) = this->bessel_basis.Polynomial_Interpolation2(il, ie, gk[ig].norm() * GlobalC::ucell.tpiba );
	return flq;
}

ModuleBase::matrix Numerical_Basis::cal_ylm(const std::vector<ModuleBase::Vector3<double>> &gk)
{
    const int total_lm = ( GlobalC::ucell.lmax + 1) * ( GlobalC::ucell.lmax + 1);
    ModuleBase::matrix ylm(total_lm, gk.size());
    ModuleBase::YlmReal::Ylm_Real(total_lm, gk.size(), gk.data(), ylm);
    return ylm;
}

std::vector<double> Numerical_Basis::cal_gpow (const std::vector<ModuleBase::Vector3<double>> &gk, const double derivative_order)
{
    constexpr double thr = 1E-12;
    std::vector<double> gpow(gk.size(), 0.0);
    for (int ig=0; ig<gpow.size(); ++ig)
    {
        if (derivative_order>=0)
        {
            gpow[ig] = std::pow(gk[ig].norm2(),derivative_order);
        }
        else
        {
            if (gk[ig].norm2() >= thr)
                gpow[ig] = std::pow(gk[ig].norm2(),derivative_order);
        }
    }
    return gpow;
}

std::vector<ModuleBase::IntArray> Numerical_Basis::init_mu_index(void)
{
	GlobalV::ofs_running << " Initialize the mu index" << std::endl;
    std::vector<ModuleBase::IntArray> mu_index_(GlobalC::ucell.ntype);

    int mu = 0;
    for (int it=0; it<GlobalC::ucell.ntype; it++)
    {
        mu_index_[it].create(
            GlobalC::ucell.atoms[it].na,
            GlobalC::ucell.atoms[it].nwl+1,
            GlobalC::ucell.nmax,
            2*(GlobalC::ucell.atoms[it].nwl+1)+1); // m ==> 2*l+1

		// mohan added 2021-01-03
		GlobalV::ofs_running << "Type " << it+1
            << " number_of_atoms " << GlobalC::ucell.atoms[it].na
            << " number_of_L " << GlobalC::ucell.atoms[it].nwl+1
            << " number_of_n " << GlobalC::ucell.nmax
            << " number_of_m " << 2*(GlobalC::ucell.atoms[it].nwl+1)+1 << std::endl;

        for (int ia=0; ia<GlobalC::ucell.atoms[it].na; ia++)
        {
            for (int l=0; l< GlobalC::ucell.atoms[it].nwl+1; l++)
            {
                for (int n=0; n< GlobalC::ucell.atoms[it].l_nchi[l]; n++)
                {
                    for (int m=0; m<2*l+1; m++)
                    {
                        mu_index_[it](ia,l,n,m) = mu;
                        mu++;
                    }
                }
            }
        }
    }
    return mu_index_;
}

void Numerical_Basis::numerical_atomic_wfc(const int& ik,
                                           const ModulePW::PW_Basis_K* wfcpw,
                                           ModuleBase::ComplexMatrix& psi,
                                           const Structure_Factor& sf)
{
    ModuleBase::TITLE("Numerical_Basis", "numerical_atomic_wfc");
    const int np = wfcpw->npwk[ik];
    std::vector<ModuleBase::Vector3<double>> gk(np);
    for (int ig = 0; ig < np; ig++)
        gk[ig] = wfcpw->getgpluskcar(ik, ig);

    const int total_lm = ( GlobalC::ucell.lmax + 1) * ( GlobalC::ucell.lmax + 1);
    ModuleBase::matrix ylm(total_lm, np);
    ModuleBase::YlmReal::Ylm_Real(total_lm, np, gk.data(), ylm);

    std::vector<double> flq(np);
    for (int it = 0; it < GlobalC::ucell.ntype; it++)
    {
        //OUT("it",it);
        for (int ia = 0; ia < GlobalC::ucell.atoms[it].na; ia++)
        {
            //OUT("ia",ia);
            std::complex<double>* sk = sf.get_sk(ik, it, ia, wfcpw);
            for (int l = 0; l < GlobalC::ucell.atoms[it].nwl+1; l++)
            {
                //OUT("l",l);
                std::complex<double> lphase = pow(ModuleBase::IMAG_UNIT, l);
                for (int ic=0; ic < GlobalC::ucell.atoms[it].l_nchi[l]; ic++)
                {
                    //OUT("ic",ic);
                    for (int ig=0; ig<np; ig++)
                    {
                        flq[ig] = this->bessel_basis.
                                  Polynomial_Interpolation(it, l, ic, gk[ig].norm() * GlobalC::ucell.tpiba );
                    }

                    for (int m=0; m<2*l+1; m++)
                    {
                        //OUT("m",m);
                        const int lm = l*l+m;
                        for (int ig=0; ig<np; ig++)
                        {
                            psi( this->mu_index[it](ia,l,ic,m), ig) =
                                lphase * sk[ig] * ylm(lm, ig) * flq[ig];
                        }
                    }
                }
            }
            delete[] sk;    sk=nullptr;
        }
    }
}

void Numerical_Basis::output_info(std::ofstream& ofs, const Bessel_Basis& bessel_basis, const K_Vectors& kv)
{
    // only print out to the information by the first processor
    if (GlobalV::MY_RANK==0)
    {
        ofs.precision(10);
        ofs << GlobalC::ucell.lat0 << std::endl;

        ofs << GlobalC::ucell.latvec.e11 << " " << GlobalC::ucell.latvec.e12 << " " << GlobalC::ucell.latvec.e13 << std::endl;
        ofs << GlobalC::ucell.latvec.e21 << " " << GlobalC::ucell.latvec.e22 << " " << GlobalC::ucell.latvec.e23 << std::endl;
        ofs << GlobalC::ucell.latvec.e31 << " " << GlobalC::ucell.latvec.e32 << " " << GlobalC::ucell.latvec.e33 << std::endl;

        ofs << GlobalC::ucell.ntype << " ntype" << std::endl;
        for (int it=0; it<GlobalC::ucell.ntype; it++)
        {
            ofs << GlobalC::ucell.atoms[it].label << " label" << std::endl; // mohan add 2009-07-23
            ofs << GlobalC::ucell.atoms[it].na << " na" << std::endl;
            for (int ia=0; ia<GlobalC::ucell.atoms[it].na; ia++)
            {
                ofs << GlobalC::ucell.atoms[it].tau[ia].x << " "
                    << GlobalC::ucell.atoms[it].tau[ia].y << " "
                    << GlobalC::ucell.atoms[it].tau[ia].z << std::endl;
            }
        }
        // ecutwfc_jlq determine the jlq corresponding to plane wave calculation.
        ofs << INPUT.ecutwfc << " ecutwfc" << std::endl; // mohan add 2009-09-08

        // this parameter determine the total number of jlq.
        ofs << bessel_basis.get_ecut() << " ecutwfc_jlq" << std::endl;//mohan modify 2009-09-08
        ofs << bessel_basis.get_rcut() << " rcut_Jlq" << std::endl;

        // mohan add 'smooth' and 'sigma' 2009-08-28
        ofs << bessel_basis.get_smooth() << " smooth" << std::endl;
        ofs << bessel_basis.get_sigma() << " sigma" << std::endl;

        ofs << bessel_basis.get_tolerence() << " tolerence" << std::endl;

        ofs << GlobalC::ucell.lmax << " lmax" << std::endl;
    }

    ofs << std::scientific;

    ofs << std::setprecision(8);
    // NOTICE: ofs_warning << "\n The precison may affect the optimize result.";

    if (GlobalV::MY_RANK==0)
    {
        ofs << kv.nkstot << " nks" << std::endl;
        ofs << GlobalV::NBANDS << " nbands" << std::endl;
        ofs << GlobalV::NLOCAL << " nwfc" << std::endl;
        ofs << bessel_basis.get_ecut_number() << " ne " << std::endl;
    }
}

void Numerical_Basis::output_k(std::ofstream& ofs, const K_Vectors& kv)
{
    // (1)
    if (GlobalV::MY_RANK==0)
    {
        ofs << "<WEIGHT_OF_KPOINTS>";
    }

    // only half of nkstot should be output in "NSPIN == 2" case, k_up and k_down has same k infomation
    int nkstot = kv.nkstot;

    // (2)
    for (int ik=0; ik<nkstot; ik++)
    {
        double kx, ky, kz, wknow;
#ifdef __MPI
        //temprary restrict kpar=1 for NSPIN=2 case for generating_orbitals
        int pool = 0;
        if(GlobalV::NSPIN !=2 ) pool = GlobalC::Pkpoints.whichpool[ik];
        const int iknow = ik - GlobalC::Pkpoints.startk_pool[GlobalV::MY_POOL];
        if (GlobalV::RANK_IN_POOL==0)
        {
            if (GlobalV::MY_POOL==0)
            {
                if (pool==0)
                {
                    kx = kv.kvec_c[ik].x;
                    ky = kv.kvec_c[ik].y;
                    kz = kv.kvec_c[ik].z;
                    wknow = kv.wk[ik];
                }
                else
                {
                    MPI_Status ierror;
                    MPI_Recv(&kx, 1, MPI_DOUBLE, GlobalC::Pkpoints.startpro_pool[pool], ik*4, MPI_COMM_WORLD,&ierror);
                    MPI_Recv(&ky, 1, MPI_DOUBLE, GlobalC::Pkpoints.startpro_pool[pool], ik*4+1, MPI_COMM_WORLD,&ierror);
                    MPI_Recv(&kz, 1, MPI_DOUBLE, GlobalC::Pkpoints.startpro_pool[pool], ik*4+2, MPI_COMM_WORLD,&ierror);
                    MPI_Recv(&wknow, 1, MPI_DOUBLE, GlobalC::Pkpoints.startpro_pool[pool], ik*4+3, MPI_COMM_WORLD,&ierror);
                }
            }
            else
            {
                if (GlobalV::MY_POOL == pool)
                {
                    MPI_Send(&kv.kvec_c[iknow].x, 1, MPI_DOUBLE, 0, ik * 4, MPI_COMM_WORLD);
                    MPI_Send(&kv.kvec_c[iknow].y, 1, MPI_DOUBLE, 0, ik * 4 + 1, MPI_COMM_WORLD);
                    MPI_Send(&kv.kvec_c[iknow].z, 1, MPI_DOUBLE, 0, ik * 4 + 2, MPI_COMM_WORLD);
                    MPI_Send(&kv.wk[iknow], 1, MPI_DOUBLE, 0, ik * 4 + 3, MPI_COMM_WORLD);
                }
            }
        }
        // this barrier is very important
        MPI_Barrier(MPI_COMM_WORLD);
#else
        if (GlobalV::MY_RANK==0)
        {
            kx = kv.kvec_c[ik].x;
            ky = kv.kvec_c[ik].y;
            kz = kv.kvec_c[ik].z;
            wknow = kv.wk[ik];
        }
#endif

        if (GlobalV::MY_RANK==0)
        {
            ofs << "\n" << kx << " " << ky << " " << kz;
            ofs << " " << wknow * 0.5;
        }
    }

    if (GlobalV::MY_RANK==0)
    {
        ofs << "\n</WEIGHT_OF_KPOINTS>" << std::endl;
    }
}

void Numerical_Basis::output_overlap_Q(std::ofstream& ofs,
                                       const std::vector<ModuleBase::ComplexArray>& overlap_Q,
                                       const K_Vectors& kv)
{
    // (3)
    if (GlobalV::MY_RANK==0)
    {
        ofs << "\n<OVERLAP_Q>";
    }

    // (4)
    /*
    if(GlobalV::MY_RANK==0)
    {
    //    	for( int i=0; i<overlap_Q1.getSize(); i++)
    //    	{
    //    		if( i%2==0 ) ofs << "\n";
    //    		ofs << " " << overlap_Q1.ptr[i] << " " << overlap_Q2.ptr[i];
    //    	}
    }
    */

    // Copy to overlap_Q_k for Pkpoints.pool_collection temporaly.
    // It's better to refactor to Pkpoints.pool_collection(overlap_Q) in the future.
    // Peize Lin comments 2021.07.25
    assert(kv.nks > 0);
    ModuleBase::ComplexArray overlap_Q_k(kv.nks,
                                         overlap_Q[0].getBound1(),
                                         overlap_Q[0].getBound2(),
                                         overlap_Q[0].getBound3());
    for (int ik = 0; ik < kv.nks; ++ik)
    {
        std::memcpy(
            overlap_Q_k.ptr + ik*overlap_Q[ik].getSize(),
            overlap_Q[ik].ptr,
            overlap_Q[ik].getSize() * sizeof(std::complex<double>));
    }

    // only half of nkstot should be output in "NSPIN == 2" case, k_up and k_down has same k infomation
    int nkstot = kv.nkstot;
    int count = 0;
    for (int ik=0; ik<nkstot; ik++)
    {
        ModuleBase::ComplexArray Qtmp(overlap_Q[ik].getBound1(), overlap_Q[ik].getBound2(), overlap_Q[ik].getBound3());
        Qtmp.zero_out();
        GlobalC::Pkpoints.pool_collection(Qtmp.ptr, overlap_Q_k, ik);
        if (GlobalV::MY_RANK==0)
        {
    //        ofs << "\n ik=" << ik;
            // begin data writing.
            const int dim = Qtmp.getSize();
            for (int i=0; i<dim; i++)
            {
                if ( count%4==0 ) ofs << "\n";
                ofs << " " << Qtmp.ptr[i].real() << " " << Qtmp.ptr[i].imag();
                ++count;
            }
            // end data writing.
        }
#ifdef __MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif
    }

    // (5)
    if (GlobalV::MY_RANK==0)
    {
        ofs << "\n</OVERLAP_Q>" << std::endl;
    }
}

void Numerical_Basis::output_overlap_Sq(const std::string& name,
                                        std::ofstream& ofs,
                                        const std::vector<ModuleBase::ComplexArray>& overlap_Sq,
                                        const K_Vectors& kv)
{
    if (GlobalV::MY_RANK==0)
    {
        ofs << "\n<OVERLAP_Sq>";
        ofs.close();
    }

    // only half of nkstot should be output in "NSPIN == 2" case, k_up and k_down has same k infomation
    int ispin = 1;
    if(GlobalV::NSPIN == 2) ispin = 2;
    int nkstot = kv.nkstot / ispin;
    int count = 0;
    for(int is = 0; is<ispin; is++)
    {
        for (int ik=0; ik<nkstot; ik++)
        {
            if ( GlobalV::MY_POOL == GlobalC::Pkpoints.whichpool[ik] )
            {
                if ( GlobalV::RANK_IN_POOL == 0)
                {
                    ofs.open(name.c_str(), std::ios::app);
                    const int ik_now = ik - GlobalC::Pkpoints.startk_pool[GlobalV::MY_POOL] + is * nkstot;

                    const int size = overlap_Sq[ik_now].getSize();
                    for (int i=0; i<size; i++)
                    {
                        if (count%2==0) ofs << "\n";
                        ofs << " " << overlap_Sq[ik_now].ptr[i].real() << " " << overlap_Sq[ik_now].ptr[i].imag();
                        ++count;
                    }

                    ofs.close();
                }
    #ifdef __MPI
                MPI_Barrier(MPI_COMM_WORLD);
    #endif
            }
            else
            {
    #ifdef __MPI
                MPI_Barrier(MPI_COMM_WORLD);
    #endif
            }

            /*
            if(MY_RANK==0)
            for(int i=0; i< Sq_real[ik].getSize(); i++)
            {
                if(i%2==0) ofs << "\n";
                ofs << " " << Sq_real[ik].ptr[i] << " " << Sq_imag[ik].ptr[i];
            }
            */

        }
    }
    if (GlobalV::MY_RANK==0)
    {
        ofs.open(name.c_str(), std::ios::app);
        ofs << "\n</OVERLAP_Sq>" << std::endl;
    }
}

// Peize Lin add 2020.04.23
void Numerical_Basis::output_overlap_V(
    std::ofstream &ofs,
	const ModuleBase::matrix &overlap_V)
{
	if (GlobalV::MY_RANK==0)
    {
        ofs << "\n<OVERLAP_V>" <<std::endl;;
		overlap_V.print(ofs);
		ofs << "</OVERLAP_V>" <<std::endl;
	}
}
