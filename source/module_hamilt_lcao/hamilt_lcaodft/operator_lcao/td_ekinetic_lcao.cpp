#include "td_ekinetic_lcao.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_hamilt_lcao/hamilt_lcaodft/center2_orb-orb11.h"
#include "module_elecstate/potentials/H_TDDFT_pw.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer_funcs.h"

#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_base/libm/libm.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_matrix.h"

namespace hamilt
{
template <typename TK, typename TR>
TDEkinetic<OperatorLCAO<TK, TR>>::TDEkinetic(LCAO_Matrix* LM_in,
                                                   hamilt::HContainer<TR>* hR_in,
                                                   std::vector<TK>* hK_in,
                                                   hamilt::HContainer<TR>* SR_in,
                                                   const K_Vectors* kv_in,
                                                   const UnitCell* ucell_in,
                                                   Grid_Driver* GridD_in)
    : SR(SR_in), kv(kv_in), OperatorLCAO<TK, TR>(LM_in, kv_in->kvec_d, hR_in, hK_in)
    {
        this->LM = LM_in;
        this->ucell = ucell_in;
        this->cal_type = lcao_tddft_velocity;
        this->Grid = GridD_in;
        this->init_td();
        // initialize HR to get adjs info.
        this->initialize_HR(Grid,this->LM->ParaV);
    }
template <typename TK, typename TR>
TDEkinetic<OperatorLCAO<TK, TR>>::~TDEkinetic()
{
    if (this->hR_tmp != nullptr)
    {
        delete this->hR_tmp;
    }
}
//term A^2*S
template <typename TK, typename TR>
void TDEkinetic<OperatorLCAO<TK, TR>>::td_ekinetic_scalar(std::complex<double>* Hloc, TR* Sloc, int nnr){
	return;
}
//term A^2*S
template <>
void TDEkinetic<OperatorLCAO<std::complex<double>, double>>::td_ekinetic_scalar(std::complex<double>* Hloc, double* Sloc, int nnr){
	std::complex<double> tmp = {cart_At.norm2()*Sloc[nnr]/4.0, 0};
	Hloc[nnr] += tmp;
	return;
}
//term A dot ∇
template <typename TK, typename TR>
void TDEkinetic<OperatorLCAO<TK, TR>>::td_ekinetic_grad(std::complex<double>* Hloc, int nnr, ModuleBase::Vector3<double> grad_overlap){
	std::complex<double> tmp= {0, grad_overlap*cart_At};
	Hloc[nnr] += tmp;
	return;
}

template <typename TK, typename TR>
void TDEkinetic<OperatorLCAO<TK, TR>>::calculate_HR()
{
    ModuleBase::TITLE("TDEkinetic", "calculate_HR");
    if(this->hR_tmp==nullptr || this->hR_tmp->size_atom_pairs()<=0)
    {
        ModuleBase::WARNING_QUIT("TDEkinetic::calculate_HR", "hR_tmp is nullptr or empty");
    }
    ModuleBase::timer::tick("TDEkinetic", "calculate_HR");

    const Parallel_Orbitals* paraV = this->hR_tmp->get_atom_pair(0).get_paraV();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int iat1 = 0; iat1 < this->ucell->nat; iat1++)
    {
        auto tau1 = ucell->get_tau(iat1);
        int T1, I1;
        ucell->iat2iait(iat1, &I1, &T1);
        AdjacentAtomInfo& adjs = this->adjs_all[iat1];
        for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
        {
            const int T2 = adjs.ntype[ad];
            const int I2 = adjs.natom[ad];
            const int iat2 = ucell->itia2iat(T2, I2);
            const ModuleBase::Vector3<int>& R_index2 = adjs.box[ad];
            ModuleBase::Vector3<double> dtau = this->ucell->cal_dtau(iat1, iat2, R_index2);

            hamilt::BaseMatrix<std::complex<double>>* tmp = this->hR_tmp->find_matrix(iat1, iat2, R_index2.x, R_index2.y, R_index2.z);
            hamilt::BaseMatrix<TR>* tmp1 = this->SR->find_matrix(iat1, iat2, R_index2.x, R_index2.y, R_index2.z);
            if (tmp != nullptr)
            {
                this->cal_HR_IJR(iat1, iat2, paraV, dtau, tmp->get_pointer(),tmp1->get_pointer());
            }
            else
            {
                ModuleBase::WARNING_QUIT("TDEkinetic::calculate_HR", "R_index not found in HR");
            }
        }
    }

    ModuleBase::timer::tick("TDEkinetic", "calculate_HR");
}

template <typename TK, typename TR>
void TDEkinetic<OperatorLCAO<TK, TR>>::cal_HR_IJR(const int& iat1,
                                                              const int& iat2,
                                                              const Parallel_Orbitals* paraV,
                                                              const ModuleBase::Vector3<double>& dtau,
                                                              std::complex<double>* data_pointer,
                                                              TR* s_pointer)
{
    const ORB_gen_tables& uot = ORB_gen_tables::get_const_instance();
    const LCAO_Orbitals& orb = LCAO_Orbitals::get_const_instance();
    // ---------------------------------------------
    // get info of orbitals of atom1 and atom2 from ucell
    // ---------------------------------------------
    int T1, I1;
    this->ucell->iat2iait(iat1, &I1, &T1);
    int T2, I2;
    this->ucell->iat2iait(iat2, &I2, &T2);
    Atom& atom1 = this->ucell->atoms[T1];
    Atom& atom2 = this->ucell->atoms[T2];

    // npol is the number of polarizations,
    // 1 for non-magnetic (one Hamiltonian matrix only has spin-up or spin-down),
    // 2 for magnetic (one Hamiltonian matrix has both spin-up and spin-down)
    const int npol = this->ucell->get_npol();

    const int* iw2l1 = atom1.iw2l;
    const int* iw2n1 = atom1.iw2n;
    const int* iw2m1 = atom1.iw2m;
    const int* iw2l2 = atom2.iw2l;
    const int* iw2n2 = atom2.iw2n;
    const int* iw2m2 = atom2.iw2m;
    // ---------------------------------------------
    // get tau1 (in cell <0,0,0>) and tau2 (in cell R)
    // in principle, only dtau is needed in this function
    // snap_psipsi should be refactored to use dtau directly
    // ---------------------------------------------
    const ModuleBase::Vector3<double>& tau1 = this->ucell->get_tau(iat1);
    const ModuleBase::Vector3<double> tau2 = tau1 + dtau;
    // ---------------------------------------------
    // calculate the Ekinetic matrix for each pair of orbitals
    // ---------------------------------------------
    double olm = 0;
    auto row_indexes = paraV->get_indexes_row(iat1);
    auto col_indexes = paraV->get_indexes_col(iat2);
    const int step_trace = col_indexes.size() + 1;
    for (int iw1l = 0; iw1l < row_indexes.size(); iw1l += npol)
    {
        const int iw1 = row_indexes[iw1l] / npol;
        const int L1 = iw2l1[iw1];
        const int N1 = iw2n1[iw1];
        const int m1 = iw2m1[iw1];
#ifdef USE_NEW_TWO_CENTER
        int M1 = (m1 % 2 == 0) ? -m1/2 : (m1+1)/2;
#endif
        for (int iw2l = 0; iw2l < col_indexes.size(); iw2l += npol)
        {
            const int iw2 = col_indexes[iw2l] / npol;
            const int L2 = iw2l2[iw2];
            const int N2 = iw2n2[iw2];
            const int m2 = iw2m2[iw2];
#ifdef USE_NEW_TWO_CENTER
            // center2_orb11_s are used to calculate <psi|∇|psi> no matter whether to use new two-center or not for now.
            ModuleBase::Vector3<double> grad_overlap = center2_orb11_s.at(T1).at(T2).at(L1).at(N1).at(L2).at(N2).cal_grad_overlap(tau1 * ucell->lat0, tau2 * ucell->lat0, m1, m2);
#else
            ModuleBase::Vector3<double> grad_overlap = center2_orb11_s.at(T1).at(T2).at(L1).at(N1).at(L2).at(N2).cal_grad_overlap(tau1 * ucell->lat0, tau2 * ucell->lat0, m1, m2);
#endif
            for (int ipol = 0; ipol < npol; ipol++)
            {
                //key change
                td_ekinetic_scalar(data_pointer, s_pointer, ipol * step_trace);
                td_ekinetic_grad(data_pointer, ipol * step_trace, grad_overlap);
            }
            data_pointer += npol;
            s_pointer += npol;
        }
        data_pointer += (npol - 1) * col_indexes.size();
        s_pointer += (npol - 1) * col_indexes.size();
    }
}
//init two center integrals and vector potential for td_ekintic term
template <typename TK, typename TR>
void TDEkinetic<OperatorLCAO<TK, TR>>::init_td(void)
{
    //calculate At in cartesian coorinates.
	double l_norm[3]={GlobalC::ucell.a1.norm() ,GlobalC::ucell.a2.norm() ,GlobalC::ucell.a3.norm()};
    double (&A)[3] = elecstate::H_TDDFT_pw::At;
	cart_At = GlobalC::ucell.a1*A[0]/l_norm[0] + GlobalC::ucell.a2*A[1]/l_norm[1] + GlobalC::ucell.a3*A[2]/l_norm[2];
    std::cout << "cart_At: " << cart_At[0] << " " <<cart_At[1]<< " " << cart_At[2] << std::endl;
	
    //init MOT,MGT
    this->MOT.allocate(
		GlobalC::ORB.get_ntype(),	// number of atom types
		GlobalC::ORB.get_lmax(),	// max L used to calculate overlap
		static_cast<int>(GlobalC::ORB.get_kmesh()) | 1,				// kpoints, for integration in k space
		GlobalC::ORB.get_Rmax(),				// max value of radial table
		GlobalC::ORB.get_dR(),								// delta R, for making radial table
		GlobalC::ORB.get_dk());											// Peize Lin change 2017-04-16
	int Lmax_used, Lmax;
	this->MOT.init_Table_Spherical_Bessel (2, 1, Lmax_used, Lmax, 1, GlobalC::ORB, GlobalC::ucell.infoNL.Beta);

	//=========================================
	// (2) init Ylm Coef
	//=========================================
	ModuleBase::Ylm::set_coefficients ();

	//=========================================
	// (3) make Gaunt coefficients table
	//=========================================
	this->MGT.init_Gaunt_CH( Lmax );
	this->MGT.init_Gaunt( Lmax );

    //init_radial table
    for( size_t TA=0; TA!=GlobalC::ORB.get_ntype(); ++TA )
		for( size_t TB=0; TB!=GlobalC::ORB.get_ntype(); ++TB )
			for( int LA=0; LA<=GlobalC::ORB.Phi[TA].getLmax(); ++LA )
				for( size_t NA=0; NA!=GlobalC::ORB.Phi[TA].getNchi(LA); ++NA )
					for( int LB=0; LB<=GlobalC::ORB.Phi[TB].getLmax(); ++LB )
						for( size_t NB=0; NB!=GlobalC::ORB.Phi[TB].getNchi(LB); ++NB )
							center2_orb11_s[TA][TB][LA][NA][LB].insert(
								std::make_pair(NB, 
                                Center2_Orb::Orb11(
									GlobalC::ORB.Phi[TA].PhiLN(LA,NA),
									GlobalC::ORB.Phi[TB].PhiLN(LB,NB),
									this->MOT, this->MGT)));
	for( auto &coA : center2_orb11_s )
		for( auto &coB : coA.second )
			for( auto &coC : coB.second )
				for( auto &coD : coC.second )
					for( auto &coE : coD.second )
						for( auto &coF : coE.second )
							coF.second.init_radial_table();
}

template <typename TK, typename TR>
void hamilt::TDEkinetic<hamilt::OperatorLCAO<TK, TR>>::set_HR_fixed(void* hR_tmp_in)
{
    this->hR_tmp = static_cast<hamilt::HContainer<std::complex<double>>*>(hR_tmp_in);
    this->allocated = false;
}
template <typename TK, typename TR>
void TDEkinetic<OperatorLCAO<TK, TR>>::initialize_HR(Grid_Driver* GridD,
                                                        const Parallel_Orbitals* paraV)
{
    if (elecstate::H_TDDFT_pw::stype != 1)
    {
        return;
    }
    ModuleBase::TITLE("TDEkinetic", "initialize_HR");
    ModuleBase::timer::tick("TDEkinetic", "initialize_HR");
    
    this->adjs_all.clear();
    this->adjs_all.reserve(this->ucell->nat);
    for (int iat1 = 0; iat1 < ucell->nat; iat1++)
    {
        auto tau1 = ucell->get_tau(iat1);
        int T1, I1;
        ucell->iat2iait(iat1, &I1, &T1);
        AdjacentAtomInfo adjs;
        GridD->Find_atom(*ucell, tau1, T1, I1, &adjs);
        std::vector<bool> is_adj(adjs.adj_num + 1, false);
        for (int ad1 = 0; ad1 < adjs.adj_num + 1; ++ad1)
        {
            const int T2 = adjs.ntype[ad1];
            const int I2 = adjs.natom[ad1];
            const int iat2 = ucell->itia2iat(T2, I2);
            if (paraV->get_row_size(iat1) <= 0 || paraV->get_col_size(iat2) <= 0)
            {
                continue;
            }
            const ModuleBase::Vector3<int>& R_index2 = adjs.box[ad1];
            // choose the real adjacent atoms
            const LCAO_Orbitals& orb = LCAO_Orbitals::get_const_instance();
            // Note: the distance of atoms should less than the cutoff radius, 
            // When equal, the theoretical value of matrix element is zero, 
            // but the calculated value is not zero due to the numerical error, which would lead to result changes.
            if (this->ucell->cal_dtau(iat1, iat2, R_index2).norm() * this->ucell->lat0
                < orb.Phi[T1].getRcut() + orb.Phi[T2].getRcut())
            {
                is_adj[ad1] = true;
            }
        }
        filter_adjs(is_adj, adjs);
        this->adjs_all.push_back(adjs);
    }
    ModuleBase::timer::tick("TDEkinetic", "initialize_HR");
}
template <typename TK, typename TR>
void TDEkinetic<OperatorLCAO<TK, TR>>::initialize_HR_tmp(const Parallel_Orbitals* paraV)
{
    if (elecstate::H_TDDFT_pw::stype != 1)
    {
        return;
    }
    ModuleBase::TITLE("TDEkinetic", "initialize_HR_tmp");
    ModuleBase::timer::tick("TDEkinetic", "initialize_HR_tmp");

    for (int i = 0; i < this->hR->size_atom_pairs(); ++i)
    {
        hamilt::AtomPair<TR>& tmp = this->hR->get_atom_pair(i);
        for(int ir = 0;ir < tmp.get_R_size(); ++ir )
        {
            const int* R_index = tmp.get_R_index(ir);
            const int iat1 = tmp.get_atom_i();
            const int iat2 = tmp.get_atom_j(); 

            hamilt::AtomPair<std::complex<double>> tmp1(iat1, iat2, R_index[0], R_index[1], R_index[2], paraV);
            this->hR_tmp->insert_pair(tmp1);
        }
    }
    this->hR_tmp->allocate(nullptr,true);

    ModuleBase::timer::tick("TDEkinetic", "initialize_HR_tmp");
}


template<typename TK, typename TR>
void TDEkinetic<OperatorLCAO<TK, TR>>::contributeHR()
{
    //const Parallel_Orbitals* paraV = this->hR->get_atom_pair(0).get_paraV();
    ModuleBase::TITLE("TDEkinetic", "contributeHR");
    ModuleBase::timer::tick("TDEkinetic", "contributeHR");
    //skip if not TDDFT velocity gauge
    if (elecstate::H_TDDFT_pw::stype != 1)
    {
        return;
    }

    if (!this->hR_tmp_done)
    {
        // if this Operator is the first node of the sub_chain, then hR_tmp is nullptr
        if (this->hR_tmp == nullptr)
        {
        this->hR_tmp = new hamilt::HContainer<std::complex<double>>(this->LM->ParaV);
        //allocate memory for hR_tmp use the same memory as hR
        this->initialize_HR_tmp(this->LM->ParaV);
        this->allocated = true;
        }
        if(this->next_sub_op != nullptr)
        {
            // pass pointer of hR_tmp to the next node
            static_cast<OperatorLCAO<TK, TR>*>(this->next_sub_op)->set_HR_fixed(this->hR_tmp);
        }
        // calculate the values in hR_tmp
        this->calculate_HR();
        this->hR_tmp_done = true;
    }

    ModuleBase::timer::tick("TDEkinetic", "contributeHR");
    return;
}

template<typename TK, typename TR>
void TDEkinetic<OperatorLCAO<TK, TR>>::contributeHk(int ik)
{
    return;
}
template<>
void TDEkinetic<OperatorLCAO<std::complex<double>, double>>::contributeHk(int ik)
{
    if (GlobalV::ESOLVER_TYPE != "tddft" || elecstate::H_TDDFT_pw::stype != 1)
    {
        return;
    }
    else{        
        ModuleBase::TITLE("TDEkinetic", "contributeHk");
        ModuleBase::timer::tick("TDEkinetic", "contributeHk");
        //folding inside HR to HK
        if(ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER())
        {
            const int nrow = this->LM->ParaV->get_row_size();
            hamilt::folding_HR(*this->hR_tmp, this->hK->data(), this->kvec_d[ik], nrow, 1);
        }
        else
        {
            const int ncol = this->LM->ParaV->get_col_size();
            hamilt::folding_HR(*this->hR_tmp, this->hK->data(), this->kvec_d[ik], ncol, 0);
        }
        
        ModuleBase::timer::tick("TDEkinetic", "contributeHk");
    }
}

template class TDEkinetic<hamilt::OperatorLCAO<double, double>>;
template class TDEkinetic<hamilt::OperatorLCAO<std::complex<double>, double>>;
template class TDEkinetic<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>>;

}