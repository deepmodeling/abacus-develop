//This file contains interfaces with libtorch,
//including loading of model and calculating gradients
//as well as subroutines that prints the results for checking

//The file contains 8 subroutines:
// cal_gvepsl : gvepsl is used for training with stress label, which is derivative of
//       descriptors wrt strain tensor, calculated by
//       d(des)/d\epsilon_{ab} = d(pdm)/d\epsilon_{ab} * d(des)/d(pdm) = gdm_epsl * gvdm
//       using einsum
// cal_gvdm : d(des)/d(pdm)
//       calculated using torch::autograd::grad
// load_model : loads model for applying V_delta
// prepare_psialpha : prepare psialpha for outputting npy file
// prepare_gevdm : prepare gevdm for outputting npy file

#ifdef __DEEPKS

#include "LCAO_deepks.h"
#include "LCAO_deepks_io.h" // mohan add 2024-07-22
#include "module_base/blas_connector.h"
#include "module_base/constants.h"
#include "module_base/libm/libm.h"
#include "module_base/parallel_reduce.h"
#include "module_hamilt_lcao/module_hcontainer/atom_pair.h"
#include "module_parameter/parameter.h"

// calculates stress of descriptors from gradient of projected density matrices
void LCAO_Deepks::cal_gvepsl(const int nat) {
    ModuleBase::TITLE("LCAO_Deepks", "cal_gvepsl");
    // preconditions
    this->cal_gvdm(nat);
    if (!gdmepsl_vector.empty()) {
        gdmepsl_vector.erase(gdmepsl_vector.begin(), gdmepsl_vector.end());
    }
    // gdmr_vector : nat(derivative) * 3 * inl(projector) * nm * nm
    if (GlobalV::MY_RANK == 0) {
        // make gdmx as tensor
        int nlmax = this->inlmax / nat;
        for (int nl = 0; nl < nlmax; ++nl) {
            std::vector<torch::Tensor> bmmv;
            // for (int ipol=0;ipol<6;++ipol)
            //{
            //     std::vector<torch::Tensor> xmmv;
            for (int i = 0; i < 6; ++i) {
                std::vector<torch::Tensor> ammv;
                for (int iat = 0; iat < nat; ++iat) {
                    int inl = iat * nlmax + nl;
                    int nm = 2 * this->inl_l[inl] + 1;
                    std::vector<double> mmv;
                    for (int m1 = 0; m1 < nm; ++m1) {
                        for (int m2 = 0; m2 < nm; ++m2) {
                            mmv.push_back(this->gdm_epsl[i][inl][m1 * nm + m2]);
                        }
                    } // nm^2
                    torch::Tensor mm
                        = torch::tensor(
                              mmv,
                              torch::TensorOptions().dtype(torch::kFloat64))
                              .reshape({nm, nm}); // nm*nm
                    ammv.push_back(mm);
				}
				torch::Tensor bmm = torch::stack(ammv, 0); // nat*nm*nm
				bmmv.push_back(bmm);
			}
            // torch::Tensor bmm = torch::stack(xmmv, 0);  //3*nat*nm*nm
            // bmmv.push_back(bmm);
            //}
            this->gdmepsl_vector.push_back(
                torch::stack(bmmv, 0)); // nbt*3*nat*nm*nm
        }
        assert(this->gdmepsl_vector.size() == nlmax);

        // einsum for each inl:
        // gdmepsl_vector : b:npol * a:inl(projector) * m:nm * n:nm
        // gevdm_vector : a:inl * v:nm (descriptor) * m:nm (pdm, dim1) * n:nm
        // (pdm, dim2) gvepsl_vector : b:npol * a:inl(projector) *
        // m:nm(descriptor)
		std::vector<torch::Tensor> gvepsl_vector;
		for (int nl = 0; nl < nlmax; ++nl) {
			gvepsl_vector.push_back(
					at::einsum("bamn, avmn->bav",
						{this->gdmepsl_vector[nl], this->gevdm_vector[nl]}));
		}

		// cat nv-> \sum_nl(nv) = \sum_nl(nm_nl)=des_per_atom
		// concatenate index a(inl) and m(nm)
		this->gvepsl_tensor = torch::cat(gvepsl_vector, -1);
		assert(this->gvepsl_tensor.size(0) == 6);
		assert(this->gvepsl_tensor.size(1) == nat);
		assert(this->gvepsl_tensor.size(2) == this->des_per_atom);
	}

	return;
}


// dDescriptor / dprojected density matrix
void LCAO_Deepks::cal_gvdm(const int nat) {
    ModuleBase::TITLE("LCAO_Deepks", "cal_gvdm");
    if (!gevdm_vector.empty()) {
        gevdm_vector.erase(gevdm_vector.begin(), gevdm_vector.end());
    }
    // cal gevdm(d(EigenValue(D))/dD)
    int nlmax = inlmax / nat;
    for (int nl = 0; nl < nlmax; ++nl) {
        std::vector<torch::Tensor> avmmv;
        for (int iat = 0; iat < nat; ++iat) {
            int inl = iat * nlmax + nl;
            int nm = 2 * this->inl_l[inl] + 1;
            // repeat each block for nm times in an additional dimension
            torch::Tensor tmp_x
                = this->pdm_tensor[inl].reshape({nm, nm}).unsqueeze(0).repeat(
                    {nm, 1, 1});
            // torch::Tensor tmp_y = std::get<0>(torch::symeig(tmp_x, true));
            torch::Tensor tmp_y = std::get<0>(torch::linalg::eigh(tmp_x, "U"));
            torch::Tensor tmp_yshell
                = torch::eye(nm, torch::TensorOptions().dtype(torch::kFloat64));
            std::vector<torch::Tensor> tmp_rpt; // repeated-pdm-tensor (x)
            std::vector<torch::Tensor> tmp_rdt; // repeated-d-tensor (y)
            std::vector<torch::Tensor> tmp_gst; // gvx-shell
            tmp_rpt.push_back(tmp_x);
            tmp_rdt.push_back(tmp_y);
            tmp_gst.push_back(tmp_yshell);
            std::vector<torch::Tensor> tmp_res;
            tmp_res
                = torch::autograd::grad(tmp_rdt,
                                        tmp_rpt,
                                        tmp_gst,
                                        false,
                                        false,
                                        /*allow_unused*/ true); // nm(v)**nm*nm
            avmmv.push_back(tmp_res[0]);
        }
        torch::Tensor avmm = torch::stack(avmmv, 0); // nat*nv**nm*nm
        this->gevdm_vector.push_back(avmm);
    }
    assert(this->gevdm_vector.size() == nlmax);
    return;
}

void LCAO_Deepks::load_model(const std::string& deepks_model) {
    ModuleBase::TITLE("LCAO_Deepks", "load_model");

    try {
        this->module = torch::jit::load(deepks_model);
    } catch (const c10::Error& e)

    {
        std::cerr << "error loading the model" << std::endl;
        return;
    }
    return;
}

// prepare_psialpha and prepare_gevdm for deepks_v_delta = 2
void LCAO_Deepks::prepare_psialpha(const int nlocal,
    const int nat,
    const UnitCell &ucell,
    const LCAO_Orbitals &orb,
    Grid_Driver &GridD)
{
    ModuleBase::TITLE("LCAO_Deepks", "prepare_psialpha");
    int nlmax = this->inlmax/nat;
    int mmax = 2*this->lmaxd+1;
    this->psialpha_tensor = torch::zeros({ nat, nlmax, 1, nlocal, mmax }, torch::TensorOptions().dtype(torch::kFloat64)); // support gamma-only

    //cutoff for alpha is same for all types of atoms
    const double Rcut_Alpha = orb.Alpha[0].getRcut();

    for (int T0 = 0; T0 < ucell.ntype; T0++)
    {
		Atom* atom0 = &ucell.atoms[T0]; 
        for (int I0 =0; I0< atom0->na; I0++)
        {
            //iat: atom index on which |alpha> is located
            const int iat = ucell.itia2iat(T0,I0);
			const ModuleBase::Vector3<double> tau0 = atom0->tau[I0];
            GridD.Find_atom(ucell, atom0->tau[I0] ,T0, I0);

            //outermost loop : find all adjacent atoms
            for (int ad=0; ad<GridD.getAdjacentNum()+1 ; ++ad)
            {
                const int T1 = GridD.getType(ad);
                const int I1 = GridD.getNatom(ad);
                const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
				const double Rcut_AO1 = orb.Phi[T1].getRcut();

                const ModuleBase::Vector3<double> tau1 = GridD.getAdjacentTau(ad);
				const Atom* atom1 = &ucell.atoms[T1];
				const int nw1_tot = atom1->nw*PARAM.globalv.npol;

				const double dist1 = (tau1-tau0).norm() * ucell.lat0;

				if (dist1 > Rcut_Alpha + Rcut_AO1)
				{
					continue;
				}

                //middle loop : all atomic basis on the adjacent atom ad
				for (int iw1=0; iw1<nw1_tot; ++iw1)
				{
					const int iw1_all = start1 + iw1;
					const int iw1_local = pv->global2local_row(iw1_all);
					const int iw2_local = pv->global2local_col(iw1_all);
					if(iw1_local < 0 || iw2_local < 0) {continue;
}
					const int iw1_0 = iw1/PARAM.globalv.npol;
					std::vector<double> nlm = this->nlm_save[iat][ad][iw1][0];
                    
                    int ib=0;
                    int nl=0;
                    for (int L0 = 0; L0 <= orb.Alpha[0].getLmax();++L0)
                    {
                        for (int N0 = 0;N0 < orb.Alpha[0].getNchi(L0);++N0)
                        {
                            const int nm = 2*L0+1;
                            
                            for (int m1=0; m1<nm; ++m1) // nm = 1 for s, 3 for p, 5 for d
                            {
                                this->psialpha_tensor[iat][nl][0][iw1_all][m1] = nlm[ib+m1];
                            }
                            ib+=nm;
                            nl++;
                        }
                    }  
				}//end iw
			}//end ad
		}//end I0
	}//end T0

#ifdef __MPI
    double msg[mmax];
    for(int iat=0; iat< nat ; iat++)
    {
        for(int nl = 0; nl < nlmax; nl++)
        {
            for(int mu = 0; mu < nlocal ; mu++)
            {
                for(int m=0;m<mmax;m++) 
				{
					msg[m] = this->psialpha_tensor[iat][nl][0][mu][m].item().toDouble();
				}
				Parallel_Reduce::reduce_all(msg,mmax);
				for(int m=0;m<mmax;m++) 
				{ 
					this->psialpha_tensor[iat][nl][0][mu][m] = msg[m];
				}
			}
        }
    }

#endif  
}

void LCAO_Deepks::prepare_psialpha_k(const int nlocal,
    const int nat,
    const int nks,
    const std::vector<ModuleBase::Vector3<double>> &kvec_d,
    const UnitCell &ucell,
    const LCAO_Orbitals &orb,
    Grid_Driver &GridD)
{
    ModuleBase::TITLE("LCAO_Deepks", "prepare_psialpha");
    int nlmax = this->inlmax/nat;
    int mmax = 2*this->lmaxd+1;
    this->psialpha_tensor = torch::zeros({ nat, nlmax, nks, nlocal, mmax }, torch::kComplexDouble); // support multi-k

    //cutoff for alpha is same for all types of atoms
    const double Rcut_Alpha = orb.Alpha[0].getRcut();

    for (int T0 = 0; T0 < ucell.ntype; T0++)
    {
		Atom* atom0 = &ucell.atoms[T0]; 
        for (int I0 =0; I0< atom0->na; I0++)
        {
            //iat: atom index on which |alpha> is located
            const int iat = ucell.itia2iat(T0,I0);
			const ModuleBase::Vector3<double> tau0 = atom0->tau[I0];
            GridD.Find_atom(ucell, atom0->tau[I0] ,T0, I0);

            //outermost loop : find all adjacent atoms
            for (int ad=0; ad<GridD.getAdjacentNum()+1 ; ++ad)
            {
                const int T1 = GridD.getType(ad);
                const int I1 = GridD.getNatom(ad);
                const int ibt = ucell.itia2iat(T1, I1);
                const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
				const double Rcut_AO1 = orb.Phi[T1].getRcut();

                const ModuleBase::Vector3<double> tau1 = GridD.getAdjacentTau(ad);
				const Atom* atom1 = &ucell.atoms[T1];
				const int nw1_tot = atom1->nw*PARAM.globalv.npol;

				const double dist1 = (tau1-tau0).norm() * ucell.lat0;

				if (dist1 > Rcut_Alpha + Rcut_AO1)
				{
					continue;
				}

                ModuleBase::Vector3<double> dR(GridD.getBox(ad).x,
                                               GridD.getBox(ad).y,
                                               GridD.getBox(ad).z);

                key_tuple key(ibt, dR.x, dR.y, dR.z);

                if (this->nlm_save_k[iat].find(key)
                    == this->nlm_save_k[iat].end()) 
                {
                    continue;
                }

                //middle loop : all atomic basis on the adjacent atom ad
				for (int iw1=0; iw1<nw1_tot; ++iw1)
				{
					const int iw1_all = start1 + iw1;
					const int iw1_local = pv->global2local_row(iw1_all);
					const int iw2_local = pv->global2local_col(iw1_all);
					if(iw1_local < 0 || iw2_local < 0) {continue;}
					const int iw1_0 = iw1/PARAM.globalv.npol;
					std::vector<double> nlm = this->nlm_save_k[iat][key][iw1][0];
                    
                    for (int ik = 0; ik <nks; ik++)
                    {
                        const double arg = kvec_d[ik] * dR * ModuleBase::TWO_PI;
                        double sinp, cosp;
                        ModuleBase::libm::sincos(arg, &sinp, &cosp);
                        const std::complex<double> kphase = std::complex<double>(cosp, sinp);
                        int ib=0;
                        int nl=0;
                        for (int L0 = 0; L0 <= orb.Alpha[0].getLmax();++L0)
                        {
                            for (int N0 = 0;N0 < orb.Alpha[0].getNchi(L0);++N0)
                            {
                                const int nm = 2*L0+1;
                                
                                for (int m1=0; m1<nm; ++m1) // nm = 1 for s, 3 for p, 5 for d
                                {
                                    std::complex<double> nlm_phase = nlm[ib + m1] * kphase;
                                    torch::Tensor nlm_tensor = torch::tensor({nlm_phase.real(), nlm_phase.imag()}, torch::kDouble);
                                    torch::Tensor complex_tensor = torch::complex(nlm_tensor[0], nlm_tensor[1]);
                                    this->psialpha_tensor[iat][nl][ik][iw1_all][m1] = complex_tensor;
                                }
                                ib+=nm;
                                nl++;
                            }
                        }  
                    }//end ik
				}//end iw
			}//end ad
		}//end I0
	}//end T0

#ifdef __MPI
    std::complex<double> msg[mmax];
    for(int iat=0; iat< nat ; iat++)
    {
        for(int nl = 0; nl < nlmax; nl++)
        {
            for(int ik = 0; ik < nks; ik++)
            {
                for(int mu = 0; mu < nlocal ; mu++)
                {
                    for(int m=0;m<mmax;m++) 
                    {
                        auto tensor_value = this->psialpha_tensor.index({iat, nl, ik, mu, m});
                        msg[m] = std::complex<double>(torch::real(tensor_value).item<double>(), torch::imag(tensor_value).item<double>());
                    }
                    Parallel_Reduce::reduce_all(msg,mmax);
                    for(int m=0;m<mmax;m++) 
                    {
                        torch::Tensor msg_tensor = torch::tensor({msg[m].real(), msg[m].imag()}, torch::kDouble);
                        torch::Tensor complex_tensor = torch::complex(msg_tensor[0], msg_tensor[1]);
                        this->psialpha_tensor[iat][nl][ik][mu][m] = complex_tensor;
                    }
                }
            }
        }
    }

#endif  
}

void LCAO_Deepks::check_vdp_psialpha(const int nat, const int nks, const int nlocal)
{
    std::ofstream ofs("vdp_psialpha.dat");
    ofs << std::setprecision(10);
    
    int nlmax = this->inlmax/nat;
    int mmax = 2*this->lmaxd+1;
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
                        ofs << this->psialpha_tensor.index({ iat,nl, iks, mu, m }).item().toDouble() << " ";
                    }
                }                
            }
            ofs << std::endl;
        }
    }
    ofs.close();
}

void LCAO_Deepks::prepare_gevdm(
    const int nat,
    const LCAO_Orbitals &orb)
{
    int nlmax = this->inlmax/nat;
    int mmax = 2*this->lmaxd+1;
    this->gevdm_tensor = torch::zeros({ nat, nlmax, mmax, mmax, mmax}, torch::TensorOptions().dtype(torch::kFloat64));

    this->cal_gvdm(nat);

    int nl=0;
    for (int L0 = 0; L0 <= orb.Alpha[0].getLmax();++L0)
    {
        for (int N0 = 0;N0 < orb.Alpha[0].getNchi(L0);++N0)
        {
            for (int iat = 0; iat < nat; iat++)
            {
                const int nm = 2*L0+1;            
                for (int v=0; v<nm; ++v) // nm = 1 for s, 3 for p, 5 for d
                {
                    for (int m=0; m<nm; ++m)
                    {
                        for (int n=0; n<nm;++n)
                        {
                            this->gevdm_tensor[iat][nl][v][m][n] = this->gevdm_vector[nl][iat][v][m][n];
                        }
                    }
                }
            }
            nl++;
        }
    }
    assert(nl == nlmax);
}

void LCAO_Deepks::check_vdp_gevdm(const int nat)
{
    std::ofstream ofs("vdp_gevdm.dat");
    ofs << std::setprecision(10);

    int nlmax = this->inlmax/nat;
    int mmax = 2*this->lmaxd+1;
    for(int iat=0; iat< nat ; iat++) 
    {
        for(int nl = 0; nl < nlmax; nl++)
        {
            for(int v=0; v<mmax; v++)
            {
                for(int m=0; m<mmax; m++)
                {
                    for(int n=0; n<mmax; n++)
                    {
                        ofs << this->gevdm_tensor.index({ iat,nl, v, m, n }).item().toDouble() << " ";
                    }
                }
            }
            ofs << std::endl;
        }
    }
    ofs.close();
}

#endif
