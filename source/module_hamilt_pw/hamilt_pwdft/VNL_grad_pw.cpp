#include "VNL_in_pw.h"
#include "module_base/math_sphbes.h"
#include "module_base/timer.h"
#include "module_base/math_ylmreal.h"
#include "module_base/math_integral.h"
#include "module_base/math_polyint.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
void pseudopot_cell_vnl::initgradq_vnl(const UnitCell &cell)
{
    const int nbrx = 10;
	const int nbrx_nc = 20;
    const int ntype = cell.ntype;
	if(GlobalV::NSPIN!=4) 
	{
		this->tab_dq.create(ntype, nbrx, GlobalV::NQX);
	}
	else 
	{
		this->tab_dq.create(ntype, nbrx_nc, GlobalV::NQX);
	}
    gradvkb.create(3, nkb, this->wfcpw->npwk_max);

    const double pref = ModuleBase::FOUR_PI / sqrt(cell.omega);
    for (int it = 0;it < ntype;it++)  
	{
		const int nbeta = cell.atoms[it].ncpp.nbeta;
		int kkbeta = cell.atoms[it].ncpp.kkbeta;
		if ( (kkbeta%2 == 0) && kkbeta>0 )
		{
			kkbeta--;
		}

		double *djl = new double[kkbeta];
		double *aux  = new double[kkbeta];

		for (int ib = 0;ib < nbeta;ib++)
		{
			const int l = cell.atoms[it].ncpp.lll[ib];
			for (int iq=0; iq<GlobalV::NQX; iq++)  
			{
				const double q = iq * GlobalV::DQ;
				ModuleBase::Sphbes::dSpherical_Bessel_dx(kkbeta, cell.atoms[it].ncpp.r, q, l, djl);

				for (int ir = 0;ir < kkbeta;ir++)
				{
					aux[ir] = cell.atoms[it].ncpp.betar(ib, ir) *
					          djl[ir] * pow(cell.atoms[it].ncpp.r[ir],2);
				} 
				double vqint;
				ModuleBase::Integral::Simpson_Integral(kkbeta, aux, cell.atoms[it].ncpp.rab, vqint);
				this->tab_dq(it, ib, iq) = vqint * pref;
			} 
		} 
		delete[] aux;
		delete[] djl;
	}

}

void pseudopot_cell_vnl::getgradq_vnl(const int ik)
{
    if(GlobalV::test_pp) ModuleBase::TITLE("pseudopot_cell_vnl","getvnl");
	ModuleBase::timer::tick("pp_cell_vnl","getvnl");

	if(lmaxkb < 0) 
	{
		return;
	}

	const int npw = this->wfcpw->npwk[ik];

    // When the internal memory is large enough, it is better to make tmpgradvkb and tmpvkb be the number of pseudopot_cell_vnl
    // We only need to initialize them once as long as the cell is unchanged.
	ModuleBase::realArray tmpgradvkb(3, nhm, npw);
    ModuleBase::matrix tmpvkb(nhm, npw);
	double *vq = new double[npw];
	double *dvq = new double[npw];

	const int x1= (lmaxkb + 1)*(lmaxkb + 1);

	ModuleBase::matrix ylm(x1, npw);
    ModuleBase::matrix *dylm = new ModuleBase::matrix[3]{ModuleBase::matrix(x1,npw),ModuleBase::matrix(x1,npw),ModuleBase::matrix(x1,npw)};

	ModuleBase::Vector3<double> *gk = new ModuleBase::Vector3<double>[npw];
	for (int ig = 0;ig < npw;ig++) 
	{
		gk[ig] = this->wfcpw->getgpluskcar(ik,ig);
	}

	ModuleBase::YlmReal::grad_Ylm_Real(x1, npw, gk, ylm, dylm[0], dylm[1], dylm[2]);

	int jkb = 0;
	for(int it = 0;it < GlobalC::ucell.ntype;it++)
	{
		// calculate beta in G-space using an interpolation table
		const int nbeta = GlobalC::ucell.atoms[it].ncpp.nbeta;
		const int nh = GlobalC::ucell.atoms[it].ncpp.nh;
        int nb0 = -1;
        for( int ih = 0; ih < nh; ++ih)
        {
            int nb = this->indv(it, ih);
            if(nb != nb0)
            {
                for (int ig = 0;ig < npw;++ig)
			    {
			    	const double gnorm = gk[ig].norm() * GlobalC::ucell.tpiba;
			    	vq [ig] = ModuleBase::PolyInt::Polynomial_Interpolation(
			    			this->tab, it, nb, GlobalV::NQX, GlobalV::DQ, gnorm );
			    	dvq[ig] =ModuleBase::PolyInt::Polynomial_Interpolation(
			    			this->tab_dq, it, nb, GlobalV::NQX, GlobalV::DQ, gnorm );
			    }
                nb0 = nb;
            }
			double lmmat[9] = {0, 0, 1, -1, 0, 0, 0, -1, 0};
            for(int id = 0; id < 3; ++id)
            {
                const int lm = static_cast<int>( nhtolm(it, ih) );
			    for (int ig = 0;ig < npw;++ig)
			    {
                    ModuleBase::Vector3<double> gg = gk[ig];
                    double ggnorm = gg.norm();
					if(ggnorm < 1e-8)
					{
						if(lm == 0 || lm > 3)
						{
							tmpgradvkb(id, ih, ig) = 0.0;
						}
						else//lm = 1,2,3;  l = 1
						{
							//q \to 0 : \nabla(f(q)Y(\hat{q})) = f(q)/q*sqrt(3/4/pi)*vec(-delta(lm,2), -delta(lm,3), delta(lm,1))
							// tmpgradvkb(id, ih, ig) = 0.0;
							tmpgradvkb(id, ih, ig) = dvq[ig] * sqrt(3.0/4.0/M_PI) * lmmat[(lm-1)*3 + id];
						}
					}
					else
					{
                        tmpgradvkb(id, ih, ig) = ylm(lm, ig) * dvq[ig] * gg[id] / ggnorm
                                                 + dylm[id](lm, ig) / this->wfcpw->tpiba
                                                       * vq[ig]; // note: dylm/d(tpiba * gx) = 1/tpiba * dylm/dgx
                    }
					tmpvkb(ih, ig) = ylm(lm,ig) * vq[ig];
			    }
            }
        }

		// vkb1 contains all betas including angular part for type nt
		// now add the structure factor and factor (-i)^l
		for (int ia=0; ia<GlobalC::ucell.atoms[it].na; ia++) 
		{
            std::complex<double> *sk = this->psf->get_sk(ik, it, ia, this->wfcpw);

            for (int ih = 0;ih < nh;++ih)
			{
				std::complex<double> pref = pow( ModuleBase::NEG_IMAG_UNIT, nhtol(it, ih));
				std::complex<double>* pvkb = &this->vkb(jkb, 0);
				for (int id = 0; id < 3 ; ++id)
            	{
					std::complex<double>* pgvkb = &this->gradvkb(id, jkb, 0);
					for (int ig = 0;ig < npw;++ig)
					{
                	    std::complex<double> skig = sk[ig];
                	    pvkb[ig] = tmpvkb(ih, ig) * skig * pref;
						// std::complex<double> dskig = ModuleBase::NEG_IMAG_UNIT * (GlobalC::ucell.atoms[it].tau[ia][id] * this->wfcpw->lat0) * skig;
						// pgvkb[ig] = tmpgradvkb(id, ih, ig) * skig * pref +  tmpvkb(ih, ig) * dskig * pref;
						// The second term will be eliminate when doing <psi|beta>Dij<beta|psi> or we can say (\nabla_q+\nabla_q')S(q'-q) = 0
						pgvkb[ig] = tmpgradvkb(id, ih, ig) * skig * pref;
					}
				} //end id
				++jkb;
			} // end ih
            
			delete [] sk;
		} // end ia
	} // enddo

	delete [] gk;
	delete [] vq;
    delete [] dvq;
    delete [] dylm;

	ModuleBase::timer::tick("pp_cell_vnl","getvnl");

	return;
}