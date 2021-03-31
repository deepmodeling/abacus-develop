#include "tools.h"
#include "global.h"
#include "potential.h"
#include "xc_functional.h"
#include "xc_gga_pw.h"
#include "efield.h"
#include "math.h"
#include "potential_libxc.h"
// new
#include "H_Hartree_pw.h"
#include "H_XC_pw.h"

Potential::Potential()
{
    vltot = new double[1];
    vr_eff1 = new double[1];
    this->out_potential = 0;
}

Potential::~Potential()
{
    delete[] vltot;
    delete[] vr_eff1;
}

void Potential::allocate(const int nrxx)
{
    TITLE("Potential","allocate");
    assert(nrxx>0);

    delete[] this->vltot;
    this->vltot = new double[nrxx];
    Memory::record("Potential","vltot",nrxx,"double");

    this->vr.create(NSPIN,nrxx);
    this->vr_eff.create(NSPIN,nrxx);
    Memory::record("Potential","vr",NSPIN*nrxx,"double");
    Memory::record("Potential","vr_eff",NSPIN*nrxx,"double");

    delete[] this->vr_eff1;
    this->vr_eff1 = new double[nrxx];
    Memory::record("Potential","vr_eff1",nrxx,"double");

    this->vnew.create(NSPIN,nrxx);
    Memory::record("Potential","vnew",NSPIN*nrxx,"double");

    return;
}

//----------------------------------------------------------
//  Initializes the self consistent potential 
//----------------------------------------------------------
void Potential::init_pot(
	const int &istep, // number of ionic steps
	ComplexMatrix &sf // structure factors
)
{
    TITLE("Potential","init_pot");
    timer::tick("Potential","init_pot");

    assert(istep>=0);

	// total potential in real space
    this->vr_eff.zero_out();

    // the vltot should and must be zero here.
    ZEROS(this->vltot, pw.nrxx);

	//-------------------------------------------------------------------
	// (1) local pseudopotential + electric field (if any) in vltot
	//-------------------------------------------------------------------
	this->set_local_pot(
		this->vltot, // 3D local pseudopotentials 
		ucell.ntype,
		pw.ngmc,
		ppcell.vloc,
		pw.ig2ngg,
		sf // structure factors		
	);

	// zhengdy-soc, pauli matrix, just index 0 has vlocal term
	int nspin0=NSPIN;

	if(NSPIN==4) 
	{
		nspin0=1;
	}

	for(int is=0; is<nspin0; ++is)
	{
		for(int ir=0; ir<pw.nrxx; ++ir)
		{
			this->vr_eff(is,ir) = this->vltot[ir];	
		}
	}

	// core correction potential.
	CHR.set_rho_core( pw.strucFac );

	//--------------------------------------------------------------------
	// (2) other effective potentials need charge density,
	// choose charge density from ionic step 0.
	//--------------------------------------------------------------------
    if(istep==0)
    {
        OUT(ofs_running,"start_pot",start_pot);

        cout << " START POTENTIAL      : " << start_pot << endl;
        if (this->start_pot == "atomic")//mohan add 2007-10-17
        {
            start_from_atomic:
            CHR.atomic_rho(NSPIN, CHR.rho);
        }
        else if (this->start_pot == "file")
        {
            ofs_running << " try to start potential from file : ";
            for(int is=0; is<NSPIN; is++)
            {
                stringstream ssc;
                ssc << global_out_dir << "SPIN" << is + 1 << "_CHG";
                ofs_running << ssc.str() << endl;
                // mohan update 2012-02-10
                if(CHR.read_rho( is, ssc.str(), CHR.rho[is] )) 
                {
                    ofs_running << " Read in the charge density: " << ssc.str() << endl;
				}
				else if(is>0 && NSPIN==4)
				{
					// read only spin (up+down)
					if(PRENSPIN == 1)
					{
						ofs_running << " Didn't read in the charge density but autoset it for spin " <<is+1<< endl;
						for(int ir=0;ir<pw.nrxx;ir++)
						{
							CHR.rho[is][ir] = 0.0;
						}
					}
					// 
					else if(PRENSPIN == 2)
					{//read up and down , then rearrange them.
						if(is==1) 
						{
							WARNING_QUIT("potential::init_pot","Incomplete charge density file!");
						}
						else if(is==2) 
						{
							ofs_running << " Didn't read in the charge density but would rearrange it later. "<< endl;
						}
						else if(is==3)
						{
							ofs_running << " rearrange charge density " << endl;
							for(int ir=0;ir<pw.nrxx;ir++)
							{
								CHR.rho[3][ir] = CHR.rho[0][ir] - CHR.rho[1][ir];
								CHR.rho[0][ir] = CHR.rho[0][ir] + CHR.rho[1][ir];
								CHR.rho[1][ir] = 0.0;
								CHR.rho[2][ir] = 0.0;
							}
						}
					}
					else
					{
						WARNING_QUIT("potential::init_pot","Incomplete charge density file!");
					}
				}
				else
                {
                    ofs_running << " Start charge density from atomic charge density." << endl;
                    goto start_from_atomic;
                }
            }
        }
        else
        {
            WARNING_QUIT("potential::init_pot","start_pot is wrong!");
        }
		
		// Peize Lin add 2020.04.04
		if(restart.info_load.load_charge && !restart.info_load.load_charge_finish)
		{
			for(int is=0; is<NSPIN; ++is)
			{
				restart.load_disk("charge", is);
			}
			restart.info_load.load_charge_finish = true;
		}
    }
    else
    {
        // the extrapolation part moves to ions.cpp.
    }

	// renormalize the charge density
    CHR.renormalize_rho();


    //----------------------------------------------------------
    // (3) compute Hartree and XC potentials saves in vr 
    //----------------------------------------------------------
    this->v_of_rho(CHR.rho, vr);

    //----------------------------------------------------------
    // (4) total potentials 
    //----------------------------------------------------------
    if(vext == 0) 
	{
		this->set_vr_eff();
	}
    else 
	{
		this->set_vrs_tddft(istep);
	}

	// plots
    //figure::picture(this->vr_eff1,pw.ncx,pw.ncy,pw.ncz);
    timer::tick("Potential","init_pot");
    return;
}


//==========================================================
// This routine computes the local potential in real space
//==========================================================
void Potential::set_local_pot(
	double* vl_pseudo, // store the local pseudopotential
	const int &ntype, // number of atom types
	const int &ngmc, // number of |g|, g is plane wave
	matrix &vloc, // local pseduopotentials
	int* ig2ngg, // ig2ngg
	ComplexMatrix &sf // structure factors	
)const
{
    TITLE("Potential","set_local_pot");
    timer::tick("Potential","set_local_pot");

    complex<double> *vg = new complex<double>[ngmc];

    ZEROS( vg, ngmc );

    for (int it=0; it<ntype; it++)
    {
        for (int ig=0; ig<ngmc; ig++)
        {
            vg[ig] += vloc(it, ig2ngg[ig]) * sf(it,ig);
        }
    }

    UFFT.ToRealSpace(vg, vl_pseudo); 

    delete[] vg;

    if(EFIELD && !DIPOLE)
    {
        Efield EFID;
        // in fact, CHR.rho is not used here.
        // if charge correction due to Efield is considered,
        // the structure here need to be updated.

        static bool first = true;
        if(first)
        {
            cout << " ADD THE EFIELD (V/A) : " << Efield::eamp*51.44 << endl;
            first = false;
        }
        EFID.add_efield(CHR.rho[0], vl_pseudo);	
    }

    //ofs_running <<" set local pseudopotential done." << endl;
    timer::tick("Potential","set_local_pot");
    return;
}

//==========================================================
// This routine computes the Hartree and Exchange and Correlation
// potential and energies which corresponds to a given charge density
// The XC potential is computed in real space, while the
// Hartree potential is computed in reciprocal space.
//==========================================================
void Potential::v_of_rho
(
    double **rho_in,
    matrix &v_in
)
{
    TITLE("Potential","v_of_rho");
    v_in.zero_out();

    timer::tick("Potential","v_of_rho",'E');

//----------------------------------------------------------
//  calculate the exchange-correlation potential
//----------------------------------------------------------
	
	#ifdef USE_LIBXC
    Potential_Libxc::v_xc(rho_in, H_XC_pw::etxc, H_XC_pw::vtxc, v_in);
	#else
    H_XC_pw::v_xc(pw.nrxx, pw.ncxyz, ucell.omega, rho_in, CHR.rho_core, v_in);
	#endif

//----------------------------------------------------------
//  calculate the Hartree potential
//----------------------------------------------------------
	H_Hartree_pw::v_hartree(ucell, pw, UFFT, NSPIN, v_in, rho_in);

    // mohan add 2011-06-20
    if(EFIELD && DIPOLE)
    {
        Efield EFID;
        for (int is = 0;is < NSPIN;is++)
        {
            EFID.add_efield(rho_in[is], &v_in.c[is*pw.nrxx]);
        }
    }
    timer::tick("Potential","v_of_rho",'E');
    return;
} //end subroutine v_of_rho



//==========================================================
// set the effective potential vr_eff on the real space grid 
// used in h_psi, adding the (spin dependent) scf (H+xc)
// part and the sum of all the local pseudopotential
// contributions.
//==========================================================
void Potential::set_vr_eff(void)
{
    TITLE("Potential","set_vr_eff");
    timer::tick("Potential","set_vr_eff");

    for (int is = 0;is < NSPIN;is++)
    {
        //=================================================================
        // define the total local potential (external + scf) for each spin
        //=================================================================
		if(NSPIN==4&&is>0)
		{
			for (int i = 0;i < pw.nrxx; i++)
			{
				this->vr_eff(is, i) = this->vr(is, i);
			}
		}
		else        
		{
			for (int i = 0;i < pw.nrxx; i++)
	        {
	            this->vr_eff(is, i) = this->vltot[i] + this->vr(is, i);
			}
		}
    }

    timer::tick("Potential","set_vr_eff");
    return;
}


// ----------------------------------------------------------------------
void Potential::newd(void)
{
    TITLE("Potential","newd");

    // distringuish non-local pseudopotential in REAL or RECIPROCAL space.
    // if in real space, call new_r
    // if in reciprocal space, call new_g

    // new g:
    //----------------------------------------------------------------------
    //  This routine computes the integral of the effective potential with
    //  the Q function and adds it to the bare ionic D term which is used
    //  to compute the non-local term in the US scheme.

	// no ultrasoft potentials: use bare coefficients for projectors
	// if( spin_orbital) ....
	// else if(noncolin) ....
	for (int iat=0; iat<ucell.nat; iat++)
	{
		const int it = ucell.iat2it[iat];
		const int nht = ucell.atoms[it].nh;
		// nht: number of beta functions per atom type
		for (int is = 0; is < NSPIN; is++)
		{
			for (int ih=0; ih<nht; ih++)
			{
				for (int jh=ih; jh<nht; jh++)
				{
					if(LSPINORB)
					{
						ppcell.deeq_nc(is , iat , ih , jh)= ppcell.dvan_so(is , it , ih , jh);
						ppcell.deeq_nc(is , iat , jh , ih)= ppcell.dvan_so(is , it , jh , ih);
					}
					else if( NSPIN==4 )
					{
						if(is==0)
						{
							ppcell.deeq_nc(is, iat, ih, jh) = ppcell.dvan(it, ih, jh);
							ppcell.deeq_nc(is, iat, jh, ih) = ppcell.dvan(it, ih, jh);
						}
						else if(is==1)
						{
							ppcell.deeq_nc(is, iat, ih, jh) = complex<double>(0.0 , 0.0);
							ppcell.deeq_nc(is, iat, jh, ih) = complex<double>(0.0 , 0.0);
						}
						else if(is==2)
						{
							ppcell.deeq_nc(is, iat, ih, jh) = complex<double>(0.0 , 0.0);
							ppcell.deeq_nc(is, iat, jh, ih) = complex<double>(0.0 , 0.0);
						}
						else if(is==3)
						{
							ppcell.deeq_nc(is, iat, ih, jh) = ppcell.dvan(it, ih, jh);
							ppcell.deeq_nc(is, iat, jh, ih) = ppcell.dvan(it, ih, jh);
						}
					}
					else
					{
						ppcell.deeq(is, iat, ih, jh) = ppcell.dvan(it, ih, jh);
						ppcell.deeq(is, iat, jh, ih) = ppcell.dvan(it, ih, jh);
					}
				}
			}
		}
	}
	return;
} // end subroutine newd
