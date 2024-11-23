#include "read_pp.h"

#include "module_parameter/parameter.h"
#include <cmath>

#include <cstring> // Peize Lin fix bug about strcpy 2016-08-02
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "module_base/math_integral.h" // for numerical integration

Pseudopot_upf::Pseudopot_upf()
{
}

Pseudopot_upf::~Pseudopot_upf()
{
}

int Pseudopot_upf::init_pseudo_reader(const std::string &fn, std::string &type, Atom_pseudo& pp)
{
    ModuleBase::TITLE("Pseudopot_upf","init");
    // First check if this pseudo-potential has spin-orbit information
    std::ifstream ifs(fn.c_str(), std::ios::in);

	// can't find the file.
	if (!ifs)
    {
        return 1;
    }

    // if(GlobalV::global_pseudo_type=="auto") //zws
	if (type == "auto")
	{
		set_pseudo_type(fn, type);
	}

	int info = -1;
	// read in the .UPF type of pseudopotentials
	// if(GlobalV::global_pseudo_type=="upf")
	if (type == "upf")
	{
		info = read_pseudo_upf(ifs, pp);
	}
	// read in the .vwr type of pseudopotentials
	// else if(GlobalV::global_pseudo_type=="vwr")
	else if (type == "vwr")
	{
		info = read_pseudo_vwr(ifs, pp);
	}
	// else if(GlobalV::global_pseudo_type=="upf201")
	else if (type == "upf201")
	{
		info = read_pseudo_upf201(ifs, pp);
	}
	// else if(GlobalV::global_pseudo_type=="blps") // sunliang added 2021.7
	else if (type == "blps")
	{
		info = read_pseudo_blps(ifs, pp);
	}
    else
    {
        return 4;
    }

	return info;
}


//----------------------------------------------------------
// setting the type of the pseudopotential file
//----------------------------------------------------------
int Pseudopot_upf::set_pseudo_type(const std::string &fn, std::string &type) //zws add
{
    std::ifstream pptype_ifs(fn.c_str(), std::ios::in);
    std::string dummy;
	std::string strversion;

	if (pptype_ifs.good())
	{
		getline(pptype_ifs,dummy);

		std::stringstream wdsstream(dummy);
		getline(wdsstream,strversion,'"');
		getline(wdsstream,strversion,'"');

		if ( trim(strversion) == "2.0.1" )
		{
			type = "upf201";
			// GlobalV::global_pseudo_type = "upf201";
		}
		else
		{
			type = "upf";
			// GlobalV::global_pseudo_type = "upf";
		}
	}
	return 0;
}

std::string& Pseudopot_upf::trim(std::string &in_str)
{
    static const std::string deltri = " \t" ; // delete tab or space
    std::string::size_type position = in_str.find_first_of(deltri, 0);
    if (position == std::string::npos)
	{
        return in_str;
	}
    return trim(in_str.erase(position, 1) );
}

std::string Pseudopot_upf::trimend(std::string &in_str)
{
    const std::string &deltri =" \t" ;
    std::string::size_type position = in_str.find_last_not_of(deltri)+1;
    std::string tmpstr=in_str.erase(position);
    return tmpstr.erase(0,tmpstr.find_first_not_of(deltri));
} //zws


int Pseudopot_upf::average_p(const double& lambda, Atom_pseudo& pp)
{
    int error = 0;
    double lambda_ = lambda;
    if(!PARAM.inp.lspinorb) { lambda_ = 0.0; }
    if (pp.has_so && pp.tvanp)
    {
        error++;
        std::cout << "------------------------------------------------------" << std::endl;
        std::cout << " FR-USPP please use lspinorb=.true." << std::endl;
        std::cout << "------------------------------------------------------" << std::endl;
        return error;
    }
    if (!pp.has_so && PARAM.inp.lspinorb)
    {
        error++;
        std::cout << "warning_quit! no soc upf used for lspinorb calculation, error!" << std::endl;
        return error;
    }
    // ModuleBase::WARNING_QUIT("average_p", "no soc upf used for lspinorb calculation, error!");

    if (!pp.has_so || (PARAM.inp.lspinorb && std::abs(lambda_ - 1.0) < 1.0e-8))
    {
        return error;
    }

    //if(std::abs(lambda_)<1.0e-8)
	if(!PARAM.inp.lspinorb)
	{
		int new_nbeta = 0; //calculate the new nbeta
		for(int nb=0; nb< pp.nbeta; nb++)
		{
			new_nbeta++;
			if(pp.lll[nb] != 0 && std::abs(pp.jjj[nb] - pp.lll[nb] - 0.5) < 1e-6) //two J = l +- 0.5 average to one
			{
				new_nbeta--;
			}
		}

		pp.nbeta = new_nbeta;
		ModuleBase::matrix dion_new;
		dion_new.create(pp.nbeta, pp.nbeta);

		int old_nbeta=-1;
		for(int nb=0; nb<pp.nbeta; nb++)
		{
			old_nbeta++;
			int l = pp.lll[old_nbeta];
			int ind=0, ind1=0;
			if(l != 0)
			{
				if(std::abs(pp.jjj[old_nbeta] - pp.lll[old_nbeta] + 0.5) < 1e-6)
				{
					if(std::abs(pp.jjj[old_nbeta+1]-pp.lll[old_nbeta+1]-0.5)>1e-6) 
					{
						error = 1;
						std::cout<<"warning_quit! error beta function 1 !" <<std::endl;
						return error;
					}
					ind = old_nbeta +1;
					ind1 = old_nbeta;
				}
				else
				{
					if(std::abs(pp.jjj[old_nbeta+1]-pp.lll[old_nbeta+1]+0.5)>1e-6)
					{
						error = 1;
						std::cout<<"warning_quit! error beta function 2 !" <<std::endl;
						return error;
					}
					ind = old_nbeta;
					ind1 = old_nbeta +1;
				}
				double vion1 = ((l+1.0) * pp.dion(ind,ind) + l * pp.dion(ind1,ind1)) / (2.0*l+1.0);
				if(std::abs(vion1)<1.0e-8) { vion1 = 0.1;
}
				//average beta (betar)
				for(int ir = 0; ir<pp.mesh;ir++)
				{
					pp.betar(nb, ir) = 1.0 / (2.0 * l + 1.0) * 
							( (l + 1.0) * sqrt(std::abs(pp.dion(ind,ind) / vion1)) *
							pp.betar(ind, ir) + 
							l * sqrt(std::abs(pp.dion(ind1,ind1) / vion1)) *
							pp.betar(ind1, ir) ) ;
				}
				//average the dion matrix
				pp.dion(nb, nb) = vion1;
				old_nbeta++;	
			}
			else
			{
				for(int ir = 0; ir<pp.mesh;ir++) {
					pp.betar(nb, ir) = pp.betar(old_nbeta, ir);
}
				pp.dion(nb, nb) = pp.dion(old_nbeta, old_nbeta);
			}
			pp.lll[nb] = pp.lll[old_nbeta]; //reset the lll index, ignore jjj index
		}

		//store the old dion and then recreate dion 
		for(int i=0;i<pp.nbeta; i++)
		{
			for(int j=0;j<pp.nbeta;j++)
			{
				dion_new(i,j) = pp.dion(i,j);
			}
		}

		pp.dion = dion_new;
	//	pp.dion.create(pp.nbeta, pp.nbeta);
	//	for(int i=0;i<pp.nbeta; i++)
	//		for(int j=0;j<pp.nbeta;j++)
	//			pp.dion(i,j) = dion_new(i,j);
		
		int new_nwfc = 0;
		for(int nb=0; nb<pp.nchi; nb++)
		{
			new_nwfc++;
			if(pp.lchi[nb] != 0 && std::abs(pp.jchi[nb] - pp.lchi[nb] - 0.5)<1e-6)
			{
				new_nwfc--;
			}
		}

		pp.nchi = new_nwfc;
		int old_nwfc=-1;
		for(int nb=0; nb<pp.nchi; nb++)
		{
			old_nwfc++;
			int l = pp.lchi[old_nwfc];
			int ind=0, ind1=0;
			if(l!=0)
			{
				if(std::abs(pp.jchi[old_nwfc] - pp.lchi[old_nwfc] + 0.5) < 1e-6)
				{
					if(std::abs(pp.jchi[old_nwfc+1]-pp.lchi[old_nwfc+1]-0.5)>1e-6) 
					{error++; std::cout<<"warning_quit! error chi function 1 !"<<std::endl; return error;}
	//					ModuleBase::WARNING_QUIT("average_p", "error chi function 1 !");
					ind = old_nwfc +1;
					ind1 = old_nwfc;
				}
				else
				{
					if(std::abs(pp.jchi[old_nwfc+1]-pp.lchi[old_nwfc+1]+0.5)>1e-6)
					{error++; std::cout<<"warning_quit! error chi function 2 !"<<std::endl; return error;}
	//					ModuleBase::WARNING_QUIT("average_p", "error chi function 2 !");
					ind = old_nwfc;
					ind1 = old_nwfc +1;
				}
				//average chi
				for(int ir = 0; ir<pp.mesh;ir++)
				{
					pp.chi(nb, ir) = 1.0 / (2.0 * l + 1.0) *
						( (l+1.0)*pp.chi(ind,ir) + (l*pp.chi(ind1,ir)) );
				}
				old_nwfc++;
			}
			else{
				for(int ir = 0; ir<pp.mesh;ir++) {
					pp.chi(nb, ir) = pp.chi(old_nwfc, ir);
}
			}
			pp.lchi[nb] = pp.lchi[old_nwfc]; //reset lchi index
		}
		pp.has_so = false;	
		return error;
	}
	else//lambda_ != 0, modulate the soc effect in pseudopotential
	{
		for(int nb=0; nb<pp.nbeta; nb++)
		{
			int l = pp.lll[nb];
			int ind=0, ind1=0;
			if(l != 0)
			{
				if(std::abs(pp.jjj[nb] - pp.lll[nb] + 0.5) < 1e-6)
				{
					if(std::abs(pp.jjj[nb+1]-pp.lll[nb+1]-0.5)>1e-6) 
					{
						error = 1;
						std::cout<<"warning_quit! error beta function 1 !" <<std::endl;
						return error;
					}
					ind = nb +1;
					ind1 = nb;
				}
				else
				{
					if(std::abs(pp.jjj[nb+1]-pp.lll[nb+1]+0.5)>1e-6)
					{
						error = 1;
						std::cout<<"warning_quit! error beta function 2 !" <<std::endl;
						return error;
					}
					ind = nb;
					ind1 = nb +1;
				}
				double vion1 = ((l+1.0) * pp.dion(ind,ind) + l * pp.dion(ind1,ind1)) / (2.0*l+1.0);
				if(std::abs(vion1)<1.0e-10) { vion1 = 0.1;
}
				//average beta (betar)
				const double sqrtDplus = sqrt(std::abs(pp.dion(ind,ind) / vion1));
				const double sqrtDminus = sqrt(std::abs(pp.dion(ind1,ind1) / vion1));
				pp.dion(ind, ind) = vion1;
				pp.dion(ind1, ind1) = vion1;
				for(int ir = 0; ir<pp.mesh;ir++)
				{
					double avera = 1.0 / (2.0 * l + 1.0) * 
							( (l + 1.0) * sqrtDplus *
							pp.betar(ind, ir) + 
							l * sqrtDminus *
							pp.betar(ind1, ir) ) ;
					double delta = 1.0 / (2.0 * l + 1.0) * 
							( sqrtDplus *
							pp.betar(ind, ir) - 
							sqrtDminus *
							pp.betar(ind1, ir) ) ;
					pp.betar(ind, ir) = (avera + l * delta * lambda_) ;
					pp.betar(ind1, ir) = (avera - (l + 1) * delta * lambda_); 
				}
				nb++;
			}
		}

		for(int nb=0; nb<pp.nchi; nb++)
		{
			int l = pp.lchi[nb];
			int ind=0, ind1=0;
			if(l!=0)
			{
				if(std::abs(pp.jchi[nb] - pp.lchi[nb] + 0.5) < 1e-6)
				{
					if(std::abs(pp.jchi[nb+1]-pp.lchi[nb+1]-0.5)>1e-6) 
					{error++; std::cout<<"warning_quit! error chi function 1 !"<<std::endl; return error;}
					ind = nb +1;
					ind1 = nb;
				}
				else
				{
					if(std::abs(pp.jchi[nb+1]-pp.lchi[nb+1]+0.5)>1e-6)
					{error++; std::cout<<"warning_quit! error chi function 2 !"<<std::endl; return error;}
					ind = nb;
					ind1 = nb +1;
				}
				//average chi
				for(int ir = 0; ir<pp.mesh;ir++)
				{
					double avera = 0.5 * 
						( pp.chi(ind,ir) + pp.chi(ind1,ir) );
					double delta = 0.5 * 
						( pp.chi(ind,ir) - pp.chi(ind1,ir) );
					pp.chi(ind, ir) = avera + delta * lambda_ ; 
					pp.chi(ind1, ir) = avera - delta * lambda_ ; 
				}
				nb++;
			}
		}
		return error;
	}
}

// Peize Lin add for bsse 2021.04.07
void Pseudopot_upf::set_empty_element(Atom_pseudo& pp)
{
	pp.zv = 0;
	for(int ir=0; ir<pp.mesh; ++ir)
	{
		pp.vloc_at[ir] = 0;
	}
	for(int i=0; i<pp.nbeta; ++i)
	{
		for(int j=0; j<pp.nbeta; ++j)
		{
			pp.dion(i,j) = 0;
		}
	}
	for(int ir=0; ir<pp.mesh; ++ir)
	{
		pp.rho_at[ir] = 0;
	}
	return;
}

/**
 * For USPP we set the augmentation charge as an l-dependent array in all
 * cases. This is already the case when upf%q_with_l is .true.
 * For vanderbilt US pseudos, where nqf and rinner are non zero, we do here
 * what otherwise would be done multiple times in many parts of the code
 * (such as in init_us_1, addusforce_r, bp_calc_btq, compute_qdipol)
 * whenever the q_l(r) were to be constructed.
 * For simple rrkj3 pseudos we duplicate the information contained in q(r)
 * for all q_l(r).
 *
 * This requires a little extra memory but unifies the treatment of q_l(r)
 * and allows further weaking with the augmentation charge.
 */
void Pseudopot_upf::set_upf_q(Atom_pseudo& pp)
{
    if (pp.tvanp && !q_with_l)
    {
        pp.qfuncl.create(pp.nqlc, pp.nbeta * (pp.nbeta + 1) / 2, pp.mesh);
        for (int nb = 0; nb < pp.nbeta; nb++)
        {
            int ln = pp.lll[nb];
            for (int mb = nb; mb < pp.nbeta; mb++)
            {
                int lm = pp.lll[mb];
                int nmb = mb * (mb + 1) / 2 + nb;

                for (int l = std::abs(ln - lm); l <= ln + lm; l += 2)
                {
                    // copy q(r) to the l-dependent grid
                    for (int ir = 0; ir < pp.mesh; ir++)
                    {
                        pp.qfuncl(l, nmb, ir) = qfunc(nmb, ir);
                    }

                    // adjust the inner values on the l-dependent grid if nqf and rinner are defined
                    if (nqf > 0 && rinner[l] > 0.0)
                    {
                        int ilast = 0;
                        for (int ir = 0; ir < pp.kkbeta; ++ir)
                        {
                            if (pp.r[ir] < rinner[l])
                            {
                                ilast = ir + 1;
                            }
                            else
                            {
                                break;
                            }
                        }
                        this->setqfnew(nqf, ilast, l, 2, &(qfcoef(nb, mb, l, 0)), pp.r.data(), &(pp.qfuncl(l, nmb, 0)));
                    }
                }
            }
        }
    }
}

void Pseudopot_upf::setqfnew(const int& nqf,
                             const int& mesh,
                             const int& l,
                             const int& n,
                             const double* qfcoef,
                             const double* r,
                             double* rho)
{
    for (int ir = 0; ir < mesh; ++ir)
    {
        double rr = r[ir] * r[ir];
        rho[ir] = qfcoef[0];
        for (int iq = 1; iq < nqf; ++iq)
        {
            rho[ir] += qfcoef[iq] * pow(rr, iq);
        }
        rho[ir] *= pow(r[ir], l + n);
    }
}

void Pseudopot_upf::skip_number(std::ifstream& ifs, bool mesh_changed)
{
	if (mesh_changed)
	{
		double temp = 0.;
		ifs >> temp;
	}
}