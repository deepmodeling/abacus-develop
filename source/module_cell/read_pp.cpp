#include "read_pp.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <sstream>
#include <cstring> // Peize Lin fix bug about strcpy 2016-08-02



Pseudopot_upf::Pseudopot_upf()
{
    functional_error = 0; // xiaohui add 2015-03-24
}

Pseudopot_upf::~Pseudopot_upf()
{
    delete[] r;
    delete[] rab;
    delete[] rho_atc;
    delete[] vloc;
    delete[] rho_at;
    delete[] lll;
    delete[] kbeta;
    delete[] els;
    delete[] els_beta;
    delete[] nchi;
    delete[] lchi;
    delete[] oc;
    delete[] epseu;
    delete[] rcut_chi;
    delete[] rcutus_chi;
    delete[] rinner;
    delete[] rcut;
    delete[] rcutus;
    delete[] nn;
    delete[] jchi;
    delete[] jjj;
}

int Pseudopot_upf::init_pseudo_reader(const std::string &fn, std::string &type)
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

	// read in the .UPF type of pseudopotentials
	// if(GlobalV::global_pseudo_type=="upf")
	if (type == "upf")
	{
		int info = read_pseudo_upf(ifs);
		return info;
	}
	// read in the .vwr type of pseudopotentials
	// else if(GlobalV::global_pseudo_type=="vwr")
	else if (type == "vwr")
	{
		int info = read_pseudo_vwr(ifs);
		return info;
	}
	// else if(GlobalV::global_pseudo_type=="upf201")
	else if (type == "upf201")
	{
		int info = read_pseudo_upf201(ifs);
		return info;
	}
	// else if(GlobalV::global_pseudo_type=="blps") // sunliang added 2021.7
	else if (type == "blps")
	{
		int info = read_pseudo_blps(ifs);
		return info;
	}

	return 0;
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


int Pseudopot_upf::average_p(const double& lambda)
{
	int error = 0;
	double lambda_ = lambda;
	if(!GlobalV::LSPINORB) lambda_ = 0.0;
	if(!this->has_so && GlobalV::LSPINORB) 
	{
		error++; 
		std::cout<<"warning_quit! no soc upf used for lspinorb calculation, error!"<<std::endl; 
		return error;
	}
	//ModuleBase::WARNING_QUIT("average_p", "no soc upf used for lspinorb calculation, error!");

	if(!this->has_so || (GlobalV::LSPINORB && std::abs(lambda_ - 1.0) < 1.0e-8) )
	{
		return error; 
	}

	//if(std::abs(lambda_)<1.0e-8)
	if(!GlobalV::LSPINORB)
	{
		int new_nbeta = 0; //calculate the new nbeta
		for(int nb=0; nb< this->nbeta; nb++)
		{
			new_nbeta++;
			if(this->lll[nb] != 0 && std::abs(this->jjj[nb] - this->lll[nb] - 0.5) < 1e-6) //two J = l +- 0.5 average to one
			{
				new_nbeta--;
			}
		}

		this->nbeta = new_nbeta;
		ModuleBase::matrix dion_new;
		dion_new.create(this->nbeta, this->nbeta);

		int old_nbeta=-1;
		for(int nb=0; nb<this->nbeta; nb++)
		{
			old_nbeta++;
			int l = this->lll[old_nbeta];
			int ind=0, ind1=0;
			if(l != 0)
			{
				if(std::abs(this->jjj[old_nbeta] - this->lll[old_nbeta] + 0.5) < 1e-6)
				{
					if(std::abs(this->jjj[old_nbeta+1]-this->lll[old_nbeta+1]-0.5)>1e-6) 
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
					if(std::abs(this->jjj[old_nbeta+1]-this->lll[old_nbeta+1]+0.5)>1e-6)
					{
						error = 1;
						std::cout<<"warning_quit! error beta function 2 !" <<std::endl;
						return error;
					}
					ind = old_nbeta;
					ind1 = old_nbeta +1;
				}
				double vion1 = ((l+1.0) * this->dion(ind,ind) + l * this->dion(ind1,ind1)) / (2.0*l+1.0);
				if(std::abs(vion1)<1.0e-8) vion1 = 0.1;
				//average beta (betar)
				for(int ir = 0; ir<this->mesh;ir++)
				{
					this->beta(nb, ir) = 1.0 / (2.0 * l + 1.0) * 
							( (l + 1.0) * sqrt(std::abs(this->dion(ind,ind) / vion1)) *
							this->beta(ind, ir) + 
							l * sqrt(std::abs(this->dion(ind1,ind1) / vion1)) *
							this->beta(ind1, ir) ) ;
				}
				//average the dion matrix
				this->dion(nb, nb) = vion1;
				old_nbeta++;	
			}
			else
			{
				for(int ir = 0; ir<this->mesh;ir++)
					this->beta(nb, ir) = this->beta(old_nbeta, ir);
				this->dion(nb, nb) = this->dion(old_nbeta, old_nbeta);
			}
			this->lll[nb] = this->lll[old_nbeta]; //reset the lll index, ignore jjj index
		}

		//store the old dion and then recreate dion 
		for(int i=0;i<this->nbeta; i++)
		{
			for(int j=0;j<this->nbeta;j++)
			{
				dion_new(i,j) = this->dion(i,j);
			}
		}

		this->dion = dion_new;
	//	this->dion.create(this->nbeta, this->nbeta);
	//	for(int i=0;i<this->nbeta; i++)
	//		for(int j=0;j<this->nbeta;j++)
	//			this->dion(i,j) = dion_new(i,j);
		
		int new_nwfc = 0;
		for(int nb=0; nb<this->nwfc; nb++)
		{
			new_nwfc++;
			if(this->lchi[nb] != 0 && std::abs(this->jchi[nb] - this->lchi[nb] - 0.5)<1e-6)
			{
				new_nwfc--;
			}
		}

		this->nwfc = new_nwfc;
		int old_nwfc=-1;
		for(int nb=0; nb<this->nwfc; nb++)
		{
			old_nwfc++;
			int l = this->lchi[old_nwfc];
			int ind=0, ind1=0;
			if(l!=0)
			{
				if(std::abs(this->jchi[old_nwfc] - this->lchi[old_nwfc] + 0.5) < 1e-6)
				{
					if(std::abs(this->jchi[old_nwfc+1]-this->lchi[old_nwfc+1]-0.5)>1e-6) 
					{error++; std::cout<<"warning_quit! error chi function 1 !"<<std::endl; return error;}
	//					ModuleBase::WARNING_QUIT("average_p", "error chi function 1 !");
					ind = old_nwfc +1;
					ind1 = old_nwfc;
				}
				else
				{
					if(std::abs(this->jchi[old_nwfc+1]-this->lchi[old_nwfc+1]+0.5)>1e-6)
					{error++; std::cout<<"warning_quit! error chi function 2 !"<<std::endl; return error;}
	//					ModuleBase::WARNING_QUIT("average_p", "error chi function 2 !");
					ind = old_nwfc;
					ind1 = old_nwfc +1;
				}
				//average chi
				for(int ir = 0; ir<this->mesh;ir++)
				{
					this->chi(nb, ir) = 1.0 / (2.0 * l + 1.0) *
						( (l+1.0)*this->chi(ind,ir) + (l*this->chi(ind1,ir)) );
				}
				old_nwfc++;
			}
			else{
				for(int ir = 0; ir<this->mesh;ir++)
					this->chi(nb, ir) = this->chi(old_nwfc, ir);
			}
			this->lchi[nb] = this->lchi[old_nwfc]; //reset lchi index
		}
		this->has_so = 0;	
		return error;
	}
	else//lambda_ != 0, modulate the soc effect in pseudopotential
	{
		for(int nb=0; nb<this->nbeta; nb++)
		{
			int l = this->lll[nb];
			int ind=0, ind1=0;
			if(l != 0)
			{
				if(std::abs(this->jjj[nb] - this->lll[nb] + 0.5) < 1e-6)
				{
					if(std::abs(this->jjj[nb+1]-this->lll[nb+1]-0.5)>1e-6) 
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
					if(std::abs(this->jjj[nb+1]-this->lll[nb+1]+0.5)>1e-6)
					{
						error = 1;
						std::cout<<"warning_quit! error beta function 2 !" <<std::endl;
						return error;
					}
					ind = nb;
					ind1 = nb +1;
				}
				double vion1 = ((l+1.0) * this->dion(ind,ind) + l * this->dion(ind1,ind1)) / (2.0*l+1.0);
				if(std::abs(vion1)<1.0e-10) vion1 = 0.1;
				//average beta (betar)
				const double sqrtDplus = sqrt(std::abs(this->dion(ind,ind) / vion1));
				const double sqrtDminus = sqrt(std::abs(this->dion(ind1,ind1) / vion1));
				this->dion(ind, ind) = vion1;
				this->dion(ind1, ind1) = vion1;
				for(int ir = 0; ir<this->mesh;ir++)
				{
					double avera = 1.0 / (2.0 * l + 1.0) * 
							( (l + 1.0) * sqrtDplus *
							this->beta(ind, ir) + 
							l * sqrtDminus *
							this->beta(ind1, ir) ) ;
					double delta = 1.0 / (2.0 * l + 1.0) * 
							( sqrtDplus *
							this->beta(ind, ir) - 
							sqrtDminus *
							this->beta(ind1, ir) ) ;
					this->beta(ind, ir) = (avera + l * delta * lambda_) ;
					this->beta(ind1, ir) = (avera - (l + 1) * delta * lambda_); 
				}
				nb++;
			}
		}

		for(int nb=0; nb<this->nwfc; nb++)
		{
			int l = this->lchi[nb];
			int ind=0, ind1=0;
			if(l!=0)
			{
				if(std::abs(this->jchi[nb] - this->lchi[nb] + 0.5) < 1e-6)
				{
					if(std::abs(this->jchi[nb+1]-this->lchi[nb+1]-0.5)>1e-6) 
					{error++; std::cout<<"warning_quit! error chi function 1 !"<<std::endl; return error;}
					ind = nb +1;
					ind1 = nb;
				}
				else
				{
					if(std::abs(this->jchi[nb+1]-this->lchi[nb+1]+0.5)>1e-6)
					{error++; std::cout<<"warning_quit! error chi function 2 !"<<std::endl; return error;}
					ind = nb;
					ind1 = nb +1;
				}
				//average chi
				for(int ir = 0; ir<this->mesh;ir++)
				{
					double avera = 0.5 * 
						( this->chi(ind,ir) + this->chi(ind1,ir) );
					double delta = 0.5 * 
						( this->chi(ind,ir) - this->chi(ind1,ir) );
					this->chi(ind, ir) = avera + delta * lambda_ ; 
					this->chi(ind1, ir) = avera - delta * lambda_ ; 
				}
				nb++;
			}
		}
		return error;
	}
}

// Peize Lin add for bsse 2021.04.07
void Pseudopot_upf::set_empty_element(void)
{
	this->zp = 0;
	for(int ir=0; ir<mesh; ++ir)
	{
		this->vloc[ir] = 0;
	}
	for(int i=0; i<nbeta; ++i)
	{
		for(int j=0; j<nbeta; ++j)
		{
			this->dion(i,j) = 0;
		}
	}
	for(int ir=0; ir<mesh; ++ir)
	{
		this->rho_at[ir] = 0;
	}
	return;
}
