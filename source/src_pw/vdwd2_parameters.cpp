//==========================================================
// AUTHOR,	 Peize Lin
// DATE ,	 2021-03-09
//==========================================================

#include "vdwd2_parameters.h"

Vdwd2_Parameters::Vdwd2_Parameters()
{
	default_C6();
	default_R0();
	C6 = C6_default;
	R0 = R0_default;
}

void Vdwd2_Parameters::initial_parameters(const Input &input)
{
	this->flag_vdwd2 = true;
	this->scaling = std::stod(input.vdw_s6);
	this->damping = input.vdw_d;
	this->C6_input(input.vdw_C6_file, input.vdw_C6_unit);
	this->R0_input(input.vdw_R0_file, input.vdw_R0_unit);
	this->model = input.vdw_cutoff_type;
	if(input.vdw_cutoff_type=="radius")
	{
		if(input.vdw_radius_unit=="Bohr")
		{
			this->radius = std::stod(input.vdw_cutoff_radius);
		}
		else
		{
			this->radius = std::stod(input.vdw_cutoff_radius) * ModuleBase::BOHR_TO_A;
		}
	}
	else if(input.vdw_cutoff_type=="period")
	{
		this->period = input.vdw_cutoff_period;
	}
}

void Vdwd2_Parameters::initset(const UnitCell_pseudo &ucell)
{
	if(model=="radius")
	{
		period.x = 2*ceil(radius/ucell.lat0/sqrt(ucell.a1.norm2())) +1;
		period.y = 2*ceil(radius/ucell.lat0/sqrt(ucell.a2.norm2())) +1;
		period.z = 2*ceil(radius/ucell.lat0/sqrt(ucell.a3.norm2())) +1;		
	}
}

void Vdwd2_Parameters::C6_input(const std::string &file, const std::string &unit)
{
	C6 = C6_default;
	if( file != "default" )
	{
		std::ifstream ifs(file);
		if(!ifs)
			ModuleBase::WARNING_QUIT("Vdwd2::C6_input", "Can not find the file "+ModuleBase::GlobalFunc::TO_STRING(file));
		std::string element;
		double value;
		while( ifs >> element >> value )
			C6[element]=value;
		ifs.close();
	}
	for(auto &c6 : C6)
	{
		if( unit == "Jnm6/mol")
			c6.second *= 1e6/(ModuleBase::ELECTRONVOLT_SI*ModuleBase::NA)/pow(ModuleBase::BOHR_TO_A,6)/ModuleBase::Ry_to_eV;
		else if( unit == "eVA6")
			c6.second /= pow(ModuleBase::BOHR_TO_A,6)/ModuleBase::Ry_to_eV;
//		else if( unit == "RyBohr6");
		else
			ModuleBase::WARNING_QUIT("Input","vdwD2_C6_unit must be Jnm6/mol or eVA6");
	}
}

void Vdwd2_Parameters::R0_input(const std::string &file, const std::string &unit)
{
	R0 = R0_default;
	if( file != "default" )
	{
		std::ifstream ifs(file.c_str());
		if(!ifs)
			ModuleBase::WARNING_QUIT("Vdwd2::R0_input", "Can not find the file "+ModuleBase::GlobalFunc::TO_STRING(file));
		std::string element;
		double value;
		while( ifs >> element >> value )
			R0[element]=value;
		ifs.close();
	}
	for(auto &r0 : R0)
	{
		if( unit == "A")
			r0.second/= ModuleBase::BOHR_TO_A;
		else if( unit == "Bohr") ;
		else
			ModuleBase::WARNING_QUIT("Input","vdwD2_R0_unit must be A or Bohr");			
	}
}

void Vdwd2_Parameters::default_C6()
{
	C6_default =
	{
		{"H",	0.14     },
		{"He",	0.08     },
		{"Li",	1.61     },
		{"Be",	1.61     },
		{"B",	3.13     },
		{"C",	1.75     },
		{"N",	1.23     },
		{"O",	0.70     },
		{"F",	0.75     },
		{"Ne",	0.63     },
		{"Na",	5.71     },
		{"Mg",	5.71     },
		{"Al",	10.79    },
		{"Si",	9.23     },
		{"P",	7.84     },
		{"S",	5.57     },
		{"Cl",	5.07     },
		{"Ar",	4.61     },
		{"K",	10.8     },
		{"Ca",	10.8     },
		{"Sc",	10.8     },
		{"Ti",	10.8     },
		{"V",	10.8     },
		{"Cr",	10.8     },
		{"Mn",	10.8     },
		{"Fe",	10.8     },
		{"Co",	10.8     },
		{"Ni",	10.8     },
		{"Cu",	10.8     },
		{"Zn",	10.8     },
		{"Ga",	16.99    },
		{"Ge",	17.10    },
		{"As",	16.37    },
		{"Se",	12.64    },
		{"Br",	12.47    },
		{"Kr",	12.01    },
		{"Rb",	24.67    },
		{"Sr",	24.67    },
		{"Y",	24.67    },
		{"Zr",	24.67    },
		{"Nb",	24.67    },
		{"Mo",	24.67    },
		{"Tc",	24.67    },
		{"Ru",	24.67    },
		{"Rh",	24.67    },
		{"Pd",	24.67    },
		{"Ag",	24.67    },
		{"Cd",	24.67    },
		{"In",	37.32    },
		{"Sn",	38.71    },
		{"Sb",	38.44    },
		{"Te",	31.74    },
		{"I",	31.50    },
		{"Xe",	29.99    },
		{"Cs",	315.275  },
		{"Ba",	226.994  },
		{"La",	176.252  },
		{"Ce",	140.68   },
		{"Pr",	140.68   },
		{"Nd",	140.68   },
		{"Pm",	140.68   },
		{"Sm",	140.68   },
		{"Eu",	140.68   },
		{"Gd",	140.68   },
		{"Tb",	140.68   },
		{"Dy",	140.68   },
		{"Ho",	140.68   },
		{"Er",	140.68   },
		{"Tm",	140.68   },
		{"Yb",	140.68   },
		{"Lu",	140.68   },
		{"Hf",	105.112  },
		{"Ta",	81.24    },
		{"W",	81.24    },
		{"Re",	81.24    },
		{"Os",	81.24    }, 
		{"Ir",	81.24    },
		{"Pt",	81.24    },
		{"Au",	81.24    },
		{"Hg",	57.364   },
		{"Tl",	57.254   },
		{"Pb",	63.162   },
		{"Bi",	63.540   },
		{"Po",	55.283   },
		{"At",	57.171   },
		{"Rn",	56.64    }
	};	
}

void Vdwd2_Parameters::default_R0()
{
	R0_default =
	{
		{"H",	1.001  },
		{"He",	1.012  },
		{"Li",	0.825  },
		{"Be",	1.408  },
		{"B",	1.485  },
		{"C",	1.452  },
		{"N",	1.397  },
		{"O",	1.342  },
		{"F",	1.287  },
		{"Ne",	1.243  },
		{"Na",	1.144  },
		{"Mg",	1.364  },
		{"Al",	1.639  },
		{"Si",	1.716  },
		{"P",	1.705  },
		{"S",	1.683  },
		{"Cl",	1.639  },
		{"Ar",	1.595  },
		{"K",	1.485  },
		{"Ca",	1.474  },
		{"Sc",	1.562  },
		{"Ti",	1.562  },
		{"V",	1.562  },
		{"Cr",	1.562  },
		{"Mn",	1.562  },
		{"Fe",	1.562  },
		{"Co",	1.562  },
		{"Ni",	1.562  },
		{"Cu",	1.562  },
		{"Zn",	1.562  },
		{"Ga",	1.65   },
		{"Ge",	1.727  },
		{"As",	1.76   },
		{"Se",	1.771  },
		{"Br",	1.749  },
		{"Kr",	1.727  },
		{"Rb",	1.628  },
		{"Sr",	1.606  },
		{"Y",	1.639  },
		{"Zr",	1.639  },
		{"Nb",	1.639  },
		{"Mo",	1.639  },
		{"Tc",	1.639  },
		{"Ru",	1.639  },
		{"Rh",	1.639  },
		{"Pd",	1.639  },
		{"Ag",	1.639  },
		{"Cd",	1.639  },
		{"In",	1.672  },
		{"Sn",	1.804  },
		{"Sb",	1.881  },
		{"Te",	1.892  },
		{"I",	1.892  },
		{"Xe",	1.881  },
		{"Cs",	1.8018 },
		{"Ba",	1.7622 },
		{"La",	1.7204 },
		{"Ce",	1.7534 },
		{"Pr",	1.7534 },
		{"Nd",	1.7534 },
		{"Pm",	1.7534 },
		{"Sm",	1.7534 },
		{"Eu",	1.7534 },
		{"Gd",	1.7534 },
		{"Tb",	1.7534 },
		{"Dy",	1.7534 },
		{"Ho",	1.7534 },
		{"Er",	1.7534 },
		{"Tm",	1.7534 },
		{"Yb",	1.7534 },
		{"Lu",	1.7534 },
		{"Hf",	1.7875 },
		{"Ta",	1.7721 },
		{"W",	1.7721 },
		{"Re",	1.7721 },
		{"Os",	1.7721 },
		{"Ir",	1.7721 },
		{"Pt",	1.7721 },
		{"Au",	1.7721 },
		{"Hg",	1.7578 },
		{"Tl",	1.9855 },
		{"Pb",	1.9437 },
		{"Bi",	1.8975 },
		{"Po",	2.0053 },
		{"At",	1.991  },
		{"Rn",	1.9239 }
	};
}