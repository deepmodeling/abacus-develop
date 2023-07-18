#include "paw_element.h"
#include "module_base/tool_title.h"
#include "module_base/tool_quit.h"

void Paw_Element::init_paw_element(const double ecutwfc_in, const double cell_factor_in)
{
    this -> ecutwfc = ecutwfc_in;
    this -> cell_factor = cell_factor_in;
}

void Paw_Element::read_paw_xml(std::string filename)
{
    ModuleBase::TITLE("Paw_Element","read_paw_xml");

    std::string line;

    std::ifstream ifs(filename.c_str());
    if(ifs.fail())
    {
        ModuleBase::WARNING_QUIT("paw_element.cpp","xml file not found!");
    }

// ============================================================
// 1. symbol, Zat, core and valence electrons
// example : <atom symbol="H" Z="1.00" core="0.00" valence="1.00"/>
// ============================================================
    line = this->scan_file(ifs, "<atom");

    this->symbol = this->extract_string(line,"symbol=");
    this->Zat    = this->extract_double(line,"Z=");
    this->core   = this->extract_double(line,"core=");
    this->val    = this->extract_double(line,"valence=");

    this->reset_buffer(ifs);

// ============================================================
// 2. cutoff radius
// example : <paw_radius rc=" 0.9949503343"/>
// ============================================================
    line = this->scan_file(ifs, "<paw_radius");

    this->rcut   = this->extract_double(line,"rc=");

    this->reset_buffer(ifs);

// ============================================================
// 3. number of projector channels and corresponding l values
// example : 
// <valence_states>
//   <state n=" 1" l="0" f=" 1.0000000E+00" rc=" 0.9949503343" e="-2.3345876E-01" id=  "H1"/>
//   <state        l="0"                    rc=" 0.9949503343" e=" 6.0000000E+00" id=  "H2"/>
//   <state        l="1"                    rc=" 0.9949503343" e=" 1.2500000E+00" id=  "H3"/>
//  </valence_states>
// ============================================================
    nstates = this->count_nstates(ifs);

    this->reset_buffer(ifs);

    lstate.resize(nstates);
    lmax = 0;
    for(int istate = 0; istate < nstates; istate ++)
    {
        line = this->scan_file(ifs, "<state");

        this->lstate[istate] = this->extract_int(line,"l=");
        lmax = std::max(lmax, lstate[istate]);
    }

    this->nstates_to_mstates();

    this->reset_buffer(ifs);

// ============================================================
// 4. the radial grid
// example :
// <radial_grid eq="r=a*(exp(d*i)-1)" a=" 6.3033848776412630E-03" d=" 6.3033848776412630E-03" istart="0" iend=" 1499" id="log1">
//  <values>
//  ...
//  </values>
//  <derivatives>
//  ...
//  </derivatives>
// ============================================================
    line = this->scan_file(ifs, "<radial_grid");

    rstep = this->extract_double(line,"a=");
    lstep = this->extract_double(line,"d=");
    int istart = this->extract_int(line,"istart=");
    int iend   = this->extract_int(line,"iend=");

    nr = iend - istart + 1;

    rr.resize(nr);
    dr.resize(nr);

    line = this->scan_file(ifs, "<values>");

    for(int ir = 0; ir < nr; ir ++)
    {
        ifs >> rr[ir];
        //double tmp = rstep * (exp(lstep * double(ir)) - 1);
        //assert(std::abs(rr[ir] - tmp) < 1e-8);
    }

    line = this->scan_file(ifs, "<derivatives>");

    for(int ir = 0; ir < nr; ir ++)
    {
        ifs >> dr[ir];
    }

    this->get_nrcut();

// ============================================================
// 5. real-space projector functions
// example: <projector_function state=  "H1" grid="log1">
// ============================================================
    ptilde_r.resize(nstates);

    for(int istate = 0; istate < nstates; istate ++)
    {
        ptilde_r[istate].resize(nr);
        line = this->scan_file(ifs, "<projector_function");

        for(int ir = 0; ir < nr; ir ++)
        {
            // Note that in ABINIT, tproj stores r*p(r), not p(r) as in here
            ifs >> ptilde_r[istate][ir];
        }
    }

// ============================================================
// 6. reciprocal-space projector functions
// performing spherical bessel transformation
// ============================================================

    this -> transform_ptilde();
}

std::string Paw_Element::scan_file(std::ifstream &ifs, std::string pattern)
{
    std::string line;

    while (!ifs.eof())
    {
        getline(ifs,line);
        if (line.find(pattern) != std::string::npos)
        {
            //replace all quotation marks by space
            //to make it easier for later operations
            std::replace( line.begin(), line.end(), '"', ' ');

            return line;
        }
    }

    ModuleBase::WARNING_QUIT("Paw_Element::scan_file","pattern not found in xml file!");
    return 0;    
}

double Paw_Element::extract_double(std::string line, std::string key)
{
    int index = line.find(key);
    if (index != std::string::npos)
    {
        std::stringstream tmp;
        double value;
        tmp << line.substr(index+key.length());
        tmp >> value;

        return value;
    }
    else
    {
        ModuleBase::WARNING_QUIT("Paw_Element::extract_double","key not found in line!");
    }

    return 0;
}

std::string Paw_Element::extract_string(std::string line, std::string key)
{
    int index = line.find(key);
    if (index != std::string::npos)
    {
        std::stringstream tmp;
        std::string value;
        tmp << line.substr(index+key.length());
        tmp >> value;

        return value;
    }
    else
    {
        ModuleBase::WARNING_QUIT("Paw_Element::extract_double","key not found in line!");
    }

    return 0;
}

int Paw_Element::extract_int(std::string line, std::string key)
{
    int index = line.find(key);
    if (index != std::string::npos)
    {
        std::stringstream tmp;
        int value;
        tmp << line.substr(index+key.length());
        tmp >> value;

        return value;
    }
    else
    {
        ModuleBase::WARNING_QUIT("Paw_Element::extract_double","key not found in line!");
    }

    return 0;
}

int Paw_Element::count_nstates(std::ifstream &ifs)
{

    std::string line;

    while (!ifs.eof())
    {
        getline(ifs,line);
        if (line.find("<valence_states>") != std::string::npos)
        {
            break;
        }
    }

    nstates = 0;
    while (!ifs.eof())
    {
        getline(ifs,line);
        if (line.find("</valence_states>") != std::string::npos)
        {
            break;
        }
        nstates ++;
    }

    return nstates;

}

void Paw_Element::reset_buffer(std::ifstream &ifs)
{
    ifs.clear();
    ifs.seekg (0, std::ios::beg);
}

// s orbital : 0
// p orbital : -1,0,1
// d orbital : -2,-1,0,1,2
// etc.
void Paw_Element::nstates_to_mstates()
{
    mstates = 0;
    for(int istate = 0; istate < nstates; istate ++)
    {
        mstates += (2 * lstate[istate] + 1);
    }

    mstate.resize(mstates);
    im_to_istate.resize(mstates);

    int index = 0;
    for(int istate = 0; istate < nstates; istate ++)
    {
        int nm = 2 * lstate[istate] + 1;
        for(int im = 0; im < nm; im ++)
        {
            mstate[index] = im - lstate[istate];
            im_to_istate[index] = istate;
            index ++;
        }
    }

    assert(index == mstates);
}

void Paw_Element::get_nrcut()
{

    if(rcut > rr[nr-1])
    {
        ModuleBase::WARNING_QUIT("Paw_Element::get_nrcut","rcut > rr[nr-1], something wrong!");
    }

    // binary search
    int i,j,jm;

    i=0;
    j=nr-1;

    while (j-i>1)
    {
        jm=(i+j)/2;
        if(rcut < rr[jm])
        {
            j=jm;
        }
        else
        {
            i=jm;
        }
    }

    nrcut = j;
}