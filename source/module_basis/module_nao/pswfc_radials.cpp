#include "module_basis/module_nao/pswfc_radials.h"
#include "module_base/tool_quit.h"
#include "module_base/math_integral.h"
#include <algorithm>
#include <cmath>

PswfcRadials& PswfcRadials::operator=(const PswfcRadials& rhs)
{
    RadialSet::operator=(rhs);
    return *this;
}

void PswfcRadials::build(const std::string& file, 
                         const int itype, 
                         const double screening_coeff,
                         const double conv_thr,
                         std::ofstream* ptr_log,
                         const int rank)
{
    // deallocates all arrays and reset variables (excluding sbt_)
    cleanup();

    std::ifstream ifs;
    bool is_open = false;

    if (rank == 0)
    {
        ifs.open(file);
        is_open = ifs.is_open();
    }

#ifdef __MPI
    Parallel_Common::bcast_bool(is_open);
#endif

    if (!is_open)
    {
        ModuleBase::WARNING_QUIT("AtomicRadials::read", "Couldn't open pseudopotential file: " + file);
    }

    itype_ = itype;
    read_upf_pswfc(ifs, screening_coeff, conv_thr, ptr_log, rank);
    set_rcut_max();
    
    if (rank == 0)
    {
        ifs.close();
    }
}

bool PswfcRadials::startswith(std::string word, std::string pattern)
{
    if(word.size() < pattern.size()) return false;
    int score = 1;
    for(int ic = 0; ic < pattern.size(); ic++)
    {
        if(word[ic] != pattern[ic])
        {
            score *= 0;
        }
        else
        {
            score *= 1;
        }
    }
    return bool(score);
}

std::string PswfcRadials::steal_from_quotes(std::string word)
{
    // first make sure there are even number of quotes in this word
    int num_quote = 0;
    for(auto letter: word)
    {
        if(letter == '\"') num_quote += 1;
    }
    assert(num_quote % 2 == 0);
    // then steal from quotes
    std::string result;
    size_t _left = word.find_first_of("\"");
    size_t _right = word.find_last_of("\"");
    result = word.substr(_left+1, _right-_left-1);
    // then remove all spaces ahead
    while(result[0] == ' ')
    {
        result.erase(0, 1);
    }
    return result;
}

std::string PswfcRadials::steal_from_quotes(std::ifstream& ifs, std::string word)
{
    // concatenate all words until the second quote, no matter how many lines and spaces between
    std::string concatenated = word.substr(
        word.find_first_of("\"")+1, 
        word.size()-word.find_first_of("\"")-1
        );
    int num_quote = 1;
    while(num_quote < 2)
    {
        std::string line;
        ifs >> line;
        for(auto letter: line)
        {
            if(letter == '\"') num_quote += 1;
            if(num_quote == 2) break;
            concatenated += letter;
        }
    }
    // then remove all spaces ahead
    while(concatenated[0] == ' ')
    {
        concatenated.erase(0, 1);
    }
    return concatenated;
}

std::string PswfcRadials::read_keyword_value(std::ifstream& ifs, std::string word)
{
    // count the number of quotes, only 1 or 2 cases are considered for pseudopotential reading
    int num_quote = 0;
    for(auto letter: word)
    {
        if(letter == '\"') num_quote += 1;
    }
    assert(num_quote == 1 || num_quote == 2);
    if(num_quote == 1) return steal_from_quotes(ifs, word);
    else return steal_from_quotes(word);
}

double PswfcRadials::radial_norm(const std::vector<double> rgrid,
                                 const std::vector<double> rvalue)
{
    std::vector<double> integrand(rvalue.size());
    for(int ir = 0; ir != rvalue.size(); ++ir)
    {
        integrand[ir] = rvalue[ir] * rvalue[ir] * rgrid[ir] * rgrid[ir];
    }
    double dr = rgrid[1] - rgrid[0];
    double norm = ModuleBase::Integral::simpson(rvalue.size(), integrand.data(), dr);
    norm = sqrt(norm);
    return norm;
}

double PswfcRadials::cut_to_convergence(const std::vector<double>& rgrid,
                                        std::vector<double>& rvalue, 
                                        const double& conv_thr)
{
    double norm = 0.0;
    int ir_ = 0;
    int ir_min_ = 0;
    int delta_ir = 5; // stepsize for radius cutoff searching, in Bohr

    int ir_max_ = rgrid.size() - 1;
    
    double dr = rgrid[1] - rgrid[0]; // radial function realspace grid stepsize, in Bohr
    printf("Norm of pseudowavefunction before cutoff: %6.4e\n", radial_norm(rgrid, rvalue));
    int istep = 1;
    double delta_norm = 1.0;
    printf("Searching for the cutoff radius for pseudowavefunction, conv_thr = %6.4e\n", conv_thr);
    printf("%10s%12s%14s%18s", "Step Nr.", "Rmax (a.u.)", "Norm", "Delta Norm\n");
    while((std::fabs(delta_norm) > conv_thr)&&(ir_ <= ir_max_))
    {
        ir_ = std::min(ir_ + delta_ir, ir_max_); // update ir_, but be careful not to exceed ir_max_
        delta_norm = norm;
        std::vector<double> rgrid_slice = std::vector<double>(rgrid.begin() + ir_min_, rgrid.begin() + ir_ + 1);
        std::vector<double> rvalue_slice = std::vector<double>(rvalue.begin() + ir_min_, rvalue.begin() + ir_ + 1);
        norm = radial_norm(rgrid_slice, rvalue_slice);
        delta_norm = norm - delta_norm;
        if(istep == 1) printf("%10d%12.2f%14.10f%18.10e\n", istep, rgrid[ir_], norm, delta_norm);
        ++istep;
    }
    printf("...\n");
    printf("%10d%12.2f%14.10f%18.10e\n", istep, rgrid[ir_], norm, delta_norm);

    rvalue = std::vector<double>(rvalue.begin() + ir_min_, rvalue.begin() + ir_ + 1);
    return rgrid[ir_max_];
}

void PswfcRadials::smooth(std::vector<double>& rgrid,
                          std::vector<double>& rvalue,
                          const double sigma)
{
    double prefactor = 1.0 / sqrt(2.0 * M_PI) / sigma;
    double rmax = rgrid.back();
    for(int ir = 0; ir != rgrid.size(); ++ir)
    {
        double delta_r = rgrid[ir] - rmax;
        double smooth = prefactor * exp(-delta_r * delta_r / 2.0 / sigma / sigma);
        rvalue[ir] *= (1 - smooth);
    }
}

std::vector<double> PswfcRadials::pswfc_prepossess(std::map<std::pair<int, int>, std::vector<double>>& lzeta_rvalues,
                                                   const double conv_thr)
{
    double nmax = 0.0;
    for(auto it = lzeta_rvalues.begin(); it != lzeta_rvalues.end(); it++)
    {
        int l = it->first.first;
        int iz = it->first.second;
        std::vector<double> rvalue = it->second;
        std::vector<double> rgrid = std::vector<double>(rvalue.size(), 0.0);
        for(int ir = 0; ir < rvalue.size(); ir++)
        {
            rgrid[ir] = ir * 0.01;
        }
        double rcut_i = cut_to_convergence(rgrid, rvalue, conv_thr);
        if(rvalue.size() > nmax) nmax = rvalue.size();
        lzeta_rvalues[it->first] = rvalue; // stores in map
    }
    // generate rgrid
    std::vector<double> rgrid = std::vector<double>(nmax, 0.0);
    for(int ir = 0; ir < nmax; ir++)
    {
        rgrid[ir] = ir * 0.01;
    }
    // zero padding on rvalue
    for(auto it = lzeta_rvalues.begin(); it != lzeta_rvalues.end(); it++)
    {
        int l = it->first.first;
        int iz = it->first.second;
        std::vector<double> rvalue = it->second;
        std::vector<double> rvalue_padded = std::vector<double>(nmax, 0.0);
        for(int ir = 0; ir < rvalue.size(); ir++)
        {
            rvalue_padded[ir] = rvalue[ir];
        }
        smooth(rgrid, rvalue_padded, 0.1); // smooth the radial function to avoid high frequency noise in FFT-spherical bessel transform
        lzeta_rvalues[it->first] = rvalue_padded; // stores in map
    }

    return rgrid;
}

void PswfcRadials::read_upf_pswfc(std::ifstream& ifs, 
                                  const double screening_coeff, 
                                  const double conv_thr,
                                  std::ofstream* ptr_log, 
                                  const int rank)
{
    int ngrid = 0;
    int nzeta = 0;
    double dr = 0.01; // in most cases, this is correct

    if(rank == 0)
    {
        std::map<std::pair<int, int>, std::vector<double>> result;

        std::string line = "";
        // read element
        while(!startswith(line, "element=")&&!ifs.eof()) ifs >> line;
        symbol_ = read_keyword_value(ifs, line);
        // read lmax
        while(!startswith(line, "l_max=")&&!ifs.eof()) ifs >> line;
        lmax_ = std::stoi(read_keyword_value(ifs, line));
        // read ngrid
        while(!startswith(line, "mesh_size=")&&!ifs.eof()) ifs >> line;
        ngrid = std::stoi(read_keyword_value(ifs, line));
        // read nzeta
        while(!startswith(line, "number_of_wfc=")&&!ifs.eof()) ifs >> line;
        nzeta = std::stoi(read_keyword_value(ifs, line));
        nchi_ = nzeta;
        // read contents of pseudowavefunction
        while(line != "<PP_PSWFC>") ifs >> line;
        nzeta_ = new int[lmax_ + 1];
        for(int il = 0; il < lmax_ + 1; il++) nzeta_[il] = 0;
        for(int iz = 0; iz < nzeta; iz++) // read chi-by-chi
        {
            // find the next <PP_CHI.> tag
            while(!startswith(line, "<PP_CHI.")&&!ifs.eof()) ifs >> line;
            // read l
            while(!startswith(line, "l=")&&!ifs.eof()) ifs >> line;
            int l = std::stoi(read_keyword_value(ifs, line));
            nzeta_[l] += 1;
            // to data
            while(line != ">"&&!ifs.eof()) ifs >> line;
            // before read data, first create container to store
            std::vector<double> rvalue = std::vector<double>(ngrid, 0.0);
            for(int ir=0; ir < ngrid; ir++)
            {
                ifs >> line;
                double screening = std::exp(-screening_coeff * ir * dr);
                rvalue[ir] = std::stod(line) * screening;
            }
            result[std::make_pair(l, nzeta_[l] - 1)] = rvalue;
            ifs >> line;
            assert(startswith(line, "</PP_CHI."));
        }
        
        if(result.size() == 0)
        {
            ModuleBase::WARNING_QUIT("AtomicRadials::read", "pseudowavefunction information is absent in pseudopotential.");
        }

        nzeta_max_ = *std::max_element(nzeta_, nzeta_ + lmax_ + 1);
        indexing(); // build index_map_

        // cut rvalue to convergence and generate rgrid
        std::vector<double> rgrid = pswfc_prepossess(result, conv_thr);
        // refresh ngird value
        ngrid = rgrid.size();

        chi_ = new NumericalRadial[nchi_];

        for(auto it = result.begin(); it != result.end(); it++)
        {
            int l = it->first.first;
            int iz = it->first.second;
            chi_[index(l, iz)].build(l, true, ngrid, rgrid.data(), it->second.data(), 0, iz, symbol_, itype_, false);       
            if(std::fabs(screening_coeff - 0.0) > 1e-6) // PHYSICAL REVIEW B 78, 245112 2008
            {
                chi_[index(l, iz)].normalize();
            }
        }
    } // rank 0 does almost everything, then broadcast one-by-one
#ifdef __MPI
    Parallel_Common::bcast_string(symbol_);
    Parallel_Common::bcast_int(lmax_);

    Parallel_Common::bcast_int(nchi_);
    Parallel_Common::bcast_int(nzeta_max_);

    Parallel_Common::bcast_int(ngrid);
    // Parallel_Common::bcast_double(dr); // we dont need to broadcast dr again because it is fixed to 0.01
#endif
    if(rank != 0)
    {
        nzeta_ = new int[lmax_ + 1];
        index_map_ = new int[(lmax_ + 1) * nzeta_max_];
    }
#ifdef __MPI
    Parallel_Common::bcast_int(nzeta_, lmax_ + 1);
    Parallel_Common::bcast_int(index_map_, (lmax_ + 1) * nzeta_max_);
#endif

    double* rvalue = new double[ngrid];
    double* rgrid = new double[ngrid];
    if(rank != 0)
    {
        for(int ir = 0; ir != ngrid; ++ir)
        {
            rgrid[ir] = ir * dr;
        }
        chi_ = new NumericalRadial[nchi_];
    }
}
