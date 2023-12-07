#include "module_basis/module_nao/pswfc_radials.h"
#include "module_base/tool_quit.h"
#include <algorithm>

void PswfcRadials::build(const std::string& file, const int itype, std::ofstream* ptr_log, const int rank)
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
    if (!is_open)
    {
        ModuleBase::WARNING_QUIT("AtomicRadials::read", "Couldn't open pseudopotential file: " + file);
    }

    itype_ = itype;
    read_upf_pswfc(ifs, ptr_log, rank);
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

void PswfcRadials::read_upf_pswfc(std::ifstream& ifs, std::ofstream* ptr_log, const int rank)
{
    int ngrid = 0;
    int nzeta = 0;
    double dr = 0.01; // in most cases, this is correct


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
            rvalue[ir] = std::stod(line);
        }
        result[std::make_pair(l, nzeta_[l] - 1)] = rvalue;
        ifs >> line;
        assert(startswith(line, "</PP_CHI."));
    }
    
    nzeta_max_ = *std::max_element(nzeta_, nzeta_ + lmax_ + 1);
    indexing(); // build index_map_

    std::vector<double> rgrid = std::vector<double>(ngrid, 0.0);
    for(int ir = 0; ir < ngrid; ir++)
    {
        rgrid[ir] = ir * dr;
    }

    chi_ = new NumericalRadial[nchi_];

    for(auto it = result.begin(); it != result.end(); it++)
    {
        int l = it->first.first;
        int iz = it->first.second;
        chi_[index(l, iz)].build(l, true, ngrid, rgrid.data(), it->second.data(), 0, iz, symbol_, itype_, false);       
        //chi_[index(l, iz)].normalize();
    }
}
