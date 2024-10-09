#ifndef ATOM_IN_H
#define ATOM_IN_H

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

class atom_in
{
  public:
    atom_in(){};
    ~atom_in(){};
    std::map<std::string, int> atom_Z
        = {{"H", 1},    {"He", 2},   {"Li", 3},   {"Be", 4},   {"B", 5},    {"C", 6},    {"N", 7},   {"O", 8},   {"F", 9},
           {"Ne", 10},  {"Na", 11},  {"Mg", 12},  {"Al", 13},  {"Si", 14},  {"P", 15},   {"S", 16},  {"Cl", 17}, {"Ar", 18},
           {"K", 19},   {"Ca", 20},  {"Sc", 21},  {"Ti", 22},  {"V", 23},   {"Cr", 24},  {"Mn", 25}, {"Fe", 26}, {"Co", 27},
           {"Ni", 28},  {"Cu", 29},  {"Zn", 30},  {"Ga", 31},  {"Ge", 32},  {"As", 33},  {"Se", 34}, {"Br", 35}, {"Kr", 36},
           {"Rb", 37},  {"Sr", 38},  {"Y", 39},   {"Zr", 40},  {"Nb", 41},  {"Mo", 42},  {"Tc", 43}, {"Ru", 44}, {"Rh", 45},
           {"Pd", 46},  {"Ag", 47},  {"Cd", 48},  {"In", 49},  {"Sn", 50},  {"Sb", 51},  {"Te", 52}, {"I", 53},  {"Xe", 54},
           {"Cs", 55},  {"Ba", 56},  {"La", 57},  {"Ce", 58},  {"Pr", 59},  {"Nd", 60},  {"Pm", 61}, {"Sm", 62}, {"Eu", 63},
           {"Gd", 64},  {"Tb", 65},  {"Dy", 66},  {"Ho", 67},  {"Er", 68},  {"Tm", 69},  {"Yb", 70}, {"Lu", 71}, {"Hf", 72},
           {"Ta", 73},  {"W", 74},   {"Re", 75},  {"Os", 76},  {"Ir", 77},  {"Pt", 78},  {"Au", 79}, {"Hg", 80}, {"Tl", 81},
           {"Pb", 82},  {"Bi", 83},  {"Po", 84},  {"At", 85},  {"Rn", 86},  
           {"Fr", 87},  {"Ra", 88},  {"Ac", 89},  {"Th", 90},  {"Pa", 91},  
           {"U", 92},   {"Np", 93},  {"Pu", 94},  {"Am", 95},  {"Cm", 96},
           {"Bk", 97},  {"Cf", 98},  {"Es", 99},  {"Fm", 100}, {"Md", 101}, {"No", 102},
           {"Lr", 103}, {"Rf", 104}, {"Db", 105}, {"Sg", 106}, {"Bh", 107}, {"Hs", 108},
           {"Mt", 109}, {"Ds", 110}, {"Rg", 111}, {"Cn", 112}, {"Nh", 113}, {"Fl", 114},
           {"Mc", 115}, {"Lv", 116}, {"Ts", 117}, {"Og", 118}};

    std::map<std::string, double> atom_RCS
        = {{"H", 0.603774}, {"He", 1.75472}, {"Li", 2.32075}, {"Be", 1.69811}, {"B", 1.54717},  {"C", 1.45283},
           {"N", 1.41509},  {"O", 1.37736},  {"F", 1.35849},  {"Ne", 1.33962}, {"Na", 2.90566}, {"Mg", 2.56604},
           {"Al", 2.22642}, {"Si", 2.09434}, {"P", 2},        {"S", 1.92453},  {"Cl", 1.86792}, {"Ar", 1.84906},
           {"K", 3.83019},  {"Ca", 3.28302}, {"Sc", 2.71698}, {"Ti", 2.49057}, {"V", 2.30189},  {"Cr", 2.22642},
           {"Mn", 2.20755}, {"Fe", 2.20755}, {"Co", 2.18868}, {"Ni", 2.16981}, {"Cu", 2.20755}, {"Zn", 2.35849},
           {"Ga", 2.37736}, {"Ge", 2.30189}, {"As", 2.26415}, {"Se", 2.18868}, {"Br", 2.15094}, {"Kr", 2.11321},
           {"Rb", 4.07547}, {"Sr", 3.60377}, {"Y", 3.0566},   {"Zr", 2.73585}, {"Nb", 2.5283},  {"Mo", 2.45283},
           {"Tc", 2.39623}, {"Ru", 2.35849}, {"Rh", 2.35849}, {"Pd", 2.41509}, {"Ag", 2.5283},  {"Cd", 2.79245},
           {"In", 2.71698}, {"Sn", 2.66038}, {"Sb", 2.64151}, {"Te", 2.56604}, {"I", 2.50943},  {"Xe", 2.4717},
           {"Cs", 4.43396}, {"Ba", 3.73585}, {"La", 3.18868}, {"Ce", 3.11321}, {"Pr", 3.11321}, {"Nd", 3.09434},
           {"Pm", 3.07547}, {"Sm", 3.0566},  {"Eu", 3.49057}, {"Gd", 3.03774}, {"Tb", 3},       {"Dy", 3},
           {"Ho", 2.98113}, {"Er", 2.96226}, {"Tm", 2.9434},  {"Yb", 3.28302}, {"Lu", 2.9434},  {"Hf", 2.71698},
           {"Ta", 2.5283},  {"W", 2.45283},  {"Re", 2.41509}, {"Os", 2.37736}, {"Ir", 2.39623}, {"Pt", 2.45283},
           {"Au", 2.5283},  {"Hg", 2.81132}, {"Tl", 2.79245}, {"Pb", 2.77358}, {"Bi", 2.75472}, {"Po", 2.75472},
           {"At", 2.73585}, {"Rn", 2.69811} };

    std::map<std::string, std::string> atom_symbol
        = { {"H", "Hydrogen"},      {"He", "Helium"},        {"Li", "Lithium"},     {"Be", "Beryllium"},   {"B", "Boron"},         {"C", "Carbon"},
            {"N", "Nitrogen"},      {"O", "Oxygen"},         {"F", "Fluorine"},     {"Ne", "Neon"},        {"Na", "Sodium"},       {"Mg", "Magnesium"},
            {"Al", "Aluminum"},     {"Si", "Silicon"},       {"P", "Phosphorus"},   {"S", "Sulfur"},       {"Cl", "Chlorine"},     {"Ar", "Argon"},
            {"K", "Potassium"},     {"Ca", "Calcium"},       {"Sc", "Scandium"},    {"Ti", "Titanium"},    {"V", "Vanadium"},      {"Cr", "Chromium"},
            {"Mn", "Manganese"},    {"Fe", "Iron"},          {"Co", "Cobalt"},      {"Ni", "Nickel"},      {"Cu", "Copper"},       {"Zn", "Zinc"},
            {"Ga", "Gallium"},      {"Ge", "Germanium"},     {"As", "Arsenic"},     {"Se", "Selenium"},    {"Br", "Bromine"},      {"Kr", "Krypton"},
            {"Rb", "Rubidium"},     {"Sr", "Strontium"},     {"Y", "Yttrium"},      {"Zr", "Zirconium"},   {"Nb", "Niobium"},      {"Mo", "Molybdenum"},
            {"Tc", "Technetium"},   {"Ru", "Ruthenium"},     {"Rh", "Rhodium"},     {"Pd", "Palladium"},   {"Ag", "Silver"},       {"Cd", "Cadmium"},
            {"In", "Indium"},       {"Sn", "Tin"},           {"Sb", "Antimony"},    {"Te", "Tellurium"},   {"I", "Iodine"},        {"Xe", "Xenon"},
            {"Cs", "Cesium"},       {"Ba", "Barium"},        {"La", "Lanthanum"},   {"Ce", "Cerium"},      {"Pr", "Praseodymium"}, {"Nd", "Neodymium"},
            {"Pm", "Promethium"},   {"Sm", "Samarium"},      {"Eu", "Europium"},    {"Gd", "Gadolinium"},  {"Tb", "Terbium"},      {"Dy", "Dysprosium"},
            {"Ho", "Holmium"},      {"Er", "Erbium"},        {"Tm", "Thulium"},     {"Yb", "Ytterbium"},   {"Lu", "Lutetium"},     {"Hf", "Hafnium"},
            {"Ta", "Tantalum"},     {"W", "Tungsten"},       {"Re", "Rhenium"},     {"Os", "Osmium"},      {"Ir", "Iridium"},      {"Pt", "Platinum"},
            {"Au", "Gold"},         {"Hg", "Mercury"},       {"Tl", "Thallium"},    {"Pb", "Lead"},        {"Bi", "Bismuth"},      {"Po", "Polonium"},
            {"At", "Astatine"},     {"Rn", "Radon"},         {"Fr", "Francium"},    {"Ra", "Radium"},      {"Ac", "Actinium"},     {"Th", "Thorium"},
            {"Pa", "Protactinium"}, {"U", "Uranium"},        {"Np", "Neptunium"},   {"Pu", "Plutonium"},   {"Am", "Americium"},    {"Cm", "Curium"},
            {"Bk", "Berkelium"},    {"Cf", "Californium"},   {"Es", "Einsteinium"}, {"Fm", "Fermium"},     {"Md", "Mendelevium"},  {"No", "Nobelium"},
            {"Lr", "Lawrencium"},   {"Rf", "Rutherfordium"}, {"Db", "Dubnium"},     {"Sg", "Seaborgium"},  {"Bh", "Bohrium"},      {"Hs", "Hassium"},
            {"Mt", "Meitnerium"},   {"Ds", "Darmstadtium"},  {"Rg", "Roentgenium"}, {"Cn", "Copernicium"}, {"Nh", "Nihonium"},     {"Fl", "Flerovium"},
            {"Mc", "Moscovium"},    {"Lv", "Livermorium"},   {"Ts", "Tennessine"},  {"Og", "Oganesson"}};
  
    std::map<std::string, int> symbol_Z
        = { {"Hydrogen", 1},      {"Helium", 2},          {"Lithium", 3},       {"Beryllium", 4},     {"Boron", 5},         {"Carbon", 6},
            {"Nitrogen", 7},      {"Oxygen", 8},          {"Fluorine", 9},      {"Neon", 10},         {"Sodium", 11},       {"Magnesium", 12},
            {"Aluminum", 13},     {"Silicon", 14},        {"Phosphorus", 15},   {"Sulfur", 16},       {"Chlorine", 17},     {"Argon", 18},
            {"Potassium", 19},    {"Calcium", 20},        {"Scandium", 21},     {"Titanium", 22},     {"Vanadium", 23},     {"Chromium", 24},
            {"Manganese", 25},    {"Iron", 26},           {"Cobalt", 27},       {"Nickel", 28},       {"Copper", 29},       {"Zinc", 30},
            {"Gallium", 31},      {"Germanium", 32},      {"Arsenic", 33},      {"Selenium", 34},     {"Bromine", 35},      {"Krypton", 36},
            {"Rubidium", 37},     {"Strontium", 38},      {"Yttrium", 39},      {"Zirconium", 40},    {"Niobium", 41},      {"Molybdenum", 42},
            {"Technetium", 43},   {"Ruthenium", 44},      {"Rhodium", 45},      {"Palladium", 46},    {"Silver", 47},       {"Cadmium", 48},
            {"Indium", 49},       {"Tin", 50},            {"Antimony", 51},     {"Tellurium", 52},    {"Iodine", 53},       {"Xenon", 54},
            {"Cesium", 55},       {"Barium", 56},         {"Lanthanum", 57},    {"Cerium", 58},       {"Praseodymium", 59}, {"Neodymium", 60},
            {"Promethium", 61},   {"Samarium", 62},       {"Europium", 63},     {"Gadolinium", 64},   {"Terbium", 65},      {"Dysprosium", 66},
            {"Holmium", 67},      {"Erbium", 68},         {"Thulium", 69},      {"Ytterbium", 70},    {"Lutetium", 71},     {"Hafnium", 72},
            {"Tantalum", 73},     {"Tungsten", 74},       {"Rhenium", 75},      {"Osmium", 76},       {"Iridium", 77},      {"Platinum", 78},
            {"Gold", 79},         {"Mercury", 80},        {"Thallium", 81},     {"Lead", 82},         {"Bismuth", 83},      {"Polonium", 84},
            {"Astatine", 85},     {"Radon", 86},          {"Francium", 87},     {"Radium", 88},       {"Actinium", 89},     {"Thorium", 90},
            {"Protactinium", 91}, {"Uranium", 92},        {"Neptunium", 93},    {"Plutonium", 94},    {"Americium", 95},    {"Curium", 96},
            {"Berkelium", 97},    {"Californium", 98},    {"Einsteinium", 99},  {"Fermium", 100},     {"Mendelevium", 101}, {"Nobelium", 102},
            {"Lawrencium", 103},  {"Rutherfordium", 104}, {"Dubnium", 105},     {"Seaborgium", 106},  {"Bohrium", 107},     {"Hassium", 108},
            {"Meitnerium", 109},  {"Darmstadtium", 110},  {"Roentgenium", 111}, {"Copernicium", 112}, {"Nihonium", 113},    {"Flerovium", 114},
            {"Moscovium", 115},   {"Livermorium", 116},   {"Tennessine", 117},  {"Oganesson", 118}
          };             
    std::map<std::string, int> principle_quantum_number
        = {
            {"H", 1},   {"He", 1},   {"Li", 2},   {"Be", 2},   {"B", 2},    {"C", 2},    {"N", 2},   {"O", 2},   {"F", 2},
            {"Ne", 1},   {"Na", 2},   {"Mg", 3},   {"Al", 3},   {"Si", 3},   {"P", 3},    {"S", 3},   {"Cl", 3},  {"Ar", 2},
            {"K", 3},    {"Ca", 4},   {"Sc", 4},   {"Ti", 4},   {"V", 4},    {"Cr", 4},   {"Mn", 4},  {"Fe", 4},  {"Co", 4},
            {"Ni", 4},   {"Cu", 4},   {"Zn", 4},   {"Ga", 4},   {"Ge", 4},   {"As", 4},   {"Se", 4},  {"Br", 4},  {"Kr", 3},
            {"Rb", 4},   {"Sr", 5},   {"Y", 5},    {"Zr", 5},   {"Nb", 5},   {"Mo", 5},   {"Tc", 5},  {"Ru", 5},  {"Rh", 5},
            {"Pd", 5},   {"Ag", 5},   {"Cd", 5},   {"In", 5},   {"Sn", 5},   {"Sb", 5},   {"Te", 5},  {"I", 5},   {"Xe", 4},
            {"Cs", 5},   {"Ba", 6},   {"La", 6},   {"Ce", 6},   {"Pr", 6},   {"Nd", 6},   {"Pm", 6},  {"Sm", 6},  {"Eu", 6},
            {"Gd", 6},   {"Tb", 6},   {"Dy", 6},   {"Ho", 6},   {"Er", 6},   {"Tm", 6},   {"Yb", 6},  {"Lu", 6},  {"Hf", 6},
            {"Ta", 6},   {"W", 6},    {"Re", 6},   {"Os", 6},   {"Ir", 6},   {"Pt", 6},   {"Au", 6},  {"Hg", 6},  {"Tl", 6},
            {"Pb", 6},   {"Bi", 6},   {"Po", 6},   {"At", 6},   {"Rn", 6},   {"Fr", 7},   {"Ra", 7},  {"Ac", 7},  {"Th", 7},
            {"Pa", 7},   {"U", 7},    {"Np", 7},   {"Pu", 7},   {"Am", 7},   {"Cm", 7},   {"Bk", 7},  {"Cf", 7},  {"Es", 7},
            {"Fm", 7},   {"Md", 7},   {"No", 7},   {"Lr", 7},   {"Rf", 7},   {"Db", 7},   {"Sg", 7},  {"Bh", 7},  {"Hs", 7},
            {"Mt", 7},   {"Ds", 7},   {"Rg", 7},   {"Cn", 7},   {"Nh", 7},   {"Fl", 7},   {"Mc", 7},  {"Lv", 7},  {"Ts", 7},
            {"Og", 7}
        };
    /// @brief ground state electron configuration, sequence of orbitals in key is in accord with the sequence of n then l
    /// @note 1s2s2p3s3p 4s 3d4p5s 4d5p6s 4f5d6p7s 5f6d7p, from NIST periodic table, 
    /// @details see: https://www.nist.gov/system/files/documents/2019/12/10/nist_periodictable_july2019_crop.pdf
    std::map<std::string, std::vector<int>> groundstate_electronconfiguration
        = {
            //      1s
            {"H", {1}}, {"He", {2}}, // 1st period
            //      1s 2s
            {"Li", {2, 1}}, {"Be", {2, 2}}, // 2nd period
            //      1s 2s 2p
            {"B",  {2, 2, 1}}, {"C", {2, 2, 2}}, {"N", {2, 2, 3}}, {"O", {2, 2, 4}}, {"F", {2, 2, 5}}, {"Ne", {2, 2, 6}}, // 2nd period
            //      1s 2s 2p 3s           1s 2s 2p 3s
            {"Na", {2, 2, 6, 1}}, {"Mg", {2, 2, 6, 2}}, // 3rd period
            //      1s 2s 2p 3s 3p           1s 2s 2p 3s 3p           1s 2s 2p 3s 3p
            {"Al", {2, 2, 6, 2, 1}}, {"Si", {2, 2, 6, 2, 2}}, {"P",  {2, 2, 6, 2, 3}}, // 3rd period
            {"S",  {2, 2, 6, 2, 4}}, {"Cl", {2, 2, 6, 2, 5}}, {"Ar", {2, 2, 6, 2, 6}}, // 3rd period
            //      1s 2s 2p 3s 3p 3d 4s           1s 2s 2p 3s 3p 3d 4s            1s 2s 2p 3s 3p 3d 4s
            {"K",  {2, 2, 6, 2, 6, 0, 1}}, {"Ca", {2, 2, 6, 2, 6, 0, 2}},  {"Sc", {2, 2, 6, 2, 6, 1, 2}}, // 4th period
            {"Ti", {2, 2, 6, 2, 6, 2, 2}}, {"V",  {2, 2, 6, 2, 6, 3, 2}},  {"Cr", {2, 2, 6, 2, 6, 4, 2}}, // 4th period
            {"Mn", {2, 2, 6, 2, 6, 5, 2}}, {"Fe", {2, 2, 6, 2, 6, 6, 2}},  {"Co", {2, 2, 6, 2, 6, 7, 2}}, // 4th period
            {"Ni", {2, 2, 6, 2, 6, 8, 2}}, {"Cu", {2, 2, 6, 2, 6, 10, 1}}, {"Zn", {2, 2, 6, 2, 6, 10, 2}}, // 4th period
            //      1s 2s 2p 3s 3p 3d  4s 4p           1s 2s 2p 3s 3p 3d  4s 4p           1s 2s 2p 3s 3p 3d  4s 4p
            {"Ga", {2, 2, 6, 2, 6, 10, 2, 1}}, {"Ge", {2, 2, 6, 2, 6, 10, 2, 2}}, {"As", {2, 2, 6, 2, 6, 10, 2, 3}}, // 4th period
            {"Se", {2, 2, 6, 2, 6, 10, 2, 4}}, {"Br", {2, 2, 6, 2, 6, 10, 2, 5}}, {"Kr", {2, 2, 6, 2, 6, 10, 2, 6}}, // 4th period
            //      1s 2s 2p 3s 3p 3d  4s 4p 4d 4f 5s            1s 2s 2p 3s 3p 3d  4s 4p 4d 4f 5s            1s 2s 2p 3s 3p 3d  4s 4p 4d 4f 5s
            {"Rb", {2, 2, 6, 2, 6, 10, 2, 6, 0, 0, 1}},  {"Sr", {2, 2, 6, 2, 6, 10, 2, 6, 0, 0, 2}},  {"Y",  {2, 2, 6, 2, 6, 10, 2, 6, 1, 0, 2}}, // 5th period
            {"Zr", {2, 2, 6, 2, 6, 10, 2, 6, 2, 0, 2}},  {"Nb", {2, 2, 6, 2, 6, 10, 2, 6, 4, 0, 1}},  {"Mo", {2, 2, 6, 2, 6, 10, 2, 6, 5, 0, 1}}, // 5th period
            {"Tc", {2, 2, 6, 2, 6, 10, 2, 6, 5, 0, 2}},  {"Ru", {2, 2, 6, 2, 6, 10, 2, 6, 7, 0, 1}},  {"Rh", {2, 2, 6, 2, 6, 10, 2, 6, 8, 0, 1}}, // 5th period
            {"Pd", {2, 2, 6, 2, 6, 10, 2, 6, 10, 0, 0}}, {"Ag", {2, 2, 6, 2, 6, 10, 2, 6, 10, 0, 1}}, {"Cd", {2, 2, 6, 2, 6, 10, 2, 6, 10, 0, 2}}, // 5th period
            //      1s 2s 2p 3s 3p 3d  4s 4p 4d  4f 5s 5p           1s 2s 2p 3s 3p 3d  4s 4p 4d  4f 5s 5p           1s 2s 2p 3s 3p 3d  4s 4p 4d  4f 5s 5p
            {"In", {2, 2, 6, 2, 6, 10, 2, 6, 10, 0, 2, 1}}, {"Sn", {2, 2, 6, 2, 6, 10, 2, 6, 10, 0, 2, 2}}, {"Sb", {2, 2, 6, 2, 6, 10, 2, 6, 10, 0, 2, 3}}, // 6th period
            {"Te", {2, 2, 6, 2, 6, 10, 2, 6, 10, 0, 2, 4}}, {"I",  {2, 2, 6, 2, 6, 10, 2, 6, 10, 0, 2, 5}}, {"Xe", {2, 2, 6, 2, 6, 10, 2, 6, 10, 0, 2, 6}}, // 6th period
            //      1s 2s 2p 3s 3p 3d  4s 4p 4d  4f 5s 5p 5d 5f 5g 6s           1s 2s 2p 3s 3p 3d  4s 4p 4d  4f 5s 5p 5d 5f 5g 6s
            {"Cs", {2, 2, 6, 2, 6, 10, 2, 6, 10, 0, 2, 6, 0, 0, 0, 1}}, {"Ba", {2, 2, 6, 2, 6, 10, 2, 6, 10, 0, 2, 6, 0, 0, 0, 2}}, // 6th period
            {"La", {2, 2, 6, 2, 6, 10, 2, 6, 10, 0, 2, 6, 1, 0, 0, 2}}, {"Ce", {2, 2, 6, 2, 6, 10, 2, 6, 10, 1, 2, 6, 1, 0, 0, 2}}, // 6th period
            {"Pr", {2, 2, 6, 2, 6, 10, 2, 6, 10, 3, 2, 6, 0, 0, 0, 2}}, {"Nd", {2, 2, 6, 2, 6, 10, 2, 6, 10, 4, 2, 6, 0, 0, 0, 2}}, // 6th period
            {"Pm", {2, 2, 6, 2, 6, 10, 2, 6, 10, 5, 2, 6, 0, 0, 0, 2}}, {"Sm", {2, 2, 6, 2, 6, 10, 2, 6, 10, 6, 2, 6, 0, 0, 0, 2}}, // 6th period
            {"Eu", {2, 2, 6, 2, 6, 10, 2, 6, 10, 7, 2, 6, 0, 0, 0, 2}}, {"Gd", {2, 2, 6, 2, 6, 10, 2, 6, 10, 7, 2, 6, 1, 0, 0, 2}}, // 6th period
            {"Tb", {2, 2, 6, 2, 6, 10, 2, 6, 10, 9, 2, 6, 0, 0, 0, 2}}, {"Dy", {2, 2, 6, 2, 6, 10, 2, 6, 10, 10, 2, 6, 0, 0, 0, 2}}, // 6th period
            {"Ho", {2, 2, 6, 2, 6, 10, 2, 6, 10, 11, 2, 6, 0, 0, 0, 2}}, {"Er", {2, 2, 6, 2, 6, 10, 2, 6, 10, 12, 2, 6, 0, 0, 0, 2}}, // 6th period
            {"Tm", {2, 2, 6, 2, 6, 10, 2, 6, 10, 13, 2, 6, 0, 0, 0, 2}}, {"Yb", {2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 0, 0, 0, 2}}, // 6th period
            {"Lu", {2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 1, 0, 0, 2}}, {"Hf", {2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 2, 0, 0, 2}}, // 6th period
            {"Ta", {2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 3, 0, 0, 2}}, {"W",  {2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 4, 0, 0, 2}}, // 6th period
            {"Re", {2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 5, 0, 0, 2}}, {"Os", {2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 6, 0, 0, 2}}, // 6th period
            {"Ir", {2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 7, 0, 0, 2}}, {"Pt", {2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 9, 0, 0, 1}}, // 6th period
            {"Au", {2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 0, 0, 1}}, {"Hg", {2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 0, 0, 2}}, // 6th period
            //      1s 2s 2p 3s 3p 3d  4s 4p 4d  4f  5s 5p 5d  5f 5g 6s 6p           1s 2s 2p 3s 3p 3d  4s 4p 4d  4f  5s 5p 5d  5f 5g 6s 6p
            {"Tl", {2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 0, 0, 2, 1}}, {"Pb", {2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 0, 0, 2, 2}}, // 6th period
            {"Bi", {2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 0, 0, 2, 3}}, {"Po", {2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 0, 0, 2, 4}}, // 6th period
            {"At", {2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 0, 0, 2, 5}}, {"Rn", {2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 0, 0, 2, 6}} // 6th period
        };
};

#endif