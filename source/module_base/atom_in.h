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
          
};

#endif