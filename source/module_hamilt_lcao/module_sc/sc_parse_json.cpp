#include "spin_constrain.h"

#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
#include <regex>

const std::map<int, std::vector<ScAtomData>>& SpinConstrain::get_ScData() const
{
    return this->ScData;
}

void SpinConstrain::clear_ScData()
{
    this->ScData.clear();
}

void SpinConstrain::Set_ScData_From_Json(const std::string& filename)
{
    ModuleBase::TITLE("SpinConstrain", "ScJsonFile");
    std::ifstream file(filename);
    if (!file.is_open()) {
        ModuleBase::WARNING_QUIT("SpinConstrain::parseScJsonFile","Error opening sc_file");
        return;
    }

    std::string line;
    int current_itype;
    std::string current_element;

    std::regex itype_regex("\"itype\": (\\d+)");
    std::regex element_regex("\"element\": \"([A-Za-z]+)\"");
    std::regex index_regex("\"index\": (\\d+)");
    std::regex lambda_regex("\"lambda\": \\[(.+?)\\]");
    std::regex sc_mag_regex("\"sc_mag\": \\[(.+?)\\]");
    std::regex sc_spin_val_regex("\"sc_spin_val\": ([0-9.]+)");
    std::regex sc_spin_angle1_regex("\"sc_spin_angle1\": ([0-9.]+)");
    std::regex sc_spin_angle2_regex("\"sc_spin_angle2\": ([0-9.]+)");

    while (getline(file, line)) {
        std::smatch match;

        if (std::regex_search(line, match, itype_regex)) {
            current_itype = std::stoi(match[1]);
        } else if (std::regex_search(line, match, element_regex)) {
            current_element = match[1];
        } else if (std::regex_search(line, match, index_regex)) {
            ScAtomData element_data;
            element_data.index = std::stoi(match[1]);

            getline(file, line); // Read the following line

            if (std::regex_search(line, match, lambda_regex)) {
                std::stringstream ss(match[1]);
                double value;
                while (ss >> value) {
                    element_data.lambda.push_back(value);
                    if (ss.peek() == ',') {
                        ss.ignore();
                    }
                }
            }

            getline(file, line); // Read the following line

            if (std::regex_search(line, match, sc_mag_regex)) {
                std::stringstream ss(match[1]);
                double value;
                while (ss >> value) {
                    element_data.sc_mag.push_back(value);
                    if (ss.peek() == ',') {
                        ss.ignore();
                    }
                }
            } else {
                if (std::regex_search(line, match, sc_spin_val_regex)) {
                    element_data.sc_spin_val = std::stod(match[1]);
                }
                getline(file, line); // Read the following line
                if (std::regex_search(line, match, sc_spin_angle1_regex)) {
                    element_data.sc_spin_angle1 = std::stod(match[1]);
                }
                getline(file, line); // Read the following line
                if (std::regex_search(line, match, sc_spin_angle2_regex)) {
                    element_data.sc_spin_angle2 = std::stod(match[1]);
                }
            }

            this->ScData[current_itype].push_back(element_data);
        }
    }
    file.close();
}