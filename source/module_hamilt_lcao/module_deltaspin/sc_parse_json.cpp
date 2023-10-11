#include "spin_constrain.h"

#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
#include <regex>

template<typename FPTYPE, typename Device>
const std::map<int, std::vector<ScAtomData>>& SpinConstrain<FPTYPE, Device>::get_ScData() const
{
    return this->ScData;
}

template<typename FPTYPE, typename Device>
void SpinConstrain<FPTYPE, Device>::clear_ScData()
{
    this->ScData.clear();
}

template<typename FPTYPE, typename Device>
void SpinConstrain<FPTYPE, Device>::Set_ScData_From_Json(const std::string& filename)
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
    std::regex target_mag_regex("\"target_mag\": \\[(.+?)\\]");
    std::regex target_mag_val_regex("\"target_mag_val\": ([0-9.]+)");
    std::regex target_mag_angle1_regex("\"target_mag_angle1\": ([0-9.]+)");
    std::regex target_mag_angle2_regex("\"target_mag_angle2\": ([0-9.]+)");
    std::regex constrain_regex("\"constrain\": \\[(.+?)\\]");

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

            if (std::regex_search(line, match, target_mag_regex)) {
                std::stringstream ss(match[1]);
                double value;
                while (ss >> value) {
                    element_data.target_mag.push_back(value);
                    if (ss.peek() == ',') {
                        ss.ignore();
                    }
                }
            } else {
                if (std::regex_search(line, match, target_mag_val_regex)) {
                    element_data.target_mag_val = std::stod(match[1]);
                }
                getline(file, line); // Read the following line
                if (std::regex_search(line, match, target_mag_angle1_regex)) {
                    element_data.target_mag_angle1 = std::stod(match[1]);
                }
                getline(file, line); // Read the following line
                if (std::regex_search(line, match, target_mag_angle2_regex)) {
                    element_data.target_mag_angle2 = std::stod(match[1]);
                }
            }

            getline(file, line); // Read the following line

            if (std::regex_search(line, match, constrain_regex)) {
                std::stringstream ss(match[1]);
                int value;
                while (ss >> value) {
                    element_data.constrain.push_back(value);
                    if (ss.peek() == ',') {
                        ss.ignore();
                    }
                }
            }

            this->ScData[current_itype].push_back(element_data);
        }
    }
    file.close();
}

template class SpinConstrain<std::complex<double>, psi::DEVICE_CPU>;