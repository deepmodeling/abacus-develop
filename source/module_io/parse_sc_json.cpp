#include "input_conv.h"

#include <sstream>
#include <regex>
#include "module_base/tool_title.h"

void Input_Conv::parseScJsonFile(const std::string& filename, std::map<std::string, std::vector<ScElementData>>& data)
{
    ModuleBase::TITLE("Input_Conv", "ScJsonFile");
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    std::string line;
    std::string current_element;

    std::regex element_regex("\"([A-Za-z]+)\": \\[");
    std::regex index_regex("\"index\": (\\d+)");
    std::regex lambda_regex("\"lambda\": \\[(.+?)\\]");
    std::regex sc_mag_regex("\"sc_mag\": \\[(.+?)\\]");
    std::regex sc_spin_val_regex("\"sc_spin_val\": ([0-9.]+)");
    std::regex sc_spin_angle1_regex("\"sc_spin_angle1\": ([0-9.]+)");
    std::regex sc_spin_angle2_regex("\"sc_spin_angle2\": ([0-9.]+)");

    while (getline(file, line)) {
        std::smatch match;

        if (std::regex_search(line, match, element_regex)) {
            current_element = match[1];
        } else if (std::regex_search(line, match, index_regex)) {
            ScElementData element_data;
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

            data[current_element].push_back(element_data);
        }
    }
}