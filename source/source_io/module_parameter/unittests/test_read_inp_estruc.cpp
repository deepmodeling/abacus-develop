#include "source_io/module_parameter/read_input.h"
#include <gtest/gtest.h>

#include <string>
#include <utility>
#include <vector>

TEST(LibxcDeviceInput, DefaultAndAllowedDevices)
{
    ModuleIO::ReadInput schema(0);
    for (const std::pair<std::string, ModuleIO::Input_Item>& entry : schema.get_input_lists())
    {
        if (entry.first == "device_xc")
        {
            const ModuleIO::InputDomain& domain = entry.second.domain;
            const Parameter defaults;
            EXPECT_EQ(domain.value(defaults).format(), "cpu");
            EXPECT_TRUE(domain.validate_explicit(entry.first, {"gpu"}).empty());
            EXPECT_FALSE(domain.validate_explicit(entry.first, {"auto"}).empty());

            const ModuleIO::InputValue gpu(std::string{"gpu"});
            const ModuleIO::InputValue cpu(std::string{"cpu"});
            ModuleIO::InputValues values = {{"device_xc", gpu}, {"device", cpu}};
            EXPECT_FALSE(domain.validate_requirements(entry.first, values).empty());
            values["device"] = gpu;
            EXPECT_TRUE(domain.validate_requirements(entry.first, values).empty());
            return;
        }
    }
    FAIL() << "device_xc is missing from the INPUT schema";
}

TEST(DmDeviceInput, DefaultSelectionAndInputRoundTrip)
{
    ModuleIO::ReadInput schema(0);
    for (const std::pair<std::string, ModuleIO::Input_Item>& entry : schema.get_input_lists())
    {
        if (entry.first != "device_dm") continue;
        ModuleIO::Input_Item item = entry.second;
        Parameter input;
        EXPECT_EQ(input.inp.device_dm, "cpu");
        EXPECT_EQ(item.default_value, "cpu");
        EXPECT_EQ(item.domain.format_choices(), "cpu, gpu");
        for (const std::string device : {"cpu", "gpu"})
        {
            item.str_values = {device};
            EXPECT_TRUE(item.domain.validate_explicit(item.label, item.str_values).empty());
            item.read_value(item, input);
            EXPECT_EQ(input.inp.device_dm, device);
            item.final_value.str("");
            item.get_final_value(item, input);
            EXPECT_EQ(item.final_value.str(), device);
        }
        for (const std::vector<std::string> tokens : {std::vector<std::string>{}, {"auto"}, {"direct"}, {"gpu", "cpu"}})
        {
            EXPECT_FALSE(item.domain.validate_explicit(item.label, tokens).empty());
        }
        return;
    }
    FAIL() << "device_dm is missing from the INPUT schema";
}

TEST(DmDeviceInput, RejectsIncompatibleExplicitGpuRequests)
{
    ModuleIO::ReadInput schema(0);
    for (const std::pair<std::string, ModuleIO::Input_Item>& entry : schema.get_input_lists())
    {
        if (entry.first != "device_dm") continue;
        const ModuleIO::InputDomain& domain = entry.second.domain;
        ModuleIO::InputValues valid;
        valid["device_dm"] = ModuleIO::InputValue(std::string("gpu"));
        valid["device"] = ModuleIO::InputValue(std::string("gpu"));
        valid["basis_type"] = ModuleIO::InputValue(std::string("lcao"));
        valid["ks_solver"] = ModuleIO::InputValue(std::string("cusolver"));
        valid["gamma_only"] = ModuleIO::InputValue(true);
        valid["nspin"] = ModuleIO::InputValue(1);
        valid["kpar"] = ModuleIO::InputValue(1);
        valid["precision"] = ModuleIO::InputValue(std::string("double"));
        valid["lcao_sk_cut"] = ModuleIO::InputValue(0.0);
        valid["esolver_type"] = ModuleIO::InputValue(std::string("ksdft"));
        valid["calculation"] = ModuleIO::InputValue(std::string("scf"));
        EXPECT_TRUE(domain.validate_requirements(entry.first, valid).empty());
        const ModuleIO::InputValues incompatible = {
            {"device", ModuleIO::InputValue(std::string("cpu"))},
            {"basis_type", ModuleIO::InputValue(std::string("pw"))},
            {"ks_solver", ModuleIO::InputValue(std::string("cusolvermp"))},
            {"gamma_only", ModuleIO::InputValue(false)},
            {"nspin", ModuleIO::InputValue(2)},
            {"kpar", ModuleIO::InputValue(2)},
            {"precision", ModuleIO::InputValue(std::string("single"))},
            {"lcao_sk_cut", ModuleIO::InputValue(1.0e-8)},
            {"esolver_type", ModuleIO::InputValue(std::string("tddft"))},
            {"calculation", ModuleIO::InputValue(std::string("nscf"))}};
        for (const std::pair<const std::string, ModuleIO::InputValue>& mismatch : incompatible)
        {
            ModuleIO::InputValues values = valid;
            values[mismatch.first] = mismatch.second;
            EXPECT_FALSE(domain.validate_requirements(entry.first, values).empty()) << mismatch.first;
            values["device_dm"] = ModuleIO::InputValue(std::string("cpu"));
            EXPECT_TRUE(domain.validate_requirements(entry.first, values).empty()) << mismatch.first;
        }
        return;
    }
    FAIL() << "device_dm is missing from the INPUT schema";
}
