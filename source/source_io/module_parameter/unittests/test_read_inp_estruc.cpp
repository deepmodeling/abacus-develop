#include "source_io/module_parameter/read_input.h"
#include <gtest/gtest.h>

#include <utility>

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
