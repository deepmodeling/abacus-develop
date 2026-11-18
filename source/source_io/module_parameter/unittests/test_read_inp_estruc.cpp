#include "source_io/module_parameter/read_input.h"
#include <gtest/gtest.h>

TEST(LibxcDeviceInput, DefaultAndAllowedDevices)
{
    ModuleIO::ReadInput schema(0);
    for (const auto& entry : schema.get_input_lists())
    {
        if (entry.first == "xc_device")
        {
            const auto& domain = entry.second.domain;
            const Parameter defaults;
            EXPECT_EQ(domain.value(defaults).format(), "cpu");
            EXPECT_TRUE(domain.validate_explicit(entry.first, {"gpu"}).empty());
            EXPECT_FALSE(domain.validate_explicit(entry.first, {"invalid"}).empty());

            const ModuleIO::InputValue gpu(std::string{"gpu"});
            const ModuleIO::InputValue cpu(std::string{"cpu"});
            ModuleIO::InputValues values = {{"xc_device", gpu}, {"device", cpu}};
            EXPECT_FALSE(domain.validate_requirements(entry.first, values).empty());
            values["device"] = gpu;
            EXPECT_TRUE(domain.validate_requirements(entry.first, values).empty());
            return;
        }
    }
    FAIL() << "xc_device is missing from the INPUT schema";
}
