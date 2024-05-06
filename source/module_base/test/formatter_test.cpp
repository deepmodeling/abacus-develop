#include "module_base/formatter.h"
#include <iostream>
#include <gtest/gtest.h>

TEST(FormatterTest, FmtCoreStaticFormat) {
  // const char*
  std::string result = FmtCore::format("Hello, %s!", "world");
  // remove the last '\0' character
  EXPECT_EQ(result, "Hello, world!");
  // std::string
  result = FmtCore::format("Hello, %s!", std::string("world"));
  EXPECT_EQ(result, "Hello, world!");
  // int
  result = FmtCore::format("Hello, %d!", 123);
  EXPECT_EQ(result, "Hello, 123!");
  // float
  result = FmtCore::format("Hello, %f!", 123.456);
  EXPECT_EQ(result, "Hello, 123.456000!");
  // char
  result = FmtCore::format("Hello, %c!", 'a');
  EXPECT_EQ(result, "Hello, a!");
  // bool
  result = FmtCore::format("Hello, %s!", true);
  EXPECT_EQ(result, "Hello, true!");
  // invalid format
  result = FmtCore::format("Hello, %z!", "world");
  EXPECT_EQ(result, "Hello, %!");
  // std::complex
  result = FmtCore::format("Hello, (%d,%d)!", std::complex<int>(1, 2));
  EXPECT_EQ(result, "Hello, (1,2)!");
  // std::complex, %12.3f + %20.10f
  result = FmtCore::format("Hello, (%12.3f, %20.10f)!", std::complex<double>(1.0, 2.0));
  EXPECT_EQ(result, "Hello, (       1.000,         2.0000000000)!");
}




int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}