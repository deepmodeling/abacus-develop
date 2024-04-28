#include "module_base/formatter.h"
#include <iostream>
#include <gtest/gtest.h>

TEST(FormatterTest, FmtCoreStaticFormatConstCharPtr) {
  std::string result = FmtCore::format("Hello, %s!", "world");
  EXPECT_EQ(result, "Hello, world!");
}

TEST(FormatterTest, FmtCoreStaticFormatStdString) {
  std::string result = FmtCore::format("Hello, %s!", std::string("world"));
  EXPECT_EQ(result, "Hello, world!");
}

TEST(FormatterTest, FmtCoreStaticFormatInt) {
  std::string result = FmtCore::format("Hello, %d!", 123);
  EXPECT_EQ(result, "Hello, 123!");
}

TEST(FormatterTest, FmtCoreStaticFormatDouble) {
  std::string result = FmtCore::format("Hello, %f!", 123.456);
  EXPECT_EQ(result, "Hello, 123.456000!");
}

TEST(FormatterTest, FmtCoreStaticFormatChar) {
  std::string result = FmtCore::format("Hello, %c!", 'a');
  EXPECT_EQ(result, "Hello, a!");
}

TEST(FormatterTest, FmtCoreStaticFormatBool) {
  std::string result = FmtCore::format("Hello, %s!", true);
  EXPECT_EQ(result, "Hello, true!");
}

TEST(FormatterTest, FmtCoreStaticFormatPointer) {
  std::string result = FmtCore::format("Hello, %p!", (void*)0x12345678);
  EXPECT_EQ(result, "Hello, 0x12345678!");
}

TEST(FormatterTest, FmtCoreStaticFormatMultiple) {
  std::string result = FmtCore::format("Hello, %s! %d %f %c %s %p", "world", 123, 123.456, 'a', true, (void*)0x12345678);
  EXPECT_EQ(result, "Hello, world! 123 123.456000 a true 0x12345678");
}

TEST(FormatterTest, FmtCoreStaticFormatEmpty) {
  std::string result = FmtCore::format("");
  EXPECT_EQ(result, "");
}

TEST(FormatterTest, FmtCoreStaticFormatNoPlaceholder) {
  std::string result = FmtCore::format("Hello, world!");
  EXPECT_EQ(result, "Hello, world!");
}

TEST(FormatterTest, FmtCoreStaticFormatNoPlaceholderWithArgs) {
  std::string result = FmtCore::format("Hello, world!", 123, 123.456, 'a', true, (void*)0x12345678);
  EXPECT_EQ(result, "Hello, world!");
}

TEST(FormatterTest, FmtCoreStaticFormatNoArgs) {
  std::string result = FmtCore::format("Hello, %s!");
  EXPECT_EQ(result, "Hello, %s!");
}

TEST(FormatterTest, FmtCoreStaticFormatInvalidPlaceholder) {
  std::string result = FmtCore::format("Hello, %z!", "world");
  EXPECT_EQ(result, "Hello, %z!");
}

TEST(FormatterTest, FmtCoreStaticFormatInvalidPlaceholderWithArgs) {
  std::string result = FmtCore::format("Hello, %z!", 123, 123.456, 'a', true, (void*)0x12345678);
  EXPECT_EQ(result, "Hello, %z!");
}




int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}