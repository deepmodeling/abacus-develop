#include "gtest/gtest.h"

#include "../paw_element.h"

class PawTest : public testing::Test
{
    protected:

    Paw_Element paw_element;
};

TEST_F(PawTest, ReadPaw)
{
    paw_element.read_paw_xml("H.LDA_PW-JTH.xml");
}
