#ifndef	CONVERSION_TEST_HPP
#define	CONVERSION_TEST_HPP

#include "../../include/conversion.hpp"

#include <string>

//////////////////////////////////////////////////////////////////////////////
// Test fixture class: compiler flags
//////////////////////////////////////////////////////////////////////////////

class ConversionTest : public ::testing::Test
{
 protected:

  ConversionTest()
  { }

  void SetUp() override
  {
    // All tests are skipped if the rank is dummy.
    if (qhipster::mpi::Environment::IsUsefulRank() == false)
      GTEST_SKIP();
  }
};

//////////////////////////////////////////////////////////////////////////////
// Test macros:

TEST_F(ConversionTest, BasicUse)
{
  ASSERT_TRUE(qhipster::toString(1)  =="1");
  ASSERT_TRUE(qhipster::toString('a')=="a");
  int var = 25;
  ASSERT_TRUE(qhipster::toString(var)=="25");
  std::string s = "I have 2 apples";
  ASSERT_EQ(qhipster::toString(s), s);
}

//////////////////////////////////////////////////////////////////////////////

#endif	// header guard CONVERSION_TEST_HPP
