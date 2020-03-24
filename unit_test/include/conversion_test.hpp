#ifndef	CONVERSION_TEST_HPP
#define	CONVERSION_TEST_HPP

#include "../../include/conversion.hpp"

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
}

//////////////////////////////////////////////////////////////////////////////

#endif	// header guard CONVERSION_TEST_HPP
