//------------------------------------------------------------------------------
// Copyright (C) 2019 Intel Corporation 
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//------------------------------------------------------------------------------

#ifndef	CONVERSION_TEST_HPP
#define	CONVERSION_TEST_HPP

#include "../../util/conversion.hpp"

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
