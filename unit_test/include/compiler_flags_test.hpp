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

#ifndef	COMPILER_FLAG_TEST_HPP
#define	COMPILER_FLAG_TEST_HPP

//////////////////////////////////////////////////////////////////////////////
// Test fixture class: compiler flags
//////////////////////////////////////////////////////////////////////////////

class CompileFlagTest : public ::testing::Test
{
 protected:

  CompileFlagTest()
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

TEST_F(CompileFlagTest, AllFlags)
{
  if (qhipster::mpi::Environment::GetPoolRank() == 0)
  {
    std::cout << "             |----------------------------|----------|\n"
              << "             |        Compiler Flag       | defined? |\n"
              << "             |----------------------------|----------|\n";
#ifdef INTELQS_HAS_MPI
    std::cout << "             |            INTELQS_HAS_MPI |    X     |\n";
#else
    std::cout << "             |            INTELQS_HAS_MPI |          |\n";
#endif
//
#ifdef USE_MKL
    std::cout << "             |                    USE_MKL |    X     |\n";
#else
    std::cout << "             |                    USE_MKL |          |\n";
#endif
//
#ifdef _OPENMP
    std::cout << "             |                    _OPENMP |    X     |\n";
#else
    std::cout << "             |                    _OPENMP |          |\n";
#endif
//
#ifdef __INTEL_COMPILER
    std::cout << "             |           __INTEL_COMPILER |    X     |\n";
#else
    std::cout << "             |           __INTEL_COMPILER |          |\n";
#endif
//
#ifdef __ICC
    std::cout << "             |                      __ICC |    X     |\n";
#else
    std::cout << "             |                      __ICC |          |\n";
#endif
//
#ifdef USE_MM_MALLOC
    std::cout << "             |              USE_MM_MALLOC |    X     |\n";
#else
    std::cout << "             |              USE_MM_MALLOC |          |\n";
#endif
//
#ifdef STANDALONE
    std::cout << "             |                 STANDALONE |    X     |\n";
#else
    std::cout << "             |                 STANDALONE |          |\n";
#endif
//
#ifdef NDEBUG
    std::cout << "             |                     NDEBUG |    X     |\n";
#else
    std::cout << "             |                     NDEBUG |          |\n";
#endif
//
#ifdef __ONLY_NORMALIZED_STATES__
    std::cout << "             | __ONLY_NORMALIZED_STATES__ |    X     |\n";
#else
    std::cout << "             | __ONLY_NORMALIZED_STATES__ |          |\n";
#endif
    std::cout << "             |----------------------------|----------|\n";
  }
  qhipster::mpi::PoolBarrier();
}

//////////////////////////////////////////////////////////////////////////////

#endif	// header guard COMPILER_FLAG_TEST_HPP
