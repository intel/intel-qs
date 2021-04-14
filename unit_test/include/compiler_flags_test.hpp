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
    if (iqs::mpi::Environment::IsUsefulRank() == false)
      GTEST_SKIP();
  }
};

//////////////////////////////////////////////////////////////////////////////
// Test macros:

TEST_F(CompileFlagTest, AllFlags)
{
  if (iqs::mpi::Environment::GetPoolRank() == 0)
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
  iqs::mpi::PoolBarrier();
}

//////////////////////////////////////////////////////////////////////////////

#endif	// header guard COMPILER_FLAG_TEST_HPP
