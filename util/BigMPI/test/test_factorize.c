#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <mpi.h>

#ifdef __x86_64__
static inline unsigned long long rdtsc(void)
{
      unsigned hi, lo;
        __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
          return ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 );
}
#endif

/* This function is copied from src/type_contiguous_x.c and must be updated
 * manually if the implementation changes.  It is an internal function so
 * it is static and therefore be cannot call the symbol from the library. */

/*
 * Synopsis
 *
 * int BigMPI_Factorize_count(MPI_Count c, int * a, int *b)
 *
 *  Input Parameter
 *
 *   c                  large count
 *
 * Output Parameters
 *
 *   a, b               integers such that c=a*b and a,b<INT_MAX
 *   rc                 returns 0 if a,b found (success), else 1 (failure)
 *
 */
static int BigMPI_Factorize_count(MPI_Count count, int * a, int *b)
{
    size_t in = (size_t)count;
    /* Not using bigmpi_int_max because we want to debug the case where the library is actually used. */
    size_t lo = in/INT_MAX+1;
    size_t hi = (size_t)floor(sqrt((double)in));
    for (size_t g=hi; g>lo; g--) {
        size_t rem = in%g;
        if (rem==0) {
            *a = (int)g;
            *b = (int)(in/g);
            return 0;
        }
    }
    *a = -1;
    *b = -1;
    return 1;
}

#define TIMING
#define DEBUG

int main(int argc, char* argv[])
{
    int rank=0, size=1;
#ifdef PARALLEL
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
    MPI_Count max = (argc>1) ? atol(argv[1]) : 1LL<<60;
    MPI_Count inc = size; /* Incremement by nproc to distribute test work. */

#ifdef TIMING
#ifdef PARALLEL
    double t0 = MPI_Wtime();
#else
    unsigned long long t0 = rdtsc();
#endif
#endif
    for (MPI_Count count=1; count<max; count+=inc) {
        int a, b;
        int rc = BigMPI_Factorize_count(count, &a, &b);
#ifdef DEBUG
        printf("factorized %zu = %d * %d (rc=%d)\n", (size_t)count, a, b, rc);
#endif
    }
#ifdef TIMING
#ifdef PARALLEL
    double t1 = MPI_Wtime();
    double dt = t1-t0;
#else
    unsigned long long t1 = rdtsc();
    double dt = (1.e-9)*(t1-t0);
#endif
    printf("factorize 1 to %zu in %lf s (%lf us per call)\n",
            (size_t)max, dt, 1.e6*dt/(double)max);
#endif

#ifdef PARALLEL
    MPI_Finalize();
#endif
    return 0;
}
