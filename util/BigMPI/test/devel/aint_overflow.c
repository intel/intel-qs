#include <stdio.h>
#include <stdint.h>
#include <limits.h>

#include <mpi.h>

int main(int argc, char * argv[])
{
    printf("INT_MAX = %d\n", INT_MAX);

    int      a = (int)INT_MAX*100;
    MPI_Aint b = (MPI_Aint)INT_MAX*100;
    int64_t  c = (int64_t)INT_MAX*100;
    uint64_t d = (uint64_t)INT_MAX*100;

    printf("a = %d\n", a);
    printf("b = %zu\n", b);
    printf("c = %lld\n", c);
    printf("d = %llu\n", d);

    int      a2 = 100;
    MPI_Aint b2 = 100;
    int64_t  c2 = 100;
    uint64_t d2 = 100;

    a2 *= INT_MAX;
    b2 *= INT_MAX;
    c2 *= INT_MAX;
    d2 *= INT_MAX;

    printf("a2 = %d\n", a2);
    printf("b2 = %zu\n", b2);
    printf("c2 = %lld\n", c2);
    printf("d2 = %llu\n", d2);

    return 0;
}
