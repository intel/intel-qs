#include <mpi.h>

#define SUM 0x1

void Bar(int a) { return; };

#define MAKE_FOO(OP) void FOO_##OP##_fn(){ Bar(MPI_##OP);\
                                           return; }

MAKE_FOO(SUM);

int main(void)
{
    return 0;
}
