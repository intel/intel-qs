/*
 * Copyright (C) 2014. See LICENSE in top-level directory.
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <mpi.h>
#include "bigmpi.h"

int main(int argc, char ** argv)
{
  MPI_Init(&argc, &argv);
  assert(0);
  MPI_Finalize();
  return 0;
}
