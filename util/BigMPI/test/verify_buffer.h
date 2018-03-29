#ifndef VERIFY_BUFFER_H
#define VERIFY_BUFFER_H

#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <mpi.h>

/* This is the generic buffer verification for when char is the
 * type and we use memset for assignment. */
static size_t verify_buffer(char *buf, MPI_Count count, int expected_value)
{
    assert(count<SIZE_MAX);

    size_t errors = 0;
    for (size_t i = 0; i < (size_t)count; i++) {
        errors += (buf[i] != (unsigned char)expected_value);
    }
    return errors;
}

static size_t verify_doubles(double *buf, MPI_Count count, double expected_value)
{
    assert(count<SIZE_MAX);

    size_t errors = 0;
    for (size_t i = 0; i < (size_t)count; i++) {
        double absdiff = fabs(buf[i] - expected_value);
        if (absdiff>1.e-4) errors++;
    }
    return errors;
}

static void set_doubles(double *buf, MPI_Count count, double value)
{
    assert(count<SIZE_MAX);

    for (size_t i = 0; i < (size_t)count; i++) {
        buf[i] = value;
    }
}

#endif // VERIFY_BUFFER_H
