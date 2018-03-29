#
# Copyright (C) 2014. See LICENSE in top-level directory.
#

ACLOCAL_AMFLAGS = -I m4
AM_CPPFLAGS = -I$(top_srcdir)/src

lib_LTLIBRARIES = libbigmpi.la

noinst_LTLIBRARIES = libbigmpii.la

libbigmpi_la_SOURCES = 	src/collectives_x.c \
			src/vcollectives_x.c \
			src/neighborhood_collectives_x.c \
			src/reductions_x.c \
			src/rma_x.c \
			src/sendrecv_x.c \
			src/fileio_x.c \
			src/type_contiguous_x.c \
			src/type_hindexed_x.c  \
			src/utils.c

libbigmpi_la_LDFLAGS = -version-info $(libbigmpi_abi_version)

libbigmpii_la_SOURCES = $(libbigmpi_la_SOURCES)
libbigmpii_la_LDFLAGS = $(libbigmpi_abi_version)

include_HEADERS = src/bigmpi.h src/bigmpi_impl.h

bin_PROGRAMS =
check_PROGRAMS = 
TESTS = 
XFAIL_TESTS = 

MPIEXEC = mpiexec -n 2
TESTS_ENVIRONMENT = $(MPIEXEC)

include test/Makefile.mk

.PHONY: checkprogs
checkprogs: $(check_PROGRAMS)
