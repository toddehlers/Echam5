export

SHELL = /bin/sh
 
ARCH  = LINUX

srcdir = .
top_srcdir = .

prefix = .
exec_prefix = ${prefix}

bindir = ${exec_prefix}/bin
sbindir = ${exec_prefix}/sbin
libexecdir = ${exec_prefix}/libexec
datadir = ${datarootdir}
datarootdir = ${prefix}/share
sysconfdir = ${prefix}/etc
libdir = ${exec_prefix}/lib
includedir = ${prefix}/include
oldincludedir = /usr/include
infodir = ${datarootdir}/info
mandir = ${datarootdir}/man

sharedstatedir = ${prefix}/com
localstatedir = ${prefix}/var

program_transform_name = s,x,x,

MPIROOT        = /opt/bwgrid/mpi/openmpi/1.4.3-intel-12.0
MPI_LIB        = -L$(MPIROOT)/lib -lmpi_f90 -lmpi_f77 -lmpi -lopen-rte -lopen-pal 
MPI_INCLUDE    = $(MODOPT)$(MPIROOT)/lib

NETCDFROOT     = /home-link/epifm01/NETCDFROOTsymbolic
NETCDF_LIB     = -L$(NETCDFROOT)/lib -lnetcdf
NETCDF_INCLUDE = -I$(NETCDFROOT)/include

LIB      = -L../lib -lsupport -L/opt/bwgrid/numlib/lapack/3.3.1/lib -llapack -lblas -L/opt/bwgrid/lib/netcdf/lib-netcdf-4.2.1.1_c-openmpi-1.4.3-gnu-4.1/lib -lnetcdf -L/opt/bwgrid/lib/netcdf/lib-netcdf-4.2_fortran-openmpi-1.4.3-gnu-4.1/lib -lnetcdff -L/opt/bwgrid/compiler/intel/cs_2011.5.220/composerxe-2011.5.220/compiler/lib/intel64/ -limf -lm
LIBS     = $(LIB) $(NETCDF_LIB) $(MPI_LIB)

MODOPT   = -I
MODULES  = ../modules

INCLUDE  = ../include
INCLUDES = $(MODOPT)$(MODULES) -I$(INCLUDE) $(NETCDF_INCLUDE) $(MPI_INCLUDE)

F90      = ifort 
FC       = ifort 
CC       = icc 
CPP      = icc  -E
AR       = ar
AS       = as 

DEFS     = -DHAVE_CONFIG_H

CFLAGS   = -I../config -O -DpgiFortran
FFLAGS   = -mieee-fp -O3 -ip -g 
F90FLAGS = $(INCLUDES) -mieee-fp -O3 -ip -g -fpp -DHAVE_LIBNETCDF64 
CPPFLAGS = 
ARFLAGS  = crv
LDFLAGS  = -mieee-fp -O3 -ip -g -fpp -DHAVE_LIBNETCDF64

SRCDIRS =  support modules src

all:
	@for DIR in $(SRCDIRS) ;\
	  do \
	    back=`pwd`; \
	    cd $$DIR ;\
	    $(MAKE) ; status=$$? ; \
	    if [ $$status != 0 ] ; then \
	      echo "Exit status from make was $$status" ; exit $$status ; \
	    fi ; \
	    cd $$back ; \
	  done 

clean:
	@for DIR in $(SRCDIRS) ;\
	  do \
	  (cd $$DIR ;\
	  $(MAKE) clean ; if [ $$? != 0 ] ; then \
	        echo "Exit status from make was $$?" ; exit 1 ; fi ;) ;\
	  done
	-rm -f config.cache
	-rm -f lib/*.a bin/echam5
	-rm -f html/[a-z]*

tar:
	@tarfile=../echam5.f90.`date +%y%m%d`.taz ; gtar zcvf $$tarfile \
	`ls */*.f90 */*.[fhc] */*inc */Makefile Makefile.in Makefile run/hjob*`

index:
	-rm -f html/[a-z]*
	util/f2html.pl -f util/fgenrc -d html \
          support modules src include
