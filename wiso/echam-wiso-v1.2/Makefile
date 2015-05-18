export

SHELL = /bin/sh
 
ARCH  = SX

srcdir = .
top_srcdir = .

prefix = .
exec_prefix = ${prefix}

bindir = ${exec_prefix}/bin
sbindir = ${exec_prefix}/sbin
libexecdir = ${exec_prefix}/libexec
datadir = ${prefix}/share
sysconfdir = ${prefix}/etc
libdir = ${exec_prefix}/lib
includedir = ${prefix}/include
oldincludedir = /usr/include
infodir = ${prefix}/info
mandir = ${prefix}/man

sharedstatedir = ${prefix}/com
localstatedir = ${prefix}/var

program_transform_name = s,x,x,

MPIROOT        = /SX/opt/mpisx/inst
MPI_LIB        = -L$(MPIROOT)/lib -lmpi
MPI_INCLUDE    = -I$(MPIROOT)/include

NETCDFROOT     = /sx8/user2/paleo/utils.sx8/netcdf
NETCDF_LIB     = -L$(NETCDFROOT)/lib -lnetcdf
NETCDF_INCLUDE = -I$(NETCDFROOT)/include

LIB      = -L../lib -lsupport -L/SX/opt/MathKeisan/lib -llapack -lblas
LIBS     = $(LIB) $(NETCDF_LIB) $(MPI_LIB)

MODOPT   = -I
MODULES  = ../modules

INCLUDE  = ../include
INCLUDES = $(MODOPT)$(MODULES) -I$(INCLUDE) $(NETCDF_INCLUDE) $(MPI_INCLUDE)

F90      = sxf90
FC       = sxf90
CC       = sxc++
CPP      = 
AR       = sxar
AS       = sxas 

DEFS     = -DHAVE_CONFIG_H

CFLAGS   = -I../config -sx8 -Chopt
FFLAGS   = -sx8 -Chopt
F90FLAGS = $(INCLUDES) -sx8 -Popenmp -Chopt -Ep -Wl,-Z 1000000 -pi line=1000 -DHAVE_LIBNETCDF64 
CPPFLAGS = 
ARFLAGS  = crv
LDFLAGS  = -sx8 -Popenmp -Chopt -Ep -Wl,-Z 1000000 -pi line=1000 -DHAVE_LIBNETCDF64


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
	@tarfile=../echam5.f90.`date +%y%m%d`.taz ; tar zcvf $$tarfile \
	`ls */*.f90 */*.[fhc] */*inc */Makefile Makefile.in Makefile run/hjob*`

index:
	-rm -f html/[a-z]*
	util/f2html.pl -f util/fgenrc -d html \
          support modules src include
