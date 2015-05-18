LIB =	support

SRCS =	util_reshape.c cdilib.c util_convert.c util_pbio.c util_sysinfo.c util_system.c util_timer.c util_backtrace.c

ifeq ($(ARCH), SX)
SRCS =	rtc_sx.o util_reshape.c cdilib.c util_convert.c util_pbio.c util_sysinfo.c util_system.c util_timer.c util_backtrace.c
endif

OBJS := $(SRCS:.c=.o)

.SUFFIXES:
.SUFFIXES: .c .o
%.o: %.c
	$(CC) $(CFLAGS) -c $<

all: $(LIB)

$(LIB): ../lib/lib$(LIB).a 

../lib/lib$(LIB).a: $(OBJS)
	$(AR) $(ARFLAGS) ../lib/lib$(LIB).a $(OBJS) 

cdilib.o: cdilib.c
	$(CC) $(CFLAGS) $(NETCDF_INCLUDE) -DHAVE_LIBNETCDF -DHAVE_CF_INTERFACE -c cdilib.c

ifeq ($(ARCH), SX)
cdilib.o: cdilib.c
	$(CC) $(CFLAGS) $(NETCDF_INCLUDE) -DHAVE_LIBNETCDF -DHAVE_CF_INTERFACE -pvctl,fullmsg,noassume,loopcnt=1000000 -Orestrict=all -Onooverlap -c cdilib.c
rtc_sx.o: rtc_sx.s
	$(AS) -c rtc_sx.s
endif

ifeq ($(ARCH), ES)
cdilib.o: cdilib.c
	$(CC) $(CFLAGS) $(NETCDF_INCLUDE) -DHAVE_LIBNETCDF -DHAVE_CF_INTERFACE -pvctl,fullmsg,noassume,loopcnt=1000000 -Orestrict=all -Onooverlap -c cdilib.c
rtc_sx.o: rtc_sx.s
	$(AS) -c rtc_sx.s
endif

clean:
	rm -f $(OBJS)
