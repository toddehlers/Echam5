LIB =	support

SRCS =	util_reshape.c gribex.c util_convert.c util_pbio.c util_sysinfo.c util_system.c util_timer.c

ifeq ($(ARCH), SX)
SRCS =	rtc_sx.o util_reshape.c gribex.c util_convert.c util_pbio.c util_sysinfo.c util_system.c util_timer.c
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

clean:
	rm -f $(OBJS)

ifeq ($(ARCH), SX)
gribex.o: gribex.c
	$(CC) $(CFLAGS) -pvctl,fullmsg,noassume,loopcnt=1000000 -Orestrict=all -Onooverlap -c gribex.c
rtc_sx.o: rtc_sx.s
	$(AS) -c rtc_sx.s
endif

ifeq ($(ARCH), ES)
gribex.o: gribex.c
	$(CC) $(CFLAGS) -pvctl,fullmsg,noassume,loopcnt=1000000 -Orestrict=all -Onooverlap -c gribex.c
rtc_sx.o: rtc_sx.s
	$(AS) -c rtc_sx.s
endif
