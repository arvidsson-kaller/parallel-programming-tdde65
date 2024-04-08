CC = icc
CFLAGS = 
LFLAGS = -lpthread -lrt -lm

all: blurc thresc

clean:
	-$(RM) *.o blurc thresc

blurc: ppmio.o gaussw.o blurfilter.o blurfilterpthreads.o blurmain.o
	$(CC) -o $@ ppmio.o gaussw.o blurfilter.o blurfilterpthreads.o blurmain.o $(LFLAGS)

thresc: thresmain.o ppmio.o thresfilter.o thresfilterpthreads.o
	$(CC) -o $@ thresmain.o ppmio.o thresfilter.o thresfilterpthreads.o $(LFLAGS)

arc:
	tar cf - *.c *.h Makefile data/* | gzip - > filters.tar.gz
