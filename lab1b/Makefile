CC = mpicc
CFLAGS = 
LFLAGS = -lpthread -lrt -lm -g

all: blurc thresc

clean:
	-$(RM) *.o blurc thresc ./out/*

blurc: ppmio.o gaussw.o blurfilter.o blurfiltermpi.o mpiblurmain.o
	$(CC) -o $@ ppmio.o gaussw.o blurfilter.o blurfiltermpi.o mpiblurmain.o $(LFLAGS)

thresc: thresmain.o ppmio.o thresfilter.o thresfilterpthreads.o
	$(CC) -o $@ thresmain.o ppmio.o thresfilter.o thresfilterpthreads.o $(LFLAGS)

arc:
	tar cf - *.c *.h Makefile data/* | gzip - > filters.tar.gz
