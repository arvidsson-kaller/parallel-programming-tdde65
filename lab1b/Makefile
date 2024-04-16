CC = mpicc
CFLAGS = 
LFLAGS = -lpthread -lrt -lm

all: blurc thresc

clean:
	-$(RM) *.o blurc thresc ./out/*

blurc: ppmio.o gaussw.o blurfiltermpi.o mpiblurmain.o
	$(CC) -o $@ ppmio.o gaussw.o blurfiltermpi.o mpiblurmain.o $(LFLAGS)

thresc: mpithresmain.o ppmio.o thresfiltermpi.o
	$(CC) -o $@ mpithresmain.o ppmio.o thresfiltermpi.o $(LFLAGS)

arc:
	tar cf - *.c *.h Makefile data/* | gzip - > filters.tar.gz
