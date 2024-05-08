
all: sequential parallel

sequential:main_sequential.cpp  physics.c
	g++ -g -o $@.o $^ -lrt -lm

parallel:main_parallel.cpp physics.c
	mpic++ -g -o $@.o $^ -lrt -lm

allclean:
	-rm *.o
