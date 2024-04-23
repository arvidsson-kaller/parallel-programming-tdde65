all: laplsolvomp laplsolv

laplsolv: laplsolv.c
	g++ -std=c11 -fopenmp $^ -o $@

laplsolvomp: laplsolvomp.c
	g++ -std=c11 -fopenmp $^ -o $@

compare:
	./laplsolvomp 123 100000 0.00001 > omp.txt
	./laplsolv 123 100000 0.00001 > default.txt
	diff ./omp.txt ./default.txt || true

clean:
	-$(RM) laplsolv laplsolvomp