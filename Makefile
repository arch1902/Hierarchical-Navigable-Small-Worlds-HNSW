compile:
	mpic++ -std=c++17 main.cpp -o wtf -fopenmp
run:
	mpirun -np 4 ./wtf 2 8

clean:
	rm wtf