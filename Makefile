compile:
	mpic++ -std=c++17 -o wtf *.cpp
run:
	mpirun -np 1 ./wtf 
clean:
	rm wtf