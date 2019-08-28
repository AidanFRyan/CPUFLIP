main: TriVec.h* CuGrid.h* CuSolver.h* Vec3.h* Particle.h* CuDriver.cpp Makefile
	g++ -Wall -std=c++11 -mavx2 CuDriver.cpp -o main -fopenmp -lpthread

nonmp: TriVec.h* CuGrid.h* CuSolver.h* Vec3.h* Particle.h* CuDriver.cpp Makefile
	g++ -Wall -std=c++11 -mavx CuDriver.cpp -o nonmp

debug: TriVec.h* CuGrid.h* CuSolver.h* Vec3.h* Particle.h* CuDriver.cpp Makefile
	g++ -g CuDriver.cpp -o main
