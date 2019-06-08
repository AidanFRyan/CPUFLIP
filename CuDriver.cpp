#include "CuSolver.h"
#define FRAMES 70
#define NUMCELLS 40
#define FPS 24
#define DX 0.025
using namespace Cu;

int main() {
	CuSolver<double> solver(NUMCELLS, NUMCELLS, NUMCELLS, FRAMES, DX, 1.0/FPS);
	solver.testInit();
	solver.advect();
	return 0;
}