#include "CuSolver.h"
#define NUMBLOCKS 1
#define NUMTHREADS 1
using namespace Cu;


template <typename T>
CuSolver<T>::CuSolver(int x, int y, int z, int frames, double dx, double t) {
	this->x = x;
	this->y = y;
	this->z = z;
	this->frames = frames;
	this->dt = t;
	lastSolved = -1;
	sourceGrid = new CuGrid<T>(x, y, z, dx, t, 997.0);
	targetGrid = new CuGrid<T>(x, y, z, dx, t, 997.0);
}

template <typename T>
CuSolver<T>::~CuSolver() {
	for(int i = 0; i < frames; i++){
		if(sourceGrid != nullptr){
			delete sourceGrid;
			sourceGrid = nullptr;
		}
		if(targetGrid != nullptr){
			delete targetGrid;
			targetGrid = nullptr;
		}
	}
}

template <typename T>
void CuSolver<T>::testInit(){
	for(int k = 0; k < z; k++){
		for(int j = 0; j < y; j++){
			for(int i = 0; i < x; i++){
				Voxel<T>* t = &sourceGrid->a.get(i, j, k);
				// if(i == 0 || i == x-1 || j == 0 || j == y-1 || k == 0 || k == z-1){
				if (i < 5 * x / 8 && i > 3 * x / 8){// && j < 5 * y / 8 && j > 3 * y / 8){
					t->t = SOLID;
				}
				// // else if(i < 5*x/8 && i > 3*x/8 && j < 5*y/8 && j > 3*y/8)
				else
				if(j < 5*y/8 && j > 3*y/8)
					t->t = FLUID;
				else{
					t->t = EMPTY;
				}
			}
		}
	}
}

template <typename T>
bool CuSolver<T>::initialValues(CuGrid<T>* initial) {
	if (!frames ||	initial->x*initial->y*initial->z != x * y*z) return false;
	sourceGrid = initial;
	lastSolved = 0;
	return true;
}

template <typename T>
bool CuSolver<T>::advect() {
	sourceGrid->initializeParticles();
	targetGrid->a.copyFrom(sourceGrid->a);
	this->printParts(0);
	double delT = sourceGrid->dt;
	for(int c = 1; c < frames; c++){
		printf("Frame %d\n", c);
		for(double tElapsed = 0; delT - tElapsed > 0.0001;){
			for(int i = 0; i < x*y*z; i++){
				targetGrid->a.a[i].aDiag = 0;
				targetGrid->a.a[i].aX = 0;
				targetGrid->a.a[i].aY = 0;
				targetGrid->a.a[i].aZ = 0;
				targetGrid->a.a[i].divU = 0;
                targetGrid->a.a[i].p = 0;
                targetGrid->a.a[i].pold = 0;
                targetGrid->a.a[i].u = Vec3<T>();
                targetGrid->a.a[i].uOld = Vec3<T>();
				targetGrid->a.pureFluidList.clear();
			}
			targetGrid->a.maxU();
			targetGrid->a.dx = sourceGrid->dx;              //copy dx (not copied in voxel copy copyFrom)
			// targetGrid->a.dt = targetGrid->dt = sourceGrid->dt = sourceGrid->dx / (sourceGrid->a.mU + sqrt(2*9.8*targetGrid->dx)); //set dt for subframe
			targetGrid->a.dt = targetGrid->dt = sourceGrid->dt = delT/4.0;
			if(targetGrid->dt + tElapsed > delT)            //correct subframe if overshooting frame duration
				targetGrid->a.dt = targetGrid->dt = sourceGrid->dt = delT-tElapsed;
			targetGrid->a.density = sourceGrid->density;    //copy density (not copied in copyFrom)
			targetGrid->a.interpU();
			for(int i = 0; i < x*y*z; i++){  //this loop needs to be turned into a function of TriVec, calculate aDiag, aX, aY, aZ for each fluid cell
				// targetGrid->a.a[i].applyBodyForces(targetGrid->dt);
				if(targetGrid->a.a[i].t == FLUID){                                          //only calc divU and apply body forces on cells which are fluid
					int myz = i/(x*y), myy = (i%(x*y))/x, myx = (i%(x*y))%x;
					targetGrid->a.a[i].applyBodyForces(targetGrid->dt);                     //apply body force acceleration (gravity)
					// printf("%f\n", targetGrid->a.a[i].u.z);
					targetGrid->a.a[i].divU = targetGrid->a.divergenceU(myx, myy, myz);     //calculate divergence of cell

					double scale = targetGrid->a.dt/targetGrid->a.density/targetGrid->a.dx/targetGrid->a.dx;             //scale for aDiag, aX, aY, aZ, only sum aDiag if that side is nonSolid, only store aX, aY, aZ if those sides are fluid
					
					
					if(targetGrid->a.get(myx-1, myy, myz).t != SOLID){
						targetGrid->a.a[i].aDiag += scale;
					}
					if(targetGrid->a.get(myx+1, myy, myz).t == FLUID){
						targetGrid->a.a[i].aDiag += scale;
						targetGrid->a.a[i].aX = -scale;
					}
					else if(targetGrid->a.get(myx+1, myy, myz).t == EMPTY){
						targetGrid->a.a[i].aDiag += scale;
					}
					
					if(targetGrid->a.get(myx, myy-1, myz).t != SOLID){
						targetGrid->a.a[i].aDiag += scale;
					}
					if(targetGrid->a.get(myx, myy+1, myz).t == FLUID){
						targetGrid->a.a[i].aDiag += scale;
						targetGrid->a.a[i].aY = -scale;
					}
					else if(targetGrid->a.get(myx, myy+1, myz).t == EMPTY){
						targetGrid->a.a[i].aDiag += scale;
					}

					if(targetGrid->a.get(myx, myy, myz-1).t != SOLID){
						targetGrid->a.a[i].aDiag += scale;
					}
					if(targetGrid->a.get(myx, myy, myz+1).t == FLUID){
						targetGrid->a.a[i].aDiag += scale;
						targetGrid->a.a[i].aZ = -scale;
					}
					else if(targetGrid->a.get(myx, myy, myz+1).t == EMPTY){
						targetGrid->a.a[i].aDiag += scale;
					}

					targetGrid->a.a[i].res = targetGrid->a.a[i].divU;
				}
			}
			targetGrid->a.maxResidual();	                                      //single threaded search for max residual value (max divU)
			double prevRes = 0;
			for(int i = 0; i < 1000 && targetGrid->a.maxRes > 0.000001; i++){    //loop, calculating new P based on parallel jacobi
                // printf("Iteration %d %f\n", i, targetGrid->a.maxRes);
				targetGrid->a.singleThreadGS();                             //Gauss Seidel for CPU
				targetGrid->a.calcResidualGS();                                //residual calc
				targetGrid->a.maxResidual();                                    //single threaded residual search
				if(abs(targetGrid->a.maxRes - prevRes) < 0.0000001){
					printf("not changing\n");
					break;
				}
				// if(i == 999)
				// 	printf("Warning: max iterations reached. Residual: %f\n", targetGrid->a.maxRes);
			}
			printf("update U\n");
			targetGrid->a.updateU(); 
			// printf("Filling gaps...\n");
			// targetGrid->a.fillGaps();                                           //update voxel U with pressure term
			printf("apply U\n");
			targetGrid->a.applyU();                                             //interp particle U from voxel U
			printf("advect particles\n");
			targetGrid->a.findPureFluids();
			targetGrid->a.advectParticles();                                   //move particles based on dt and velocity
			// targetGrid->a.fillGaps();
			tElapsed += targetGrid->dt;       
		}
		sourceGrid->a.copyFrom(targetGrid->a);
		this->printParts(c);
	}

	return true;
}

template <typename T>
bool CuSolver<T>::solve() {
	if (frames <= 0 || lastSolved == -1)
		return false;
	for (int i = lastSolved; i < frames-1; i++) {
		if (!advect(i))	return false;
	}
	return true;
}

template <typename T>
void CuSolver<T>::resize(int x, int y, int z){
	this->x = x;
	this->y = y;
	this->z = z;
}

template <typename T>
void CuSolver<T>::resizeTime(int t) {
	frames = t;
}

template <typename T>
void CuSolver<T>::printGrid(int f){
	printf("\n\nGrid Number %d:\n\n", f);
	sourceGrid->print();
}

template <typename T>
void CuSolver<T>::printParts(int f){
	targetGrid->printParts(to_string(f));
}
