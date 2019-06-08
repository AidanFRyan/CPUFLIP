#ifndef KERNELS_H
#define KERNELS_H
#include "CuGrid.h"

template <typename T>
__global__ void initializeGrid(CuGrid<T>* sourceGrid){
	
	int index = threadIdx.x + blockIdx.x*blockDim.x;
	int offset = gridDim.x*blockDim.x;
	int x = sourceGrid->x;
	int y = sourceGrid->y;
	int z = sourceGrid->z;
	int maxGrid = x*y*z;
	curandState_t *state = new curandState_t;
	int parts = 8;
	double dx = sourceGrid->dx;
	curand_init(clock64(), index, 0, state);
	for(int i = index; i < maxGrid; i+=offset){
		if(sourceGrid->a.a[i].t == FLUID){
			int myz = i/(x*y), myy = (i%(x*y))/x, myx = (i%(x*y))%x;
			// curList = sourceGrid->a.a[i].particles = new Particle<T>[parts];
			sourceGrid->a.a[i].numParticles = parts;
			for(int j = 0; j < parts; j++){
				// printf("%d\n", j);
				sourceGrid->a.a[i].particles[j].p.x = myx*dx + dx*curand_uniform(state);
				sourceGrid->a.a[i].particles[j].p.y = myy*dx + dx*curand_uniform(state);
				sourceGrid->a.a[i].particles[j].p.z = myz*dx + dx*curand_uniform(state);
				sourceGrid->a.a[i].particles[j].v.x = 0;//curand_uniform(state);
				sourceGrid->a.a[i].particles[j].v.y = 0;//curand_uniform(state);
				sourceGrid->a.a[i].particles[j].v.z = 0;//curand_uniform(state);
			}
		}
		else sourceGrid->a.a[i].numParticles = 0;
	}
	delete state;
}

template <typename T>
__global__ void solveForFrame(CuGrid<T>* sourceGrid, CuGrid<T>* targetGrid, CuGrid<T>** gridUpdate, int frames) {
    cg::grid_group g = cg::this_grid();
    g.sync();
	int x = sourceGrid->x;
	int y = sourceGrid->y;
	int z = sourceGrid->z;
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	int offset = gridDim.x*blockDim.x;
	int maxGridIndex = x*y*z;
	double delT = sourceGrid->dt;
    g.sync();
	targetGrid->a.copyFrom(sourceGrid->a);                      //initial copy from source to target
    g.sync();
	for(int curFrame = 0; curFrame < frames; ++curFrame){       //for each frame
        if(index == 0)	printf("Frame %d\n", curFrame);
		for(double tElapsed = 0; delT - tElapsed > 0.0001;){    //loop subframes until we reach dt of the frame

            for(int i = index; i < maxGridIndex; i+=offset){
                targetGrid->a.a[i].aDiag = 0;
				targetGrid->a.a[i].aX = 0;
				targetGrid->a.a[i].aY = 0;
				targetGrid->a.a[i].aZ = 0;
				targetGrid->a.a[i].divU = 0;
                targetGrid->a.a[i].p = 0;
                targetGrid->a.a[i].pold = 0;
                targetGrid->a.a[i].u.x = targetGrid->a.a[i].u.y = targetGrid->a.a[i].u.z = 0;
                targetGrid->a.a[i].uOld = Vec3<T>();
            }

			if(index == 0){
				targetGrid->a.dx = sourceGrid->dx;              //copy dx (not copied in voxel copy copyFrom)
				targetGrid->a.dt = targetGrid->dt = sourceGrid->dt = delT/10.0; //set dt for subframe
				if(targetGrid->dt + tElapsed > delT)            //correct subframe if overshooting frame duration
					targetGrid->a.dt = targetGrid->dt = sourceGrid->dt = delT-tElapsed;
				targetGrid->a.density = sourceGrid->density;    //copy density (not copied in copyFrom)
			}

            g.sync();
			targetGrid->a.interpU();                            //interpolate from particle velocity to grid
            g.sync();
			
			for(int i = index; i < maxGridIndex; i += offset){  //this loop needs to be turned into a function of TriVec, calculate aDiag, aX, aY, aZ for each fluid cell
                g.sync();
                                                                //Zero out all coefficients, all velocities, all pressures, and divergences

				if(targetGrid->a.a[i].t == FLUID){                                          //only calc divU and apply body forces on cells which are fluid
                    int myz = i/(x*y), myy = (i%(x*y))/x, myx = (i%(x*y))%x;
                    targetGrid->a.a[i].applyBodyForces(targetGrid->dt);                     //apply body force acceleration (gravity)
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
            g.sync();
			targetGrid->a.maxResidual();                                        //single threaded search for max residual value (max divU)
            g.sync();
			for(int i = 0; i < 1000 && targetGrid->a.maxRes > 0.00001; i++){    //loop, calculating new P based on parallel jacobi
                targetGrid->a.multiThreadJacobi(g);                             //parallel jacobi
				targetGrid->a.calcResidualGS(g);                                //parallel residual calc
				targetGrid->a.maxResidual();                                    //single threaded residual search
                g.sync();
                if(index == 0)
                    if(i == 999)
                        printf("Warning: max iterations reached. Residual: %f\n", targetGrid->a.maxRes);
			}
			if(index == 0) printf("update U\n");
			targetGrid->a.updateU();                                            //update voxel U with pressure term
            g.sync();
			if(index == 0) printf("apply U\n");
			targetGrid->a.applyU();                                             //interp particle U from voxel U
            g.sync();
			if(index == 0) printf("advect particles\n");
			targetGrid->a.advectParticles(g);                                   //move particles based on dt and velocity
			tElapsed += targetGrid->dt;                                         //add dt to total time elapsed
		}
        if(index == 0) printf("copying to source\n");
        sourceGrid->a.copyFrom(targetGrid->a);                                  //copy from target to source TriVec
        g.sync();
        if(index == 0) printf("finished copying\n");			
		while(*gridUpdate != nullptr){g.sync();}	                            //hold while cpu copies
        g.sync();
		if(index == 0){
			*gridUpdate = sourceGrid;
		}
        g.sync();
	}
}
#endif