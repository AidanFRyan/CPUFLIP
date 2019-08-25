#include "TriVec.h"
using namespace std;

// using namespace Cu;

template <typename T>
TriVec<T>::TriVec() {
	size = 0;
	a = nullptr;
	x = y = z = 0;
	invalid = Voxel<T>(true);
	maxRes = 0;
	density = 0;
	dist = std::uniform_real_distribution<double>(0.0, 1.0);
}

template <typename T>
TriVec<T>::TriVec(unsigned long long x, unsigned long long y, unsigned long long z, double dx, double dt, double density, bool dev) : TriVec() {
	this->size = x*y*z;
	if(!dev){
		this->a = new Voxel<T>[size];
	}
	else{
		this->a = nullptr;
	}
	this->dx = dx;
	this->x = x;
	this->y = y;
	this->z = z;
	invalid = Voxel<T>(true);
	this->density = density;
	this->dt = dt;
}

template <typename T>
TriVec<T>::~TriVec() {
	if(a != nullptr){
		delete[] a;
		a = nullptr;
	}
}

template <typename T>
const TriVec<T>& TriVec<T>::operator=(const TriVec<T>& in){
	size = in.size;
	a = in.a;

	x = in.x;
	y = in.y;
	z = in.z;

	dx = in.dx;

	return *this;
}

template <typename T>
const TriVec<T>& TriVec<T>::operator=(TriVec<T>&& in){	//move assignment op, useful for copying when you need to preserve x,y,z and dont actually need a copy, especially when instantiating member TriVecs
	size = in.size;
	a = in.a;
	in.a = nullptr;
	
	x = in.x;
	y = in.y;
	z = in.z;

	dx = in.dx;

	return *this;
}

template <typename T>
T TriVec<T>::pd(const T& l, const T& v, const T& r){	//linterp cell bounds for now, eventually going for cubic interp
	T l_2 = (l + v)/2, r_2 = (r + v)/2; 
	return (r_2-l_2)/dx;
}

template <typename T>
Voxel<T>& TriVec<T>::get(int x, int y, int z){
	if(x >= 0 && x < this->x && y >= 0 && y < this->y && z >= 0 && z < this->z)
		return a[z*this->x*this->y + y*this->x + x];
	// printf("invalid\n");
	return invalid;
}

template <typename T>
 Vec3<T> TriVec<T>::backsolveU(T x, T y, T z, T dt){	//Remnant from long ago
	Vec3<T> t = get(x, y, z).u;
	t = interpU(x+0.5f*dt*t.x, y+0.5f*dt*t.y, z+0.5f*dt*t.z);
	return t;
}

template <typename T>
T TriVec<T>::kWeight(Vec3<T> x){
	T tx = abs(x.x/dx), ty = abs(x.y/dx), tz = abs(x.z/dx);

	// printf("%f %f %f\n", tx,ty,tz);
	if(tx >= 0 && tx < 0.5)
		tx = .75-tx*tx;
	else if(tx >= 0.5 && tx < 1.5)
		tx = 0.5*(1.5-tx)*(1.5-tx);
	else tx = 0;

	if(ty >= 0 && ty < 0.5)
		ty = .75-ty*ty;
	else if(ty >= 0.5 && ty < 1.5)
		ty = 0.5*(1.5-ty)*(1.5-ty);
	else ty = 0;

	if(tz >= 0 && tz < 0.5)
		tz = .75-tz*tz;
	else if(tz >= 0.5 && tz < 1.5)
		tz = 0.5*(1.5-tz)*(1.5-tz);
	else tz = 0;

	return tx*ty*tz;
}

template <typename T>
void TriVec<T>::interpU(){	//interpolate from particle velocities to grid
	for(int l = 0; l < size; ++l){
		if(a[l].t != SOLID){
			int k = l/(x*y), j = (l%(x*y))/x, i = (l%(x*y))%x;
			T offset = dx/2;
			T t, x = i*dx + offset, y = j*dx + offset, z = k*dx + offset;
			Vec3<T> sumk, sumuk, curPosx(x+offset, y, z), curPosy(x, y+offset, z), curPosz(x, y, z+offset);
			for(int p = 0; p < a[l].numParticles; p++){
				a[l].particles[p].vOld = a[l].particles[p].v;
			}
			for(int ioff = -1; ioff < 2; ioff++){ //loop on x-1, x, x+1 cells
				for(int joff = -1; joff < 2; joff++){	//loop on y-1, y, and y+1 cells
					for(int koff = -1; koff < 2; koff++){	//loop on z-1, z, and z+1 cells
							for(int p = 0; p < get(i+ioff, j+joff, k+koff).numParticles; p++){	//for each particle in cell
								t = kWeight(get(i+ioff, j+joff, k+koff).particles[p].p - curPosx);
								sumk.x += t;
								sumuk.x += get(i+ioff, j+joff, k+koff).particles[p].v.x*t;
								
								t = kWeight(get(i+ioff, j+joff, k+koff).particles[p].p - curPosy);
								sumk.y += t;
								sumuk.y += get(i+ioff, j+joff, k+koff).particles[p].v.y*t;

								t = kWeight(get(i+ioff, j+joff, k+koff).particles[p].p - curPosz);
								sumk.z += t;
								sumuk.z += get(i+ioff, j+joff, k+koff).particles[p].v.z*t;
								// printf("%d %d %d | %f %f %f\n", ioff, joff, koff, p.v.x, p.v.y, p.v.z);
							}
					}
				}
			}
			sumuk.x /= (sumk.x + 0.0000000000000001);
			sumuk.y /= (sumk.y + 0.0000000000000001);
			sumuk.z /= (sumk.z + 0.0000000000000001);
			// printf("%d | %f %f %f\n", l, sumuk.x, sumuk.y, sumuk.z);
			if(get(i+1, j, k).t == SOLID){
				sumuk.x = 0;
			}
			if(get(i, j+1, k).t == SOLID){
				sumuk.y = 0;
			}
			if(get(i, j, k+1).t == SOLID){
				sumuk.z = 0;
			}
			a[l].u = sumuk;
		
		}
		a[l].uOld = a[l].u;
	}
}

template <typename T>
T linterp(T a, T p0, T p1){
	return (1-a) * p0 + a*p1;
}

template <typename T>
T trilinterp(T aX, T aY, T aZ, T p000, T p001,  T p010, T p011, T p100, T p101, T p110, T p111){
	return linterp(aZ, linterp(aY, linterp(aX, p000, p001), linterp(aX, p010, p011)), linterp(aY, linterp(aX, p100, p101), linterp(aX, p110, p111)));
}

#define ALPHA 0.97										//defines amount of FLIP to use, 1 is all FLIP, 0 is no FLIP
template <typename T>
void TriVec<T>::interpUtoP(Particle<T>& in){		//interpolate surrounding grid velocities to particles

	Vec3<T> newU, oldU;
	int tx = in.p.x/dx, ty = in.p.y/dx, tz = in.p.z/dx;		//get xyz of particle's voxel
	if(in.p.y - ty*dx < dx/2){								//x component of velocity is stored on upper x bound and halfway point in y and z, need to adjust whether y and z are from y/z and y+1/z+1 or y-1/z-1 and y/z for trilin interp
		--ty;
	}
	if(in.p.z - tz*dx < dx/2){
		--tz;
	}
	//trilinear interp of velocities of 8 cells whose x-velocities surround particle, and position scale for the linterps
	//x scale is distance from particle to negative x face of containing voxel
	//y scale is distance from particle to y component of midpoint of "left" voxel
	//z scale is distance from particle to z component of midpoint of "left" voxel
		
	newU.x = trilinterp((in.p.x - tx*dx) / dx, (in.p.y - (ty*dx + dx/2))/dx, (in.p.z - (tz*dx + dx/2))/dx, get(tx-1, ty, tz).u.x, get(tx, ty, tz).u.x, get(tx-1, ty+1, tz).u.x, get(tx, ty+1, tz).u.x, get(tx-1, ty, tz+1).u.x, get(tx, ty, tz+1).u.x, get(tx-1, ty+1, tz+1).u.x, get(tx, ty+1, tz+1).u.x);
	oldU.x = trilinterp((in.p.x - tx*dx) / dx, (in.p.y - (ty*dx + dx/2))/dx, (in.p.z - (tz*dx + dx/2))/dx, get(tx-1, ty, tz).uOld.x, get(tx, ty, tz).uOld.x, get(tx-1, ty+1, tz).uOld.x, get(tx, ty+1, tz).uOld.x, get(tx-1, ty, tz+1).uOld.x, get(tx, ty, tz+1).uOld.x, get(tx-1, ty+1, tz+1).uOld.x, get(tx, ty+1, tz+1).uOld.x);
	
	tx = in.p.x/dx, ty = in.p.y/dx, tz = in.p.z/dx;
	if(in.p.x - tx*dx < dx/2){
		--tx;
	}
	if(in.p.z - tz*dx < dx/2){
		--tz;
	}

	newU.y = trilinterp((in.p.x - (tx*dx + dx/2))/dx, (in.p.y - ty*dx)/dx, (in.p.z - (tz*dx + dx/2))/dx, get(tx, ty-1, tz).u.y, get(tx+1, ty-1, tz).u.y, get(tx, ty, tz).u.y, get(tx+1, ty, tz).u.y, get(tx, ty-1, tz+1).u.y, get(tx+1, ty-1, tz+1).u.y, get(tx, ty, tz+1).u.y, get(tx+1, ty, tz+1).u.y);
	oldU.y = trilinterp((in.p.x - (tx*dx + dx/2))/dx, (in.p.y - ty*dx)/dx, (in.p.z - (tz*dx + dx/2))/dx, get(tx, ty-1, tz).uOld.y, get(tx+1, ty-1, tz).uOld.y, get(tx, ty, tz).uOld.y, get(tx+1, ty, tz).uOld.y, get(tx, ty-1, tz+1).uOld.y, get(tx+1, ty-1, tz+1).uOld.y, get(tx, ty, tz+1).uOld.y, get(tx+1, ty, tz+1).uOld.y);

	tx = in.p.x/dx, ty = in.p.y/dx, tz = in.p.z/dx;
	if(in.p.x - tx*dx < dx/2){
		--tx;
	}
	if(in.p.y - ty*dx < dx/2){
		--ty;
	}
	newU.z = trilinterp((in.p.x - (tx*dx + dx/2))/dx, (in.p.y - (ty*dx + dx/2))/dx, (in.p.z - tz*dx)/dx, get(tx, ty, tz-1).u.z, get(tx+1, ty, tz-1).u.z, get(tx, ty+1, tz-1).u.z, get(tx+1, ty+1, tz-1).u.z, get(tx, ty, tz).u.z, get(tx+1, ty, tz).u.z, get(tx, ty+1, tz).u.z, get(tx+1, ty+1, tz).u.z);
	oldU.z = trilinterp((in.p.x - (tx*dx + dx/2))/dx, (in.p.y - (ty*dx + dx/2))/dx, (in.p.z - tz*dx)/dx, get(tx, ty, tz-1).uOld.z, get(tx+1, ty, tz-1).uOld.z, get(tx, ty+1, tz-1).uOld.z, get(tx+1, ty+1, tz-1).uOld.z, get(tx, ty, tz).uOld.z, get(tx+1, ty, tz).uOld.z, get(tx, ty+1, tz).uOld.z, get(tx+1, ty+1, tz).uOld.z);
	
	// tx = in.p.x/dx, ty = in.p.y/dx, tz = in.p.z/dx;

	// printf("old %d %d %d: %f %f %f\n", tx, ty, tz, in.v.x, in.v.y, in.v.z);
	// printf("flip %f %f %f | old %f %f %f | new %f %f %f\n", in.vOld.x + (newU.x - oldU.x), in.vOld.y + (newU.y - oldU.y), in.vOld.z + (newU.z - oldU.z), oldU.x, oldU.y, oldU.z, newU.x, newU.y, newU.z);
	in.v = (1-ALPHA)*newU + ALPHA*(in.v + (newU - oldU));
	
	// if(in.v.z == 0){
	// 	tx = in.p.x/dx, ty = in.p.y/dx, tz = in.p.z/dx;
	// 		if(in.p.y - ty*dx < dx/2){								//x component of velocity is stored on upper x bound and halfway point in y and z, need to adjust whether y and z are from y/z and y+1/z+1 or y-1/z-1 and y/z for trilin interp
	// 		--ty;
	// 	}
	// 	if(in.p.z - tz*dx < dx/2){
	// 		--tz;
	// 	}
	// 	// printf("%d %d %d | ", tx, ty, tz);
	// 	tx = in.p.x/dx, ty = in.p.y/dx, tz = in.p.z/dx;
	// 	if(in.p.x - tx*dx < dx/2){
	// 		--tx;
	// 	}
	// 	if(in.p.z - tz*dx < dx/2){
	// 		--tz;
	// 	}
	// 	// printf("%d %d %d | ", tx, ty, tz);
	// 	tx = in.p.x/dx, ty = in.p.y/dx, tz = in.p.z/dx;
	// 	if(in.p.x - tx*dx < dx/2){
	// 		--tx;
	// 	}
	// 	if(in.p.y - ty*dx < dx/2){
	// 		--ty;
	// 	}
	// 	// printf("%d %d %d\n", tx, ty, tz);
		
	// }

	// printf("new %d %d %d: %f %f %f\n", tx, ty, tz, in.v.x, in.v.y, in.v.z);
}

template <typename T>
void TriVec<T>::applyU(){											//apply grid U to particles
	int offset = 1;											//for each fluid voxel, run interpUtoP for each particle inside
	for(int l = 0; l < size; l+=offset){
		if(a[l].t == FLUID){			
			for(int p = 0; p < a[l].numParticles; p++){
				interpUtoP(a[l].particles[p]);
			}
		}
	}
}

template <typename T>
T TriVec<T>::divergenceU(int x, int y, int z){		//calc right hand side, modified to only allow for non-solid surfaces
	Vec3<T> t = get(x, y, z).u;
	double scale = -1/dx, div = 0;						//updated term in case units for pressure were wrong
	if(get(x-1, y, z).t != SOLID)
		div -= get(x-1, y, z).u.x;
	if(get(x+1, y, z).t != SOLID)
		div += get(x, y, z).u.x;
	if(get(x, y-1, z).t != SOLID)
		div -= get(x, y-1, z).u.y;
	if(get(x, y+1, z).t != SOLID)
		div += get(x, y, z).u.y;
	if(get(x, y, z-1).t != SOLID)
		div -= get(x, y, z-1).u.z;
	if(get(x, y, z+1).t != SOLID)
		div += get(x, y, z).u.z;

	div *= scale;
	return div;
}

template <typename T>
void TriVec<T>::updateU(){		//update U of fluid cells with pressure difference between it and it's positive neighbour (u,v,w on positive face of voxel)
	int offset = 1;
	T scale = dt / (dx*density);
	for(int i = 0; i < size; i += offset){
		if(a[i].t != SOLID){
		//if(a[i].t == FLUID){
			int tz = i/(x*y), ty = (i%(x*y))/x, tx = (i%(x*y))%x;
			T dpx, dpy, dpz;
			if(get(tx+1, ty, tz).t != SOLID)
			//if(get(tx+1, ty, tz).t == FLUID)
				dpx = get(tx+1, ty, tz).p - get(tx, ty, tz).p;
			else dpx = 0;

			if(get(tx, ty+1, tz).t != SOLID)
			//if (get(tx, ty+1, tz).t == FLUID)
				dpy = get(tx, ty+1, tz).p - get(tx, ty, tz).p;
			else dpy = 0;

			if(get(tx, ty, tz+1).t != SOLID)
			//if (get(tx, ty, tz+1).t == FLUID)
				dpz = get(tx, ty, tz+1).p - get(tx, ty, tz).p;
			else dpz = 0;

			Vec3<T> t(dpx, dpy, dpz);	//calculate pressure diffs in x,y,z
			t = t*scale;	
			a[i].u = a[i].u - t;																				//subtract pressure-derived term per equation 43
		}
	}
}

template <typename T>
void TriVec<T>::advectParticles(){													//advect particles!
	int offset = 1;
	for(int i = 0; i < size; i += offset){									//loop through each block, should change this to search for fluid so each thread hits at least one fluid
		if(a[i].t == FLUID){
			for(int p = 0; p < a[i].numParticles; p++){															//for each particle in fluid voxel
				// a[i].particles[p].p += a[i].particles[p].v * dt;												//update position by dt*v

				//Bridson's Interpretation of Ralston's advection:
				Particle<T> k2 = a[i].particles[p], k3;
				k2.p += k2.v*(dt/2.0);
				interpUtoP(k2);
				k3 = k2;
				k3.p += k3.v*dt*(3.0/4.0);
				interpUtoP(k3);
				a[i].particles[p].p += a[i].particles[p].v*dt*(2.0/9.0) + k2.v*dt*(3.0/9.0) + k3.v*dt*(4.0/9.0);
			}
		}
	}
														//after updating all positions, search for particles switching voxels
	for(int i = 0; i < size; i++){
		if (a[i].t == FLUID) {
			int tz = i / (x*y), ty = (i % (x*y)) / x, tx = (i % (x*y)) % x;											//get xyz of voxel
			for (int p = 0; p < a[i].numParticles; p++) {														//for all particles in this voxel
				if ((int)floor(a[i].particles[p].p.x / dx) != tx || (int)floor(a[i].particles[p].p.y / dx) != ty || (int)floor(a[i].particles[p].p.z / dx) != tz) {	//if particle is outside of this voxel
					int nx = floor(a[i].particles[p].p.x / dx), ny = floor(a[i].particles[p].p.y / dx), nz = floor(a[i].particles[p].p.z / dx);	//get voxel xyz of where particle should be
					if (get(nx, ny, nz).t == SOLID) {															//if this voxel is solid, reverse particle direction (need to implement LSG to project it out)
						// a[i].particles[p].p = a[i].particles[p].p - dt*a[i].particles[p].v;
						reinsertToFluid(a[i].particles[p]);
						--a[i].numParticles;
						for (int l = p; l < a[i].numParticles; l++) {
							a[i].particles[l] = a[i].particles[l + 1];
						}
						a[i].particles.pop_back();
					}
					else {
						if (get(nx, ny, nz).t == EMPTY)													//if voxel is empty make it a fluid
							get(nx, ny, nz).t = FLUID;
						get(nx, ny, nz).particles.push_back(a[i].particles[p]);
						++get(nx, ny, nz).numParticles;

						--a[i].numParticles;																//if there isn't space, still remove the particle from this voxel
						for (int l = p; l < a[i].numParticles; l++) {
							a[i].particles[l] = a[i].particles[l + 1];										//and shift the remaining particles back in the "stack"
						}
						a[i].particles.pop_back();
					}
				}
			}

			if (a[i].numParticles == 0)																		//if current voxel has no particles (must have been a fluid at start) then set voxel to empty
				a[i].t = EMPTY;
		}
	}			
}

template <typename T>
void TriVec<T>::findPureFluids(){
	for(int i = 0; i < size; ++i){
		if(a[i].t == FLUID){
			int tz = i/(x*y), ty = (i%(x*y))/x, tx = (i%(x*y))%x;
			if(get(tx+1, ty, tz).t == FLUID && get(tx-1, ty, tz).t == FLUID)
				if(get(tx, ty+1, tz).t == FLUID && get(tx, ty-1, tz).t == FLUID)
					if(get(tx, ty, tz+1).t == FLUID && get(tx, ty, tz-1).t == FLUID)
						pureFluidList.push_back(i);
		}
	}
	if(!pureFluidList.size()){
		for(int i = 0; i < size; ++i){
			if(a[i].t == FLUID){
				pureFluidList.push_back(i);
			}
		}
	}
}

template <typename T>
void TriVec<T>::reinsertToFluid(Particle<T>& in){
	int i = pureFluidList[dist(gen)*pureFluidList.size()];
	int tz = i/(x*y), ty = (i%(x*y))/x, tx = (i%(x*y))%x;
	in.p.x = tx+dist(gen)*dx;
	in.p.y = ty+dist(gen)*dx;
	in.p.z = tz+dist(gen)*dx;
	interpUtoP(in);	
	a[i].particles.push_back(in);
	++a[i].numParticles;
}

template <typename T>
Vec3<T> TriVec<T>::negA(int x, int y, int z){
	return Vec3<T>(get(x-1, y, z).aX, get(x, y-1, z).aY, get(x, y, z-1).aZ);
}

template <typename T>
void TriVec<T>::singleThreadGS(){	//domain decomp needed for larger than 10x10x10 grid
	for(int i = 0; i < size; i++){
		if(a[i].t == FLUID){
			int tz = i/(x*y), ty = (i%(x*y))/x, tx = (i%(x*y))%x;
			Vec3<T> t = negA(tx, ty, tz);
			a[i].p = (a[i].divU - (a[i].aX*get(tx+1, ty, tz).pold + a[i].aY*get(tx, ty+1, tz).pold + a[i].aZ*get(tx, ty, tz+1).pold + t.x*get(tx-1, ty, tz).p + t.y*get(tx, ty-1, tz).p + t.z*get(tx, ty, tz-1).p))/(a[i].aDiag+0.000000000000001);
		}
		else a[i].p = 0;
	}
}

template <typename T>
void TriVec<T>::calcResidualGS(){
	int offset = 1;
	for(int i = 0; i < size; i += offset){
		if(a[i].t == FLUID){
			int tz = i/(x*y), ty = (i%(x*y))/x, tx = (i%(x*y))%x;						//get x, y, z coords of voxel. Calculate residual of pressure step for each fluid Voxel
			
			Vec3<T> t = negA(tx, ty, tz);
			a[i].res = a[i].divU - (a[i].aDiag*a[i].p + (a[i].aX*get(tx+1, ty, tz).p + a[i].aY*get(tx, ty+1, tz).p + a[i].aZ*get(tx, ty, tz+1).p + t.x*get(tx-1, ty, tz).p + t.y*get(tx, ty-1, tz).p + t.z*get(tx, ty, tz-1).p));
			a[i].pold = a[i].p;
		}
		else a[i].res = 0;
	}
}

//Jacobi method is only truly parallelizable iterative solver, GS and SOR need to be decomposed/use multigrid to be used in parallel

template <typename T>
void TriVec<T>::multiThreadJacobi(){
	// long long numFluid = 0, index = threadIdx.x + blockDim.x*blockIdx.x;
	int offset = 1;
	for(long long i = 0; i < size; i += offset){									//each thread searches for a fluid cell so all run into a fluid unless there are more threads than fluid cells
		if(a[i].t == FLUID){
			// ++numFluid;
			// if(numFluid % index == 0){
				int tz = i/(x*y), ty = (i%(x*y))/x, tx = (i%(x*y))%x;
				Vec3<T> t = negA(tx, ty, tz);
				a[i].p = (a[i].divU - (a[i].aX*get(tx+1, ty, tz).pold + a[i].aY*get(tx, ty+1, tz).pold + a[i].aZ*get(tx, ty, tz+1).pold + t.x*get(tx-1, ty, tz).pold + t.y*get(tx, ty-1, tz).pold + t.z*get(tx, ty, tz-1).pold))/(a[i].aDiag+0.000000000000001);
				// a[i].p = 0;
		}
	}
}

template <typename T>
void TriVec<T>::maxResidual(){
	maxRes = 0;
	for(int i = 0; i < size; i++){
		if(abs(this->a[i].res) > maxRes){
			maxRes = abs(this->a[i].res);
		}
	}
}

template <typename T>
void TriVec<T>::maxU(){
	T curMax = 0;
	for(int i = 0; i < size; i++){
		if(a[i].t == FLUID){
			for(int p = 0; p < a[i].numParticles; p++){
				if(abs(a[i].particles[p].v.x) > curMax){
					curMax = abs(a[i].particles[p].v.x);
				}
				if(abs(a[i].particles[p].v.y) > curMax){
					curMax = abs(a[i].particles[p].v.y);
				}
				if(abs(a[i].particles[p].v.z) > curMax){
					curMax = abs(a[i].particles[p].v.z);
				}
			}
		}
	}
	mU = curMax;
}

template <typename T>
void TriVec<T>::copyFrom(TriVec<T>& in){
	int offset = 1;
	for(int i = 0; i < in.size; i+=offset){
		a[i] = in.a[i];
	}
}