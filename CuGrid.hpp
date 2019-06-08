#include "CuGrid.h"

using namespace Cu;

template <typename T>
CuGrid<T>::CuGrid() {
	x = 0;
	y = 0;
	z = 0;
}

template <typename T>
CuGrid<T>::CuGrid(int x, int y, int z, double dx, double dt, double density) {
	this->x = x;
	this->y = y;
	this->z = z;
	printf("allocing %lu bytes\n", sizeof(Voxel<T>)*x*y*z);
	this->a = std::move(TriVec<T>(x, y, z, dx, dt, density));

	// this->d_a = std::move(TriVec<T>(x, y, z, dx, dt, density, true));

	this->dx = dx;
	this->dt = dt;
	this->density = density;
}

// template <typename T>
// void CuGrid<T>::removeDeviceCopies() {
// 	if(d_a.a != nullptr){
// 		cudaFree(d_a.a);
// 		d_a.a = nullptr;
// 	}
// }

template <typename T>
void CuGrid<T>::setTstep(double t){
	this->dt = t;
	this->a.dt = t;
}

template <typename T>
void CuGrid<T>::nullifyHostTriVecs(){
	a.a = nullptr;
}

template <typename T>
void CuGrid<T>::nullifyDeviceTriVecs(){
	// d_a.a = nullptr;
}

template <typename T>
void CuGrid<T>::printParts(string filename){
	ofstream f;
	f.open(filename);
	for(int i = 0; i < a.size; i++){
		// if(a.a[i].t == FLUID)
			for(int p = 0; p < a.a[i].numParticles; p++){
				f<<a.a[i].particles[p].p.x<<' '<<a.a[i].particles[p].p.y<<' '<<a.a[i].particles[p].p.z<<endl;
			}
	}
	f.close();
}

template <typename T>
CuGrid<T>::CuGrid(const CuGrid<T>& in){
	this->x = in.x;
	this->y = in.y;
	this->z = in.z;
	this->dx = in.dx;
	this->dt = in.dt;
	this->density = in.density;

	this->a = std::move(TriVec<T>(x, y, z, dx, dt, density));

	// this->d_a = std::move(TriVec<T>(x, y, z, dx, dt, density, true));
}

template <typename T>
const CuGrid<T>& CuGrid<T>::operator=(const CuGrid<T>& in){
	a = in.a;
	// d_a = in.d_a;
	x = in.x;	
	y = in.y;
	z = in.z;
	dx = in.dx;
	dt = in.dt;
	density = in.density;
	return *this;
}

// template <typename T>
// CuGrid<T>* CuGrid<T>::toDeviceCopy() {
// 	CuGrid<T> *d_temp, temp = *this;

// 	cudaMemcpy(d_a.a, a.a, sizeof(Voxel<T>)*x*y*z, cudaMemcpyHostToDevice);

// 	temp.a = d_a;

// 	cudaMalloc((void**)&d_temp, sizeof(CuGrid<T>));	//need to cudaFree the return value
// 	cudaMemcpy(d_temp, &temp, sizeof(CuGrid<T>), cudaMemcpyHostToDevice);

// 	temp.nullifyDeviceTriVecs();
// 	temp.nullifyHostTriVecs();
	
// 	return d_temp;
// }

template <typename T>
CuGrid<T>::~CuGrid() {
	// removeDeviceCopies();
}

// template <typename T>
// bool CuGrid<T>::copyFromDevice(CuGrid<T>* d_grid) {
// 	CuGrid<T>* grid = new CuGrid<T>;
// 	cudaMemcpy(grid, d_grid, sizeof(CuGrid<T>), cudaMemcpyDeviceToHost);

// 	if(grid->x != x || grid->y != y || grid->z != z){
// 		delete grid;
// 		return false;
// 	}

// 	cudaMemcpy(this->a.a, grid->a.a, sizeof(Voxel<T>)*x*y*z, cudaMemcpyDeviceToHost);

// 	grid->nullifyHostTriVecs();

// 	delete grid;
// 	return true;
// }

// template <typename T>
// bool CuGrid<T>::copyFromDeviceAsync(CuGrid<T>* d_grid, cudaStream_t& stream) {
// 	CuGrid<T>* grid = new CuGrid<T>;
// 	cudaMemcpyAsync(grid, d_grid, sizeof(CuGrid<T>), cudaMemcpyDeviceToHost, stream);

// 	if(grid->x != x || grid->y != y || grid->z != z){
// 		delete grid;
// 		return false;
// 	}

// 	cudaMemcpyAsync(this->a.a, grid->a.a, sizeof(Voxel<T>)*x*y*z, cudaMemcpyDeviceToHost, stream);

// 	grid->nullifyHostTriVecs();

// 	delete grid;
// 	return true;
// }

// template <typename T>
// bool CuGrid<T>::allocateOnDevice() {
// 	cudaDeviceSynchronize();

// 	cudaMalloc((void**)&d_a.a, sizeof(Voxel<T>)*x*y*z);

// 	cudaDeviceSynchronize();

// 	return true;
// }

// template <typename T>
// CuGrid<T>* CuGrid<T>::toDeviceEmpty() {
// 	CuGrid<T> temp, *d_temp;
// 	temp.a = d_a;
// 	temp.x = x;
// 	temp.y = y;
// 	temp.z = z;

// 	cudaMalloc((void**)&d_temp, sizeof(CuGrid<T>));
// 	cudaMemcpy(d_temp, &temp, sizeof(CuGrid<T>), cudaMemcpyHostToDevice);

// 	temp.nullifyHostTriVecs();

// 	return d_temp;
// }

template <typename T>
void CuGrid<T>::print(){
	for(int k = z-1; k >= 0; k--){
		printf("\nZ = %d\n", k);
		for(int j = y-1; j >= 0; j--){
			for(int i = x-1; i >= 0; i--){
				Voxel<T> *t = &a.get(i, j, k);
				// if(t->t == FLUID)
				// 	printf(" F ");
				// else if(t->t == SOLID)
				// 	printf(" S ");
				// else printf("   ");
				printf("| %f |", t->p);
			}
			printf("\n");
		}
	}
}

template <typename T>
void CuGrid<T>::initializeParticles(){
	srand(time(NULL));
	for(int i = 0; i < x*y*z; i++){
		if(a.a[i].t == FLUID){
			int tz = i/(x*y), ty = (i%(x*y))/x, tx = (i%(x*y))%x;
			a.a[i].numParticles = 8;
			for(int j = 0; j < 8; j++){
				a.a[i].particles[j].p.x = dx*tx + dx*(double)rand()/RAND_MAX;
				a.a[i].particles[j].p.y = dx*ty + dx*(double)rand()/RAND_MAX;
				a.a[i].particles[j].p.z = dx*tz + dx*(double)rand()/RAND_MAX;
				// printf("%d %d\n", i, j);
			}
		}
	}
}