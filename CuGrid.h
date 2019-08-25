#ifndef CUGRID_H
#define CUGRID_H
#include <fstream>
#include <string>
#include "TriVec.h"
using namespace Cu;

namespace Cu {

	template <typename T>
	class CuGrid {
	public:
		CuGrid();
		CuGrid(int x, int y, int z, double dx, double dt, double density);
		CuGrid(const CuGrid&);
		~CuGrid();
		// CuGrid* toDeviceCopy();
		// CuGrid* toDeviceEmpty();
		// bool copyFromDevice(CuGrid* d_grid);
		// bool copyFromDeviceAsync(CuGrid* d_grid, cudaStream_t& stream);
		// bool allocateOnDevice();
		// void removeDeviceCopies();
		void initializeParticles();
		void setTstep(double t);
		void print();
		void printParts(string filename);
		const CuGrid& operator=(const CuGrid&);
		TriVec<T> a;
		unsigned long long x, y, z;
		double dx, dt, density;
		
	private:
		void nullifyHostTriVecs();
		void nullifyDeviceTriVecs();
	};

#include "CuGrid.hpp"
}

#endif