#ifndef VOXEL_H
#define VOXEL_H
#include "Particle.h"
using namespace Cu;
namespace Cu{
    enum VoxType {SOLID, FLUID, EMPTY};
	template <typename T>
	class Voxel{
	public:
		Voxel();
		Voxel(bool in);
		Voxel(const Voxel& in);
		const Voxel& operator=(const Voxel& in);
		void applyBodyForces(float dt);
		Vec3<T> u, uOld, f;	//note that f is actually acceleration vector in this version of the code
		T p, pold, divU, res, z, s;	//pressure, divergence of u term (for pressure solver), residual, z vector quantities (negative residual in CG, preconditioned value in PCG)
		VoxType t;
		bool invalid;
		std::vector< Particle<T> > particles;
		// Particle<T> particles[24];
		short numParticles;
		T aDiag, aX, aY, aZ;//, dotTemp;
	};
    #include "Voxel.hpp"
}
#endif