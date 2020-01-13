#ifndef _PARTICLE_H_
#define _PARTICLE_H_

#include "bunnykiller.h"

#include "LinearAlgebra/VectorN.h"
#include "RayTracing/Intersection.h"

/** Struct used to store the particles in the ParticleTracing class. 
	It stores all the information of the hit point, allowing  to know 
	the origin of the particle and its hit point, plus the energy and 
	the path travelled. 
	This allows using the particle for all particle-tracing-based algorithms,
	such as Photon Mapping [Jensen01], Instant Radiosity [Keller97], or 
	Photon Beams [Jarosz11].
	In order to include heterogeneous participating media for the Photon 
	Beams algorithm, some modifications must be performed to include the
	information of the travelled medium. */
template<int D>
struct Particle
{
	unsigned char m_level;
	Spectrum m_power;
	Real m_time;
	Intersection<D> m_hit_point;

	Particle(const Intersection<D> &it, const Spectrum &power, Real time, int level)
		:m_hit_point(it), m_power(power), m_level(level), m_time(time){}
}; //Particle

#endif //_PARTICLE_H_