#ifndef _PHOTON_H_
#define _PHOTON_H_

#include "Utils/globals.h"
#include "LinearAlgebra/VectorN.h"
#include "Particle.h"

template<int D>
class Photon
{
public:
	unsigned char m_level;
	VectorN<D> m_position;

#ifdef _PHOTON_MAP_OPTIMIZED_
	unsigned char m_phi, m_thita;
	char m_power[4];
	unsigned char m_time;

	Photon():m_level(0),m_phi(0),m_thita(0),m_time(0){m_power[0]=m_power[1]=m_power[2]=m_power[3]=0;}
#else
	VectorN<D> m_direction;
	Spectrum m_power;
	Real m_time;
	Photon():m_level(0),m_time(0){}
	Photon(const VectorN<D> &position, const VectorN<D> &direction, const Spectrum &power, const Real time, const unsigned int level)
		:m_position(position), m_direction(direction), m_power(power), m_time(time), m_level(level){}
	Photon(const Particle<D> &particle)
		:m_position(particle.m_hit_point.get_position()), m_direction(particle.m_hit_point.get_ray().get_direction()), 
		m_power(particle.m_power), m_time(particle.m_time), m_level(particle.m_level) {}
#endif

}; //Photon

#endif //_PHOTON_H_