#ifndef _VOLUME_INTEGRATOR_H_
#define _VOLUME_INTEGRATOR_H_

#include "Integrator.h"

/**	Base virtual class for volumetric integrators, 
	used to compute light transport within a 
	participating media. */


template<int D>
class VolumeIntegrator: public Integrator<D>
{
protected:
	Scattering::Level m_components;
public:
	VolumeIntegrator(World<D> *w):Integrator<D>(w),m_components(Scattering::ALL){}

	void set_scattering_components(const Scattering::Level &components){ m_components = components; }

	virtual Spectrum operator()(const Ray<D> &r) const = 0;
	virtual Spectrum operator()(const Intersection<D> &it) const = 0;
	virtual void operator()(const Intersection<D> &it, std::list<RadianceSample> &samples, const Real delta_time = 0., const unsigned int nb_time_samples = 1) const = 0;
	virtual void operator()(const Ray<D> &r, std::list<RadianceSample> &samples, const Real delta_time = 0., const unsigned int nb_time_samples = 1) const = 0;
	virtual void operator()(const VectorN<D> &p, std::list<RadianceSample> &samples, const Real delta_time = 0., const unsigned int nb_time_samples = 1 ) const {};

}; //VolumeIntegrator


#endif //_VOLUME_INTEGRATOR_H_