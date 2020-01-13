#ifndef _SURFACE_INTEGRATOR_H_
#define _SURFACE_INTEGRATOR_H_

#include "Integrator.h"

template<int D>
class SurfaceIntegrator: public Integrator<D>
{

public:
	SurfaceIntegrator(World<D> *w, unsigned int incoming_samples = 1):Integrator<D>(w,incoming_samples){}
	~SurfaceIntegrator(){}

	virtual Spectrum operator()(const Intersection<D> &it) const = 0;
	virtual void operator()(const Intersection<D> &it, std::list<RadianceSample> &samples, const Real delta_time = 0., const unsigned int nb_time_samples = 1 ) const = 0;
	virtual void operator()(const VectorN<D> &p, std::list<RadianceSample> &samples, const Real delta_time = 0., const unsigned int nb_time_samples = 1 ) const {};
};


#endif //_SURFACE_INTEGRATOR_H_