#ifndef _ALBEDO_INTEGRATOR_H_
#define _ALBEDO_INTEGRATOR_H_

#include "Integrator/SurfaceIntegrator.h"
#include "RayTracing/World.h"

template<class TermCriteria, int D>
class Albedo: public SurfaceIntegrator<D>
{

public:
	Albedo(World<D> *w, const TermCriteria &termination,  unsigned int incoming_samples = 1);
	~Albedo();

	virtual Spectrum operator()(const Intersection<D> &it) const;
	virtual void operator()(const Intersection<D> &it, std::list<RadianceSample> &samples, const Real delta_time = 0., const unsigned int nb_time_samples = 1 ) const;
};


template<class TermCriteria, int D>
Albedo<TermCriteria,D>::Albedo(World<D> *w, const TermCriteria &_termCriteria, unsigned int incoming_samples)
	:SurfaceIntegrator<D>(w,incoming_samples)
{}

template<class TermCriteria, int D>
Albedo<TermCriteria,D>::~Albedo()
{}

template<class TermCriteria, int D>
Spectrum Albedo<TermCriteria,D>::operator()(const Intersection<D> &it) const
{
	if( it.material() )
		return it.material()->get_absorption(it.get_uv());
	else
		return Spectrum();
		
}

template<class TermCriteria, int D>
void Albedo<TermCriteria,D>::operator()(const Intersection<D> &it, std::list<RadianceSample> &samples, const Real delta_time, const unsigned int nb_time_samples) const
{
}


#endif //_PATH_TRACING_INTEGRATOR_H_