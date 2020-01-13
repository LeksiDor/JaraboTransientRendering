/* 
 * Copyright (C) 2015, Adrian Jarabo (http://giga.cps.unizar.es/~ajarabo/)
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom
 * the Software is furnished to do so, subject to the following conditions:

 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.

 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
 * OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef _INTEGRATOR_H_
#define _INTEGRATOR_H_

#include "bunnykiller.h"

#include <vector>

#include "RayTracing/Intersection.h"
#include "Integrator/RadianceSample.h"
#include "Color/Spectrum.h"
#include "LinearAlgebra/VectorN.h"
#include "Film/Film.h"

template<unsigned D, class Radiance>
class World;

namespace Scattering
{
	enum Level
	{
		NONE = 0,
		SINGLE = 1,
		MULTIPLE = 2,
		ALL = SINGLE | MULTIPLE
	};
};

/**	Base virtual class for all integrators. */
template<unsigned D, class Radiance, class RadianceAttenuation>
class Integrator
{
public:
	using RadianceSampleR = RadianceSample<D, Radiance>;
	using RadSampleList = std::vector<RadianceSampleR>;
	using FilmR = Film<D, Radiance>;
public:
	World<D, Radiance>& m_world;
	unsigned int m_incoming_samples;
	Real m_inv_incoming_samples;

	Scattering::Level m_components;

	FilmR *m_film;
	FILE *f_log;
public:
	Integrator(World<D, Radiance>& w, unsigned int incoming_samples = 1,
			FILE *_f_log = nullptr, FilmR *_m_film = nullptr) :
		m_world(w),
		m_incoming_samples(incoming_samples),
		m_inv_incoming_samples(1. / Real(m_incoming_samples)),
		m_components(Scattering::ALL),
		m_film(_m_film),
		f_log(_f_log)
	{}

	virtual ~Integrator()
	{}

	void set_scattering_components(const Scattering::Level &components)
	{
		m_components = components;
	}

	/*
	 * Preprocess function for the integrator. The default behaviour
	 * is doing nothing. Other integrators (e.g. these based on
	 * particle tracing) would re-implement this function.
	 */
	virtual void preprocess()
	{}
	
	/**	Operator that integrates the radiance incoming to the
	 Intersection it. Note that these functions are not thought to
	 model volumetric scattering; only surface radiance integration.
	 NOTE: Legacy code; this function is deprecated. */
	virtual Radiance operator()(const Intersection<D> &it) const
	{
		return (*this)(it.get_ray());
	}

	/**	Operator that integrates the radiance incoming to the
	 Intersection it, similar to operator(it), but taking
	 into account different time samples. This more adequate
	 for transient rendering.
	 NOTE: Legacy code; this function is deprecated. */
	virtual void operator()(const Intersection<D> &it, RadSampleList &samples,
			const Real delta_time = 0., const unsigned int nb_time_samples = 1,
			const Real offset_time = 0.) const
	{
		(*this)(it.get_ray(), samples, delta_time, nb_time_samples);
	}

	/**	Operator that integrates the radiance incoming from Ray r.
	 This function includes the light incoming from the environment,
	 and also the light scattered by the medium through the ray. */
	virtual Radiance operator()(const Ray<D> &r) const
	{
		return Radiance(0.);
	}
	
	/**	Operator that integrates the radiance incoming to
	 r.origin(), from the direction r.direction(). This is similar
	 to operator(r), but taking into account different time samples.
	 This more adequate for transient rendering. */
	virtual void operator()(const Ray<D> &r, RadSampleList &samples, const Real delta_time = 0.,
			const unsigned int nb_time_samples = 1) const
	{}
}; // Integrator

#endif //_INTEGRATOR_H_
