/* 
 * Copyright (C) 2015, Adrian Jarabo (http://giga.cps.unizar.es/~ajarabo/)
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom
 * the Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
 * OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef _HOMOGENEOUS_MEDIUM_H_
#define _HOMOGENEOUS_MEDIUM_H_

#include "Media/Medium.h"
#include "Media/PhaseFunction/PhaseFunction.h"
#include "Media/PhaseFunction/IsotropicPhaseFunction.h"

/** Class that implements an isotropic medium. 
	The templatized parameter PF defines the phase function. */
template<unsigned D, class PF>
class HomogeneousMedium: public Medium<D>
{
	Spectrum m_absorption, m_scattering,
			 m_extinction, m_albedo;
	Spectrum m_exp_extinction;
	Real m_mean_extinction;
	PF m_phase_function;
public:
	HomogeneousMedium() :
		m_mean_extinction(0.), m_phase_function()
	{}

	HomogeneousMedium(const Spectrum &absorption, const Spectrum &scattering, PF&& pf) :
		m_absorption(absorption), m_scattering(scattering),
		m_extinction(absorption+scattering), m_albedo(scattering/(absorption+scattering)),
		m_exp_extinction(exp(-m_extinction)), m_mean_extinction(m_extinction.avg()),
		m_phase_function(pf)
	{}

	HomogeneousMedium(const Spectrum &absorption, const Spectrum &scattering, PF& pf) :
		m_absorption(absorption), m_scattering(scattering),
		m_extinction(absorption+scattering), m_albedo(scattering/(absorption+scattering)),
		m_exp_extinction(exp(-m_extinction)), m_mean_extinction(m_extinction.avg()),
		m_phase_function(pf)
	{}

	Real get_mean_free_path(const Ray<D> &r)const
	{
		return m_mean_extinction;
	}

	Real sample_mean_free_path(const Ray<D> &r, Real &pdf ) const
	{
		Real t = -log(Random::StdRNG.next_real()) / m_mean_extinction;
		pdf = exp(-m_mean_extinction*t);

		return t;
	}

	Real p(const Ray<D> &r, const Real d) const
	{
		return exp(-m_mean_extinction*d);
	}

	inline Spectrum get_transmittance(const Ray<D> &r) const
	{
		return exp(-m_extinction*r.get_parameter());
	}

	inline Spectrum get_extinction(const VectorN<D> &p) const
	{
		return m_extinction;
	}

	inline Spectrum get_max_extinction() const
	{
		return m_extinction;
	}

	inline Spectrum get_albedo(const VectorN<D> &p) const
	{
		return m_albedo;
	}

	inline Spectrum get_absorption(const VectorN<D> &p) const
	{
		return m_absorption;
	}

	inline Spectrum get_scattering(const VectorN<D> &p) const
	{
		return m_scattering;
	}
	
	// Probability of sampling a specific direction 
	inline Real p(const VectorN<D> &p, const VectorN<D> &wi, const VectorN<D> &wo) const
	{
		return m_phase_function.p(wi, wo);
	}

	inline Real p(const VectorN<D> &p, const VectorN<D> &wi, const VectorN<D> &wo, const Real delta_time) const
	{
		return m_phase_function.p(wi, wo, delta_time);
	}

	// Scattering functions based on unpolarized spectrum
	inline void f(const VectorN<D>& p, const VectorN<D>& wi, const VectorN<D>& wo, Spectrum &R) const
	{
		m_phase_function(wi, wo, R);
	}

	inline void f(const VectorN<D>& p, const VectorN<D>& wi, const VectorN<D>& wo, Spectrum &R, Real &delta_time, Real &pdf) const
	{
		m_phase_function(wi, wo, R, delta_time, pdf);
	}

	// Scattering functions based on polarized spectrum (Muller Matrices)
	inline void f(const VectorN<D>& p, const VectorN<D>& wi, const VectorN<D>& wo, PolarizedAttenuation<D> &R) const
	{
		m_phase_function(wi, wo, R);
	}

	inline void f(const VectorN<D>& p, const VectorN<D>& wi, const VectorN<D>& wo, PolarizedAttenuation<D> &R, Real &delta_time, Real &pdf) const
	{
		m_phase_function(wi, wo, R, delta_time, pdf);
	}

	inline void f(const VectorN<D>& p, const VectorN<D>& wi, const VectorN<D>& wo, FluorescentAttenuation<D> &R) const
	{
		m_phase_function(wi, wo, R);
	}

	inline void f(const VectorN<D>& p, const VectorN<D>& wi, const VectorN<D>& wo, FluorescentAttenuation<D> &R, Real &delta_time, Real &pdf) const
	{
		m_phase_function(wi, wo, R, delta_time, pdf);
	}

	inline void f(const VectorN<D>& p, const VectorN<D>& wi, const VectorN<D>& wo, FluorescentPolarizedAttenuation<D> &R) const
	{
		m_phase_function(wi, wo, R);
	}

	inline void f(const VectorN<D>& p, const VectorN<D>& wi, const VectorN<D>& wo, FluorescentPolarizedAttenuation<D> &R, Real &delta_time, Real &pdf) const
	{
		m_phase_function(wi, wo, R, delta_time, pdf);
	}

	inline void sample_direction(const VectorN<D> &p, const VectorN<D>& wi, VectorN<D> &wo, Spectrum &R, Real &pdf) const
	{
		m_phase_function.sample_direction(wi, wo, R, pdf); 
	}
    
	inline void sample_direction(const VectorN<D> &p, const VectorN<D>& wi, VectorN<D> &wo, Spectrum &R, Real &delta_time, Real &pdf) const
	{
		m_phase_function.sample_direction(wi, wo, R, delta_time, pdf); 
	}

	inline void sample_direction(const VectorN<D> &p, const VectorN<D>& wi, VectorN<D> &wo, FluorescentAttenuation<D> &R, Real &pdf) const
	{
		m_phase_function.sample_direction(wi, wo, R, pdf);
	}

	inline void sample_direction(const VectorN<D> &p, const VectorN<D>& wi, VectorN<D> &wo, FluorescentAttenuation<D> &R, Real &delta_time, Real &pdf) const
	{
		m_phase_function.sample_direction(wi, wo, R, delta_time, pdf);
	}

	inline void sample_direction(const VectorN<D> &p, const VectorN<D>& wi, VectorN<D> &wo, PolarizedAttenuation<D> &R, Real &pdf) const
	{
		m_phase_function.sample_direction(wi, wo, R, pdf);
	}

	inline void sample_direction(const VectorN<D> &p, const VectorN<D>& wi, VectorN<D> &wo, PolarizedAttenuation<D> &R, Real &delta_time, Real &pdf) const
	{
		m_phase_function.sample_direction(wi, wo, R, delta_time, pdf);
	}

	inline void sample_direction(const VectorN<D> &p, const VectorN<D>& wi, VectorN<D> &wo, FluorescentPolarizedAttenuation<D> &R, Real &pdf) const
	{
		m_phase_function.sample_direction(wi, wo, R, pdf);
	}

	inline void sample_direction(const VectorN<D> &p, const VectorN<D>& wi, VectorN<D> &wo, FluorescentPolarizedAttenuation<D> &R, Real &delta_time, Real &pdf) const
	{
		m_phase_function.sample_direction(wi, wo, R, delta_time, pdf);
	}
}; // HomogeneousMedium

typedef HomogeneousMedium<3, PhaseFunction::Isotropic<3> > HomogeneousIsotropicMedium3D;
typedef HomogeneousMedium<2, PhaseFunction::Isotropic<2> > HomogeneousIsotropicMedium2D;

#endif //_HOMOGENEOUS_MEDIUM_H_
