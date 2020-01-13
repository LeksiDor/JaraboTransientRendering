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

#ifndef _TRANSPARENT_H_
#define _TRANSPARENT_H_

#include "Material/Material.h"
#include "RayTracing/Ray.h"
#include "RayTracing/Intersection.h"

template<unsigned D>
class Transparent : public Material<D>
{
protected:
	typedef Material<D> MTR;

protected:
	// Absorption coefficient
	Spectrum m_absorption;

	// Medium inside and outside the material
	Medium<D> *m_medium_0;
	Medium<D> *m_medium_1;

	//Reflectance::DELTA_TRANSMISSION
public:
	Transparent(const Spectrum &absorption = Spectrum(0.f)) :
		Material<D>(Reflectance::DELTA_TRANSMISSION), m_absorption(absorption),
		m_medium_0(MTR::m_default_medium), m_medium_1(nullptr)
	{}

	Transparent(Medium<D> *m_in, Medium<D> *m_out = nullptr) :
		Material<D>(Reflectance::DELTA_TRANSMISSION), m_absorption(0.),
		m_medium_0(m_out ? m_out : MTR::m_default_medium), m_medium_1(m_in)
	{}

	Transparent(const Spectrum &absorption, Medium<D> *m_in, Medium<D> *m_out = nullptr) :
		Material<D>(Reflectance::DELTA_TRANSMISSION), m_absorption(absorption),
		m_medium_0(m_out ? m_out : MTR::m_default_medium), m_medium_1(m_in)
	{}

	virtual ~Transparent() {}

	/*---------------------------------------------------------------------------*/
	/** PDF */
	Real p(const VectorN<D> &omega_i, const VectorN<D> &omega_o, const VectorN<D> &normal, const Vector2 &uv) const
	{
		return 0.;
	}

	/*---------------------------------------------------------------------------*/
	/** Spectrum-based attenuation */
	void f(const VectorN<D> &omega_i, const VectorN<D> &omega_o, const VectorN<D> &normal, const Vector2 &uv, Spectrum &R) const
	{
		R = Spectrum(0.);
	}

	void sample_direction(const VectorN<D> &omega_i, VectorN<D> &omega_o, const VectorN<D> &normal, const Vector2 &uv, Spectrum &R, Real &pdf) const
	{
		omega_o = omega_i;
		pdf = 1.;
		R = (1. - m_absorption);
	}

	void sample_outgoing_ray(const Intersection<D> &it, Ray<D> &new_ray, Spectrum &R, Real &pdf) const
	{
		// Test if the ray enters the media or if it is leaving
		new_ray = Ray<D>(it.get_position(), it.get_ray().get_direction(), true,
			it.get_ray().get_level() + 1, it.get_ray().get_ior(),
			(dot(it.get_ray().get_direction(), it.get_normal()) > 0.) ? m_medium_0 : m_medium_1);

		R = (1. - m_absorption);
		pdf = 1.;
	}

	/*---------------------------------------------------------------------------*/
	/** Polarized spectrum-based attenuation */
	virtual void f(const VectorN<D> &, const VectorN<D> &, const VectorN<D> &, const Vector2 &, PolarizedAttenuation<D> &R) const
	{
		R = PolarizedAttenuation<D>(0.);
	}

	/** Samples a new outgoing direction */
	virtual void sample_direction(const VectorN<D> &omega_i, VectorN<D> &omega_o, const VectorN<D> &, const Vector2 &, PolarizedAttenuation<D> &R, Real &pdf) const
	{
		omega_o = omega_i;
		pdf = 1.;
		R = PolarizedAttenuation<D>(1. - m_absorption);
	}

	/** Samples a new outgoing ray in the intersection */
	virtual void sample_outgoing_ray(const Intersection<D> &it, Ray<D> &new_ray, PolarizedAttenuation<D> &R, Real &pdf) const
	{
		//Test if the ray enters the media or if it is leaving
		new_ray = Ray<D>(it.get_position(), it.get_ray().get_direction(), true,
			it.get_ray().get_level() + 1., it.get_ray().get_ior(),
			(dot(it.get_ray().get_direction(), it.get_normal()) > 0.) ? m_medium_0 : m_medium_1);

		R = PolarizedAttenuation<D>(1. - m_absorption);
		pdf = 1.;
	}

	/*---------------------------------------------------------------------------*/
	Spectrum get_absorption(const Intersection<D> &) const
	{
		return m_absorption;
	}
}; // class Transparent

typedef Transparent<3> Transparent3D;
typedef Transparent<2> Transparent2D;

#endif // _TRANSPARENT_H_
