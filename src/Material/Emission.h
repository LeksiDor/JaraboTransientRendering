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

#ifndef _EMISSION_H_
#define _EMISSION_H_

#include <Image/Image.h>
#include "Material/Material.h"
#include "RayTracing/Ray.h"
#include "RayTracing/Intersection.h"
#include <vector>

/**	Class defining the simplest emissing
	material. 
	TODO: Parametrize the light distribution,
		now it's only gaussian. */
template<unsigned D, class Radiance>
class Emission: public Material<D>
{
	Radiance m_power;
	Imaging::Image *m_emission_map;

	Real m_delay;
public:
	Emission() : Material<D>(Reflectance::EMISSIVE), m_power(1.), m_delay(0.), m_emission_map(0){}
	Emission(const char* name) :Material<D>(Reflectance::EMISSIVE), m_power(1.), m_delay(0.), m_emission_map(0)
	{
		m_emission_map = new Imaging::Image(name);
	}
	Emission(Radiance power, const char* name) : Material<D>(Reflectance::EMISSIVE), m_power(power), m_delay(0.), m_emission_map(0)
	{
		m_emission_map = new Imaging::Image(name);
	}
	~Emission()
	{	if(m_emission_map) delete m_emission_map;	}

	
	void f(const VectorN<D> &omega_i, const VectorN<D> &omega_o, const VectorN<D> &normal, const Vector2 &uv, Spectrum &R) const
	{	
		R = 0.;
	}

	/** Computes the pdf */
	Real p(const VectorN<D> &omega_i, const VectorN<D> &omega_o, const VectorN<D> &normal, const Vector2 &uv )const
	{	return 1.; }

	void sample_direction(const VectorN<D> &omega_i, VectorN<D> &omega_o, const VectorN<D> &normal, const Vector2 &uv, Spectrum &R, Real &pdf) const
	{
		pdf = 1.;
		// Need for importance sampling...	
		R = Spectrum(0.);
	}
	void sample_outgoing_ray(const Intersection<D> &it, Ray<D> &new_ray, Spectrum &R, Real &pdf) const
	{
		pdf = 1.;
		R = Spectrum(0.);
	}
	
	Spectrum emit_light( const VectorN<D> &omega_o, const VectorN<D> &normal, const Vector2 &uv )const
	{
		Real di = dot_clamped(normal, omega_o);
		Spectrum p = m_power * di;

		if (m_emission_map) {
			Imaging::RGBColor rgb = (*m_emission_map)(m_emission_map->width()*uv[0], m_emission_map->height()*uv[1]);
			//printf("RGB Emission @[%f,%f]: %f, %f, %f\n", m_emission_map->width()*uv[0],m_emission_map->height()*uv[1], rgb.r(), rgb.g(), rgb.b());
			p *= Spectrum(rgb.r(), rgb.g(), rgb.b());
		}
		
		return p;	
	}

	Spectrum emit_light( const VectorN<D> &omega_o, const VectorN<D> &normal, const Vector2 &uv, Real &delta_time, Real &pdf )const
	{
		pdf = 1; delta_time = m_delay; return emit_light(omega_o, normal, uv);
	}

	Spectrum get_absorption(const Intersection<D> &it ) const
	{
		return Spectrum(1.);
	}
}; // Emission

typedef Emission<3, Spectrum> Emission3D;
typedef Emission<2, Spectrum> Emission2D;

#endif //_EMISSION_H_
