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

#ifndef _MATERIAL_H_
#define _MATERIAL_H_

#include "bunnykiller.h"
#include "LinearAlgebra/Vector3.h"
#include "LinearAlgebra/Vector2.h"
#include "LightSource/LightSample.h"
#include "Color/Spectrum.h"
#include "Color/PolarizedAttenuation.h"
#include "Color/FluorescentAttenuation.h"
#include "Color/FluorescentPolarizedAttenuation.h"
#include "Reflectance/Reflectance.h"
#include "RayTracing/Ray.h"

template<unsigned D> class Intersection;
template<unsigned D> class Medium;

/**	Base class for materials. */
template<unsigned D>
class Material
{
protected:
	// Default world media
	static Medium<D> *m_default_medium;
	// Default world index of refraction
	static Real m_default_n;
	// Type of material
	Reflectance::Type m_type;
public:
	Material(Reflectance::Type type = Reflectance::NOTAG):m_type(type) {}

	virtual ~Material() {}
	
	/*---------------------------------------------------------------------------*/
	/** PDFs */

	/** Virtual function for the PDF */
	virtual Real p(const VectorN<D> &omega_i, const VectorN<D> &omega_o,
			const VectorN<D> &normal, const Vector2 &uv) const = 0;

	/** Virtual function for the transient PDF */
	virtual Real p(const VectorN<D> &omega_i, const VectorN<D> &omega_o,
			const VectorN<D> &normal, const Vector2 &uv, Real delta_time) const
	{
		return p(omega_i, omega_o, normal, uv);
	}

	/*---------------------------------------------------------------------------*/
	/** Spectrum-based attenuation */
	/** Computes the local direct reflected (or refracted) light. 
		omega_i is the view (incoming) direction, while omega_o is light (outgoing)*/
	virtual void f(const VectorN<D> &omega_i, const VectorN<D> &omega_o,
			const VectorN<D> &normal, const Vector2 &uv, Spectrum &R) const = 0;

	virtual void f(const VectorN<D> &omega_i, const VectorN<D> &omega_o,
			const VectorN<D> &normal, const Vector2 &uv, Spectrum &R,
			Real &delta_time, Real &pdf) const
	{
		delta_time = 0.; pdf = 1.f;
		f(omega_i, omega_o, normal, uv, R);
	}

	/** Samples a new outgoing direction */
	virtual void sample_direction(const VectorN<D> &omega_i, VectorN<D> &omega_o,
			const VectorN<D> &normal, const Vector2 &uv, Spectrum &R, Real &pdf) const = 0;

	virtual void sample_direction(const VectorN<D> &omega_i, VectorN<D> &omega_o,
			const VectorN<D> &normal, const Vector2 &uv, Spectrum &R, Real &delta_time,
			Real &pdf) const
	{
		delta_time = 0.;
		sample_direction(omega_i, omega_o, normal, uv, R, pdf);
	}

	/** Samples a new outgoing ray in the intersection */	
	virtual void sample_outgoing_ray(const Intersection<D> &it, Ray<D> &new_ray,
			Spectrum &R, Real &pdf) const = 0;

	virtual void sample_outgoing_ray(const Intersection<D> &it, Ray<D> &new_ray,
			Spectrum &R, Real &delta_time, Real &pdf) const
	{
		delta_time = 0.;
		sample_outgoing_ray(it, new_ray, R, pdf);
	}

	/** Returns the emited light in direction omega_o */
	virtual void emit_light(const VectorN<D> &omega_o, const VectorN<D> &normal,
			const Vector2 &uv, Spectrum &R) const
	{
		R = Spectrum(0.);
	}

	virtual void emit_light(const VectorN<D> &omega_o, const VectorN<D> &normal,
			const Vector2 &uv, Spectrum &R, Real &delta_time, Real &pdf)const
	{
		delta_time = 0; pdf = 1;
		emit_light(omega_o, normal, uv, R);
	}

	/*---------------------------------------------------------------------------*/
	/** Polarized spectrum-based attenuation */
	/** Computes the local direct reflected (or refracted) light.
	    omega_i is the view (incoming) direction, while omega_o is light (outgoing)
	    direction */
	virtual void f(const VectorN<D> &omega_i, const VectorN<D> &omega_o,
			const VectorN<D> &normal, const Vector2 &uv, PolarizedAttenuation<D> &R) const
	{
		Spectrum r;
		f(omega_i, omega_o, normal, uv, r);
		R = PolarizedAttenuation<D>(r, PolarizationFrame<D>(omega_o, omega_i),
			    					PolarizationFrame<D>(-omega_i, omega_o));
	}

	virtual void f(const VectorN<D> &omega_i, const VectorN<D> &omega_o,
			const VectorN<D> &normal, const Vector2 &uv, PolarizedAttenuation<D> &R,
			Real &delta_time, Real &pdf) const
	{
		delta_time = 0.; pdf = 1.f;
		f(omega_i, omega_o, normal, uv, R);
	}

	/** Samples a new outgoing direction */
	virtual void sample_direction(const VectorN<D> &omega_i, VectorN<D> &omega_o,
			const VectorN<D> &normal, const Vector2 &uv, PolarizedAttenuation<D> &R,
			Real &pdf) const
	{
		Spectrum r;
		sample_direction(omega_i, omega_o, normal, uv, r, pdf);
		R = PolarizedAttenuation<D>(r, PolarizationFrame<D>(omega_o, omega_i),
									PolarizationFrame<D>(-omega_i, omega_o));
	}

	virtual void sample_direction(const VectorN<D> &omega_i, VectorN<D> &omega_o,
			const VectorN<D> &normal, const Vector2 &uv, PolarizedAttenuation<D> &R,
			Real &delta_time, Real &pdf) const
	{
		delta_time = 0.;
		sample_direction(omega_i, omega_o, normal, uv, R, pdf);
	}

	/** Samples a new outgoing ray in the intersection */
	virtual void sample_outgoing_ray(const Intersection<D> &it, Ray<D> &new_ray,
			PolarizedAttenuation<D> &R, Real &pdf) const
	{
		Spectrum r;
		sample_outgoing_ray(it, new_ray, r, pdf);
		R = PolarizedAttenuation<D>(r,
				                    PolarizationFrame<D>(new_ray.get_direction(),
									                     it.get_ray().get_direction()),
									PolarizationFrame<D>(-it.get_ray().get_direction(),
									                     new_ray.get_direction()));
	}

	virtual void sample_outgoing_ray(const Intersection<D> &it, Ray<D> &new_ray,
			PolarizedAttenuation<D> &R, Real &delta_time, Real &pdf) const
	{
		delta_time = 0.;
		sample_outgoing_ray(it, new_ray, R, pdf);
	}

	/** Returns the emited light in direction omega_o */
	virtual void emit_light(const VectorN<D> &omega_o, const VectorN<D> &normal,
			const Vector2 &uv, PolarizedAttenuation<D> &R)const
	{
		R = PolarizedAttenuation<D>(0);
	}

	virtual void emit_light(const VectorN<D> &omega_o, const VectorN<D> &normal,
			const Vector2 &uv, PolarizedAttenuation<D> &R, Real &delta_time,
			Real &pdf) const
	{
		delta_time = 0; pdf = 1;
		emit_light(omega_o, normal, uv, R);
	}

	/*---------------------------------------------------------------------------*/
	/** Fluorescence spectrum-based attenuation (i hate copy-pasting code)*/
	virtual void f(const VectorN<D> &omega_i, const VectorN<D> &omega_o,
			const VectorN<D> &normal, const Vector2 &uv, FluorescentAttenuation<D> &R) const
	{
		Spectrum r;
		f(omega_i, omega_o, normal, uv, r);

		R = FluorescentAttenuation<D>(r);
	}

	virtual void f(const VectorN<D> &omega_i, const VectorN<D> &omega_o,
			const VectorN<D> &normal, const Vector2 &uv, FluorescentAttenuation<D> &R,
			Real &delta_time, Real &pdf) const
	{
		delta_time = 1.; pdf = 1.f;
		f(omega_i, omega_o, normal, uv, R);
	}

	/** Samples a new outgoing direction */
	virtual void sample_direction(const VectorN<D> &omega_i, VectorN<D> &omega_o,
			const VectorN<D> &normal, const Vector2 &uv, FluorescentAttenuation<D> &R,
			Real &pdf) const
	{
		Spectrum r;
		sample_direction(omega_i, omega_o, normal, uv, r, pdf);
		R = FluorescentAttenuation<D>(r);
	}

	virtual void sample_direction(const VectorN<D> &omega_i, VectorN<D> &omega_o,
			const VectorN<D> &normal, const Vector2 &uv, FluorescentAttenuation<D> &R,
			Real &delta_time, Real &pdf) const
	{
		delta_time = 0.;
		sample_direction(omega_i, omega_o, normal, uv, R, pdf);
	}

	/** Samples a new outgoing ray in the intersection */
	virtual void sample_outgoing_ray(const Intersection<D> &it, Ray<D> &new_ray,
			FluorescentAttenuation<D> &R, Real &pdf) const
	{
		Spectrum r;
		sample_outgoing_ray(it, new_ray, r, pdf);
		R = FluorescentAttenuation<D>(r);
	}

	virtual void sample_outgoing_ray(const Intersection<D> &it, Ray<D> &new_ray,
			FluorescentAttenuation<D> &R, Real &delta_time, Real &pdf) const
	{
		delta_time = 0.; sample_outgoing_ray(it, new_ray, R, pdf);
	}

	/** Returns the emited light in direction omega_o */
	virtual void emit_light(const VectorN<D> &omega_o, const VectorN<D> &normal,
			const Vector2 &uv, FluorescentAttenuation<D> &R) const
	{
		R = FluorescentAttenuation<D>(0.0);
	}

	virtual void emit_light(const VectorN<D> &omega_o, const VectorN<D> &normal,
			const Vector2 &uv, FluorescentAttenuation<D> &R, Real &delta_time,
			Real &pdf)const
	{
		delta_time = 0; pdf = 1;
		emit_light(omega_o, normal, uv, R);
	}

	/*-------------------------------------------------------------------------- - */
	/** Polarized fluorescent spectrum-based attenuation */
	virtual void f(const VectorN<D> &omega_i, const VectorN<D> &omega_o,
			const VectorN<D> &normal, const Vector2 &uv,
			FluorescentPolarizedAttenuation<D> &R) const
	{
		Spectrum r;
		f(omega_i, omega_o, normal, uv, r);

		R = FluorescentPolarizedAttenuation<D>(r,
											   PolarizationFrame<D>(omega_o, omega_i),
					                           PolarizationFrame<D>(-omega_i, omega_o));
	}

	virtual void f(const VectorN<D> &omega_i, const VectorN<D> &omega_o,
			const VectorN<D> &normal, const Vector2 &uv, FluorescentPolarizedAttenuation<D> &R,
			Real &delta_time, Real &pdf) const
	{
		delta_time = 0.; pdf = 1.f;
		f(omega_i, omega_o, normal, uv, R);
	}

	/** Samples a new outgoing direction */
	virtual void sample_direction(const VectorN<D> &omega_i, VectorN<D> &omega_o,
			const VectorN<D> &normal, const Vector2 &uv, FluorescentPolarizedAttenuation<D> &R,
			Real &pdf) const
	{
		Spectrum r;
		sample_direction(omega_i, omega_o, normal, uv, r, pdf);
		R = FluorescentPolarizedAttenuation<D>(r,
				                               PolarizationFrame<D>(omega_o, omega_i),
					                           PolarizationFrame<D>(-omega_i, omega_o));
	}

	virtual void sample_direction(const VectorN<D> &omega_i, VectorN<D> &omega_o,
			const VectorN<D> &normal, const Vector2 &uv, FluorescentPolarizedAttenuation<D> &R,
			Real &delta_time, Real &pdf) const
	{
		delta_time = 0.;
		sample_direction(omega_i, omega_o, normal, uv, R, pdf);
	}

	/** Samples a new outgoing ray in the intersection */
	virtual void sample_outgoing_ray(const Intersection<D> &it, Ray<D> &new_ray,
			FluorescentPolarizedAttenuation<D> &R, Real &pdf) const
	{
		Spectrum r;
		sample_outgoing_ray(it, new_ray, r, pdf);
		R = FluorescentPolarizedAttenuation<D>(r,
				                               PolarizationFrame<D>(new_ray.get_direction(),
				                                                    it.get_ray().get_direction()),
					                           PolarizationFrame<D>(-it.get_ray().get_direction(),
					                                                new_ray.get_direction()));
	}

	virtual void sample_outgoing_ray(const Intersection<D> &it, Ray<D> &new_ray,
			FluorescentPolarizedAttenuation<D> &R, Real &delta_time, Real &pdf) const
	{
		delta_time = 0.;
		sample_outgoing_ray(it, new_ray, R, pdf);
	}

	/** Returns the emited light in direction omega_o */
	virtual void emit_light(const VectorN<D> &omega_o, const VectorN<D> &normal,
			const Vector2 &uv, FluorescentPolarizedAttenuation<D> &R) const
	{
		R = FluorescentPolarizedAttenuation<D>(0.);
	}

	virtual void emit_light(const VectorN<D> &omega_o, const VectorN<D> &normal,
			const Vector2 &uv, FluorescentPolarizedAttenuation<D> &R, Real &delta_time,
			Real &pdf) const
	{
		delta_time = 0; pdf = 1;
		emit_light(omega_o, normal, uv, R);
	}

	/*---------------------------------------------------------------------------*/

	/** Get the absorption of the material */
	virtual Spectrum get_absorption(const Intersection<D> &it) const = 0;

	/** Returns material's type */
	Reflectance::Type get_type() const
	{
		return m_type;
	}

	/** Returns material's type */
	bool is_type(Reflectance::Type type) const
	{
		return (type & m_type);
	}
	
	static void set_default_medium(Medium<D> *m)
	{
		m_default_medium = m;
	}

	static void set_default_refraction_index(Real n)
	{
		m_default_n = n;
	}
}; // Material

typedef Material<3> Material3D;
typedef Material<2> Material2D;

#endif
