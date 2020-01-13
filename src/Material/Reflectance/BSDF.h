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

#ifndef _BSDF_H_
#define _BSDF_H_

#include "bunnykiller.h"
#include "Material/Reflectance/Reflectance.h"
#include "Color/Spectrum.h"
#include "Color/PolarizedAttenuation.h"
#include "Color/FluorescentAttenuation.h"
#include "LinearAlgebra/Vector3.h"
#include "LinearAlgebra/Vector2.h"
#include "RayTracing/Ray.h"
#include "RayTracing/Intersection.h"

//template<int D> class Ray;
//template<int D> class Intersection;

/**	Base virtual class for BSDFs. 
 These model both reflectance (BRDF) and transmitance (BTDF) functions. */
template<unsigned D>
class BSDF
{
protected:
	Reflectance::Type m_type;
public:
	BSDF(Reflectance::Type type = Reflectance::NOTAG) :
			m_type(type)
	{
	}

	virtual ~BSDF()
	{
	}
	
	/*---------------------------------------------------------------------------*/
	/** Virtual function for the PDF */
	virtual Real p(const VectorN<D> &omega_i, const VectorN<D> &omega_o, const VectorN<D> &normal,
		const Vector2 &uv) const = 0;

	virtual Real p(const VectorN<D> &omega_i, const VectorN<D> &omega_o, const VectorN<D> &normal,
		const Vector2 &uv, Real delta_time) const
	{
		return p(omega_i, omega_o, normal, uv);
	}

	/*---------------------------------------------------------------------------*/
	/** Spectrum-based attenuation */

	/** Virtual function for the BRDF */
	virtual void f(const VectorN<D> &omega_i, const VectorN<D> &omega_o, const VectorN<D> &normal,
		const Vector2 &uv, Spectrum &R) const = 0;

	/** Virtual function for sampling an outgoing direction, where */
	virtual void sample_direction(const VectorN<D> &omega_i, VectorN<D> &omega_o,
		const VectorN<D> &normal, const Vector2 &uv, Spectrum &R, Real &pdf) const = 0;

	/** Virtual function for sampling an outgoing ray */
	virtual void sample_outgoing_ray(const Intersection<D> &it, Ray<D> &new_ray, Spectrum &R,
		Real &pdf) const = 0;

	// TIME RESOLVED FUNCTIONS (see Material/Material.h)

	/** Virtual function for the BRDF */
	virtual void f(const VectorN<D> &omega_i, const VectorN<D> &omega_o, const VectorN<D> &normal,
		const Vector2 &uv, Spectrum &R, Real &delta_time, Real &pdf) const
	{
		delta_time = 0.;
		pdf = 1.f;
		f(omega_i, omega_o, normal, uv, R);
	}

	/** Virtual function for sampling an outgoing direction, where */
	virtual void sample_direction(const VectorN<D> &omega_i, VectorN<D> &omega_o,
		const VectorN<D> &normal, const Vector2 &uv, Spectrum &R, Real &delta_time, Real &pdf) const
	{
		delta_time = 0.;
		sample_direction(omega_i, omega_o, normal, uv, R, pdf);
	}

	/** Virtual function for sampling an outgoing ray */
	virtual void sample_outgoing_ray(const Intersection<D> &it, Ray<D> &new_ray, Spectrum &R,
		Real &delta_time, Real &pdf) const
	{
		delta_time = 0.;
		sample_outgoing_ray(it, new_ray, R, pdf);
	}

	/*---------------------------------------------------------------------------*/
	/** Fluorescence-based attenuation */
	/** Did i said i hate copy pasting code? */
	virtual void f(const VectorN<D> &omega_i, const VectorN<D> &omega_o, const VectorN<D> &normal,
		const Vector2 &uv, FluorescentAttenuation<D> &R) const
	{
		Spectrum r;
		f(omega_i, omega_o, normal, uv, r);
		R = FluorescentAttenuation<D>(r);
	}

	/** Virtual function for sampling an outgoing direction, where */
	virtual void sample_direction(const VectorN<D> &omega_i, VectorN<D> &omega_o,
		const VectorN<D> &normal, const Vector2 &uv, FluorescentAttenuation<D> &R, Real &pdf) const
	{
		Spectrum r;
		sample_direction(omega_i, omega_o, normal, uv, r, pdf);
		R = FluorescentAttenuation<D>(r);
	}

	/** Virtual function for sampling an outgoing ray */
	virtual void sample_outgoing_ray(const Intersection<D> &it, Ray<D> &new_ray,
		FluorescentAttenuation<D> &R, Real &pdf) const
	{
		Spectrum r;
		sample_outgoing_ray(it, new_ray, r, pdf);
		R = FluorescentAttenuation<D>(r);
	}

	// TIME RESOLVED FUNCTIONS
	/** Virtual function for the BRDF */
	virtual void f(const VectorN<D> &omega_i, const VectorN<D> &omega_o, const VectorN<D> &normal,
		const Vector2 &uv, FluorescentAttenuation<D> &R, Real &delta_time, Real &pdf) const
	{
		delta_time = 0.;
		pdf = 1.f;
		f(omega_i, omega_o, normal, uv, R);
	}

	/** Virtual function for sampling an outgoing direction, where */
	virtual void sample_direction(const VectorN<D> &omega_i, VectorN<D> &omega_o,
		const VectorN<D> &normal, const Vector2 &uv, FluorescentAttenuation<D> &R, Real &delta_time,
		Real &pdf) const
	{
		delta_time = 0.;
		sample_direction(omega_i, omega_o, normal, uv, R, pdf);
	}

	/** Virtual function for sampling an outgoing ray */
	virtual void sample_outgoing_ray(const Intersection<D> &it, Ray<D> &new_ray,
		FluorescentAttenuation<D> &R, Real &delta_time, Real &pdf) const
	{
		delta_time = 0.;
		sample_outgoing_ray(it, new_ray, R, pdf);
	}

	/*---------------------------------------------------------------------------*/
	/** Polarized spectrum-based attenuation */
	/** Computes the local direct reflected (or refracted) light.
	 omega_i is the view (incoming) direction, while omega_o is light (outgoing)*/
	virtual void f(const VectorN<D> &omega_i, const VectorN<D> &omega_o, const VectorN<D> &normal,
		const Vector2 &uv, PolarizedAttenuation<D> &R) const
	{
		Spectrum r;
		f(omega_i, omega_o, normal, uv, r);
		R = PolarizedAttenuation<D>(r, PolarizationFrame<D>(omega_o, omega_i),
				PolarizationFrame<D>(-omega_i, omega_o));
	}

	virtual void f(const VectorN<D> &omega_i, const VectorN<D> &omega_o, const VectorN<D> &normal,
		const Vector2 &uv, PolarizedAttenuation<D> &R, Real &delta_time, Real &pdf) const
	{
		delta_time = 0.;
		pdf = 1.f;
		f(omega_i, omega_o, normal, uv, R);
	}

	/** Samples a new outgoing direction */
	virtual void sample_direction(const VectorN<D> &omega_i, VectorN<D> &omega_o,
		const VectorN<D> &normal, const Vector2 &uv, PolarizedAttenuation<D> &R, Real &pdf) const
	{
		Spectrum r;
		sample_direction(omega_i, omega_o, normal, uv, r, pdf);
		R = PolarizedAttenuation<D>(r, PolarizationFrame<D>(omega_o, omega_i),
				PolarizationFrame<D>(-omega_i, omega_o));
	}

	virtual void sample_direction(const VectorN<D> &omega_i, VectorN<D> &omega_o,
		const VectorN<D> &normal, const Vector2 &uv, PolarizedAttenuation<D> &R, Real &delta_time,
		Real &pdf) const
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
				PolarizationFrame<D>(new_ray.get_direction(), it.get_ray().get_direction()),
				PolarizationFrame<D>(-it.get_ray().get_direction(), new_ray.get_direction()));
	}

	virtual void sample_outgoing_ray(const Intersection<D> &it, Ray<D> &new_ray,
		PolarizedAttenuation<D> &R, Real &delta_time, Real &pdf) const
	{
		delta_time = 0.;
		sample_outgoing_ray(it, new_ray, R, pdf);
	}

	/*---------------------------------------------------------------------------*/
	/** Fluorescent Polarized spectrum-based attenuation */
	/** Computes the local direct reflected (or refracted) light.
	 omega_i is the view (incoming) direction, while omega_o is light (outgoing)*/
	virtual void f(const VectorN<D> &omega_i, const VectorN<D> &omega_o, const VectorN<D> &normal,
		const Vector2 &uv, FluorescentPolarizedAttenuation<D> &R) const
	{
		Spectrum r;
		f(omega_i, omega_o, normal, uv, r);

		R = FluorescentPolarizedAttenuation<D>(r, PolarizationFrame<D>(omega_o, omega_i),
				PolarizationFrame<D>(-omega_i, omega_o));
	}

	virtual void f(const VectorN<D> &omega_i, const VectorN<D> &omega_o, const VectorN<D> &normal,
		const Vector2 &uv, FluorescentPolarizedAttenuation<D> &R, Real &delta_time, Real &pdf) const
	{
		delta_time = 0.;
		pdf = 1.f;
		f(omega_i, omega_o, normal, uv, R);
	}

	/** Samples a new outgoing direction */
	virtual void sample_direction(const VectorN<D> &omega_i, VectorN<D> &omega_o,
		const VectorN<D> &normal, const Vector2 &uv, FluorescentPolarizedAttenuation<D> &R,
		Real &pdf) const
	{
		Spectrum r;
		sample_direction(omega_i, omega_o, normal, uv, r, pdf);

		R = FluorescentPolarizedAttenuation<D>(r, PolarizationFrame<D>(omega_o, omega_i),
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
				PolarizationFrame<D>(new_ray.get_direction(), it.get_ray().get_direction()),
				PolarizationFrame<D>(-it.get_ray().get_direction(), new_ray.get_direction()));
	}

	virtual void sample_outgoing_ray(const Intersection<D> &it, Ray<D> &new_ray,
		FluorescentPolarizedAttenuation<D> &R, Real &delta_time, Real &pdf) const
	{
		delta_time = 0.;
		sample_outgoing_ray(it, new_ray, R, pdf);
	}

	/*---------------------------------------------------------------------------*/
	/** Get the absorption of the BSDF */
	virtual Spectrum get_absorption(const Vector2 &uv) const = 0;

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
};
//BSDF

typedef BSDF<3> BSDF3D;
typedef BSDF<2> BSDF2D;

#endif //_BSDF_H_
