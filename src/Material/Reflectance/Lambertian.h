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

#ifndef _LAMBERTIAN_H_
#define _LAMBERTIAN_H_

#include "bunnykiller.h"

#include <cmath>

#include "Material/Reflectance/BSDF.h"
#include "RayTracing/Ray.h"
#include "RayTracing/Intersection.h"
#include "Sampling/Sampling.h"

/**	Simple Lambertian BRDF. */
template<unsigned D>
class Lambertian : public BSDF<D>
{
	// Diffuse coefficient
	Spectrum Kd;
public:
	Lambertian(Spectrum _Kd = Spectrum(.67)) :
			BSDF<D>(Reflectance::DIFFUSE),
			Kd(_Kd)
	{
	}
	
	void f(const VectorN<D> &omega_i, const VectorN<D> &omega_o, const VectorN<D> &normal,
		const Vector2 &uv, Spectrum &R) const;

	Real p(const VectorN<D> &omega_i, const VectorN<D> &omega_o, const VectorN<D> &normal,
		const Vector2 &uv) const;

	void sample_direction(const VectorN<D> &omega_i, VectorN<D> &omega_o, const VectorN<D> &normal,
		const Vector2 &uv, Spectrum &R, Real &pdf) const;

	void sample_outgoing_ray(const Intersection<D> &it, Ray<D> &new_ray, Spectrum &R,
		Real &pdf) const;

	Spectrum get_absorption(const Vector2 &) const
	{
		return Spectrum(1.) - Kd;
	}
}; // Lambertian

//------------------------------------------
//Specializations for 3D
template<unsigned D>
void Lambertian<D>::f(const VectorN<D> &omega_i, const VectorN<D> &omega_o,
	const VectorN<D> &normal, const Vector2 &uv, Spectrum &R) const
{
	Real norm(0.);
	
	if (D == 2)
		norm = static_cast<Real>(.5);
	if (D == 3)
		norm = static_cast<Real>(M_1_PI);

	R = Kd * norm;
}

template<unsigned D>
Real Lambertian<D>::p(const VectorN<D> &omega_i, const VectorN<D> &omega_o,
	const VectorN<D> &normal, const Vector2 &uv) const
{
	Real norm(0.);

	if (D == 2)
		norm = static_cast<Real>(.5);
	if (D == 3)
		norm = static_cast<Real>(M_1_PI);

	return dot_clamped(omega_o, normal) * norm;
}

template<unsigned D>
void Lambertian<D>::sample_direction(const VectorN<D> &omega_i, VectorN<D> &omega_o,
	const VectorN<D> &normal, const Vector2 &uv, Spectrum &R, Real &pdf) const
{
	Sampling.cosine_weighted(omega_o, pdf);
	omega_o.transform_matrix_to(VectorN<D>(0., 1.), normal);

	f(-omega_i, -omega_o, normal, uv, R);
}

template<unsigned D>
void Lambertian<D>::sample_outgoing_ray(const Intersection<D> &it, Ray<D> &new_ray, Spectrum &R,
	Real &pdf) const
{
	VectorN<D> new_omega_o;
	Sampling.cosine_weighted(new_omega_o, pdf);
	
	new_omega_o = new_omega_o[0] * it.get_tangent_x() + new_omega_o[1] * it.get_normal()
			+ new_omega_o[2] * it.get_tangent_y();

	const Ray<D>& old_ray = it.get_ray();
	
	f(-old_ray.get_direction(), -new_omega_o, it.get_normal(), it.get_uv(), R);

	new_ray = Ray<D>(it.get_position(), new_omega_o, true, old_ray.get_level() + 1,
			old_ray.get_ior(), old_ray.get_medium());
}

typedef Lambertian<3> Lambertian3D;
typedef Lambertian<2> Lambertian2D;

#endif //_LAMBERTIAN_H_
