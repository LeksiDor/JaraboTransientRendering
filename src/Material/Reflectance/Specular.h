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

#ifndef _SPECULAR_H_
#define _SPECULAR_H_

#include "Material/Reflectance/BSDF.h"
#include "LinearAlgebra/Vector3.h"
#include "RayTracing/Ray.h"
#include "RayTracing/Intersection.h"

/**	Simple Delta Specular BRDF. */
template<unsigned D>
class Specular: public BSDF<D>
{
	//Absorption coefficient
	Spectrum absorption;
public:
	Specular<D>():absorption(0.f),BSDF<D>(Reflectance::DELTA){}
	Specular<D>(const Spectrum &_absorption):absorption(_absorption),BSDF<D>(Reflectance::DELTA){}
	
	Spectrum f(const VectorN<D> &omega_i, const VectorN<D> &omega_o, const VectorN<D> &normal, const Vector2 &uv) const
	{
		if (omega_o.reflect(normal).dot(omega_i) > 0.999)
			return (1. - absorption)/dot_abs(normal, omega_o);

		return Spectrum();
	}

	Spectrum sample_direction(const VectorN<D> &omega_i, VectorN<D> &omega_o, const VectorN<D> &normal, const Vector2 &uv, Real &pdf) const
	{
		omega_o = omega_i.reflect(normal);
		pdf = 1;//./dot(normal, omega_o); //Must take care with absorption somehow...
		return (1.-absorption)/dot_abs(normal, omega_o);
	}

	Spectrum sample_outgoing_ray(const Intersection<D> &it, Ray<D> &new_ray, Real &pdf) const
	{
		new_ray = Ray<D>(it.get_position(), it.get_ray().get_direction().reflect(it.get_normal()), true,
			it.get_ray().get_level()+1, it.get_ray().get_refraction_index(), it.get_ray().get_medium());
		new_ray.shift();
		pdf = 1.;

		return (1. - absorption)/dot_abs(it.get_normal(), new_ray.get_direction());
	}

	Spectrum get_absorption( const Vector2 &uv ) const
	{
		return absorption;
	}
}; //Specular

typedef Specular<3> Specular3D;
typedef Specular<2> Specular2D;

#endif //_SPECULAR_H_
