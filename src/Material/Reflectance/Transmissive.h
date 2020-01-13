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

#ifndef _TRANSMISSIVE_H_
#define _TRANSMISSIVE_H_

#include "Material/Reflectance/BSDF.h"
#include "RayTracing/Ray.h"
#include "RayTracing/Intersection.h"

/**	Simple Delta Transmissive BSDF. 
	Note that it does not use real fresnel transmission*/
template<unsigned D>
class Transmissive: public BSDF<D>
{
	//Absorption coefficient
	Spectrum absorption;
	//Index of refraction
	Real refraction_index;
public:
	Transmissive():absorption(0.f), refraction_index(1.f), BSDF<D>( Reflectance::DELTA_TRANSMISSION ){}
	Transmissive(const Spectrum &_absorption, const Real _refraction)
		:absorption(_absorption),refraction_index(_refraction),BSDF<D>(Reflectance::DELTA_TRANSMISSION){}
	Transmissive(const Real _refraction)
		:absorption(0.f),refraction_index(_refraction),BSDF<D>(Reflectance::DELTA_TRANSMISSION){}
	
	Spectrum f(const VectorN<D> &omega_i, const VectorN<D> &omega_o, const VectorN<D> &normal, const Vector2 &uv) const
	{
		if( dot(omega_i.refract(normal, refraction_index), omega_o) > 0.99 )
			return (1.-absorption)/dot_abs(normal, omega_o);

		return Spectrum(0.);
	}
	Spectrum sample_direction(const VectorN<D> &omega_i, VectorN<D> &omega_o, const VectorN<D> &normal, const Vector2 &uv, Real &pdf) const
	{
		pdf = 1.; //Must take care with absorption somehow...
		omega_o = omega_i.refract(normal, refraction_index);
		return (1.-absorption)/dot_abs(normal, omega_o);
	}
	Spectrum sample_outgoing_ray(const Intersection<D> &it, Ray<D> &new_ray, Real &pdf) const
	{
		VectorN<D> omega=it.get_ray().get_direction().refract(it.get_normal(), refraction_index);
		//Test if the ray enters the media or if it is leaving
		new_ray = Ray<D>(it.get_position(), omega, true,
						 it.get_ray().get_level()+1, (dot(omega, it.get_normal()) < 0.)?refraction_index:DEFAULT_REFRACTION_INDEX, 
						 it.get_ray().get_medium());
		//printf("Transmissive new_ray.medium = %x\n", it.get_ray().get_medium());
		new_ray.shift();
		pdf = 1;
		return (1.-absorption)/dot_abs(it.get_normal(), new_ray.get_direction());
	}

	Spectrum get_absorption( const Vector2 &uv ) const{ return absorption; }
}; //Transmissive

typedef Transmissive<3> Transmissive3D;
typedef Transmissive<2> Transmissive2D;

#endif //_TRANSMISSIVE_H_
