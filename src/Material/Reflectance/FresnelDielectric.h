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

#ifndef _FRESNEL_DIELECTRIC_H_
#define _FRESNEL_DIELECTRIC_H_

#include "Material/Reflectance/BSDF.h"
#include "Material/Reflectance/Fresnel.h"
#include "RayTracing/Ray.h"
#include "RayTracing/Intersection.h"
#include "Utils/RandomNumbers.h"

template<unsigned D>
class FresnelDielectric : public BSDF<D>
{
	//Absorption coefficient
	Spectrum m_absorption;

	//Index of refraction
	Real m_n0;
	Real m_n1;

	void get_eta(const VectorN<D> &omega, const VectorN<D> &normal, Real &n, Real &in) const
	{
		n = (dot(omega, normal) > 0.) ? (m_n1 / m_n0) : (m_n0 / m_n1);
		in = (dot(omega, normal) < 0.) ? (m_n1 / m_n0) : (m_n0 / m_n1);
	}

//Reflectance::DELTA_TRANSMISSION
public:
	FresnelDielectric() :m_absorption(0.f), m_n1(DEFAULT_REFRACTION_INDEX), m_n0(DEFAULT_REFRACTION_INDEX), BSDF<D>(Reflectance::DELTA){}
	FresnelDielectric(const Spectrum &absorption, const Real n_in, const Real n_out = DEFAULT_REFRACTION_INDEX)
		:m_absorption(absorption), m_n1(n_in), m_n0(n_out), BSDF<D>(Reflectance::DELTA){}
	FresnelDielectric(const Real n_in, const Real n_out = DEFAULT_REFRACTION_INDEX)
		:m_absorption(0.), m_n1(n_in), m_n0(n_out), BSDF<D>(Reflectance::DELTA){}

	virtual ~FresnelDielectric() {}

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
		/*if( dot(omega_i.refract(normal, m_n), omega_o) > 0.99 )
			R = Spectrum(1. - m_absorption) / dot_abs(normal, omega_o);
		else*/
		R = Spectrum(0.);
	}

	void sample_direction(const VectorN<D> &omega_i, VectorN<D> &omega_o, const VectorN<D> &normal, const Vector2 &uv, Spectrum &R, Real &pdf) const
	{
		Real eta, ieta;
		get_eta(omega_i, normal, eta, ieta);

		Real r = Fresnel::fresnel_reflectance(omega_i, normal, ieta, 0.);

		// Russian roulette
		Real epsilon1 = RNG::StdRandom::get_real();

		if (epsilon1 < r ) // Fresnel reflection
		{
			omega_o = omega_i.reflect(normal);
			pdf = r;
			R = r;
		}
		else // Refraction
		{
			pdf = 1. - r;
			R = 1. - r;
			omega_o = omega_i.refract(normal, eta).normalized();
		}

		R *= (1 - m_absorption);
	}
	
	void sample_outgoing_ray(const Intersection<D> &it, Ray<D> &new_ray, Spectrum &R, Real &pdf) const
	{
		VectorN<D> omega;
		Real eta, ieta;
		get_eta(it.get_ray().get_direction(), it.get_normal(), eta, ieta);

		Real r = Fresnel::fresnel_reflectance(it.get_ray().get_direction(), it.get_normal(), ieta, 0.);
		Real ior;
		bool reflected = 1;

		// Russian roulette
		Real epsilon1 = RNG::StdRandom::get_real();

		if (epsilon1 < r ) // Fresnel reflection
		{
			omega = it.get_ray().get_direction().reflect(it.get_normal());
			pdf = r;
			R = r;
			ior = it.get_ray().get_refraction_index();
		}
		else // Refraction
		{
			omega = it.get_ray().get_direction().refract(it.get_normal(), eta).normalized();
			pdf =  1 - r;
			R = 1 - r;
			ior = (dot(omega, it.get_normal()) < 0.) ? m_n0 : m_n1;
		}

		//Test if the ray enters the media or if it is leaving
		new_ray = Ray<D>(it.get_position(), omega, true,
						 it.get_ray().get_level()+1, ior, 
						 NULL);
		
		new_ray.shift();

		R *= (1 - m_absorption);
	}

	/*---------------------------------------------------------------------------*/
	/** Polarized spectrum-based attenuation */
	virtual void f(const VectorN<D> &omega_i, const VectorN<D> &omega_o, const VectorN<D> &normal, const Vector2 &uv, PolarizedLightAttenuation<D> &R) const
	{
		R = PolarizedLightAttenuation<D>(0.,
			PolarizationFrame<D>(omega_o, omega_i),
			PolarizationFrame<D>(-omega_i, omega_o));
	}

	/** Samples a new outgoing direction */
	virtual void sample_direction(const VectorN<D> &omega_i, VectorN<D> &omega_o, const VectorN<D> &normal, const Vector2 &uv, PolarizedLightAttenuation<D> &R, Real &pdf) const
	{
		Real eta, ieta;
		get_eta(omega_i, normal, eta, ieta);

		Real Rp, Rs, Pp, Ps;
		Fresnel::fresnel_reflectance(omega_i, normal, ieta, 0., Rp, Rs, Pp, Ps);

		Real r = (Rs + Rp)*.5;

		// Russian roulette
		Real epsilon1 = RNG::StdRandom::get_real();

		if (epsilon1 > r) // Fresnel transmission
		{
			omega_o = omega_i.refract(normal, eta).normalized();
			pdf = 1 - r;
			R = PolarizedLightAttenuation<D>(1 - r, 
					PolarizationFrame<D>(-omega_o, omega_i),
					PolarizationFrame<D>(omega_i, omega_o));
		}
		else
		{
			Real sqRsRp = sqrt(Rs*Rp);
			Real A = r,
				B = (Rs - Rp)*.5,
				C = cosf(Ps - Pp)*sqRsRp,
				S = sinf(Ps - Pp)*sqRsRp;

			omega_o = omega_i.reflect(normal);

			R = PolarizedLightAttenuation<D>(A, B, 0, 0, B, A, 0, 0, 0, 0, C, S, 0, 0, -S, C, 
					PolarizationFrame<D>(omega_o, omega_i),
					PolarizationFrame<D>(-omega_i, omega_o));

			pdf = r;
		}
		R *= (1 - m_absorption);
	}

	/** Samples a new outgoing ray in the intersection */
	virtual void sample_outgoing_ray(const Intersection<D> &it, Ray<D> &new_ray, PolarizedLightAttenuation<D> &R, Real &pdf) const
	{
		VectorN<D> omega;
		Real ior;

		Real eta, ieta;
		get_eta(it.get_ray().get_direction(), it.get_normal(), eta, ieta);

		Real Rp, Rs, Pp, Ps;
		Fresnel::fresnel_reflectance(it.get_ray().get_direction(), it.get_normal(), ieta, 0., Rp, Rs, Pp, Ps);

		Real r = (Rs + Rp)*.5;

		// Russian roulette
		Real epsilon1 = RNG::StdRandom::get_real();

		if (epsilon1 > r) {
			// Fresnel transmission
			omega = it.get_ray().get_direction().refract(it.get_normal(), eta).normalized();
			pdf = 1 - r;
			R = PolarizedLightAttenuation<D>(1 - r, 
				PolarizationFrame<D>(omega, it.get_ray().get_direction()),
				PolarizationFrame<D>(-it.get_ray().get_direction(), omega));
			ior = (dot(omega, it.get_normal()) < 0.) ? m_n0 : m_n1;
		} else {
			Real sqRsRp = sqrt(Rs*Rp);
			Real A = r,
				B = (Rs - Rp)*.5,
				C = cosf(Ps - Pp)*sqRsRp,
				S = sinf(Ps - Pp)*sqRsRp;

			omega = it.get_ray().get_direction().reflect(it.get_normal());
			ior = it.get_ray().get_refraction_index();

			R = PolarizedLightAttenuation<D>(A, B, 0, 0, B, A, 0, 0, 0, 0, C, S, 0, 0, -S, C,
				PolarizationFrame<D>(omega, it.get_ray().get_direction()),
				PolarizationFrame<D>(-it.get_ray().get_direction(), omega));

			pdf = 1;//r;
		}
		new_ray = Ray<D>(it.get_position(), omega, true,
			it.get_ray().get_level() + 1, ior,
			NULL);

		new_ray.shift();

		R *= (1 - m_absorption);
	}

	/*---------------------------------------------------------------------------*/
	Spectrum get_absorption(const Vector2 &uv) const{ return m_absorption; }

}; //FresnelDielectric

typedef FresnelDielectric<3> FresnelDielectric3D;
typedef FresnelDielectric<2> FresnelDielectric2D;

#endif //_TRANSMISSIVE_H_
