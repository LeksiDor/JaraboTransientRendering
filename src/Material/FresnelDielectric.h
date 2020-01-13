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

#ifndef _FRESNEL_DIELECTRIC_H_
#define _FRESNEL_DIELECTRIC_H_

#include "bunnykiller.h"

#include "Material/Material.h"
#include "Material/Reflectance/Fresnel.h"
#include "RayTracing/Ray.h"
#include "RayTracing/Intersection.h"


template<unsigned D>
class FresnelDielectric : public Material<D>
{
protected:
	typedef Material<D> MTR;

protected:
	// Absorption coefficient
	// Note that this absorption has no real physical sense, and complex IOR or an absorbing media
	// should be used instead. However, it allows for easily coloring the interactions.
	Spectrum m_absorption;

	// Medium inside and outside the material
	Medium<D>* m_medium_0;
	Medium<D>* m_medium_1;

	//Index of refraction
	Real m_n0;
	Real m_n1;

	void get_eta(const VectorN<D> &omega, const VectorN<D> &normal, Real &n, Real &in) const
	{
		n = (dot(omega, normal) > 0.) ? (m_n1 / m_n0) : (m_n0 / m_n1);
		in = 1/n;
	}

//Reflectance::DELTA_TRANSMISSION
public:
	FresnelDielectric() :
		Material<D>(Reflectance::DELTA), m_absorption(0.f),
		m_medium_0(MTR::m_default_medium), m_medium_1(nullptr),
		m_n0(MTR::m_default_n), m_n1(MTR::m_default_n)
	{}

	FresnelDielectric(const Spectrum &absorption, const Real n_in, const Real n_out = 0.) :
		Material<D>(Reflectance::DELTA), m_absorption(absorption),
		m_medium_0(MTR::m_default_medium), m_medium_1(nullptr),
		m_n0(n_out ? n_out : MTR::m_default_n), m_n1(n_in)
	{}

	FresnelDielectric(const Real n_in, const Real n_out = 0) :
		Material<D>(Reflectance::DELTA), m_absorption(0.),
		m_medium_0(MTR::m_default_medium), m_medium_1(nullptr),
		m_n0(n_out ? n_out : MTR::m_default_n), m_n1(n_in)
	{}

	FresnelDielectric(const Real n_in, const Real n_out, Medium<D> *m_in, Medium<D> *m_out = nullptr) :
		Material<D>(Reflectance::DELTA), m_absorption(0.),
		m_medium_0(m_out ? m_out : MTR::m_default_medium), m_medium_1(m_in),
		m_n0(n_out ? n_out : MTR::m_default_n), m_n1(n_in)
	{}

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

	void sample_direction(const VectorN<D> &omega_i, VectorN<D> &omega_o, const VectorN<D> &normal, const Vector2 &, Spectrum &R, Real &pdf) const
	{
		Real eta, ieta;
		get_eta(omega_i, normal, eta, ieta);

		Real r = Fresnel::fresnel_reflectance(omega_i, normal, ieta, 0.);

		// Russian roulette
		Real epsilon1 = Random::StdRNG.next_real();
		
		if (epsilon1 < r) {
			// Fresnel reflection
			omega_o = omega_i.reflect(normal);
			pdf = r;
			R = r;
		} else {
			// Refraction
			pdf = 1-r;
			R = 1 - r;
			omega_o = omega_i.refract(normal, eta); omega_o.normalize();
		}

		R *= (1. - m_absorption);
	}
	
	void sample_outgoing_ray(const Intersection<D> &it, Ray<D> &new_ray, Spectrum &R, Real &pdf) const
	{
		//printf("Fresnel Dielectric - ");
		VectorN<D> omega;
		Real eta, ieta;
		get_eta(it.get_ray().get_direction(), it.get_normal(), eta, ieta);


		Real r = Fresnel::fresnel_reflectance(it.get_ray().get_direction(), it.get_normal(), ieta);
		Real ior;
		Medium<D> *medium_n;

		// Russian roulette
		Real epsilon1 = Random::StdRNG.next_real();
		
		/* Test if the ray enters the media or if it is leaving */
		if (epsilon1 < r) { /* Reflection */
			omega = it.get_ray().get_direction().reflect(it.get_normal());
			pdf = r;
			R = r;
			ior = it.get_ray().get_ior();
			medium_n = it.get_ray().get_medium();
		} else { /* Refraction */
			omega = it.get_ray().get_direction().refract(it.get_normal(), eta);
			omega.normalize();
			pdf =  1 - r;
			R = 1 - r;
			ior = (dot(omega, it.get_normal()) > 0.) ? m_n0 : m_n1;
			medium_n = (dot(omega, it.get_normal()) > 0.) ? m_medium_0 : m_medium_1;
		}

		new_ray = Ray<D>(it.get_position(), omega, true, it.get_ray().get_level() + 1,
				ior, medium_n);

		R *= (1. - m_absorption);
	}

	/*---------------------------------------------------------------------------*/
	/** Polarized spectrum-based attenuation */
	virtual void f(const VectorN<D> &omega_i, const VectorN<D> &omega_o, const VectorN<D> &,
			const Vector2 &, PolarizedAttenuation<D> &R) const
	{
		R = PolarizedAttenuation<D>(0.,
			PolarizationFrame<D>(omega_o, omega_i),
			PolarizationFrame<D>(-omega_i, omega_o));
	}

	/** Samples a new outgoing direction */
	virtual void sample_direction(const VectorN<D> &omega_i, VectorN<D> &omega_o,
			const VectorN<D> &normal, const Vector2 &, PolarizedAttenuation<D> &R, Real &pdf) const
	{
		Real eta, ieta;
		get_eta(omega_i, normal, eta, ieta);

		Real Rp, Rs, Pp, Ps;
		Fresnel::fresnel_reflectance(omega_i, normal, ieta, 0., Rp, Rs, Pp, Ps);

		Real r = (Rs + Rp)*.5;

		// Russian roulette
		Real epsilon1 = Random::StdRNG.next_real();

		if (epsilon1 > r) { // Fresnel transmission
			omega_o = omega_i.refract(normal, eta);
			omega_o.normalize();

			Real RsRp = Rs*Rp;
			Real A = 1 - r,
				B = (Rp - Rs)*.5,
				C = (RsRp < 0) ? Real(0.) : sqrt(RsRp),
				S = 0.;

			pdf = 1 - r;
			R = PolarizedAttenuation<D>(
					A, B,  0, 0,
					B, A,  0, 0,
					0, 0,  C, S,
					0, 0, -S, C,
					PolarizationFrame<D>(omega_o, omega_i),
					PolarizationFrame<D>(-omega_i, omega_o));
		} else { // Fresnel reflectance
			Real RsRp = Rs*Rp;
			Real sqRsRp = (RsRp < 0) ? Real(0.) : sqrt(RsRp);
			Real A = r,
				B = (Rs - Rp)*.5,
				C = cosf(Ps - Pp)*sqRsRp,
				S = sinf(Ps - Pp)*sqRsRp;

			omega_o = omega_i.reflect(normal);

			R = PolarizedAttenuation<D>(
					A, B,  0, 0,
					B, A,  0, 0,
					0, 0,  C, S,
					0, 0, -S, C,
					PolarizationFrame<D>(omega_o, omega_i),
					PolarizationFrame<D>(-omega_i, omega_o));

			pdf = r;
		}
		R *= (1. - m_absorption);
	}

	/** Samples a new outgoing ray in the intersection */
	virtual void sample_outgoing_ray(const Intersection<D> &it, Ray<D> &new_ray,
			PolarizedAttenuation<D> &R, Real &pdf) const
	{
		VectorN<D> omega;
		Real ior;

		Real eta, ieta;
		Medium<D> *medium_n;

		get_eta(it.get_ray().get_direction(), it.get_normal(), eta, ieta);

		Real Rp, Rs, Pp, Ps;
		Fresnel::fresnel_reflectance(it.get_ray().get_direction(), it.get_normal(), ieta, Rp, Rs, Pp, Ps);

		Real r = (Rs + Rp)*.5;

		// Russian roulette
		Real epsilon1 = Random::StdRNG.next_real();

		if (epsilon1 > r) { // Fresnel transmission
			omega = it.get_ray().get_direction().refract(it.get_normal(), eta);
			omega.normalize();

			Real RsRp = Rs*Rp;
			Real A = 1-r,
				 B = (Rp - Rs)*.5,
				 C = (RsRp < 0) ? Real(0.) : sqrt(RsRp),
				 S = 0.;
			
			pdf = 1 - r;
			R = PolarizedAttenuation<D>(
					A, B,  0, 0,
					B, A,  0, 0,
					0, 0,  C, S,
					0, 0, -S, C,
				    PolarizationFrame<D>(omega, it.get_ray().get_direction()),
				    PolarizationFrame<D>(-it.get_ray().get_direction(), omega));

			ior = (dot(omega, it.get_normal()) > 0.) ? m_n0 : m_n1;
			medium_n = (dot(omega, it.get_normal()) > 0.) ? m_medium_0 : m_medium_1;
		} else { // Fresnel reflectance
			Real RsRp = Rs*Rp;
			Real sqRsRp = (RsRp < 0) ? Real(0.) : sqrt(RsRp);
			Real A = r,
				 B = (Rs - Rp)*.5,
				 C = cos(Ps - Pp)*sqRsRp,
				 S = sin(Ps - Pp)*sqRsRp;

			omega = it.get_ray().get_direction().reflect(it.get_normal());
			omega.normalize();

			ior = it.get_ray().get_ior();
			medium_n = it.get_ray().get_medium();


			R = PolarizedAttenuation<D>(
					A, B,  0, 0,
					B, A,  0, 0,
					0, 0,  C, S,
					0, 0, -S, C,
					PolarizationFrame<D>(omega, it.get_ray().get_direction()),
				    PolarizationFrame<D>(-it.get_ray().get_direction(), omega));

			pdf = r;
		}
		new_ray = Ray<D>(it.get_position(), omega, true,
			it.get_ray().get_level() + 1, ior,
			medium_n);

		R *= (1. - m_absorption);
	}

	/*POLARIZED FLUORESCENCE I-DONT-LIKE-COPYPASTING-CODE-SQUAD*/
	/*---------------------------------------------------------------------------*/
	/** Polarized spectrum-based attenuation */
	virtual void f(const VectorN<D> &omega_i, const VectorN<D> &omega_o, const VectorN<D> &, const Vector2 &, FluorescentPolarizedAttenuation<D> &R) const
	{
		R = FluorescentPolarizedAttenuation<D>(0.,
			PolarizationFrame<D>(omega_o, omega_i),
			PolarizationFrame<D>(-omega_i, omega_o));
	}

	/** Samples a new outgoing direction */
	virtual void sample_direction(const VectorN<D> &omega_i, VectorN<D> &omega_o, const VectorN<D> &normal, const Vector2 &, FluorescentPolarizedAttenuation<D> &R, Real &pdf) const
	{
		Real eta, ieta;
		get_eta(omega_i, normal, eta, ieta);

		Real Rp, Rs, Pp, Ps;
		Fresnel::fresnel_reflectance(omega_i, normal, ieta, 0., Rp, Rs, Pp, Ps);

		Real r = (Rs + Rp)*.5;

		// Russian roulette
		Real epsilon1 = Random::StdRNG.next_real();

		if (epsilon1 > r) {
			// Fresnel transmission
			omega_o = omega_i.refract(normal, eta);
			omega_o.normalize();

			Real RsRp = Rs*Rp;
			Real sqRsRp = (RsRp < 0) ? Real(0.) : sqrt(RsRp);
			Real A = 1 - r,
				B = (Rp - Rs)*.5,
				C = sqRsRp,
				S = 0.;

			pdf = 1 - r;
			R = FluorescentPolarizedAttenuation<D>(
					A, B,  0, 0,
					B, A,  0, 0,
					0, 0,  C, S,
					0, 0, -S, C,
				    PolarizationFrame<D>(omega_o, omega_i),
				    PolarizationFrame<D>(-omega_i, omega_o));
		} else {
			Real RsRp = Rs*Rp;
			Real sqRsRp = (RsRp < 0) ? Real(0.) : sqrt(RsRp);
			Real A = r,
				 B = (Rs - Rp)*.5,
				 C = cos(Ps - Pp)*sqRsRp,
				 S = sin(Ps - Pp)*sqRsRp;

			omega_o = omega_i.reflect(normal);

			R = FluorescentPolarizedAttenuation<D>(
					A, B,  0, 0,
					B, A,  0, 0,
					0, 0,  C, S,
					0, 0, -S, C,
				    PolarizationFrame<D>(omega_o, omega_i),
				    PolarizationFrame<D>(-omega_i, omega_o));

			pdf = r;
		}

		R *= Spectrum(1. - m_absorption);
	}

	/** Samples a new outgoing ray in the intersection */
	virtual void sample_outgoing_ray(const Intersection<D> &it, Ray<D> &new_ray, FluorescentPolarizedAttenuation<D> &R, Real &pdf) const
	{
		VectorN<D> omega;
		Real ior;

		Real eta, ieta;
		Medium<D> *medium_n;

		get_eta(it.get_ray().get_direction(), it.get_normal(), eta, ieta);

		Real Rp, Rs, Pp, Ps;
		Fresnel::fresnel_reflectance(it.get_ray().get_direction(), it.get_normal(), ieta, Rp, Rs, Pp, Ps);

		Real r = (Rs + Rp)*.5;

		// Russian roulette
		Real epsilon1 = Random::StdRNG.next_real();

		if (epsilon1 > r) {
			// Fresnel transmission
			omega = it.get_ray().get_direction().refract(it.get_normal(), eta);
			omega.normalize();

			Real RsRp = Rs*Rp;
			Real A = 1 - r,
				B = (Rp - Rs)*.5,
				C = (RsRp < 0) ? Real(0.) : sqrt(RsRp),
				S = 0.;

			pdf = 1 - r;
			R = FluorescentPolarizedAttenuation<D>(
					A, B,  0, 0,
					B, A,  0, 0,
					0, 0,  C, S,
					0, 0, -S, C,
				    PolarizationFrame<D>(omega, it.get_ray().get_direction()),
				    PolarizationFrame<D>(-it.get_ray().get_direction(), omega));

			ior = (dot(omega, it.get_normal()) > 0.) ? m_n0 : m_n1;
			medium_n = (dot(omega, it.get_normal()) > 0.) ? m_medium_0 : m_medium_1;

		} else {
			Real RsRp = Rs*Rp;
			Real sqRsRp = (RsRp < 0.) ? Real(0.) : sqrt(RsRp);
			Real A = r,
				B = (Rs - Rp)*.5,
				C = cosf(Ps - Pp)*sqRsRp,
				S = sinf(Ps - Pp)*sqRsRp;

			omega = it.get_ray().get_direction().reflect(it.get_normal());
			omega.normalize();

			ior = it.get_ray().get_ior();
			medium_n = it.get_ray().get_medium();


			R = FluorescentPolarizedAttenuation<D>(
					A, B,  0, 0,
					B, A,  0, 0,
					0, 0,  C, S,
					0, 0, -S, C,
				    PolarizationFrame<D>(omega, it.get_ray().get_direction()),
				    PolarizationFrame<D>(-it.get_ray().get_direction(), omega));

			pdf = r;
		}

		new_ray = Ray<D>(it.get_position(), omega, true, it.get_ray().get_level() + 1, ior, medium_n);

		R *= Spectrum(1. - m_absorption);
	}

	/*---------------------------------------------------------------------------*/
	Spectrum get_absorption(const Intersection<D> &) const
	{
		return m_absorption;
	}
}; //FresnelDielectric

typedef FresnelDielectric<3> FresnelDielectric3D;
typedef FresnelDielectric<2> FresnelDielectric2D;

#endif //_TRANSMISSIVE_H_
