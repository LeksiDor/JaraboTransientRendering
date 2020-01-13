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

#ifndef _FRESNEL_CONDUCTOR_H_
#define _FRESNEL_CONDUCTOR_H_

#include "bunnykiller.h"

#include "Material/Material.h"
#include "Material/Reflectance/Fresnel.h"
#include "RayTracing/Ray.h"
#include "RayTracing/Intersection.h"

template<unsigned D>
class FresnelConductor : public Material<D>
{
protected:
	typedef Material<D> MTR;

protected:
	/* Index of refraction */
	Real m_n0, m_n1, m_k1;

	void get_eta(const VectorN<D> &omega, const VectorN<D> &normal, Real &n, Real &in) const
	{
		n = (dot(omega, normal) > 0.) ? (m_n1 / m_n0) : (m_n0 / m_n1);
		in = 1. / n;
	}

public:
	FresnelConductor() :
			Material<D>(Reflectance::DELTA),
			m_n0(MTR::m_default_n),
			m_n1(MTR::m_default_n),
			m_k1(0.)
	{
	}
	FresnelConductor(const Real n_in, const Real k_in, const Real n_out = 0.) :
			Material<D>(Reflectance::DELTA),
			m_n0(n_out ? n_out : MTR::m_default_n),
			m_n1(n_in),
			m_k1(k_in)
	{
	}

	virtual ~FresnelConductor()
	{
	}

	/** PDF */
	Real p(const VectorN<D>&, const VectorN<D>&, const VectorN<D>&, const Vector2&) const
	{
		return 0.;
	}

	/** Spectrum-based attenuation */
	void f(const VectorN<D>&, const VectorN<D>&, const VectorN<D>&, const Vector2&,
		Spectrum &R) const
	{
		R = Spectrum(0.);
	}

	void sample_direction(const VectorN<D>& omega_i, VectorN<D>& omega_o, const VectorN<D>& normal,
		const Vector2&, Spectrum& R, Real& pdf) const
	{
		Real eta, ieta;
		get_eta(omega_i, normal, eta, ieta);

		Real r = Fresnel::fresnel_reflectance(omega_i, normal, ieta, m_k1);

		omega_o = omega_i.reflect(normal);
		pdf = 1.;
		R = r;
	}

	void sample_outgoing_ray(const Intersection<D> &it, Ray<D> &new_ray, Spectrum &R,
		Real &pdf) const
	{
		const Ray<D>& old_ray = it.get_ray();

		Real eta, ieta;
		get_eta(old_ray.get_direction(), it.get_normal(), eta, ieta);

		Real r = Fresnel::fresnel_reflectance(old_ray.get_direction(), it.get_normal(), ieta, m_k1);

		VectorN<D> omega = old_ray.get_direction().reflect(it.get_normal());
		pdf = 1.;
		R = r;

		// Test if the ray enters the media or is leaving
		new_ray = Ray<D>(it.get_position(), omega, true, old_ray.get_level() + 1, old_ray.get_ior(),
				old_ray.get_medium());
	}

	/** Polarized spectrum-based attenuation */
	virtual void f(const VectorN<D> &omega_i, const VectorN<D> &omega_o, const VectorN<D> &,
		const Vector2 &, PolarizedAttenuation<D> &R) const
	{
		R = PolarizedAttenuation<D>(0., PolarizationFrame<D>(omega_o, omega_i),
				PolarizationFrame<D>(-omega_i, omega_o));
	}

	/** Samples a new outgoing direction */
	virtual void sample_direction(const VectorN<D> &omega_i, VectorN<D> &omega_o,
		const VectorN<D> &normal, const Vector2 &, PolarizedAttenuation<D> &R, Real &pdf) const
	{
		Real eta, ieta;
		get_eta(omega_i, normal, eta, ieta);

		Real Rp, Rs, Pp, Ps;
		Fresnel::fresnel_reflectance(omega_i, normal, ieta, m_k1, Rp, Rs, Pp, Ps);

		Real sqRsRp = sqrt(Rs * Rp);
		Real A = (Rs + Rp) * .5, B = (Rs - Rp) * .5, C = cos(Ps - Pp) * sqRsRp, S = sin(Ps - Pp)
				* sqRsRp;

		omega_o = omega_i.reflect(normal);

		R = PolarizedAttenuation<D>(A, B, 0., 0., B, A, 0., 0., 0., 0., C, S, 0., 0., -S, C,
				PolarizationFrame<D>(omega_o, omega_i), PolarizationFrame<D>(-omega_i, omega_o));
		pdf = 1.;
	}

	/** Samples a new outgoing ray in the intersection */
	virtual void sample_outgoing_ray(const Intersection<D> &it, Ray<D> &new_ray,
		PolarizedAttenuation<D> &R, Real &pdf) const
	{
		const Ray<D>& old_ray = it.get_ray();

		Real eta, ieta;
		get_eta(old_ray.get_direction(), it.get_normal(), eta, ieta);

		Real Rp, Rs, Pp, Ps;
		Fresnel::fresnel_reflectance(old_ray.get_direction(), it.get_normal(), ieta, m_k1, Rp, Rs,
				Pp, Ps);

		Real sqRsRp = sqrt(Rs * Rp);
		Real A = (Rs + Rp) * .5, B = (Rs - Rp) * .5, C = cos(Ps - Pp) * sqRsRp, S = sin(Ps - Pp)
				* sqRsRp;
		
		VectorN<D> omega = old_ray.get_direction().reflect(it.get_normal());

		R = PolarizedAttenuation<D>(A, B, 0., 0., B, A, 0., 0., 0., 0., C, S, 0., 0., -S, C,
				PolarizationFrame<D>(omega, old_ray.get_direction()),
				PolarizationFrame<D>(-old_ray.get_direction(), omega));
		pdf = 1.;

		new_ray = Ray<D>(it.get_position(), omega, true, old_ray.get_level() + 1, old_ray.get_ior(),
				old_ray.get_medium());
	}

	Spectrum get_absorption(const Intersection<D> &it) const
	{
		Real eta, ieta;
		get_eta(it.get_ray().get_direction(), it.get_normal(), eta, ieta);

		Real Rp, Rs, Pp, Ps;
		Fresnel::fresnel_reflectance(it.get_ray().get_direction(), it.get_normal(), ieta, m_k1, Rp,
				Rs, Pp, Ps);

		return Spectrum(1. - (Rs + Rp) * .5);
	}
}; // FresnelConductor

typedef FresnelConductor<3> FresnelConductor3D;
typedef FresnelConductor<2> FresnelConductor2D;

#endif //_FRESNEL_CONDUCTOR_H_
