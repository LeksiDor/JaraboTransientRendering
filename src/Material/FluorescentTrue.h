/*
 * Copyright (C) 2017, Victor Arellano (http://giga.cps.unizar.es/~varella/)
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

#ifndef _FLUORESCENT_H_
#define _FLUORESCENT_H_

#include "Material/Material.h"
#include "Material/Reflectance/Lambertian.h"
#include "Color/FluorescentMatrix.h"
#include "Utils/RandomNumbers.h"

/**
 * Fluorescent material.
 */
template<unsigned D>
class FluorescentTrue : public Material<D>
{
	typedef Material<D> MTR;

	/* Let's hardcode! */
	FluorescentMatrix direct;
	FluorescentMatrix reemision;

	/* Pre-baked linear polarized reemision (cool!) */
	FluorescentMatrix linearPolarizedReemision[16];

	Real reemision_time;

	Lambertian<D> m_def_mat;
public:
	FluorescentTrue(FluorescentMatrix direct,FluorescentMatrix reemision,Real reemisionTime)
		: MTR(Reflectance::DIFFUSE), direct(direct), reemision(reemision), reemision_time(reemisionTime), m_def_mat(1.)
	{
		for (size_t i = 0; i < 16; i++)
			linearPolarizedReemision[i] = FluorescentMatrix();

		linearPolarizedReemision[0] = linearPolarizedReemision[1] = linearPolarizedReemision[4] = linearPolarizedReemision[5] = reemision;
	}

	virtual ~FluorescentTrue() {}

	virtual void f(const VectorN<D> &omega_i, const VectorN<D> &omega_o, const VectorN<D> &normal, const Vector2 &uv, Spectrum &R)const
	{
		m_def_mat.f(omega_i, omega_o, normal, uv, R);
	}

	virtual Real p(const VectorN<D> &omega_i, const VectorN<D> &omega_o, const VectorN<D> &normal, const Vector2 &uv)const
	{
		return m_def_mat.p(omega_i,omega_o,normal,uv);
	}

	virtual void sample_direction(const VectorN<D> &omega_i, VectorN<D> &omega_o, const VectorN<D> &normal, const Vector2 &uv, Spectrum &R, Real &pdf) const
	{
		m_def_mat.sample_direction(omega_i, omega_o, normal, uv, R, pdf);
	}

	virtual void sample_outgoing_ray(const Intersection<D> &it, Ray<D> &new_ray, Spectrum &R, Real &pdf) const
	{
		m_def_mat.sample_outgoing_ray(it, new_ray, R, pdf);
	}

	virtual Spectrum get_absorption(const Intersection<D> &it) const
	{
		return m_def_mat.get_absorption(it.get_uv());
	}

	virtual void f(const VectorN<D> &omega_i, const VectorN<D> &omega_o, const VectorN<D> &normal, const Vector2 &uv, PolarizedAttenuation<D> &R, Real &delta_time, Real &pdf) const
	{
		/* Direct reflection */
		Spectrum r;
		delta_time = reemision_time;
		f(omega_i, omega_o, normal, uv, r);

		R = PolarizedAttenuation<D>(
			r, r, 0, 0,
			r, r, 0, 0,
			0, 0, 0, 0,
			0, 0, 0, 0,
			PolarizationFrame<D>(omega_o, omega_i),
			PolarizationFrame<D>(-omega_i, omega_o));
		R *= r;

		pdf = 1;
	}

	/** Samples a new outgoing direction */
	virtual void sample_direction(const VectorN<D> &omega_i, VectorN<D> &omega_o, const VectorN<D> &normal, const Vector2 &uv, PolarizedAttenuation<D> &R, Real &pdf) const
	{
		Spectrum r;
		sample_direction(omega_i, omega_o, normal, uv, r, pdf);
		R = PolarizedAttenuation<D>(r);
	}
	virtual void sample_direction(const VectorN<D> &omega_i, VectorN<D> &omega_o, const VectorN<D> &normal, const Vector2 &uv, PolarizedAttenuation<D> &R, Real &delta_time, Real &pdf) const
	{
		delta_time = 0.;
		sample_direction(omega_i, omega_o, normal, uv, R, pdf);
	}

	/** Samples a new outgoing ray in the intersection */
	virtual void sample_outgoing_ray(const Intersection<D> &it, Ray<D> &new_ray, PolarizedAttenuation<D> &R, Real &delta_time, Real &pdf) const
	{
		/* Direct reflection */
		Spectrum r;
		delta_time = reemision_time;
		sample_outgoing_ray(it, new_ray, r, pdf);

		R = PolarizedAttenuation<D>(
			r, r, 0, 0,
			r, r, 0, 0,
			0, 0, 0, 0,
			0, 0, 0, 0,
			PolarizationFrame<D>(new_ray.get_direction(), it.get_ray().get_direction()),
			PolarizationFrame<D>(-it.get_ray().get_direction(), new_ray.get_direction()));
		R *= r;

		pdf = 1.;
	}

	virtual void f(const VectorN<D> &omega_i, const VectorN<D> &omega_o, const VectorN<D> &normal, const Vector2 &uv, FluorescentPolarizedAttenuation<D> &R, Real &delta_time, Real &pdf) const
	{
		Spectrum r;
		
		/* Russian roulette to determine if light reflected directly
		 * or reemited as fluorescense.
		 */
		Real epsilon1 = Random::StdRNG.next_real();

		Real pdf_time = .5;
		if (epsilon1 < pdf_time) {
			/* Direct reflection */
			delta_time = 0.;
			f(omega_i, omega_o, normal, uv, r);
			R = FluorescentPolarizedAttenuation<D>(direct*r, PolarizationFrame<D>(omega_o, omega_i),
				PolarizationFrame<D>(-omega_i, omega_o));
		} else {
			/* Direct reflection */
			delta_time = reemision_time;
			f(omega_i, omega_o, normal, uv, r);

			R = FluorescentPolarizedAttenuation<D>(
				linearPolarizedReemision[0]*r,  linearPolarizedReemision[1]*r,  linearPolarizedReemision[2]*r,  linearPolarizedReemision[3]*r,
				linearPolarizedReemision[4]*r,  linearPolarizedReemision[5]*r,  linearPolarizedReemision[6]*r,  linearPolarizedReemision[7]*r,
				linearPolarizedReemision[8]*r,  linearPolarizedReemision[9]*r,  linearPolarizedReemision[10]*r, linearPolarizedReemision[11]*r,
				linearPolarizedReemision[12]*r, linearPolarizedReemision[13]*r, linearPolarizedReemision[14]*r, linearPolarizedReemision[15]*r,
				PolarizationFrame<D>(omega_o, omega_i), PolarizationFrame<D>(-omega_i, omega_o));
		}
		
		pdf = pdf_time;
	}

	/** Samples a new outgoing direction */
	virtual void sample_direction(const VectorN<D> &omega_i, VectorN<D> &omega_o, const VectorN<D> &normal, const Vector2 &uv, FluorescentPolarizedAttenuation<D> &R, Real &pdf) const
	{
		Spectrum r;
		sample_direction(omega_i, omega_o, normal, uv, r, pdf);
		R = FluorescentPolarizedAttenuation<D>(r);
	}
	virtual void sample_direction(const VectorN<D> &omega_i, VectorN<D> &omega_o, const VectorN<D> &normal, const Vector2 &uv, FluorescentPolarizedAttenuation<D> &R, Real &delta_time, Real &pdf) const
	{
		delta_time = 0.;
		sample_direction(omega_i, omega_o, normal, uv, R, pdf);
	}

	/** Samples a new outgoing ray in the intersection */
	virtual void sample_outgoing_ray(const Intersection<D> &it, Ray<D> &new_ray, FluorescentPolarizedAttenuation<D> &R, Real &delta_time, Real &pdf) const
	{
		Spectrum r;
		/* Russian roulette to determine if light reflected directly
		 * or reemited as fluorescense.
		 */
		Real epsilon1 = Random::StdRNG.next_real();

		Real pdf_time = .5;
		if (epsilon1 < pdf_time) {
			/* Direct reflection */
			delta_time = 0.;
			sample_outgoing_ray(it, new_ray, r, pdf);
			R = FluorescentPolarizedAttenuation<D>(direct*r, PolarizationFrame<D>(new_ray.get_direction(), it.get_ray().get_direction()),
				PolarizationFrame<D>(-it.get_ray().get_direction(), new_ray.get_direction()));
		} else {
			/* Direct reflection */
			delta_time = reemision_time;
			sample_outgoing_ray(it, new_ray, r, pdf);

			R = FluorescentPolarizedAttenuation<D>(
				linearPolarizedReemision[0]*r,  linearPolarizedReemision[1]*r,  linearPolarizedReemision[2]*r,  linearPolarizedReemision[3]*r,
				linearPolarizedReemision[4]*r,  linearPolarizedReemision[5]*r,  linearPolarizedReemision[6]*r,  linearPolarizedReemision[7]*r,
				linearPolarizedReemision[8]*r,  linearPolarizedReemision[9]*r,  linearPolarizedReemision[10]*r, linearPolarizedReemision[11]*r,
				linearPolarizedReemision[12]*r, linearPolarizedReemision[13]*r, linearPolarizedReemision[14]*r, linearPolarizedReemision[15]*r,
				PolarizationFrame<D>(new_ray.get_direction(), it.get_ray().get_direction()),
				PolarizationFrame<D>(-it.get_ray().get_direction(), new_ray.get_direction()));
		}
		
		pdf = pdf_time;
	}
	
	virtual void f(const VectorN<D> &omega_i, const VectorN<D> &omega_o, const VectorN<D> &normal, const Vector2 &uv, FluorescentAttenuation<D> &R, Real &delta_time, Real &pdf) const
	{
		Spectrum r;

		/* Russian roulette to determine if light reflected directly
		 * or reemited as fluorescense.
		 */
		Real epsilon1 = Random::StdRNG.next_real();

		Real pdf_time = .5;
		if (epsilon1 < pdf_time) {
			/* Direct reflection */
			delta_time = 0.;
			f(omega_i, omega_o, normal, uv, r);
			R = FluorescentAttenuation<D>(direct);
			R *= r;
		} else {
			/* Direct reflection */
			delta_time = reemision_time;
			f(omega_i, omega_o, normal, uv, r);

			R = FluorescentAttenuation<D>(reemision);
			R *= r;
		}

		pdf = pdf_time;
	}

	/** Samples a new outgoing direction */
	virtual void sample_direction(const VectorN<D> &omega_i, VectorN<D> &omega_o, const VectorN<D> &normal, const Vector2 &uv, FluorescentAttenuation<D> &R, Real &pdf) const
	{
		Spectrum r;
		sample_direction(omega_i, omega_o, normal, uv, r, pdf);
		R = FluorescentAttenuation<D>(r);
	}
	virtual void sample_direction(const VectorN<D> &omega_i, VectorN<D> &omega_o, const VectorN<D> &normal, const Vector2 &uv, FluorescentAttenuation<D> &R, Real &delta_time, Real &pdf) const
	{
		delta_time = 0.;
		sample_direction(omega_i, omega_o, normal, uv, R, pdf);
	}

	/** Samples a new outgoing ray in the intersection */
	virtual void sample_outgoing_ray(const Intersection<D> &it, Ray<D> &new_ray, FluorescentAttenuation<D> &R, Real &delta_time, Real &pdf) const
	{
		Spectrum r;
		/* Russian roulette to determine if light reflected directly
		 * or reemited as fluorescense.
		 */
		Real epsilon1 = Random::StdRNG.next_real();

		Real pdf_time = .5;
		if (epsilon1 < pdf_time) {
			/* Direct reflection */
			delta_time = 0.;
			sample_outgoing_ray(it, new_ray, r, pdf);
			R = FluorescentAttenuation<D>(direct);
			R *= r;
		} else {
			/* Direct reflection */
			delta_time = reemision_time;
			sample_outgoing_ray(it, new_ray, r, pdf);

			R = FluorescentAttenuation<D>(reemision);
			R *= r;
		}

		pdf = pdf_time;
	}
}; /* Fluorescent */

typedef FluorescentTrue<3> Fluorescent3D;
typedef FluorescentTrue<2> Fluorescent2D;

#endif //_FLUORESCENT_H_
