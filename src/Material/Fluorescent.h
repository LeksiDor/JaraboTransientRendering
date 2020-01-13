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
#include "Utils/RandomNumbers.h"

// Fluorescent material.
template<unsigned D>
class Fluorescent: public Material<D>
{
	// Let's hardcode!
	Spectrum m_absorption;
	Spectrum m_scattering;
	Spectrum m_albedo;
	Lambertian<D> m_reflected;
	Real m_C;
	Real m_delta_time;
	Real m_quantum_field;

	Spectrum reemit(const Spectrum &s0)const
	{
		// Rough approximation of Eq.17 in [Gutierrez 2005]
		Spectrum val_lambda = m_absorption*s0*m_quantum_field*Spectrum(680,550,450)/680.;
		return Spectrum(val_lambda.sum(),0.,0.);
	}
public:
	Fluorescent(Real C = 100, Real delta_time = 100000):Material(Reflectance::DIFFUSE),m_reflected(1.), 
														m_C(C),m_delta_time(delta_time),m_quantum_field(.05)
	{
		// Really rough approximation, R=680; G=550; B=450
		m_absorption = m_C*Spectrum(0.015,0.007,0.035); //Values obtained from Tb.1 in [Gutierrez 088]
		m_scattering = 550.*.3*powf(m_C,0.62f)/Spectrum(680,550,450); //Eq.11 [Gutierrez 088]
		m_albedo = m_scattering/(m_scattering+m_absorption);
	}

	virtual ~Fluorescent() {}

	/** Computes the local direct reflected (or refracted) light. 
		omega_o is the view (incoming) direction, while omega_i is light (outgoing)*/
	virtual Spectrum f( const VectorN<D> &omega_o, const VectorN<D> &omega_i, const VectorN<D> &normal, const Vector2 &uv ) const
	{
		return m_reflected.f(omega_o,omega_i,normal,uv)*m_albedo;
	}
	/** Computes the pdf */
	virtual Real p(const VectorN<D> &omega_o, const VectorN<D> &omega_i, const VectorN<D> &normal, const Vector2 &uv )const
	{
		return m_reflected.p(omega_o,omega_i,normal,uv);
	}
	/** Samples a new outgoing direction */
	virtual Spectrum sample_direction( const VectorN<D> &omega_o, VectorN<D> &omega_i, const VectorN<D> &normal, const Vector2 &uv, Real &pdf ) const
	{
		return m_reflected.sample_direction(omega_o, omega_i, normal, uv, pdf)*m_albedo;
	}
	/** Samples a new outgoing ray in the intersection */	
	virtual Spectrum sample_outgoing_ray( const Intersection<D> &it, Ray<D> &new_ray, Real &pdf ) const
	{
		return m_reflected.sample_outgoing_ray(it, new_ray, pdf)*m_albedo;
	}
	
	// TIME RESOLVED MATERIAL FUNCTIONS
	/** Computes the local direct reflected (or refracted) light. 
		omega_o is the view (incoming) direction, while omega_i is light (outgoing)*/
	virtual Spectrum f( const Spectrum &L_i, const VectorN<D> &omega_o, const VectorN<D> &omega_i, const VectorN<D> &normal, const Vector2 &uv, Real &delta_time, Real &pdf  ) const 
	{	
		// Russian roulette to determine if light reflected directly
		// or reemited as fluorescense.
		Real epsilon1 = RNG::StdRandom::get_real();
		while( epsilon1 <= 0. )
			epsilon1 = RNG::StdRandom::get_real();

		Real pdf_time = .5;
		if (epsilon1 < pdf_time) {
			//Direct reflection
			delta_time = 0.;
			pdf = pdf_time;

			return m_reflected.f(omega_o,omega_i,normal,uv)*m_albedo*L_i;
		} else {
			//Direct reflection
			delta_time = m_delta_time;
			pdf = pdf_time;
			return m_reflected.f(omega_o,omega_i,normal,uv)*reemit(L_i);
		}
	}

	/** Samples a new outgoing direction */
	virtual Spectrum sample_direction( const VectorN<D> &omega_o, VectorN<D> &omega_i, const VectorN<D> &normal, const Vector2 &uv, Real &delta_time, Real &pdf ) const
	{	
		delta_time = 0.;
		return sample_direction(omega_o,omega_i,normal,uv,pdf);
	}

	/** Samples a new outgoing ray in the intersection */	
	virtual Spectrum sample_outgoing_ray( const Intersection<D> &it, Ray<D> &new_ray, Real &delta_time, Real &pdf ) const
	{
		delta_time = 0.;
		return sample_outgoing_ray(it, new_ray, pdf);
	}

	/** Get the absorption of the material */
	virtual Spectrum get_absorption( const Vector2 &uv ) const
	{
		return m_absorption;
	}
}; //Fluorescent

typedef Fluorescent<3> Fluorescent3D;
typedef Fluorescent<2> Fluorescent2D;

#endif //_FLUORESCENT_H_
