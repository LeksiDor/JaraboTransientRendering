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

#ifndef _BIREFRINGENCE_H_
#define _BIREFRINGENCE_H_

#include "Material/Material.h"
#include "RayTracing/Intersection.h"
#include "RayTracing/Ray.h"
#include "Material/Reflectance/FresnelDielectric.h"

// Birefringent material. We assume a fully polarized light.
// NOTE: this is just partially based on actual Birefringency, and its
// an almost pure approximation (at best).
template<unsigned int D>
class Birefringence: public Material<D>
{
	// Let's hardcode!
	Spectrum m_absorption;
	VectorN<D> m_vector_field;
	VectorN<D> m_up_vector; // Note that this is only needed for D=3

	Real m_no, m_ne;
	FresnelDielectric<D> m_ordinary_refraction;

	// Following Weildich and Wilkie (Eq.2)
	Real compute_extraordinary_ior( const Real no, const Real ne, const Real dot_thita ) const
	{
		Real cos_thita2 = dot_thita*dot_thita;
		Real sin_thita2 = 1.-cos_thita2;

		return no*ne/sqrtf(no*no*sin_thita2 + ne*ne*cos_thita2);
	}

	// Crappy approximation, oh well... It'll look like good... Sort of, at least...
	VectorN<D> compute_normal( const VectorN<D> &normal0 ) const
	{
		VectorN<D> normal1 = m_up_vector + (m_no-m_ne)*m_vector_field;

		// Transform the world to match the reference system...
		return normal1.transform_matrix_to( m_up_vector, normal0 );		
	}
public:
	Birefringence(const Spectrum &absorption, const Real ior_o, const Real ior_e, const VectorN<D> vector_field )
		: Material(Reflectance::DELTA_TRANSMISSION), m_absorption(absorption), m_vector_field(vector_field),
	   	m_up_vector((D==3)?VectorN<D>(0.,0.,1.):VectorN<D>(0.)), m_no(ior_o), m_ne(ior_e),m_ordinary_refraction(ior_o)
	{
		m_vector_field.normalize();
	}

	virtual ~Birefringence() {}

	/** Computes the local direct reflected (or refracted) light. 
		omega_o is the view (incoming) direction, while omega_i is light (outgoing)*/
	virtual Spectrum f( const VectorN<D> &omega_o, const VectorN<D> &omega_i, const VectorN<D> &normal, const Vector2 &uv ) const
	{	
		if( dot(omega_o.refract(normal, m_no), omega_i) > 0.99 )
			return (1.-m_absorption)/dot_abs(normal, omega_i);

		return Spectrum(0.);
	}
	/** Computes the pdf */
	virtual Real p(const VectorN<D> &omega_o, const VectorN<D> &omega_i, const VectorN<D> &normal, const Vector2 &uv )const
	{
		return 0.;
	}

	/** Samples a new outgoing direction */
	virtual Spectrum sample_direction( const VectorN<D> &omega_o, VectorN<D> &omega_i, const VectorN<D> &normal, const Vector2 &uv, Real &pdf ) const
	{	
		//printf("Sample Direction\n");
		Real pdf_phasor = 0.5;
		Spectrum fvalue(0.);

		if( abs(normal[2]) < (1.-1.e-2) )
			pdf_phasor = 1;

		if( dot(normal,omega_o)>0. )
			pdf_phasor = 1;

		// Russian roulette to choose between ordinary and extraordinary ray
		Real epsilon1 = RNG::StdRandom::get_real();
		while( epsilon1 <= 0. )
			epsilon1 = RNG::StdRandom::get_real();
		
		if ( epsilon1 <= pdf_phasor ) {
			// Ordinary transmission
			fvalue = m_ordinary_refraction.sample_direction( omega_o, omega_i, normal, uv, pdf );
		} else {
			// Extraordinary transmission
			Real n2 = compute_extraordinary_ior(m_no, m_ne, dot_abs(-omega_o,normal));
			
			// And here is where the crap begins!
			VectorN<D> normal_e = compute_normal(normal);
			FresnelDielectric<D> extraordinary_refraction(n2);

			fvalue = extraordinary_refraction.sample_direction( omega_o, omega_i, normal_e, uv, pdf );
		}

		return (1.-m_absorption)*fvalue;
	}

	/** Samples a new outgoing ray in the intersection */	
	virtual Spectrum sample_outgoing_ray( const Intersection<D> &it, Ray<D> &new_ray, Real &pdf ) const
	{	
		//printf("Sample Outgoing Ray\n");
		Real pdf_phasor = .5;
		Spectrum fvalue(0.);

		if( it.get_normal()[2] < (1.-1.e-2) )
			pdf_phasor = 1;

		// Russian roulette to choose between ordinary and extraordinary ray
		Real epsilon1 = RNG::StdRandom::get_real();
		while( epsilon1 <= 0. )
			epsilon1 = RNG::StdRandom::get_real();
		
		if( epsilon1 < pdf_phasor ) {
			// Ordinary transmission
			fvalue = m_ordinary_refraction.sample_outgoing_ray( it, new_ray, pdf );
		} else {
			// Extraordinary transmission
			Real n2 = compute_extraordinary_ior(m_no, m_ne, dot_abs(it.get_ray().get_direction(), it.get_normal()));
			
			// And here is where the crap begins!
			VectorN<D> normal_e = compute_normal(it.get_normal());
			FresnelDielectric<D> extraordinary_refraction(n2);

			Intersection<D> it1(it.get_ray(),it.intersected(),normal_e,it.get_uv());

			fvalue = extraordinary_refraction.sample_outgoing_ray( it1, new_ray, pdf );
		}

		return (1.-m_absorption)*fvalue;
	}
	
	// TIME RESOLVED MATERIAL FUNCTIONS
	// For now we assume no Fresnel retardance...

	/** Get the absorption of the material */
	virtual Spectrum get_absorption( const Vector2 &uv ) const
	{
		return m_absorption;
	}

}; //Birefringence

typedef Birefringence<3> Birefringence3D;
typedef Birefringence<2> Birefringence2D;

#endif //_FLUORESCENT_H_
