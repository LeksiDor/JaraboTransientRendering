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

#ifndef _HENYEY_GREENSTEIN_PHASE_FUNCTION_H_
#define _HENYEY_GREENSTEIN_PHASE_FUNCTION_H_

#include "bunnykiller.h"

#include <cmath>

#include "Media/PhaseFunction/PhaseFunction.h"

namespace PhaseFunction
{
	//-------------------------------------------------------------------------------------------------
	// Henyey-Greenstein phase function. 
    template<unsigned D>
	class HenyeyGreenstein : public PhaseFunction<D>
	{
	protected:
		typedef PhaseFunction<D> PFTR;
	protected:
		Spectrum m_g;
		Real m_gY, m_inv_2gY;
	public:
		HenyeyGreenstein(const Real _g, const Real mean_time) : 
	        PFTR(mean_time), m_g(_g), m_gY(Color::rgb2xyz(Vector3(_g))[1]), m_inv_2gY(0.5/m_gY)
        {}

		HenyeyGreenstein(const Spectrum _g, const Real mean_time) :
			PFTR(mean_time), m_g(_g), m_gY(Color::rgb2xyz(Vector3(_g[0], _g[1], _g[2]))[1]), m_inv_2gY(0.5/m_gY)
		{}

		//HenyeyGreenstein(const PhaseFunction<D> &pf): g(0)
		//{ /*TODO: Integrate phase function numerically on the sphere to obtain g*/}

		void operator()(const VectorN<D>& wi, const VectorN<D>& wo, Spectrum &R) const;
		inline void operator()(const VectorN<D>& wi, const VectorN<D>& wo,
			Spectrum &R, Real &delta_time, Real &pdf) const
		{
			PFTR::operator()(wi, wo, R, delta_time, pdf);
		}

		Real p(const VectorN<D>& wi, const VectorN<D> &wo) const;
		inline Real p(const VectorN<D>& wi, const VectorN<D> &wo, const Real delta_time) const
		{
			return p(wi, wo);
		}

		void sample_direction(const VectorN<D> &wi, VectorN<D> &wo, Spectrum &R, Real &pdf) const;
		inline void sample_direction(const VectorN<D>& wi, VectorN<D> &wo,
			Spectrum &R, Real &delta_time, Real &pdf) const
		{
			PFTR::sample_direction(wi, wo, R, delta_time, pdf);
		}

		inline void sample_direction(const VectorN<D>& wi, VectorN<D> &wo,
			FluorescentAttenuation<D> &R, Real &pdf)const
		{
			Spectrum r;
			sample_direction(wi, wo, r, pdf);

			R.set(r);
		}

		inline void sample_direction(const VectorN<D>& wi, VectorN<D> &wo,
			FluorescentAttenuation<D>& R, Real& delta_time, Real& pdf) const
		{
			PFTR::sample_direction(wi, wo, R, delta_time, pdf);
		}

		inline void operator()(const VectorN<D>& wi, const VectorN<D>& wo,
			FluorescentAttenuation<D>& R) const
		{
			Spectrum r;
			(*this)(wi, wo, r);

			R.set(r);
		}

		inline void operator()(const VectorN<D>& wi, const VectorN<D>& wo,
			FluorescentAttenuation<D>& R, Real& delta_time, Real& pdf) const
		{
			PFTR::operator()(wi, wo, R, delta_time, pdf);
		}

		inline void sample_direction(const VectorN<D>& wi, VectorN<D> &wo, PolarizedAttenuation<D> &R, Real &pdf)const
		{
			Spectrum r;
			sample_direction(wi, wo, r, pdf);
			
			R.set(r, 0, 0, 0,
				  0, r, 0, 0,
				  0, 0, r, 0,
				  0, 0, 0, r,
				PolarizationFrame<D>(wo, wi),
				PolarizationFrame<D>(-wi, wo));
		}
		
		inline void sample_direction(const VectorN<D>& wi, VectorN<D>& wo,
			PolarizedAttenuation<D>& R, Real& delta_time, Real& pdf) const
		{
			PFTR::sample_direction(wi, wo, R, delta_time, pdf);
		}

		inline void operator()(const VectorN<D>& wi, const VectorN<D>& wo, PolarizedAttenuation<D>& R) const
		{
			Spectrum r;
			(*this)(wi, wo, r);

			R.set(r, 0, 0, 0,
				  0, r, 0, 0,
				  0, 0, r, 0,
				  0, 0, 0, r,
				PolarizationFrame<D>(wo, wi),
				PolarizationFrame<D>(-wi, wo));
		}

		inline void operator()(const VectorN<D>& wi, const VectorN<D>& wo,
			PolarizedAttenuation<D> &R, Real &delta_time, Real &pdf) const
		{
			PFTR::operator()(wi, wo, R, delta_time, pdf);
		}

		inline void sample_direction(const VectorN<D>& wi, VectorN<D>& wo, FluorescentPolarizedAttenuation<D>& R, Real& pdf)const
		{
			Spectrum r;
			sample_direction(wi, wo, r, pdf);

			R.set(r, 0., 0., 0.,
				  0., r, 0., 0.,
				  0., 0., r, 0.,
				  0., 0., 0., r,
				PolarizationFrame<D>(wo, wi),
				PolarizationFrame<D>(-wi, wo));
		}

		inline void sample_direction(const VectorN<D>& wi, VectorN<D>& wo,
			FluorescentPolarizedAttenuation<D> &R, Real &delta_time, Real &pdf) const
		{
			PFTR::sample_direction(wi, wo, R, delta_time, pdf);
		}

		inline void operator()(const VectorN<D>& wi, const VectorN<D>& wo, FluorescentPolarizedAttenuation<D>& R) const
		{
			Spectrum r;
			(*this)(wi, wo, r);

			R.set(r, 0., 0., 0.,
				  0., r, 0., 0.,
				  0., 0., r, 0.,
				  0., 0., 0., r,
				PolarizationFrame<D>(wo, wi),
				PolarizationFrame<D>(-wi, wo));
		}

		inline void operator()(const VectorN<D>& wi, const VectorN<D>& wo,
			FluorescentPolarizedAttenuation<D> &R, Real &delta_time, Real &pdf) const
		{
			PFTR::operator()(wi, wo, R, delta_time, pdf);
		}

		void set_directionality(const Real _g)
		{
			m_g = Spectrum(_g);
			m_gY = Color::rgb2xyz(Vector3(_g))[1];
			m_inv_2gY = 0.5/m_gY;
		}

		void set_directionality(const Spectrum _g)
		{
			m_g = _g; 
			m_gY = Color::rgb2xyz(Vector3(_g[0], _g[1], _g[2]))[1];
			m_inv_2gY = 0.5/m_gY;
		}
	}; // class HenyeyGreenstein

	template<>
	void HenyeyGreenstein<2>::sample_direction(const Vector2& wi, Vector2& wo, Spectrum& R, Real& pdf) const
	{
		Real epsilon1 = Random::StdRNG.next_real();

		Real cos_theta = std::cos(std::min(2.*atan((1. - m_gY)*std::tan(M_PI_2*(1. - 2.*epsilon1))/(1. + m_gY)), 1.));
		Real sin_theta = ((Random::StdRNG.next_int() % 2) == 1 ? 1. : -1.)*std::sqrt(std::max(0., 1. - cos_theta*cos_theta));
		
		Vector2 dir(sin_theta, cos_theta);
		wo = dir.transform_matrix_to(Vector2(0., 1.), wi);
		
		R = Spectrum(
			(1. - m_g[0]*m_g[0])/(2.*M_PI*(1. + m_g[0]*m_g[0] - 2.*m_g[0]*cos_theta)),
			(1. - m_g[1]*m_g[1])/(2.*M_PI*(1. + m_g[1]*m_g[1] - 2.*m_g[1]*cos_theta)),
			(1. - m_g[2]*m_g[2])/(2.*M_PI*(1. + m_g[2]*m_g[2] - 2.*m_g[2]*cos_theta))
		);

		pdf = (1. - m_gY*m_gY)/(2.*M_PI*(1. + m_gY*m_gY - 2.*m_gY*cos_theta));
	}

	template<>
	void HenyeyGreenstein<3>::sample_direction(const Vector3& wi, Vector3& wo, Spectrum& R, Real& pdf) const
	{
		Real epsilon1 = Random::StdRNG.next_real();
		Real epsilon2 = Random::StdRNG.next_real();
		
		Real phi = 2.*M_PI*epsilon1;

		Real aux = (1. - m_gY*m_gY)/(1. - m_gY + 2.*m_gY*epsilon2);

		Real cos_theta = Utils::clamp(m_inv_2gY*(1. + m_gY*m_gY - aux*aux), -1., 1.);
		Real sin_theta = std::sqrt(1. - cos_theta*cos_theta);
		
		Vector3 dir(std::cos(phi)*sin_theta, cos_theta, std::sin(phi)*sin_theta);
		wo = dir.transform_matrix_to(Vector3(0., 1., 0.), wi);

		pdf = (1. - m_gY*m_gY)/(4.*M_PI*pow(1. + m_gY*m_gY - 2.*m_gY*cos_theta, Real(1.5)));

		R = Spectrum(
			(1. - m_g[0]*m_g[0])/(4.*M_PI*pow(1. + m_g[0]*m_g[0] - 2.*m_g[0]*cos_theta, Real(1.5))),
			(1. - m_g[1]*m_g[1])/(4.*M_PI*pow(1. + m_g[1]*m_g[1] - 2.*m_g[1]*cos_theta, Real(1.5))),
			(1. - m_g[2]*m_g[2])/(4.*M_PI*pow(1. + m_g[2]*m_g[2] - 2.*m_g[2]*cos_theta, Real(1.5)))
		);
	}

	template<>
	void HenyeyGreenstein<2>::operator()(const Vector2& wi, const Vector2& wo, Spectrum &R) const
	{
		Real cos_theta = dot(wi, wo);
		R = Spectrum(
			(1. - m_g[0]*m_g[0])/(2.*M_PI*(1. + m_g[0]*m_g[0] - 2.*m_g[0]*cos_theta)),
			(1. - m_g[1]*m_g[1])/(2.*M_PI*(1. + m_g[1]*m_g[1] - 2.*m_g[1]*cos_theta)),
			(1. - m_g[2]*m_g[2])/(2.*M_PI*(1. + m_g[2]*m_g[2] - 2.*m_g[2]*cos_theta))
		);
	}

	template<>
	void HenyeyGreenstein<3>::operator()(const Vector3& wi, const Vector3& wo, Spectrum &R) const
	{
		Real cos_theta = dot(wi, wo);
		R = Spectrum(
			(1. - m_g[0]*m_g[0])/(4.*M_PI*pow(1. + m_g[0]*m_g[0] - 2.*m_g[0]*cos_theta, Real(1.5))),
			(1. - m_g[1]*m_g[1])/(4.*M_PI*pow(1. + m_g[1]*m_g[1] - 2.*m_g[1]*cos_theta, Real(1.5))),
			(1. - m_g[2]*m_g[2])/(4.*M_PI*pow(1. + m_g[2]*m_g[2] - 2.*m_g[2]*cos_theta, Real(1.5)))
		);
	}
    
	template<>
	Real HenyeyGreenstein<2>::p(const Vector2& wi, const Vector2& wo) const
	{
		Real cos_theta = dot(wi, wo);
		return (1. - m_gY*m_gY)/(2.*M_PI*(1. + m_gY*m_gY - 2.*m_gY*cos_theta));
	}

	template<>
	Real HenyeyGreenstein<3>::p(const Vector3& wi, const Vector3& wo) const
	{
		Real cos_theta = dot(wi, wo);
		return (1. - m_gY*m_gY)/(4.*M_PI*pow(1. + m_gY*m_gY - 2.*m_gY*cos_theta, Real(1.5)));
	}
}; // namespace PhaseFunction

#endif // _HENYEY_GREENSTEIN_PHASE_FUNCTION_H_
