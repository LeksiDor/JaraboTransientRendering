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

#ifndef _RAYLEIGH_PHASE_FUNCTION_H_
#define _RAYLEIGH_PHASE_FUNCTION_H_

#include "bunnykiller.h"

#include "Utils/Utils.h"
#include "Media/PhaseFunction/PhaseFunction.h"

namespace PhaseFunction
{
	/**
	 *	Rayleigh phase function.
	 *	Computes the rayleigh scattering phase function, based on a simplified approach.
	 *
	 *	Unpolarized sampling and evaluation following Frisvad technique [1].
	 *	Complete Mueller matrices following Bickel formulae [2]
	 *
	 *	[1] J. R. Frisvad 2011
	 *	    Importance sampling the Rayleigh phase function
	 *		JOSA: Optics, Image Science, and Vision, 28(12)
	 *	[2] W. S. Bickel 1985
	 *	    The Mueller scattering matrix elements for Rayleigh spheres
	 *		Applications of Circularly Polarized Radiation Using Synchrotron and Ordinary Sources.
	 *		Springer US, 1985. 69-76.
	 */
	template<unsigned D>
	class Rayleigh : public PhaseFunction<D>
	{
	protected:
		using PFTR = PhaseFunction<D>;
		/* Rayleigh function is undefined for D != 3 */
		static_assert(D == 3, "Invalid dimensions: Rayleigh scattering only accepts D=3");
	private:
		void sample_importance(const Vector3 &wi, Vector3 &wo, Real &pdf) const
		{
			Real epsilon1 = Random::StdRNG.next_real();
			Real epsilon2 = Random::StdRNG.next_real();

			Real phi = 2. * M_PI * epsilon1;

			Real two_eps_minus = 2. * epsilon2 - 1.;
			Real u = -std::cbrt(
					2. * two_eps_minus + std::sqrt(4. * two_eps_minus * two_eps_minus + 1.));

			Real cos_theta = Utils::clamp(u - 1. / u, -1., 1.);
			Real sin_theta = std::sqrt(1. - cos_theta * cos_theta);

			pdf = 3. * (1. + cos_theta * cos_theta) / (16. * M_PI);

			Vector3 dir(std::cos(phi) * sin_theta, cos_theta, std::sin(phi) * sin_theta);
			wo = dir.transform_matrix_to(Vector3(0., 1., 0.), wi);
		}
	public:
		Rayleigh() :
				PFTR(0.)
		{
		}

		void sample_direction(const VectorN<D> &wi, VectorN<D> &wo, Spectrum &R, Real &pdf) const;
		void sample_direction(const VectorN<D>& wi, VectorN<D> &wo, Spectrum &R, Real &delta_time,
			Real &pdf) const
		{
			delta_time = 0.;
			pdf = 1.;
			sample_direction(wi, wo, R, pdf);
		}

		void operator()(const VectorN<D>& wi, const VectorN<D>& wo, Spectrum &R) const;
		void operator()(const VectorN<D>& wi, const VectorN<D>& wo, Spectrum &R, Real &delta_time,
			Real &pdf) const
		{
			delta_time = 0.;
			pdf = 1.;
			(*this)(wi, wo, R);
		}

		Real p(const VectorN<D>& wi, const VectorN<D> &wo) const;
		Real p(const VectorN<D>& wi, const VectorN<D> &wo, const Real delta_time) const
		{
			return p(wi, wo);
		}

		void sample_direction(const VectorN<D> &wi, VectorN<D> &wo, FluorescentAttenuation<D> &R,
			Real &pdf) const;
		void sample_direction(const VectorN<D> &wi, VectorN<D> &wo, FluorescentAttenuation<D> &R,
			Real &delta_time, Real &pdf) const
		{
			delta_time = 0.;
			pdf = 1.;
			sample_direction(wi, wo, R, pdf);
		}

		void operator()(const VectorN<D>& wi, const VectorN<D>& wo,
			FluorescentAttenuation<D> &R) const;
		void operator()(const VectorN<D> &wi, const VectorN<D> &wo, FluorescentAttenuation<D> &R,
			Real &delta_time, Real &pdf) const
		{
			delta_time = 0.;
			pdf = 1.;
			(*this)(wi, wo, R);
		}

		void sample_direction(const VectorN<D> &wi, VectorN<D> &wo, PolarizedAttenuation<D> &R,
			Real &pdf) const;
		void sample_direction(const VectorN<D> &wi, VectorN<D> &wo, PolarizedAttenuation<D> &R,
			Real &delta_time, Real &pdf) const
		{
			delta_time = 0.;
			pdf = 1.;
			sample_direction(wi, wo, R, pdf);
		}

		void operator()(const VectorN<D>& wi, const VectorN<D>& wo,
			PolarizedAttenuation<D> &R) const;
		void operator()(const VectorN<D> &wi, const VectorN<D> &wo, PolarizedAttenuation<D> &R,
			Real &delta_time, Real &pdf) const
		{
			delta_time = 0.;
			pdf = 1.;
			(*this)(wi, wo, R);
		}

		void sample_direction(const VectorN<D> &wi, VectorN<D> &wo,
			FluorescentPolarizedAttenuation<D> &R, Real &pdf) const;
		void sample_direction(const VectorN<D> &wi, VectorN<D> &wo,
			FluorescentPolarizedAttenuation<D> &R, Real &delta_time, Real &pdf) const
		{
			delta_time = 0.;
			pdf = 1.;
			sample_direction(wi, wo, R, pdf);
		}

		void operator()(const VectorN<D>& wi, const VectorN<D>& wo,
			FluorescentPolarizedAttenuation<D> &R) const;
		void operator()(const VectorN<D> &wi, const VectorN<D> &wo,
			FluorescentPolarizedAttenuation<D> &R, Real &delta_time, Real &pdf) const
		{
			delta_time = 0.;
			pdf = 1.;
			(*this)(wi, wo, R);
		}
	}; /* class Rayleigh */

	template<>
	void Rayleigh<3>::operator()(const Vector3& wi, const Vector3& wo, Spectrum& R) const
	{
		/* Simplified version based on [4]:
		 *
		 * [4] http://www.oceanopticsbook.info/view/radiative_transfer_theory/level_2/the_vector_radiative_transfer_equation
		 */
		Real cos_theta = dot(wi, wo);
		R = 3. * (1. + cos_theta * cos_theta) / (16. * M_PI);
	}

	template<>
	void Rayleigh<3>::sample_direction(const Vector3& wi, Vector3& wo, Spectrum& R, Real& pdf) const
	{
		sample_importance(wi, wo, pdf);
		(*this)(wi, wo, R);
	}

	template<>
	void Rayleigh<3>::operator()(const Vector3& wi, const Vector3& wo,
		FluorescentAttenuation<3> &R) const
	{
		/* Simplified version based on [4]:
		 *
		 * [4] http://www.oceanopticsbook.info/view/radiative_transfer_theory/level_2/the_vector_radiative_transfer_equation
		 */
		Real cos_theta = dot(wi, wo);
		R = 3. * (1. + cos_theta * cos_theta) / (16. * M_PI);
	}

	template<>
	void Rayleigh<3>::sample_direction(const Vector3& wi, Vector3 &wo, FluorescentAttenuation<3> &R,
		Real &pdf) const
	{
		sample_importance(wi, wo, pdf);
		(*this)(wi, wo, R);
	}

	template<>
	Real Rayleigh<3>::p(const Vector3& wi, const Vector3& wo) const
	{
		Real cos_theta = dot(wi, wo);
		Real pdf = 3. * (1. + cos_theta * cos_theta) / (16. * M_PI);

		return pdf;
	}

	template<>
	void Rayleigh<3>::operator()(const Vector3& wi, const Vector3& wo,
		PolarizedAttenuation<3> &R) const
	{
		Real cos_theta = Utils::clamp(dot(wi, wo), -1., 1.);

		Real cos_theta2 = cos_theta * cos_theta;
		Real sin_theta2 = 1. - cos_theta2;

		// Norm factor
		Real r = 3. / (16. * M_PI);

		Spectrum S11(r * (1. + cos_theta2));
		Spectrum S12(-r * (sin_theta2));
		Spectrum S33(r * (2. * cos_theta));

		R = PolarizedAttenuation<3>(
						S11, S12,   0,   0,
						S12, S11,   0,   0,
						  0,   0, S33,   0,
						  0,   0,   0, S33,
						PolarizationFrame<3>( wo, wi),
						PolarizationFrame<3>(-wi, wo));
	}

	template<>
	void Rayleigh<3>::sample_direction(const Vector3& wi, Vector3& wo, PolarizedAttenuation<3> &R,
		Real &pdf) const
	{
		sample_importance(wi, wo, pdf);
		(*this)(wi, wo, R);
	}

	template<>
	void Rayleigh<3>::operator()(const Vector3& wi, const Vector3& wo,
		FluorescentPolarizedAttenuation<3> &R) const
	{
		Real cos_theta = Utils::clamp(dot(wi, wo), -1., 1.);

		Real cos_theta2 = cos_theta * cos_theta;
		Real sin_theta2 = 1. - cos_theta2;

		/* Norm factor */
		Real r = 3. / (16. * M_PI);

		Spectrum S11(r * (1. + cos_theta2));
		Spectrum S12(-r * (sin_theta2));
		Spectrum S33(r * (2. * cos_theta));

		R = FluorescentPolarizedAttenuation<3>(
						S11, S12,   0,   0,
						S12, S11,   0,   0,
						  0,   0, S33,   0,
						  0,   0,   0, S33,
						PolarizationFrame<3>( wo, wi),
						PolarizationFrame<3>(-wi, wo));
	}

	template<>
	void Rayleigh<3>::sample_direction(const Vector3& wi, Vector3& wo,
		FluorescentPolarizedAttenuation<3> &R, Real &pdf) const
	{
		sample_importance(wi, wo, pdf);
		(*this)(wi, wo, R);
	}
}; /* namespace PhaseFunction */

#endif /* _RAYLEIGH_PHASE_FUNCTION_H_ */
