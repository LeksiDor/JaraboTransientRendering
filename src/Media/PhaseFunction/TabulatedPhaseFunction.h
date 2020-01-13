/*
 * Copyright (C) 2017, Ibon Guillen (http://giga.cps.unizar.es/~ibon/)
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

#ifndef _TABULATED_PHASE_FUNCTION_H_
#define _TABULATED_PHASE_FUNCTION_H_

#include "bunnykiller.h"

#include <complex>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "Media/PhaseFunction/PhaseFunction.h"
#include "Sampling/Distributions/TabulatedDistribution.h"
#include "Sampling/Distributions/ConstantDistribution.h"
#include "Sampling/Distributions/LinearDistribution.h"
#include "Utils/RandomNumbers.h"

#define _CONST_DISTR
//#define _LINEAR_DISTR

namespace PhaseFunction
{
	//-------------------------------------------------------------------------------------------------
	// Tabulated phase function.
    template<unsigned D>
	class Tabulated : public PhaseFunction<D>
	{
	protected:
		typedef PhaseFunction<D> PFTR;

	protected:
		struct MieParams {
			Spectrum S11, S12, S33, S34;

			MieParams() :
				S11(0.), S12(0.), S33(0.), S34(0.)
			{}

			MieParams(const Spectrum _S11, const Spectrum _S12,
					  const Spectrum _S33, const Spectrum _S34) :
				S11(_S11), S12(_S12), S33(_S33), S34(_S34)
			{}
		};

		// Mie params
		std::vector<MieParams> m_S;

		Real m_delta, m_inv_delta;

		// Importance sampling
#if defined(_CONST_DISTR)
		ConstantDistribution m_D;
#elif defined(_LINEAR_DISTR)
		LinearDistribution m_D;
#endif
	private:
		inline MieParams interp_params(Real cos_theta) const {
			cos_theta = cos_theta + 1.;

			size_t idx = size_t(std::floor(cos_theta * m_inv_delta));

			const MieParams& S0 = m_S[idx];
			const MieParams& S1 = m_S[idx+1];

			Real cos0 = Real(idx)*m_delta;

			Real t1 = (cos_theta - cos0)*m_inv_delta;
			Real t0 = 1. - t1;

			return MieParams(
				t0*S0.S11 + t1*S1.S11,
				t0*S0.S12 + t1*S1.S12,
				t0*S0.S33 + t1*S1.S33,
				t0*S0.S34 + t1*S1.S34
			);
		}

		const MieParams sample_dist(const Vector3 &wi, Vector3 &wo, Real &pdf) const {
#if defined(_CONST_DISTR) || defined(_LINEAR_DISTR)
			Real epsilon1 = Random::StdRNG.next_real();
			Real epsilon2 = Random::StdRNG.next_real();

			Real phi = 2.*M_PI*epsilon1;

			Real cos_theta = m_D.sample(epsilon2, pdf);
			Real sin_theta = sqrt(1. - cos_theta*cos_theta);

			pdf *= (.5*M_1_PI);

			Vector3 dir(cos(phi)*sin_theta, cos_theta, sin(phi)*sin_theta);
			wo = dir.transform_matrix_to(Vector3(0., 1., 0.), wi);

#else
			Sampling::spherical(wo, pdf);
			wo = wo.transform_matrix_to(Vector3(0., 1., 0.), wi);

			Real cos_theta = dot(wi, wo);
#endif
			return interp_params(cos_theta);
		}
	public:
		Tabulated(std::ifstream& source, const Real mean_time) :
			PFTR(mean_time)
		{
			// Read tabulated pf
			std::string line;

			size_t N;
			std::getline(source, line);
			{
				std::istringstream in(line);
				in >> N;
			}
			m_S.reserve(N);

			std::vector<Real> avg;
			avg.reserve(N);

			std::getline(source, line);
			{
				std::istringstream in(line);
				for (size_t i = 0; i < N; i++) {
					float r, g, b;

					// Recover Muller matrix coefficients
					in >> r >> g >> b;
					Spectrum S11 = (.25*M_1_PI) * Spectrum(r, g, b);

					in >> r >> g >> b;
					Spectrum S12 = (.25*M_1_PI) * Spectrum(r, g, b);

					in >> r >> g >> b;
					Spectrum S33 = (.25*M_1_PI) * Spectrum(r, g, b);

					in >> r >> g >> b;
					Spectrum S34 = (.25*M_1_PI) * Spectrum(r, g, b);

					m_S.push_back(MieParams(
						S11, S12, S33, S34
					));

					// Store importance
					avg.push_back(S11.avg());
				}
			}
#if defined(_CONST_DISTR)
			m_D = ConstantDistribution(avg.data(), N, -1., 1.);
#elif defined(_LINEAR_DISTR)
			m_D = LinearDistribution(avg.data(), N, -1., 1.);
#endif

			m_delta = 2. / Real(N-1);
			m_inv_delta = 0.5 * Real(N-1);
		}


#if defined(_CONST_DISTR)
		Tabulated(Tabulated& t) :
			PFTR(t.m_mean_time), m_S(std::move(t.m_S)), m_delta(t.m_delta),  m_inv_delta(t.m_inv_delta),
			m_D(t.m_D)
		{};
#else
		Tabulated(Tabulated& t) :
			PFTR(t.m_mean_time), m_S(std::move(t.m_S)), m_delta(t.m_delta),  m_inv_delta(t.m_inv_delta)
		{};
#endif


		virtual ~Tabulated() {}
	public:
		void sample_direction(const VectorN<D> &wi, VectorN<D> &wo, Spectrum &R, Real &pdf) const;
		inline void sample_direction(const VectorN<D>& wi, VectorN<D> &wo,
			Spectrum &R, Real &delta_time, Real &pdf) const
		{
			delta_time = 0.; pdf = 1.;
			sample_direction(wi, wo, R, pdf);
		}

		void operator()(const VectorN<D>& wi, const VectorN<D>& wo, Spectrum &R) const;
		inline void operator()(const VectorN<D>& wi, const VectorN<D>& wo,
			Spectrum &R, Real &delta_time, Real &pdf) const
		{
			delta_time = 0.; pdf = 1.;
			(*this)(wi, wo, R);
		}

		Real p(const VectorN<D>& wi, const VectorN<D> &wo) const;
		inline Real p(const VectorN<D>& wi, const VectorN<D> &wo, const Real delta_time)const
		{
			return p(wi, wo);
		}

		void sample_direction(const VectorN<D> &wi, VectorN<D> &wo,
			FluorescentAttenuation<D> &R, Real &pdf) const;
		void sample_direction(const VectorN<D> &wi, VectorN<D> &wo,
			FluorescentAttenuation<D> &R, Real &delta_time, Real &pdf) const
		{
			delta_time = 0; pdf = 1;
			sample_direction(wi, wo, R, pdf);
		}

		void operator()(const VectorN<D>& wi, const VectorN<D>& wo,
			FluorescentAttenuation<D> &R) const;
		void operator()(const VectorN<D> &wi, const VectorN<D> &wo,
			FluorescentAttenuation<D> &R, Real &delta_time, Real &pdf) const
		{
			delta_time = 0; pdf = 1;
			(*this)(wi, wo, R);
		}

		void sample_direction(const VectorN<D> &wi, VectorN<D> &wo,
			PolarizedAttenuation<D> &R, Real &pdf) const;
		void sample_direction(const VectorN<D> &wi, VectorN<D> &wo,
			PolarizedAttenuation<D> &R, Real &delta_time, Real &pdf) const
		{
			delta_time = 0; pdf = 1;
			sample_direction(wi, wo, R, pdf);
		}

		void operator()(const VectorN<D>& wi, const VectorN<D>& wo,
			PolarizedAttenuation<D> &R) const;
		void operator()(const VectorN<D> &wi, const VectorN<D> &wo,
			PolarizedAttenuation<D> &R, Real &delta_time, Real &pdf) const
		{
			delta_time = 0; pdf = 1;
			(*this)(wi, wo, R);
		}

		void sample_direction(const VectorN<D> &wi, VectorN<D> &wo,
			FluorescentPolarizedAttenuation<D> &R, Real &pdf) const;
		void sample_direction(const VectorN<D> &wi, VectorN<D> &wo,
			FluorescentPolarizedAttenuation<D> &R, Real &delta_time, Real &pdf) const
		{
			delta_time = 0; pdf = 1;
			sample_direction(wi, wo, R, pdf);
		}

		void operator()(const VectorN<D>& wi, const VectorN<D>& wo,
			FluorescentPolarizedAttenuation<D> &R) const;
		void operator()(const VectorN<D> &wi, const VectorN<D> &wo,
			FluorescentPolarizedAttenuation<D> &R, Real &delta_time, Real &pdf) const
		{
			delta_time = 0; pdf = 1;
			(*this)(wi, wo, R);
		}
	}; // class Tabulated

    template<unsigned int D>
	void Tabulated<D>::operator()(const VectorN<D>& wi, const VectorN<D>& wo, Spectrum& R) const
	{
		throw std::runtime_error("Invalid dimensions: Tabulated pf scattering only accepts N=3\n");
	}

	template<>
	void Tabulated<3>::operator()(const Vector3& wi, const Vector3& wo, Spectrum& R) const
	{
		MieParams S = interp_params(dot(wi, wo));

		R = S.S11;
	}

	template<unsigned int D>
	void Tabulated<D>::sample_direction(const VectorN<D>& wi, VectorN<D>& wo, Spectrum& R, Real& pdf) const
	{
		throw std::runtime_error("Invalid dimensions: Tabulated pf scattering only accepts N=3\n");
	}

	template<>
	void Tabulated<3>::sample_direction(const Vector3& wi, Vector3& wo, Spectrum& R, Real& pdf) const
	{
		MieParams S = sample_dist(wi, wo, pdf);

		R = S.S11;
	}

	template<unsigned int D>
	void Tabulated<D>::operator()(const VectorN<D>& wi, const VectorN<D>& wo,
		FluorescentAttenuation<D>& R) const
	{
		throw std::runtime_error("Invalid dimensions: Tabulated pf scattering only accepts N=3\n");
	}

	template<>
	void Tabulated<3>::operator()(const Vector3& wi, const Vector3& wo,
		FluorescentAttenuation<3>& R) const
	{
		MieParams S = interp_params(dot(wi, wo));

		R = S.S11;
	}

	template<unsigned int D>
	void Tabulated<D>::sample_direction(const VectorN<D>& wi, VectorN<D>& wo,
		FluorescentAttenuation<D>& R, Real& pdf) const
	{
		throw std::runtime_error("Invalid dimensions: Tabulated pf scattering only accepts N=3\n");
	}

	template<>
	void Tabulated<3>::sample_direction(const Vector3& wi, Vector3& wo,
		FluorescentAttenuation<3>& R, Real& pdf) const
	{
		MieParams S = sample_dist(wi, wo, pdf);

		R = S.S11;
	}

	template<unsigned int D>
	void Tabulated<D>::operator()(const VectorN<D>& wi, const VectorN<D>& wo,
		PolarizedAttenuation<D> &R) const
	{
		throw("Invalid dimensions: Tabulated pf scattering only accepts N=3\n");
	}

	template<>
	void Tabulated<3>::operator()(const VectorN<3>& wi, const VectorN<3>& wo,
		PolarizedAttenuation<3> &R) const
	{
		MieParams S = interp_params(dot(wi, wo));

		R = PolarizedAttenuation<3>(
				S.S11, S.S12,      0,     0,
				S.S12, S.S11,      0,     0,
					0,     0,  S.S33, S.S34,
					0,     0, -S.S34, S.S33,
				PolarizationFrame<3>( wo, wi),
				PolarizationFrame<3>(-wi, wo));
	}

	template<unsigned int D>
	void Tabulated<D>::sample_direction(const VectorN<D>& wi, VectorN<D> &wo,
		PolarizedAttenuation<D> &R, Real &pdf) const
	{
		throw std::runtime_error("Invalid dimensions: Tabulated pf scattering only accepts N=3\n");
	}

	template<>
	void Tabulated<3>::sample_direction(const VectorN<3>& wi, VectorN<3> &wo,
		PolarizedAttenuation<3> &R, Real &pdf) const
	{
		MieParams S = sample_dist(wi, wo, pdf);

		R = PolarizedAttenuation<3>(
				S.S11, S.S12,      0,     0,
				S.S12, S.S11,      0,     0,
					0,     0,  S.S33, S.S34,
					0,     0, -S.S34, S.S33,
				PolarizationFrame<3>( wo, wi),
				PolarizationFrame<3>(-wi, wo));
	}

	template<unsigned int D>
	void Tabulated<D>::operator()(const VectorN<D>& wi, const VectorN<D>& wo,
		FluorescentPolarizedAttenuation<D> &R) const
	{
		throw std::runtime_error("Invalid dimensions: Tabulated pf scattering only accepts N=3\n");
	}

	template<>
	void Tabulated<3>::operator()(const VectorN<3>& wi, const VectorN<3>& wo,
		FluorescentPolarizedAttenuation<3> &R) const
	{
		MieParams S = interp_params(dot(wi, wo));

		R = FluorescentPolarizedAttenuation<3>(
				S.S11, S.S12,     0.,    0.,
				S.S12, S.S11,     0.,    0.,
				   0.,    0.,  S.S33, S.S34,
				   0.,    0., -S.S34, S.S33,
				PolarizationFrame<3>( wo, wi),
				PolarizationFrame<3>(-wi, wo));
	}

	template<unsigned int D>
	void Tabulated<D>::sample_direction(const VectorN<D>& wi, VectorN<D> &wo,
		FluorescentPolarizedAttenuation<D> &R, Real &pdf) const
	{
		throw std::runtime_error("Invalid dimensions: Tabulated pf scattering only accepts N=3\n");
	}

	template<>
	void Tabulated<3>::sample_direction(const VectorN<3>& wi, VectorN<3> &wo,
		FluorescentPolarizedAttenuation<3> &R, Real &pdf) const
	{
		MieParams S = sample_dist(wi, wo, pdf);

		R = FluorescentPolarizedAttenuation<3>(
				S.S11, S.S12,     0.,    0.,
				S.S12, S.S11,     0.,    0.,
				   0.,    0.,  S.S33, S.S34,
				   0.,    0., -S.S34, S.S33,
				PolarizationFrame<3>( wo, wi),
				PolarizationFrame<3>(-wi, wo));
	}

	template<>
	Real Tabulated<3>::p(const VectorN<3>& wi, const VectorN<3> &wo) const
	{
#if defined(_CONST_DISTR) || defined(_LINEAR_DISTR)
		return m_D.pdf(dot(wi, wo))*(.5*M_1_PI);
#else
		return (.25*M_1_PI);
#endif
	}
}; // namespace PhaseFunction

#endif // _TABULATED_PHASE_FUNCTION_H_
