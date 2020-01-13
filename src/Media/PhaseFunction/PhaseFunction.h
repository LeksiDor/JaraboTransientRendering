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

#ifndef _PHASE_FUNCTION_H_
#define _PHASE_FUNCTION_H_

#include "bunnykiller.h"

#include <cmath>

#include "Sampling/Sampling.h"
#include "Color/Spectrum.h"
#include "Color/PolarizedAttenuation.h"
#include "Color/FluorescentAttenuation.h"
#include "LinearAlgebra/Vector3.h"
#include "LinearAlgebra/Vector2.h"

namespace PhaseFunction
{
	//-------------------------------------------------------------------------------------------------
	// Base virtual phase function. 
    template <unsigned D>
	class PhaseFunction
	{
	protected:
		Real m_mean_time; //this should be changed for some more intuitive parameter
		Real m_inv_mean_time;
	protected:
		// Sample the time delay
		inline Real sample_time_delay(Real &pdf) const
		{	
			Real epsilon1;
            if (m_mean_time == 0.) {
                pdf = 1.;
                return 0.;
            } else {
                epsilon1 = Random::StdRNG.next_real();
                Real delta_time = -std::log(epsilon1)/m_inv_mean_time;
                pdf = exp(-m_inv_mean_time*delta_time)*m_inv_mean_time;
				return delta_time;
            }
		}

		inline Real pdf_time_delay(const Real delta_time) const
		{
			if (m_mean_time == 0.) {
				return 1.;
			}

			return exp(-m_inv_mean_time*delta_time)*m_inv_mean_time;
		}
	public:
		PhaseFunction(const Real mean_time): m_mean_time(mean_time), m_inv_mean_time(mean_time > 0 ? 1./mean_time : 0.) {}

		virtual ~PhaseFunction() {}

		virtual Real p(const VectorN<D>& wi, const VectorN<D> &wo) const = 0;

		virtual Real p(const VectorN<D>& wi, const VectorN<D> &wo, const Real delta_time) const
		{
			return pdf_time_delay(delta_time)*this->p(wi, wo);
		}

		virtual void operator()(const VectorN<D>& wi, const VectorN<D>& wo, Spectrum &R) const = 0;

		virtual void operator()(const VectorN<D>& wi, const VectorN<D>& wo,
			Spectrum &R, Real &delta_time, Real &pdf) const
		{
			delta_time = sample_time_delay(pdf);

			(*this)(wi, wo, R);
			R *= pdf;
		}

        virtual void sample_direction(const VectorN<D>& wi, VectorN<D> &wo, Spectrum &R, Real &pdf) const = 0;
		virtual void sample_direction(const VectorN<D>& wi, VectorN<D> &wo, Spectrum &R, Real &delta_time, Real &pdf) const
		{        
			this->sample_direction(wi, wo, R, pdf);
			
			Real pdf_time;
			delta_time = sample_time_delay(pdf_time);
			pdf *= pdf_time;

			R *= pdf_time;
		}

		/*
		* By default, the fluorescent attenuation is just an spectrum (no fluorescency)
		*/
		virtual void operator()(const VectorN<D>& wi, const VectorN<D>& wo, FluorescentAttenuation<D> &R) const
		{
			Spectrum r;
			(*this)(wi, wo, r);
			R.set(r);
		}

		virtual void operator()(const VectorN<D>& wi, const VectorN<D>& wo,
			FluorescentAttenuation<D> &R, Real &delta_time, Real &pdf) const
		{
			delta_time = sample_time_delay(pdf);

			(*this)(wi, wo, R);
			R *= pdf;
		}

		virtual void sample_direction(const VectorN<D>& wi, VectorN<D> &wo, FluorescentAttenuation<D> &R, Real &pdf) const
		{
			Spectrum r;
			sample_direction(wi, wo, r, pdf);
			R.set(r);
		}

		virtual void sample_direction(const VectorN<D>& wi, VectorN<D> &wo,
			FluorescentAttenuation<D> &R, Real &delta_time, Real &pdf) const
		{
			sample_direction(wi, wo, R, pdf);

			Real pdf_time;
			delta_time = sample_time_delay(pdf_time);
			pdf *= pdf_time;

			R *= pdf_time;
		}

		virtual void operator()(const VectorN<D>& wi, const VectorN<D>& wo, PolarizedAttenuation<D> &R) const
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
		virtual void operator()(const VectorN<D>& wi, const VectorN<D>& wo,
			PolarizedAttenuation<D> &R, Real &delta_time, Real &pdf) const
		{
			delta_time = sample_time_delay(pdf);

			(*this)(wi, wo, R);
			R *= pdf;
		};

		virtual void sample_direction(const VectorN<D>& wi, VectorN<D> &wo, PolarizedAttenuation<D> &R, Real &pdf) const
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

		virtual void sample_direction(const VectorN<D>& wi, VectorN<D> &wo,
			PolarizedAttenuation<D> &R, Real &delta_time, Real &pdf) const
		{
			sample_direction(wi, wo, R, pdf);

			Real pdf_time;
			delta_time = sample_time_delay(pdf_time);
			pdf *= pdf_time;

			R *= pdf_time;
		}

		virtual void operator()(const VectorN<D>& wi, const VectorN<D>& wo, FluorescentPolarizedAttenuation<D> &R) const
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

		virtual void operator()(const VectorN<D>& wi, const VectorN<D>& wo,
			FluorescentPolarizedAttenuation<D> &R, Real &delta_time, Real &pdf) const
		{
			delta_time = sample_time_delay(pdf);

			(*this)(wi, wo, R);
			R *= pdf;
		};

		virtual void sample_direction(const VectorN<D>& wi, VectorN<D> &wo, FluorescentPolarizedAttenuation<D> &R, Real &pdf) const
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

		virtual void sample_direction(const VectorN<D>& wi, VectorN<D> &wo,
			FluorescentPolarizedAttenuation<D> &R, Real &delta_time, Real &pdf) const
		{
			sample_direction(wi, wo, R, pdf);

			Real pdf_time;
			delta_time = sample_time_delay(pdf_time);
			pdf *= pdf_time;

			R *= pdf_time;
		}
	}; // class PhaseFunction
}; // namespace PhaseFunction

#endif //_PHASE_FUNCTION_H_
