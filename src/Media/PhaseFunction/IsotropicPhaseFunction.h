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

#ifndef _ISOTROPIC_PHASE_FUNCTION_H_
#define _ISOTROPIC_PHASE_FUNCTION_H_

#include "Media/PhaseFunction/PhaseFunction.h"

namespace PhaseFunction
{
	//-------------------------------------------------------------------------------------------------
	// Isotropic phase function. 
	template<unsigned D>
	class Isotropic: public PhaseFunction<D>
	{
	protected:
		using PFTR = PhaseFunction<D>;
	public:
        Isotropic(const Real mean_time) : PFTR(mean_time)
        {}
	
		inline Real p(const VectorN<D>& wi, const VectorN<D> &wo) const
		{
			if (D == 2)
				return (.5*M_1_PI);
			else
				return (.25*M_1_PI);
		}

		inline Real p(const VectorN<D>& wi, const VectorN<D> &wo, const Real delta_time)const
		{
			return PFTR::p(wi, wo, delta_time);
		}

		inline void sample_direction(const VectorN<D>& wi, VectorN<D> &wo, Spectrum &R, Real &pdf)const
		{
			Sampling.direction_uniformly(wo, pdf);
			wo = wo.transform_matrix_to(Vector3(0., 1., 0.), wi);

			(*this)(wi, wo, R);
		}

		inline void sample_direction(const VectorN<D>& wi, VectorN<D> &wo,
			Spectrum &R, Real &delta_time, Real &pdf) const
		{
			PFTR::sample_direction(wi, wo, R, delta_time, pdf);
		}

		inline void operator()(const VectorN<D>& wi, const VectorN<D>& wo, Spectrum &R) const
		{ 
			if (D == 2)
				R = Spectrum(.5*M_1_PI);
			else
				R = Spectrum(.25*M_1_PI);
		}

		inline void operator()(const VectorN<D>& wi, const VectorN<D>& wo,
			Spectrum &R, Real &delta_time, Real &pdf) const
		{
			PFTR::operator()(wi, wo, R, delta_time, pdf);
		}

		inline void sample_direction(const VectorN<D>& wi, VectorN<D> &wo, FluorescentAttenuation<D> &R, Real &pdf)const
		{
			PFTR::sample_direction(wi, wo, R, pdf);
		}

		inline void sample_direction(const VectorN<D>& wi, VectorN<D> &wo,
			FluorescentAttenuation<D> &R, Real &delta_time, Real &pdf) const
		{
			PFTR::sample_direction(wi, wo, R, delta_time, pdf);
		}

		inline void operator()(const VectorN<D>& wi, const VectorN<D>& wo, FluorescentAttenuation<D> &R) const
		{
			PFTR::operator()(wi, wo, R);
		}

		inline void operator()(const VectorN<D>& wi, const VectorN<D>& wo,
			FluorescentAttenuation<D> &R, Real &delta_time, Real &pdf) const
		{
			PFTR::operator()(wi, wo, R, delta_time, pdf);
		}

		inline void sample_direction(const VectorN<D>& wi, VectorN<D> &wo, PolarizedAttenuation<D> &R, Real &pdf)const
		{
			PFTR::sample_direction(wi, wo, R, pdf);
		}

		inline void sample_direction(const VectorN<D>& wi, VectorN<D> &wo,
			PolarizedAttenuation<D> &R, Real &delta_time, Real &pdf) const
		{
			PFTR::sample_direction(wi, wo, R, delta_time, pdf);
		}

		inline void operator()(const VectorN<D>& wi, const VectorN<D>& wo, PolarizedAttenuation<D> &R) const
		{
			PFTR::operator()(wi, wo, R);
		}

		inline void operator()(const VectorN<D>& wi, const VectorN<D>& wo,
			PolarizedAttenuation<D> &R, Real &delta_time, Real &pdf) const
		{
			PFTR::operator()(wi, wo, R, delta_time, pdf);
		}

		inline void sample_direction(const VectorN<D>& wi, VectorN<D> &wo, FluorescentPolarizedAttenuation<D> &R, Real &pdf)const
		{
			PFTR::sample_direction(wi, wo, R, pdf);
		}

		inline void sample_direction(const VectorN<D>& wi, VectorN<D> &wo,
			FluorescentPolarizedAttenuation<D> &R, Real &delta_time, Real &pdf) const
		{
			PFTR::sample_direction(wi, wo, R, delta_time, pdf);
		}

		inline void operator()(const VectorN<D>& wi, const VectorN<D>& wo, FluorescentPolarizedAttenuation<D> &R) const
		{
			PFTR::operator()(wi, wo, R);
		}

		inline void operator()(const VectorN<D>& wi, const VectorN<D>& wo,
			FluorescentPolarizedAttenuation<D> &R, Real &delta_time, Real &pdf) const
		{
			PFTR::operator()(wi, wo, R, delta_time, pdf);
		}
	}; // class Isotropic
}; // namespace PhaseFunction

#endif //_ISOTROPIC_PHASE_FUNCTION_H_
