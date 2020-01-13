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

#ifndef _WOODCOCK_TRACKING_H_
#define _WOODCOCK_TRACKING_H_

#include <cmath>
#include "RayTracing/Ray.h"
#include "Utils/RandomNumbers.h"

namespace WoodcockTracking
{
	// Implements Woodcock tracking [1] to efficiently and unbiasedly sampling 
	// the mean free path in heterogeneous media, as described in [2].
	//
	// [1]	E. R. Woodcock, T. Murphy, P. J. Hemmings, T. C. Longworth. 1965. 
	//		Techniques used in the GEM code for Monte Carlo neutronics calculations
	//		in reactors and other systems of complex geometry. ANL-7050, Argonne Na-
	//		tional Laboratory
	// [2]	M. Raab, D. Seibert, A. Keller. 2006. Unbiased global illumination with 
	//		participating media. In Monte Carlo and Quasi-Monte Carlo Methods
	//		http://www.uni-ulm.de/fileadmin/website_uni_ulm/iui.inst.100/institut/Papers/ugiwpm.pdf
	template<int D, class M>
	Real sample_mean_free_path( const Ray<D> &r, const M &medium, Real &pdf)
	{
		Real sigma_t = medium.get_max_extinction().avg();
		Real total_sigma_t = sigma_t;

		Real t = -log(RNG::Random::get_real())/sigma_t;
		
		while( medium.get_extinction(r+t).avg()/sigma_t < RNG::Random::get_real()) {
			t -= log(RNG::Random::get_real())/sigma_t; total_sigma_t += sigma_t;
		}

		pdf = medium.get_extinction(r+t).avg()*exp(-total_sigma_t*t);
		return t;
	}

	template<int D, class M>
	Real p(const Ray<D> &r, const M &medium, const Real t )
	{
		Real sigma_t = medium.get_max_extinction().avg();
		Real int_sigma_t = static_cast<int>(t*sigma_t);

		return exp(-(int_sigma_t*sigma_t) - sigma_t*(t-int_sigma_t));
	}

	template<int D, class M>
	Spectrum compute_transmitance(const Ray<D> &r, const M &medium, const Real t )
	{
		if( t == std::numeric_limits<Real>::infinity() )
			return Spectrum(0.);

		Real sigma_t = medium.get_max_extinction().get_max();
		Spectrum int_sigma_t;

		Real step = t/sigma_t;

		for( Real s = 0; s<t; s=step ) {
			int_sigma_t += medium.get_extinction(r+s)/sigma_t;
		}

		return exp(-int_sigma_t*sigma_t);
	}
};


#endif
