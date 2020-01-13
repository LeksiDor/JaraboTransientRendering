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

#ifndef _SAMPLING_TECHNIQUES_H_
#define _SAMPLING_TECHNIQUES_H_

#include "bunnykiller.h"

#include <cmath>

#include "Utils/Utils.h"
#include "LinearAlgebra/VectorN.h"

const Real LIGHT_SPEED = 299792458.;

namespace Probability
{
	template<unsigned D>
	class Equiangular
	{
		Real minimum_distance_r, distance_ray_light;
		Real u_max, u_min;
	public:
		Real ru(Real u) const
		{
			return tan(u)*distance_ray_light + minimum_distance_r;
		}

		Real ur(Real r) const
		{
			return atan((r - minimum_distance_r)/distance_ray_light);
		}

		Real dudr(Real r) const
		{
			return distance_ray_light/
					  (r*r - 2.*minimum_distance_r*r + minimum_distance_r*minimum_distance_r + distance_ray_light*distance_ray_light);
		}

		Equiangular(const VectorN<D>& w, const VectorN<D>& l) :
			minimum_distance_r(dot(normalize(w), l)),
			distance_ray_light(VectorN<D>(l - normalize(w)*minimum_distance_r).length()),
			u_max(M_PI/2.), u_min(atan(-minimum_distance_r/distance_ray_light))
		{}

		Real pdf(Real r) const
		{
			return dudr(r)/(u_max - u_min);
		}

		Real cdf(Real r) const
		{
			return (ur(r) - u_min)/(u_max - u_min);
		}

		Real cdf_inv(Real e) const
		{
			return ru(e*(u_max - u_min) + u_min);
		}

		Real dominion_min() const
		{
			return Real(0);
		}

		Real dominion_max() const
		{
			return std::numeric_limits<Real>::max();
		}
	};
	
	template<unsigned D>
	class TimeLineToPoint
	{
		Real t_a, t_b;
		VectorN<D> w, l;	
		Real ec, ec2;

	public:
		Real t_max() const
		{
			return t_b;
		}

		Real t_min() const
		{
			return t_a;
		}

		Real rt(Real t) const
		{
			return (t*t - ec2*dot(l,l)) / (2.*ec*t - 2.*ec2*dot(l,w));
		}

		Real tr(Real r) const
		{
			return ec*(r + std::sqrt(r*r - 2.*r*dot(l,w) + dot(l,l)));
		}

		Real dtdr(Real r) const
		{
			return ec*(1. + (r - dot(l,w)) / std::sqrt(r*r - 2.*r*dot(l,w) + dot(l,l)));
		}

		//Sets t_a and t_b to possible outcomes given the geometry.
		void set_temporal_boundaries()
		{
			t_a = std::max(t_a, (l.length())*ec);
		}

		TimeLineToPoint(Real _t_a, Real _t_b, const VectorN<D> &_w, const VectorN<D> &_l, Real _eta = 1.0, Real _c = LIGHT_SPEED)
				: t_a(_t_a), t_b(_t_b), w(_w), l(_l), ec(_eta/_c), ec2(ec*ec)
		{
			set_temporal_boundaries();
		}

		TimeLineToPoint(Real _t_a, Real _t_b, const VectorN<D> &_w, const VectorN<D> &_xi, const VectorN<D> &_xo, 
							Real _eta = 1.0, Real _c = LIGHT_SPEED)
				: t_a(_t_a), t_b(_t_b), w(_w), l(_xi-_xo), ec(_eta/_c), ec2(ec*ec)
		{
			set_temporal_boundaries();
		}

		Real pdf(Real r) const
		{
			return dtdr(r)/(t_b - t_a);
		}

		Real cdf(Real r) const
		{
			return (tr(r) - t_a)/(t_b - t_a);
		}

		Real cdf_inv(Real e) const
		{
			return rt(e*(t_b - t_a) + t_a);
		}

		Real dominion_min() const
		{
			return rt(t_min());
		}

		Real dominion_max() const
		{
			return rt(t_max());
		}
	};
	
	template<unsigned D>
	class TimeAngular 
	{
		Real t_a, t_b;
		Real r, modl;
		Real ec;
	public:
		Real t_max() const
		{
			return t_b;
		}

		Real t_min() const
		{
			return t_a;
		}

		Real theta_t(Real t) const
		{
			// Watch out at the boundaries of this stuff...
			Real cos_t = (modl*modl + 2.*t*r/ec - t*t/(ec*ec)) / (2.*r*modl);
			return acos(Utils::clamp(cos_t, -1., 1.));
		}


		Real dominion_min() const
		{
			return theta_t(t_min());
		}

		Real dominion_max() const
		{
			return theta_t(t_max());
		}

		Real t_theta(Real theta) const
		{
			return ec*(r + sqrt(r*r + modl*modl - 2.0*r*modl*cos(theta)));
	//		return ec*(r*r + r + modl*modl - 2.0*r*modl*cos(theta));
		}

		Real dtdtheta(Real theta) const
		{
			return (ec*modl*r*sin(theta))/sqrt(r*r + modl*modl - 2.0*r*modl*cos(theta));
	//		return ec*(2.0*r*modl*sin(theta));
		}

		// Sets t_a and t_b to possible outcomes given the geometry.
		void set_temporal_boundaries()
		{
			if (r > modl) t_a = std::max(t_a,(static_cast<Real>(2.0)*r - modl)*ec);
			else          t_a = std::max(t_a, modl*ec);
			t_b = std::min(t_b,(static_cast<Real>(2.0)*r + modl)*ec);
		}

		TimeAngular( Real _t_a, Real _t_b, Real _r, Real _modl, Real _eta = 1.0, Real _c = LIGHT_SPEED) :
			t_a(_t_a), t_b(_t_b), r(_r), modl(_modl), ec(_eta/_c)
		{
			set_temporal_boundaries();
		}

		TimeAngular( Real _t_a, Real _t_b, Real _r, const VectorN<D>& _l, Real _eta = 1.0, Real _c = LIGHT_SPEED) :
			t_a(_t_a), t_b(_t_b), r(_r), modl(_l.length()), ec(_eta/_c)
		{
			set_temporal_boundaries();
		}

		Real pdf(Real theta) const
		{
			return dtdtheta(theta)/(t_b - t_a);
		}

		Real cdf(Real theta) const
		{
			return (t_theta(theta) - t_a)/(t_b - t_a);
		}

		Real cdf_inv(Real e) const
		{
			return theta_t(e*(t_b - t_a) + t_a);
		}
	};
}; //SamplingTechniques

#endif //_SAMPLING_TECHNIQUES_H_
