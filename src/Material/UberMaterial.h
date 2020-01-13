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

#ifndef _UBER_MATERIAL_H_
#define _UBER_MATERIAL_H_

#include "bunnykiller.h"

#include <vector>

#include "Material/Material.h"
#include "RayTracing/Ray.h"
#include "RayTracing/Intersection.h"
#include "Material/Reflectance/BSDF.h"
#include "Utils/RandomNumbers.h"

template<unsigned D> class BSDF;

template<unsigned D>
class UberMaterial: public Material<D>
{
protected:
	typedef Material<D> MTR;

protected:
	std::vector<BSDF<D>*> bsdfs;
	std::vector<Real> weights;
	Real total_weight, inv_total_weight;
public:
	UberMaterial() : total_weight(0.), inv_total_weight(0.) {}

	void add_bsdf(BSDF<D>* bsdf, const Real weight = 1.)
	{
		bsdfs.push_back(bsdf);
		weights.push_back(weight);
		total_weight += weight;
		inv_total_weight = 1./total_weight;

		MTR::m_type = (Reflectance::Type)(MTR::m_type | bsdf->get_type());
	}

	Real p(const VectorN<D> &omega_i, const VectorN<D> &omega_o,
			const VectorN<D> &normal, const Vector2 &uv)const
	{
		Real pdf(0.);
		for (size_t i = 0; i<bsdfs.size(); ++i)
			pdf += bsdfs[i]->p(omega_i, omega_o, normal, uv)*weights[i];

		return pdf * inv_total_weight;
	}

	void f(const VectorN<D> &omega_i, const VectorN<D> &omega_o,
			const VectorN<D> &normal, const Vector2 &uv, Spectrum &R) const
	{
		R = Spectrum(0.);
		for (size_t i = 0; i < bsdfs.size(); ++i) {
			Spectrum Ri;
			bsdfs[i]->f(omega_i, omega_o, normal, uv, Ri);
			R += Ri*weights[i];
		}

		R *= inv_total_weight;
	}
	
	void sample_direction(const VectorN<D> &omega_i, VectorN<D> &omega_o,
			const VectorN<D> &normal, const Vector2 &uv, Spectrum &R, Real &pdf) const
	{
		// Need for importance sampling...
		int idx = (int)(Random::StdRNG.next_real()*bsdfs.size());
		
		bsdfs[idx]->sample_direction(omega_i, omega_o, normal, uv, R, pdf);
		R *= weights[idx] / total_weight;
		pdf /= (Real)(bsdfs.size());
	}
	
	void sample_outgoing_ray(const Intersection<D> &it, Ray<D> &new_ray,
			Spectrum &R, Real &pdf) const
	{
		// Need for importance sampling...
		int idx = (int)(Random::StdRNG.next_real()*bsdfs.size());

		bsdfs[idx]->sample_outgoing_ray(it, new_ray, R, pdf); 
		R *= weights[idx] / total_weight;
		pdf /= (Real)(bsdfs.size());
	}

	void f(const VectorN<D> &omega_i, const VectorN<D> &omega_o, const VectorN<D> &normal,
			const Vector2 &uv, PolarizedAttenuation<D> &R) const
	{
//		R = PolarizedAttenuation<D>(0.);
//		for (int i = 0; i < bsdfs.size(); ++i)
//		{
//			PolarizedAttenuation<D> Ri;
//			bsdfs[i]->f(omega_i, omega_o, normal, uv, Ri);
//			R += Ri*weights[i];
//		}
//
//		R *= inv_total_weight;

		bsdfs[0]->f(omega_i, omega_o, normal, uv, R);
	}

	void sample_direction(const VectorN<D> &omega_i, VectorN<D> &omega_o,
			const VectorN<D> &normal, const Vector2 &uv, PolarizedAttenuation<D> &R,
			Real &pdf) const
	{
		// Need for importance sampling...
		int idx = Random::StdRNG.next_real()*bsdfs.size();

		bsdfs[idx]->sample_direction(omega_i, omega_o, normal, uv, R, pdf);
		R *= weights[idx] / total_weight;
		pdf /= (Real)(bsdfs.size());
	}

	void sample_outgoing_ray(const Intersection<D> &it, Ray<D> &new_ray,
			PolarizedAttenuation<D> &R, Real &pdf) const
	{
		// Need for importance sampling...
		int idx = Random::StdRNG.next_real()*bsdfs.size();

		bsdfs[idx]->sample_outgoing_ray(it, new_ray, R, pdf);
		R *= weights[idx] / total_weight;
		pdf /= (Real)(bsdfs.size());
	}

	void f(const VectorN<D> &omega_i, const VectorN<D> &omega_o, const VectorN<D> &normal,
			const Vector2 &uv, FluorescentPolarizedAttenuation<D> &R) const
	{
		/*R = PolarizedLightAttenuation<D>(0);
		for (int i = 0; i < bsdfs.size(); ++i)
		{
		PolarizedLightAttenuation<D> Ri;
		bsdfs[i]->f(omega_i, omega_o, normal, uv, Ri);
		R += Ri*weights[i];
		}

		R *= inv_total_weight;*/

		bsdfs[0]->f(omega_i, omega_o, normal, uv, R);
	}

	void sample_direction(const VectorN<D> &omega_i, VectorN<D> &omega_o,
			const VectorN<D> &normal, const Vector2 &uv, FluorescentPolarizedAttenuation<D> &R,
			Real &pdf) const
	{
		// Need for importance sampling...
		int idx = Random::StdRNG.next_real()*bsdfs.size();

		bsdfs[idx]->sample_direction(omega_i, omega_o, normal, uv, R, pdf);
		R *= weights[idx] / total_weight;
		pdf /= (Real)(bsdfs.size());
	}

	void sample_outgoing_ray(const Intersection<D> &it, Ray<D> &new_ray,
			FluorescentPolarizedAttenuation<D> &R, Real &pdf) const
	{
		// Need for importance sampling...
		int idx = Random::StdRNG.next_real()*bsdfs.size();

		bsdfs[idx]->sample_outgoing_ray(it, new_ray, R, pdf);
		R *= weights[idx] / total_weight;
		pdf /= (Real)(bsdfs.size());
	}

	void f(const VectorN<D> &omega_i, const VectorN<D> &omega_o, const VectorN<D> &normal,
			const Vector2 &uv, FluorescentAttenuation<D> &R) const
	{
		bsdfs[0]->f(omega_i, omega_o, normal, uv, R);
	}

	void sample_direction(const VectorN<D> &omega_i, VectorN<D> &omega_o,
			const VectorN<D> &normal, const Vector2 &uv, FluorescentAttenuation<D> &R, Real &pdf) const
	{
		// Need for importance sampling...
		int idx = Random::StdRNG.next_real()*bsdfs.size();

		bsdfs[idx]->sample_direction(omega_i, omega_o, normal, uv, R, pdf);
		R *= weights[idx] / total_weight;
		pdf /= (Real)(bsdfs.size());
	}

	void sample_outgoing_ray(const Intersection<D> &it, Ray<D> &new_ray,
			FluorescentAttenuation<D> &R, Real &pdf) const
	{
		// Need for importance sampling...
		int idx = Random::StdRNG.next_real()*bsdfs.size();

		bsdfs[idx]->sample_outgoing_ray(it, new_ray, R, pdf);
		R *= weights[idx] / total_weight;
		pdf /= (Real)(bsdfs.size());
	}

	Spectrum get_absorption(const Intersection<D> &it) const
	{
		Spectrum result;
		for (size_t i=0; i<bsdfs.size(); ++i)
			result += bsdfs[i]->get_absorption(it.get_uv())*weights[i];

		return result * inv_total_weight;
	}
}; //UberMaterial

typedef UberMaterial<3> UberMaterial3D;
typedef UberMaterial<2> UberMaterial2D;

#endif //_UBER_MATERIAL_H_
