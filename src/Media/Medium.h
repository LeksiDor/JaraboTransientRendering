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

#ifndef _MEDIUM_H_
#define _MEDIUM_H_

#include "bunnykiller.h"

#include "Color/Spectrum.h"
#include "Color/PolarizedAttenuation.h"
#include "Color/FluorescentAttenuation.h"
#include "Color/FluorescentPolarizedAttenuation.h"
#include "LinearAlgebra/VectorN.h"
#include "RayTracing/Ray.h"

/** Base virtual class for participating media. 
	It includes all functions describing the medium, including albedo, transmittance, mean free path
	and evaluation of the phase function.  
	The templatized parameter D defines the number of dimensions of the medium*/
template<unsigned D>
class Medium
{
public:
	virtual ~Medium() {}

	virtual Real get_mean_free_path(const Ray<D> &r)const = 0;
	virtual Real sample_mean_free_path(const Ray<D> &r, Real &pdf ) const = 0;
	virtual Real p(const Ray<D> &r, const Real d ) const = 0;
	virtual Spectrum get_transmittance(const Ray<D> &r)const = 0;
	virtual Spectrum get_extinction(const VectorN<D> &p)const = 0;
	virtual Spectrum get_max_extinction()const = 0;
	virtual Spectrum get_albedo(const VectorN<D> &p)const = 0;
	virtual Spectrum get_absorption(const VectorN<D> &p)const = 0;
	virtual Spectrum get_scattering(const VectorN<D> &p)const = 0;
	
	// Probability of sampling a specific direction 
	virtual Real p(const VectorN<D> &p, const VectorN<D> &wi, const VectorN<D> &wo) const = 0;
	virtual Real p(const VectorN<D> &p, const VectorN<D> &wi, const VectorN<D> &wo, const Real delta_time) const = 0;

	// Scattering functions based on unpolarized spectrum
	virtual void f(const VectorN<D> &p, const VectorN<D> &wi, const VectorN<D> &wo, Spectrum &R) const = 0;
	virtual void f(const VectorN<D> &p, const VectorN<D> &wi, const VectorN<D> &wo, Spectrum &R, Real &delta_time, Real &pdf)const = 0;

	virtual void sample_direction(const VectorN<D> &p, const VectorN<D> &wi, VectorN<D> &wo, Spectrum &R, Real &pdf) const = 0;
	virtual void sample_direction(const VectorN<D> &p, const VectorN<D> &wi, VectorN<D> &wo, Spectrum &R, Real &delta_time, Real &pdf)const = 0;

	// Scattering functions based on polarized spectrum (Muller Matrices)
	virtual void f(const VectorN<D> &p, const VectorN<D> &wi, const VectorN<D> &wo, PolarizedAttenuation<D> &R) const = 0;
	virtual void f(const VectorN<D> &p, const VectorN<D> &wi, const VectorN<D> &wo, PolarizedAttenuation<D> &R, Real &delta_time, Real &pdf)const = 0;

	virtual void sample_direction(const VectorN<D> &p, const VectorN<D> &wi, VectorN<D> &wo, PolarizedAttenuation<D> &R, Real &pdf)const = 0;
	virtual void sample_direction(const VectorN<D> &p, const VectorN<D> &wi, VectorN<D> &wo, PolarizedAttenuation<D> &R, Real &delta_time, Real &pdf)const = 0;

	// Scattering functions based on radiance
	virtual void f(const VectorN<D> &p, const VectorN<D> &wi, const VectorN<D> &wo, FluorescentAttenuation<D> &R) const = 0;
	virtual void f(const VectorN<D> &p, const VectorN<D> &wi, const VectorN<D> &wo, FluorescentAttenuation<D> &R, Real &delta_time, Real &pdf)const = 0;

	virtual void sample_direction(const VectorN<D> &p, const VectorN<D> &wi, VectorN<D> &wo, FluorescentAttenuation<D> &R, Real &pdf)const = 0;
	virtual void sample_direction(const VectorN<D> &p, const VectorN<D> &wi, VectorN<D> &wo, FluorescentAttenuation<D> &R, Real &delta_time, Real &pdf)const = 0;

	// Scattering functions based on polarized fluorescent spectrum (Muller Matrices)
	virtual void f(const VectorN<D> &p, const VectorN<D> &wi, const VectorN<D> &wo, FluorescentPolarizedAttenuation<D> &R) const = 0;
	virtual void f(const VectorN<D> &p, const VectorN<D> &wi, const VectorN<D> &wo, FluorescentPolarizedAttenuation<D> &R, Real &delta_time, Real &pdf) const = 0;

	virtual void sample_direction(const VectorN<D> &p, const VectorN<D> &wi, VectorN<D> &wo, FluorescentPolarizedAttenuation<D> &R, Real &pdf) const = 0;
	virtual void sample_direction(const VectorN<D> &p, const VectorN<D> &wi, VectorN<D> &wo, FluorescentPolarizedAttenuation<D> &R, Real &delta_time, Real &pdf) const = 0;
}; // Medium


template<unsigned D>
Medium<D>* Material<D>::m_default_medium = nullptr;

template<unsigned D>
Real Material<D>::m_default_n = 1.;

#endif //_MEDIUM_H_
