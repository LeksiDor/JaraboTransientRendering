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

#ifndef _DIRECTIONAL_LIGHT_SOURCE_H_
#define _DIRECTIONAL_LIGHT_SOURCE_H_

#include "LightSource/LightSource.h"
#include "LinearAlgebra/Vector3.h"
#include "RayTracing/World.h"
#include "RayTracing/Ray.h"
#include "RayTracing/Intersection.h"

/** DirectionalLightSource class represents a directional light sending light equally
 in a single direction without attenuation. */
template<int D, class Radiance>
class DirectionalLightSource : public LightSource<D, Radiance>
{
protected:
	typedef LightSource<D, Radiance> LS;

	VectorN<D> direction;
	
	VectorN<D> m_sampling_position;
	VectorN<D> m_sampling_frame[D-1];
	Real m_sampling_radius, m_sampling_pdf;
	bool m_precomputed;

public:
	DirectionalLightSource(World<D, Radiance>* _world, const VectorN<D>& dir, const Radiance& ints, Real time = 0.) :
		LightSource<D, Radiance>(_world, ints, time), direction(normalize(dir)), m_precomputed(false)
	{}
	
	virtual ~DirectionalLightSource() {}

	VectorN<D> get_direction() const { return direction; }
	VectorN<D> get_incoming_direction(const VectorN<D> &point_lighted) const
	{
		return direction;
	}
	bool is_visible(const VectorN<D> &point_lighted) const
	{
		VectorN<D> dir = direction;
		Ray<D> ray(point_lighted, -dir, false);
		
		return( !this->world->intersects(ray)) ;
	}
	virtual Radiance get_incoming_light(const VectorN<D> &point_lighted) const;
	
	virtual bool sample(const Intersection<D> &it, LightSample<D, Radiance> &light_sample, Real &pdf)const;
	virtual bool sample(const VectorN<D> &p, LightSample<D, Radiance> &light_sample, Real &pdf)const;
	
	virtual void sample(LightSample<D, Radiance> &light_sample, Real &pdf)const;


	void setup()
	{
		m_precomputed = true;

		VectorN<D> mid = LS::world->get_bounding_volume().get_center();
		Real r = LS::world->get_bounding_volume().get_diagonal();

		m_sampling_position = mid - direction*r * 2;

		if (D == 3)
		{
			m_sampling_frame[1] = cross(direction, VectorN<D>(0, 1));
			m_sampling_frame[0] = cross(m_sampling_frame[1], direction);
			m_sampling_pdf = 1 / (4*r*r);
		}

		m_sampling_radius = 2*r;
	}

}; //DirectionalLightSource



template<int D, class Radiance>
Radiance DirectionalLightSource<D, Radiance>::get_incoming_light(const VectorN<D> &point_lighted) const
{
	return is_visible(point_lighted)  ? LS::intensities : Spectrum(0.);
}

template<int D, class Radiance>
bool DirectionalLightSource<D, Radiance>::sample(const Intersection<D> &it, LightSample<D, Radiance> &light_sample, Real &pdf)const
{
	light_sample.pos = VectorN<D>(0.);
	light_sample.dir = direction;
	light_sample.dist = 0.;
	light_sample.irradiance = 0.;
	light_sample.instant = LS::time;
	
	pdf = 1.;
	
	if (it.get_normal().dot(-light_sample.dir) > 0. && is_visible(it.get_position()))
	{
		light_sample.irradiance = LS::intensities;
		return true;
	}
	
	return false;
}

template<int D, class Radiance>
bool DirectionalLightSource<D, Radiance>::sample(const VectorN<D> &p, LightSample<D, Radiance> &light_sample, Real &pdf)const
{
	light_sample.pos = VectorN<D>(0.);
	light_sample.dir = direction;
	light_sample.dist = 0.;
	light_sample.irradiance = 0.;
	light_sample.instant = LS::time;
	
	pdf = 1.;
	
	if (is_visible(p))
	{
		light_sample.irradiance = LS::intensities;
		return true;
	}
	
	return false;
}

template<int D, class Radiance>
void DirectionalLightSource<D, Radiance>::sample(LightSample<D, Radiance> &light_sample, Real &pdf)const
{
	if (!m_precomputed)
		throw("Cannot sample directional light without precomputation");
	
	
	Real u, v;
	Sampling.jittered(u, v, pdf, -m_sampling_radius, m_sampling_radius);

	pdf = m_sampling_pdf;

	light_sample.pos = m_sampling_position + m_sampling_frame[0] * u + m_sampling_frame[1] * v;
	light_sample.dist = 0.;
	light_sample.irradiance = LS::intensities;
	light_sample.instant = LS::time;
	light_sample.dir = direction;

	//printf("Pos: %f, %f, %f - Dir: %f, %f, %f\n", light_sample.pos[0], light_sample.pos[1], light_sample.pos[2], light_sample.dir[0], light_sample.dir[1], light_sample.dir[2]);
}



#endif //_DIRECTIONAL_LIGHT_SOURCE_H_
