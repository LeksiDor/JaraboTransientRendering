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

#ifndef _SPOT_LIGHT_SOURCE_H
#define _SPOT_LIGHT_SOURCE_H

#include <math.h>

#include "bunnykiller.h"
#include "LightSource/PointLightSource.h"
#include "RayTracing/World.h"

/** SpotLightSource class represents a spot point light emiting. 
	TODO: Comment this properly...*/
template<int D, class Radiance>
class SpotLightSource : public PointLightSource<D, Radiance>
{
protected:
	typedef PointLightSource<D, Radiance> PLSTR;

protected:
	Real m_angle;
	Real m_cos_angle;
	Real m_p;
	VectorN<D> direction;
public:
	SpotLightSource(World<D, Radiance>* _world, const VectorN<D>& pos, const VectorN<D>& target, const Real &angle, const Radiance& ints, Real time = 0.) :
		PointLightSource<D, Radiance>(_world, pos, ints, time), m_angle(angle*M_PI / 360), m_cos_angle((angle>180.) ? 0 : cos(m_angle)),
		  m_p((D==2) ? (m_angle) : (2*M_PI*(1 - m_cos_angle))), direction(normalize(target - pos)) {}
  
	Radiance get_incoming_light(const VectorN<D> &point_lighted) const;
	bool sample(const Intersection<D> &it, LightSample<D, Radiance> &light_sample, Real &pdf)const;
	void sample(LightSample<D, Radiance> &light_sample, Real &pdf)const;
	bool sample(const VectorN<D> &p, LightSample<D, Radiance> &light_sample, Real &pdf)const;

}; //SpotLightSource

template<int D, class Radiance>
Radiance SpotLightSource<D, Radiance>::get_incoming_light(const VectorN<D> &point_lighted) const
{
	VectorN<D> dir = point_lighted - PLSTR::position;
	Real inv_distance2 = 1. / dir.length2();
	Real inv_distance = sqrt(inv_distance2);

	dir *= inv_distance;

	Real att(0.);
	if (D == 2)
		att = inv_distance;

	if (D == 3)
		att = inv_distance2;

	return PLSTR::intensities*att*((dot(direction, dir)<m_cos_angle) ? 0. : 1.);
}	

template<int D, class Radiance>
bool SpotLightSource<D, Radiance>::sample(const Intersection<D> &it, LightSample<D, Radiance> &light_sample, Real &pdf)const
{
	light_sample.pos = PLSTR::position;
	light_sample.dir = it.get_position() - PLSTR::position;
	Real distance2 = light_sample.dir.length2();
	light_sample.dist = sqrt(distance2);
	light_sample.dir /= light_sample.dist;
	light_sample.irradiance = 0.;
	light_sample.instant = PLSTR::time;

	pdf = 1.;

	Real dot_light = dot(light_sample.dir, direction);

	if (it.get_normal().dot(-light_sample.dir) > 0. && dot_light >= m_cos_angle && PLSTR::is_visible(it.get_position()))
	{
		Real att(0.);
		if (D == 2)
			att = 1. / light_sample.dist;

		if (D == 3)
			att = 1. / distance2;

		light_sample.irradiance = PLSTR::intensities * att;
		return true;
	}

	return false;

}

template<int D, class Radiance>
bool SpotLightSource<D, Radiance>::sample(const VectorN<D> &p, LightSample<D, Radiance> &light_sample, Real &pdf)const
{
	light_sample.pos = PLSTR::position;
	light_sample.dir = p - PLSTR::position;
	Real distance2 = light_sample.dir.length2();
	light_sample.dist = sqrt(distance2);
	light_sample.dir /= light_sample.dist;
	light_sample.irradiance = 0.;
	light_sample.instant = PLSTR::time;

	pdf = 1.;

	Real dot_light = dot(light_sample.dir, direction);

	if ( dot_light >= m_cos_angle && PLSTR::is_visible(p))
	{
		Real att(0.);
		if (D == 2)
			att = 1. / light_sample.dist;

		if (D == 3)
			att = 1. / distance2;

		light_sample.irradiance = PLSTR::intensities * att;
		return true;
	}

	return false;
}

template<int D, class Radiance>
void SpotLightSource<D, Radiance>::sample(LightSample<D, Radiance> &light_sample, Real &pdf)const
{
	light_sample.pos = PLSTR::position;
	light_sample.dist = 0.;
	light_sample.irradiance = PLSTR::intensities;
	light_sample.instant = PLSTR::time;

	Sampling.cone(m_angle, m_cos_angle, light_sample.dir, pdf);

	light_sample.dir.transform_matrix_to(VectorN<D>(0, 1), direction);
}

#endif //_SPOT_LIGHT_SOURCE_H


