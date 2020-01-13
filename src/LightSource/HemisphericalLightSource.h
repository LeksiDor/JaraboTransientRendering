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

#ifndef _HEMISPHERICAL_LIGHT_SOURCE_H
#define _HEMISPHERICAL_LIGHT_SOURCE_H

#include <math.h>

#include "bunnykiller.h"
#include "LightSource/PointLightSource.h"
#include "RayTracing/World.h"

/** HemisphericalLightSource class represents a point light emiting
light with a cosine distribution. */
template<unsigned int D, class Radiance>
class HemisphericalLightSource : public PointLightSource<D, Radiance>
{
protected:
	typedef PointLightSource<D, Radiance> PLS;

	VectorN<D> direction;
public:
	HemisphericalLightSource(World<D, Radiance>* _world, const VectorN<D>& pos,
							 const VectorN<D>& dir, const Radiance& ints, Real time = 0.) :
		PointLightSource<D, Radiance>(_world, pos, ints, time), direction(dir)
	{}

	virtual ~HemisphericalLightSource() {}

	Radiance get_incoming_light(const VectorN<D> &point_lighted) const;
	bool sample(const Intersection<D> &it, LightSample<D, Radiance> &light_sample, Real &pdf)const;
	void sample(LightSample<D, Radiance> &light_sample, Real &pdf)const;
	bool sample(const VectorN<D> &p, LightSample<D, Radiance> &light_sample, Real &pdf)const;

	virtual void print(FILE *_f_log)const;
}; //HemisphericalLightSource

template<unsigned int D, class Radiance>
Radiance HemisphericalLightSource<D, Radiance>::get_incoming_light(const VectorN<D> &point_lighted) const
{
	VectorN<D> dir = point_lighted - PLS::position;
	Real inv_distance2 = 1. / dir.length2();
	Real inv_distance = sqrt(inv_distance2);

	dir *= inv_distance;

	if (dot(direction, dir) <= 0)
		return Radiance(0);

	Real att(0.);
	if (D == 2)
		att = inv_distance;

	if (D == 3)
		att = inv_distance2;

	return PLS::intensities*att;
}

template<unsigned int D, class Radiance>
bool HemisphericalLightSource<D, Radiance>::sample(const Intersection<D> &it, LightSample<D, Radiance> &light_sample, Real &pdf)const
{
	light_sample.pos = PLS::position;
	light_sample.dir = it.get_position() - PLS::position;
	Real distance2 = light_sample.dir.length2();
	light_sample.dist = sqrt(distance2);
	light_sample.dir /= light_sample.dist;
	light_sample.irradiance = 0.;
	light_sample.instant = PLS::time;

	pdf = 1.;

	if (it.get_normal().dot(-light_sample.dir) > 0. && dot(light_sample.dir, direction)>0. && PLS::is_visible(it.get_position()))
	{
		Real att(0.);
		if (D == 2)
			att = 1. / light_sample.dist;

		if (D == 3)
			att = 1. / distance2;

		light_sample.irradiance = PLS::intensities * att;
		return true;
	}

	return false;
}

template<unsigned int D, class Radiance>
bool HemisphericalLightSource<D, Radiance>::sample(const VectorN<D> &p, LightSample<D, Radiance> &light_sample, Real &pdf)const
{
	light_sample.pos = PLS::position;
	light_sample.dir = p - PLS::position;
	Real distance2 = light_sample.dir.length2();
	light_sample.dist = sqrt(distance2);
	light_sample.dir /= light_sample.dist;
	light_sample.irradiance = 0.;
	light_sample.instant = PLS::time;

	pdf = 1.;

	if (dot(light_sample.dir, direction)>0. && PLS::is_visible(p))
	{
		Real att(0.);
		if (D == 2)
			att = 1. / light_sample.dist;

		if (D == 3)
			att = 1. / distance2;

		light_sample.irradiance = PLS::intensities * att;
		return true;
	}

	return false;
}

template<unsigned int D, class Radiance>
void HemisphericalLightSource<D, Radiance>::sample(LightSample<D, Radiance> &light_sample, Real &pdf)const
{
	light_sample.pos = PLS::position;
	light_sample.dist = 0.;
	light_sample.irradiance = PLS::intensities;
	light_sample.instant = PLS::time;

	// Sample semicircle
	Sampling.hemispherical(light_sample.dir, pdf);
	light_sample.dir.transform_matrix_to(VectorN<D>(0, 1), direction);
}

template<unsigned int D, class Radiance>
void HemisphericalLightSource<D, Radiance>::print(FILE *_f_log)const
{
	fprintf(_f_log, "Hemispherical Light Source:\n-Position [%f %f %f], Normal [%f %f %f]\n\n", PLS::position[0], PLS::position[1], PLS::position[2], direction[0], direction[1], direction[2]);
}

#endif //_HEMISPHERICAL_LIGHT_SOURCE_H
