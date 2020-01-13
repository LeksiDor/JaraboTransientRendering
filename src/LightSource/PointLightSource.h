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

#ifndef _POINT_LIGHT_SOURCE_H_
#define _POINT_LIGHT_SOURCE_H_

#include "bunnykiller.h"

#include "Color/Spectrum.h"
#include "LightSource/LightSource.h"
#include "LinearAlgebra/Vector3.h"
#include "RayTracing/World.h"
#include "RayTracing/Ray.h"
#include "RayTracing/Intersection.h"
#include "Utils/RandomNumbers.h"

/** PointLightSource class represents a single point sending light equally
		in all directions. */
template<unsigned D, class Radiance>
class PointLightSource : public LightSource<D, Radiance>
{
protected:
	typedef LightSource<D, Radiance> LS;

	VectorN<D> position;
public:
	PointLightSource(World<D, Radiance>* _world, const VectorN<D>& pos, const Radiance& ints, Real time = 0.) :
		LightSource<D, Radiance>(_world, ints, time), position(pos) {}

	virtual ~PointLightSource() {}
  
	VectorN<D> get_position() const
	{
		return position;
	}

	VectorN<D> get_incoming_direction(const VectorN<D> &point_lighted) const
	{
		return normalize(point_lighted-position);
	}

	bool is_visible(const VectorN<D> &point_lighted) const
	{
		VectorN<D> dir = point_lighted - position;
		Real l = dir.length();
		dir /= l;

		Ray<D> ray(position, dir, false);

		return !this->world->intersects(ray, l);
	}

	virtual Radiance get_incoming_light(const VectorN<D> &point_lighted) const;

	virtual bool sample(const Intersection<D> &it, LightSample<D, Radiance> &light_sample, Real &pdf) const;

	virtual void sample(LightSample<D, Radiance> &light_sample, Real &pdf) const;
	
	virtual bool sample(const VectorN<D> &p, LightSample<D, Radiance> &light_sample, Real &pdf) const;

	virtual void print(FILE* f_log) const;
}; // PointLightSource

// Specializations for 2D and 3D
template<unsigned D, class Radiance>
Radiance PointLightSource<D, Radiance>::get_incoming_light(const VectorN<D> &point_lighted) const
{
	Real inv_distance2 = 1./(point_lighted - position).length2();
	Real inv_distance = sqrt(inv_distance2);

	Real att(0.);
	if (D == 2)
		att = inv_distance;

	if (D == 3)
		att = inv_distance2;

	Radiance light = this->intensities*att;
	VectorN<D> dir = (point_lighted - position)*inv_distance;
	set_direction(dir, light);

	return light;
}

template<unsigned D, class Radiance>
bool PointLightSource<D, Radiance>::sample(const Intersection<D> &it, LightSample<D, Radiance> &light_sample, Real &pdf) const
{
	light_sample.irradiance = 0.;
	light_sample.pos = position;
	light_sample.dir = it.get_position() - position;

	/* Check surface orientation */
	if (it.get_normal().dot(-light_sample.dir) <= 0.) {
		return false;
	}

	Real distance2 = light_sample.dir.length2();
	light_sample.dist = std::sqrt(distance2);
	light_sample.dir /= light_sample.dist;

	/* Check visibility */
	Ray<D> ray(light_sample.pos, light_sample.dir, false);
	if (this->world->intersects(ray, light_sample.dist)) {
		return false;
	}

	Real att(0.);
	if (D == 2)
		att = 1. / light_sample.dist;
	if (D == 3)
		att = 1. / distance2;

	light_sample.irradiance = LS::intensities * att;
	set_direction(light_sample.dir, light_sample.irradiance);
	light_sample.instant = LS::time;

	pdf = 1.;

	return true;
}

template<unsigned D, class Radiance>
bool PointLightSource<D, Radiance>::sample(const VectorN<D> &p, LightSample<D, Radiance> &light_sample, Real &pdf) const
{
	light_sample.irradiance = 0.;
	light_sample.pos = position;
	light_sample.dir = p - position;

	Real distance2 = light_sample.dir.length2();
	light_sample.dist = std::sqrt(distance2);
	light_sample.dir /= light_sample.dist;

	/* Check visibility */
	Ray<D> ray(light_sample.pos, light_sample.dir, false);
	if (this->world->intersects(ray, light_sample.dist)) {
		return false;
	}

	Real att(0.);
	if (D == 2)
		att = 1. / light_sample.dist;
	if (D == 3)
		att = 1. / distance2;

	light_sample.irradiance = LS::intensities * att;
	set_direction(light_sample.dir, light_sample.irradiance);
	light_sample.instant = LS::time;

	pdf = 1.;

	return true;
}

template<unsigned D, class Radiance>
void PointLightSource<D, Radiance>::sample(LightSample<D, Radiance> &light_sample, Real &pdf) const
{
	light_sample.pos = position;
	light_sample.dist = 0.;
	light_sample.irradiance = LS::intensities;
	light_sample.instant = LS::time;

	Sampling.direction_uniformly(light_sample.dir, pdf);
	set_direction(light_sample.dir, light_sample.irradiance);
}

template<unsigned D, class Radiance>
void PointLightSource<D, Radiance>::print(FILE* f_log) const
{
	fprintf(f_log, "Point Light Source:\n");
	fprintf(f_log, "- Position: [%f %f %f]\n", position[0], position[1], position[2]);
	fprintf(f_log, "- Time: %f\n", LS::time);
	LS::intensities.print(f_log);
}

#endif //_POINT_LIGHT_SOURCE_H_
