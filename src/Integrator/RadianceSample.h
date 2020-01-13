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

#ifndef _RADIANCE_SAMPLE_H_
#define _RADIANCE_SAMPLE_H_

#include "bunnykiller.h"

#include <vector>

#include "Color/Spectrum.h"
#include "LinearAlgebra/VectorN.h"

template<unsigned D, class Radiance>
struct RadianceSample
{
	Radiance radiance;
	Real time; // meters
#ifdef _USE_PM_
	Real camera_range; //camera spatial time range of radiance sample in 2D blur
	Real beam_range; //beam time range of radiance sample for beams 2D blur
#endif
	unsigned int bounce;

	constexpr RadianceSample() :
		radiance(0.),
		time(0.),
#ifdef _USE_PM_
		camera_range(0.),
		beam_range(0.),
#endif
		bounce(0)
	{}

	RadianceSample(const Radiance &_radiance, Real _time = 0., int _bounce = 0) :
		radiance(_radiance),
		time(_time),
#ifdef _USE_PM_
		camera_range(0.),
		beam_range(0.),
#endif
		bounce(_bounce)
	{}

	RadianceSample(const Radiance &_radiance, int _bounce) :
		radiance(_radiance),
		time(0.),
#ifdef _USE_PM_
		camera_range(0.),
		beam_range(0.),
#endif
		bounce(_bounce)
	{}
}; // RadianceSample

template<unsigned D, class Radiance>
struct RadianceSampleRecord
{
	RadianceSample<D, Radiance> sample;
	Real distance;
	VectorN<D> pos;
	VectorN<D> normal;

	RadianceSampleRecord(RadianceSample<D, Radiance> _sample, Real _distance = 0.,
			VectorN<D> _pos = {0., 0., 0.}, VectorN<D> _norm = {0., 0., 0.}) :
		sample(_sample),
		distance(_distance),
		pos(_pos),
		normal(_norm)
	{}
}; // RadianceSampleRecord

template<unsigned D, class Radiance>
struct RadianceSampleRecordVector
{
	std::vector<RadianceSample<D, Radiance>> samples;
	Real distance;
	VectorN<D> pos;
	VectorN<D> normal;

	RadianceSampleRecordVector(Real _distance = 0., VectorN<D> _pos = {0., 0., 0.},
			VectorN<D> _norm = {0., 0., 0.}) :
		samples(),
		distance(_distance),
		pos(_pos),
		normal(_norm)
	{}

	RadianceSampleRecordVector(std::vector<RadianceSample<D, Radiance>> _samples,
			Real _distance = 0., VectorN<D> _pos = {0., 0., 0.}, VectorN<D> _norm = {0., 0., 0.}) :
		samples(_samples),
		distance(_distance),
		pos(_pos),
		normal(_norm)
	{}

	void add_sample(const RadianceSample<D, Radiance>& sample)
	{
		samples.push_back(sample);
	}

	bool has_samples() const
	{
		return !samples.empty();
	}
}; // RadianceSampleRecord

#endif // _RADIANCE_SAMPLE_H_
