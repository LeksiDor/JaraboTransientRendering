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

#ifndef _SAMPLING_H_
#define _SAMPLING_H_

#include "bunnykiller.h"

#include <cmath>

#include "Utils/Utils.h"
#include "Utils/RandomNumbers.h"

template<class RNG>
class TSampling
{
private:
	RNG& rng;
public:
	constexpr TSampling(RNG& _rng) :
			rng(_rng)
	{
	}
public:
	void jittered(Real &x, Real &pdf, Real min, Real max)
	{
		x = rng.next_real();
		x *= (max - min);
		x += min;
		pdf = 1. / (max - min);
	}

	void jittered(Real &x, Real &y, Real &pdf, Real min, Real max)
	{
		Real pdfx;
		jittered(x, pdfx, min, max);
		jittered(y, pdf, min, max);
		pdf *= pdfx;
	}

	void jittered(Real &x, Real &y, Real &z, Real &pdf, Real min, Real max)
	{
		Real pdfx, pdfy;
		jittered(x, pdfx, min, max);
		jittered(y, pdfy, min, max);
		jittered(z, pdf, min, max);
		pdf *= pdfx * pdfy;
	}
	
	void circular(Vector2 &v, Real &pdf)
	{
		Real epsilon1 = rng.next_real();
		Real phi = 2. * M_PI * epsilon1;

		pdf = (.5 * M_1_PI);
		
		v[0] = std::cos(phi);
		v[1] = std::sin(phi);
	}

	void spherical(Vector3 &v, Real &pdf)
	{
		Real epsilon1 = rng.next_real();
		Real epsilon2 = rng.next_real();
		
		Real phi = 2. * M_PI * epsilon1;
		Real cos_theta = 1. - 2. * epsilon2;
		Real sin_theta = std::sqrt(1. - cos_theta * cos_theta);

		pdf = (.25 * M_1_PI);

		v[0] = std::cos(phi) * sin_theta;
		v[1] = cos_theta;
		v[2] = std::sin(phi) * sin_theta;
	}

	inline void direction_uniformly(Vector2 &v, Real &pdf)
	{
		circular(v, pdf);
	}

	inline void direction_uniformly(Vector3 &v, Real &pdf)
	{
		spherical(v, pdf);
	}
	
	void hemispherical(Vector3 &v, Real &pdf)
	{
		Real epsilon1 = rng.next_real();
		Real epsilon2 = rng.next_real();
		
		Real phi = 2. * M_PI * epsilon1;
		Real sin_theta = std::sqrt(1. - epsilon2 * epsilon2);

		pdf = (.5 * M_1_PI);

		v = Vector3(std::cos(phi) * sin_theta, epsilon2, std::sin(phi) * sin_theta);
	}

	void cosine_weighted(Vector2 &v, Real &pdf)
	{
		Real epsilon = rng.next_real();

		Real cos_theta = std::sqrt(epsilon);
		Real sin_theta = std::sqrt(1. - cos_theta * cos_theta);

		v = Vector2(sin_theta, cos_theta);

		pdf = (cos_theta * .5);
	}

	void cosine_weighted(Vector3 &v, Real &pdf)
	{
		Real epsilon1 = rng.next_real();
		Real epsilon2 = rng.next_real();

		Real phi = 2. * M_PI * epsilon1;
		Real cos_theta = std::sqrt(epsilon2);
		Real sin_theta = std::sqrt(1. - cos_theta * cos_theta);

		pdf = (cos_theta * M_1_PI);
		
		v = Vector3(std::cos(phi) * sin_theta, cos_theta, std::sin(phi) * sin_theta);
	}

	void cone(Real angle, Real cos_angle, Vector2 &v, Real &pdf)
	{
		Real epsilon1 = rng.next_real();
		Real phi = (epsilon1 - .5) * angle;

		pdf = Real(1.) / angle;

		v[0] = std::cos(phi);
		v[1] = std::sin(phi);
	}

	void cone(Real angle, Real cos_angle, Vector3 &v, Real &pdf)
	{
		Real epsilon1 = rng.next_real();
		Real epsilon2 = rng.next_real();

		Real cos_theta = 1. - epsilon1 * (1. - cos_angle);
		Real sin_theta = std::sqrt(1. - cos_theta * cos_theta);

		Real phi = 2. * M_PI * epsilon2;

		v = Vector3(std::cos(phi) * sin_theta, cos_theta, std::sin(phi) * sin_theta);
		pdf = 1. / (2. * M_PI * (1. - cos_angle));
	}

}; // Sampling

TSampling<Random::RNG> Sampling = TSampling<Random::RNG>(Random::StdRNG);

#endif // _SAMPLING_H_
