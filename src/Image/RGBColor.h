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

#ifndef _IMAGE_RGBCOLOR_H_
#define _IMAGE_RGBCOLOR_H_

#include "bunnykiller.h"

#include <vector>

namespace Imaging
{
class RGBColor
{
	Real data[3];

  public:
	RGBColor() :
		data {0., 0., 0.}
	{}

	explicit RGBColor(Real c) :
		data {c, c, c}
	{}

	explicit RGBColor(Real r, Real g, Real b) :
		data {r, g, b}
	{}

	explicit RGBColor(const std::vector<Real> &that) :
		data {that[0], that[1], that[2]}
	{}

	explicit RGBColor(const Real* that) :
		data {that[0], that[1], that[2]}
	{}

	const RGBColor& operator=(const std::vector<Real> &that)
	{
		for (size_t i = 0; i < 3; i++)
			(*this)[i] = that[i];
		return (*this);
	}

	Real& operator[](const unsigned i)
	{
		return data[i];
	}

	const Real& operator[](const unsigned i) const
	{
		return data[i];
	}

	Real r() const
	{
		return (*this)[0];
	}

	Real g() const
	{
		return (*this)[1];
	}

	Real b() const
	{
		return (*this)[2];
	}

	const Real* pointer() const
	{
		return &((*this)[0]);
	}

	RGBColor &operator+=(const RGBColor &s)
	{
		for (size_t i = 0; i < 3; i++)
			(*this)[i] += s[i];
		return (*this);
	}

	RGBColor operator+(const RGBColor &s) const
	{
		RGBColor sol = (*this);
		return sol += s;
	}

	RGBColor& operator-=(const RGBColor &s)
	{
		for (size_t i = 0; i < 3; i++)
			(*this)[i] -= s[i];
		return (*this);
	}

	RGBColor operator-(const RGBColor &s) const
	{
		RGBColor sol = (*this);
		return sol -= s;
	}

	RGBColor& operator*=(const RGBColor &s)
	{
		for (size_t i = 0; i < 3; i++)
			(*this)[i] *= s[i];
		return (*this);
	}

	RGBColor operator*(const RGBColor &s) const
	{
		RGBColor sol = (*this);
		return sol *= s;
	}

	RGBColor& operator*=(float f)
	{
		for (size_t i = 0; i < 3; i++)
			(*this)[i] *= f;
		return (*this);
	}

	RGBColor operator*(float f) const
	{
		RGBColor sol = (*this);
		return sol *= f;
	}

	RGBColor& operator*=(double f)
	{
		for (size_t i = 0; i < 3; i++)
			(*this)[i] *= f;
		return (*this);
	}

	RGBColor operator*(double f) const
	{
		RGBColor sol = (*this);
		return sol *= f;
	}

	RGBColor& operator/=(float f)
	{
		for (size_t i = 0; i < 3; i++)
			(*this)[i] /= f;
		return (*this);
	}

	RGBColor operator/(float f) const
	{
		RGBColor sol = (*this);
		return sol /= f;
	}

	RGBColor& operator/=(double f)
	{
		for (size_t i = 0; i < 3; i++)
			(*this)[i] /= f;
		return (*this);
	}

	RGBColor operator/(double f) const
	{
		RGBColor sol = (*this);
		return sol /= f;
	}
};

template <typename T>
RGBColor operator*(const T &t, const RGBColor &s)
{
	return s * t;
}

template <typename Stream>
Stream &operator<<(Stream &os, const RGBColor &s)
{
	os.output(s);
	return os;
}

}; // namespace Imaging

#endif // _IMAGE_RGBCOLOR_H_
