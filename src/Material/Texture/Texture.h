/*
 * Copyright (C) 2018, Ibon Guillen (http://giga.cps.unizar.es/~ibon/)
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

#ifndef _TEXTURE_H_
#define _TEXTURE_H_

#include "bunnykiller.h"

#include <array>
#include <cmath>

#include "Image/RGBColor.h"
#include "LinearAlgebra/Vector2.h"
#include "Utils/Utils.h"

class Texture
{
public:
	enum Filtering
	{
		NearestNeighbor, Bilinear
	};

	enum WrapMode
	{
		Black, Clamp, Repeat, Mirror
	};
protected:
	Imaging::Image<Real> m_image;
	WrapMode m_wrap_u, m_wrap_v;
	Filtering m_filtering;
public:
	Texture() :
			m_image(),
			m_wrap_u(Repeat),
			m_wrap_v(Repeat),
			m_filtering(NearestNeighbor)
	{
	}
	
	Texture(Imaging::Image<Real> image, WrapMode wrap_u = Repeat, WrapMode wrap_v = Repeat,
		Filtering filtering = NearestNeighbor) :
			m_image(image),
			m_wrap_u(wrap_u),
			m_wrap_v(wrap_v),
			m_filtering(filtering)
	{
		m_image.normalize(1.);
	}
	
	virtual ~Texture()
	{
	}

	unsigned width() const
	{
		return m_image.width();
	}

	unsigned height() const
	{
		return m_image.height();
	}

	void set_filtering(Filtering filtering)
	{
		m_filtering = filtering;
	}

	void set_wrapping(WrapMode wrap_u, WrapMode wrap_v)
	{
		m_wrap_u = wrap_u;
		m_wrap_v = wrap_v;
	}

	/* Access data */
	virtual const Spectrum operator()(float u, float v) const
	{
		/* Vertices without texture coordinates have negative uvs,
		 reject them */
		if (u < 0. || v < 0.) {
			return Spectrum(0.);
		}

		/* Check if uv coordinates are within limits */
		if (1. <= u) {
			switch (m_wrap_u) {
				case Black:
					return Spectrum(0.);
				case Clamp:
					u = Utils::clamp(u, 0., std::nexttoward(1., 0.));
					break;
				case Repeat:
					u = std::fmod(u, 1.0);
					break;
				case Mirror:
					u = 1. - std::fmod(u, 1.0);
					break;
			}
		}

		if (1. <= v) {
			switch (m_wrap_v) {
				case Black:
					return Spectrum(0.);
				case Clamp:
					v = Utils::clamp(v, 0., std::nexttoward(1., 0.));
					break;
				case Repeat:
					v = std::fmod(v, 1.0);
					break;
				case Mirror:
					v = 1. - std::fmod(v, 1.0);
					break;
			}
		}

		/* Access texture */
		Imaging::RGBColor rgb;

		switch (m_filtering) {
			default:
			case NearestNeighbor: {
				unsigned x = m_image.width() * u;
				unsigned y = m_image.height() * v;

				rgb = m_image(x, y);
			}
			break;
			case Bilinear: {
				float x = (m_image.width() - 1) * u;
				float y = (m_image.height() - 1) * v;

				unsigned x_0 = std::floor(x);
				unsigned x_1 = x_0 + 1;
				unsigned y_0 = std::floor(y);
				unsigned y_1 = y_0 + 1;

				float x_p = x - x_0;
				float y_p = y - y_0;

				Imaging::RGBColor rgb00 = m_image(x_0, y_0);
				Imaging::RGBColor rgb01 = m_image(x_0, y_1);
				Imaging::RGBColor rgb10 = m_image(x_1, y_0);
				Imaging::RGBColor rgb11 = m_image(x_1, y_1);

				rgb = (rgb00 * (1. - x_p) + rgb10 * x_p) * (1. - y_p)
						+ (rgb01 * (1. - x_p) + rgb11 * x_p) * y_p;
			}
			break;
		};

		return Spectrum(rgb.r(), rgb.g(), rgb.b());
	}

	inline const Spectrum operator()(const Vector2 &uv) const
	{
		return (*this)((float) uv[0], (float) uv[1]);
	}

	/* Direct access to image data (ignores wrapping and filtering) */
	virtual const Spectrum operator()(unsigned x, unsigned y) const
	{
		Imaging::RGBColor rgb = m_image(x, y);
		return Spectrum(rgb.r(), rgb.g(), rgb.b());
	}
}; // Texture

#endif //_TEXTURE_H_
