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

#ifndef _POLARIZATION_FRAME_H_
#define _POLARIZATION_FRAME_H_

#include "bunnykiller.h"

#include <array>
#include <cstdio>

#include "LinearAlgebra/VectorN.h"
#include "Utils/Utils.h"

template<unsigned D>
struct PolarizationFrame
{
	std::array<VectorN<D>, D> m_frame;

	PolarizationFrame();

	PolarizationFrame(const VectorN<D> &up, const VectorN<D> &towards);

	PolarizationFrame(const PolarizationFrame<D> &f)
	{
		m_frame = f.m_frame;
	}

	const PolarizationFrame<D>& operator=(const PolarizationFrame<D> &f)
	{
		m_frame = f.m_frame;
		return *this;
	}

	inline VectorN<D> operator[](const int dim) const
	{
		return m_frame[dim];
	}

	Real get_angle(const PolarizationFrame<D> &f) const;

	void get_cos_sin(const PolarizationFrame<D> &f, Real& cos, Real& sin) const;

	void print(FILE* stream)const
	{
		fprintf(stream, "- Frame: [%f, %f, %f] - [%f, %f, %f]\n",
		   m_frame[1][0], m_frame[1][1], m_frame[1][2],
		   m_frame[0][0], m_frame[0][1], m_frame[0][2]);
	}
};

template<>
PolarizationFrame<3>::PolarizationFrame()
{
	m_frame[0] = VectorN<3>(0., 1., 0.);
	m_frame[1] = VectorN<3>(0., 0., 1.);
	m_frame[2] = VectorN<3>(1., 0., 0.);
}

template<>
PolarizationFrame<2>::PolarizationFrame()
{
	m_frame[0] = VectorN<2>(0., 1.);
	m_frame[1] = VectorN<2>(1., 0.);
}

template<>
PolarizationFrame<3>::PolarizationFrame(const VectorN<3> &up, const VectorN<3> &towards)
{
	m_frame[1] = towards.normalized();
	m_frame[2] = cross(towards, up).normalized();
	m_frame[0] = cross(m_frame[2], m_frame[1]);
}

template<>
PolarizationFrame<2>::PolarizationFrame(const VectorN<2> &up, const VectorN<2> &towards)
{
	m_frame[0] = up;
	m_frame[1] = towards;
}

template<>
Real PolarizationFrame<3>::get_angle(const PolarizationFrame<3> &f) const
{
	if (dot_abs(m_frame[1], f.m_frame[1]) < (1. - 1.e-5)) {
		fprintf(stderr, "\nDOT: %f vs %f\n", dot(m_frame[1], f.m_frame[1]), (1. - 1.e-5));
		this->print(stderr);
		f.print(stderr);
		throw std::runtime_error("Frames non-collinear: Cannot compute the rotation of the polarization frames");
	}

	// Avoid precision problems
	Real cos_alpha = Utils::clamp(dot(m_frame[0], f.m_frame[0]), -1., 1.);
	Real alpha = std::acos(cos_alpha);

	return (dot(m_frame[0], f.m_frame[2]) < 0.) ? alpha : -alpha;
}

template<>
void PolarizationFrame<3>::get_cos_sin(const PolarizationFrame<3> &f, Real& cos, Real& sin) const
{
	if (dot_abs(m_frame[1], f.m_frame[1]) < (1. - 1.e-5)) {
		fprintf(stderr, "\nDOT: %f vs %f\n", dot(m_frame[1], f.m_frame[1]), (1. - 1.e-5));
		this->print(stderr);
		f.print(stderr);
		throw std::runtime_error("Frames non-collinear: Cannot compute the rotation of the polarization frames");
	}

	// Avoid precision problems
	cos = Utils::clamp(dot(m_frame[0], f.m_frame[0]), -1., 1.);
	sin = Utils::clamp(-dot(m_frame[0], f.m_frame[2]), -1., 1.);
}

template<>
Real PolarizationFrame<2>::get_angle(const PolarizationFrame<2> &v) const
{
	return 0.0;
}

template<>
void PolarizationFrame<2>::get_cos_sin(const PolarizationFrame<2> &, Real&, Real&) const
{
}

template<unsigned D>
PolarizationFrame<D> transform_frame(const PolarizationFrame<D> &input, const VectorN<D> &d)
{
	if (D != 3)
		throw std::runtime_error("Unsuported dimensions");

	return PolarizationFrame<D>(input.m_frame[0], d);
}

#endif //_POLARIZATION_FRAME_H_
