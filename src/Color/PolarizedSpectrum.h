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

#ifndef _POLARIZED_SPECTRUM_H_
#define _POLARIZED_SPECTRUM_H_

#include "bunnykiller.h"

#include <array>
#include <cmath>

#include "Color/Spectrum.h"
#include "Color/PolarizationFrame.h"

/**
 * Implemented following [1]
 *
 * [1] A. Wilkie, A. Weidlich 2012
 * 	   Polarised light in computer graphics
 *     SIGGRAPH Asia 2012 Courses
 */
template<unsigned D, class Radiance>
class tPolarizedLight
{
public:
	static constexpr unsigned spectral_samples = Radiance::spectral_samples;
	static constexpr unsigned components = 4;
protected:
	std::array<Radiance, 4> m_stokes_vector;
	PolarizationFrame<D> m_frame;
	bool m_polarized;
public:
	// ------------------------------------------------------------------------------------
	// Constructors
	constexpr tPolarizedLight() :
		m_stokes_vector({ 0., 0., 0., 0. }),
		m_frame(),
		m_polarized(false)
	{}

	template<class T>
	constexpr tPolarizedLight(const T& s0) :
		m_stokes_vector({ Radiance(s0), 0., 0., 0. }),
		m_frame(),
		m_polarized(false)
	{}

	template<class T>
	constexpr tPolarizedLight(const T& s0, const PolarizationFrame<D>& frame) :
		m_stokes_vector({ Radiance(s0), 0., 0., 0. }),
		m_frame(frame),
		m_polarized(true)
	{}

	template<class T0, class T1>
	constexpr tPolarizedLight(const T0& s0, const T1& s1,
			const PolarizationFrame<D>& frame) :
		m_stokes_vector({ Radiance(s0), Radiance(s1), 0., 0. }),
		m_frame(frame),
		m_polarized(true)
	{}

	template<class T0, class T1, class T2, class T3>
	constexpr tPolarizedLight(const T0& s0, const T1& s1, const T2& s2, const T3& s3,
			const PolarizationFrame<D>& frame) :
		m_stokes_vector({ Radiance(s0), Radiance(s1), Radiance(s2), Radiance(s3) }),
		m_frame(frame),
		m_polarized(true)
	{}

	tPolarizedLight(const tPolarizedLight<D, Radiance>& r) :
		m_stokes_vector(r.m_stokes_vector),
		m_frame(r.m_frame),
		m_polarized(r.m_polarized)
	{}
	
	// ------------------------------------------------------------------------------------
	// Assign operators
	template<class T>
	const tPolarizedLight<D, Radiance> &operator=(const T& f)
	{
		m_stokes_vector = { Radiance(f), 0., 0., 0. };
		m_frame = PolarizationFrame<D>();
		m_polarized = false;

		return *this;
	}

	const tPolarizedLight<D, Radiance>& operator=(const tPolarizedLight<D, Radiance>& r)
	{
		m_stokes_vector = r.m_stokes_vector;
		m_frame = r.m_frame;
		m_polarized = r.m_polarized;

		return *this;
	}

	// ------------------------------------------------------------------------------------
	// Access functions
	Real get_max() const
	{
		return m_stokes_vector[0].get_max();
	}

	Real get_min() const
	{
		return m_stokes_vector[0].get_min();
	}
	
	void get_lambda(Real *lambda) const
	{

		m_stokes_vector[0].get_lambda(lambda);
	}

	void extract_values(const std::vector<Real> &lambda, const Real *values)
	{
		m_stokes_vector[0].extract_values(lambda, values);
	}

	// ------------------------------------------------------------------------------------
	// Access operations
	const Radiance& operator[](const int i) const
	{
		return m_stokes_vector[i];
	}

	Radiance operator[](const int i)
	{
		return m_stokes_vector[i];
	}

	inline bool polarized() const
	{
		return m_polarized;
	}

	// ------------------------------------------------------------------------------------
	// Comparison 
	bool operator==(const tPolarizedLight<D, Radiance> &l) const
	{
		if (m_polarized != l.m_polarized)
			return false;

		for (size_t i = 0; i < 4; ++i)
			if (!(m_stokes_vector[i] == l.m_stokes_vector[i]))
				return false;

		return true;
	}

	bool operator==(const Radiance &v) const
	{
		if (m_polarized)
			return false;

		return m_stokes_vector[0] == v;
	}

	// ------------------------------------------------------------------------------------
	// Unary Operations
	Real sum() const
	{
		return m_stokes_vector[0].sum();
	}

	Real avg() const
	{
		return m_stokes_vector[0].avg();
	}

	Real maxChannel() const
	{
		return m_stokes_vector[0].maxChannel();
	}

	bool is_zero() const
	{
		return m_stokes_vector[0].is_zero();
	}

	bool is_negative() const
	{
		return m_stokes_vector[0].is_negative();
	}

	bool is_valid() const
	{
		return m_stokes_vector[0].is_valid() && m_stokes_vector[1].is_valid()
				&& m_stokes_vector[2].is_valid() && m_stokes_vector[3].is_valid();
	}

	bool is_finite() const
	{
		return m_stokes_vector[0].is_finite() && m_stokes_vector[1].is_finite()
				&& m_stokes_vector[2].is_finite() && m_stokes_vector[3].is_finite();
	}

	Real length() const
	{
		return m_stokes_vector[0].length();
	}

	Real length2() const
	{
		return m_stokes_vector[0].length();
	}

	tPolarizedLight<D, Radiance> exp() const
	{
		tPolarizedLight<D, Radiance> nlight(*this);

		nlight.m_stokes_vector[0] = nlight.m_stokes_vector[0].exp();
		nlight.m_stokes_vector[1] = nlight.m_stokes_vector[1].exp();
		nlight.m_stokes_vector[2] = nlight.m_stokes_vector[2].exp();
		nlight.m_stokes_vector[3] = nlight.m_stokes_vector[3].exp();

		return nlight;
	}

	tPolarizedLight<D, Radiance> operator-() const
	{
		tPolarizedLight<D, Radiance> nlight(*this);

		nlight.m_stokes_vector[0] = -nlight.m_stokes_vector[0];
		nlight.m_stokes_vector[1] = -nlight.m_stokes_vector[1];
		nlight.m_stokes_vector[2] = -nlight.m_stokes_vector[2];
		nlight.m_stokes_vector[3] = -nlight.m_stokes_vector[3];

		return nlight;
	}

	// ------------------------------------------------------------------------------------
	// Operations with scalars or Spectrum
	template<class T>
	tPolarizedLight<D, Radiance> operator+(const T& f) const
	{
		tPolarizedLight<D, Radiance> nlight(*this);
		nlight.m_stokes_vector[0] += f;

		return nlight;
	}
	
	template<class T>
	const tPolarizedLight<D, Radiance> &operator+=(const T& f) const
	{
		m_stokes_vector[0] += f;

		return *this;
	}

	template<class T>
	tPolarizedLight<D, Radiance> operator-(const T& f) const
	{
		tPolarizedLight<D, Radiance> nlight(*this);
		nlight.m_stokes_vector[0] -= f;

		return nlight;
	}

	template<class T>
	const tPolarizedLight<D, Radiance>& operator-=(const T& f) const
	{
		m_stokes_vector[0] -= f;

		return *this;
	}

	template<class T>
	tPolarizedLight<D, Radiance> operator*(const T& f) const
	{
		if (!m_polarized)
			return tPolarizedLight<D, Radiance>(m_stokes_vector[0] * f);

		tPolarizedLight<D, Radiance> nlight(*this);
		nlight.m_stokes_vector[0] = nlight.m_stokes_vector[0] * f;
		nlight.m_stokes_vector[1] = nlight.m_stokes_vector[1] * f;
		nlight.m_stokes_vector[2] = nlight.m_stokes_vector[2] * f;
		nlight.m_stokes_vector[3] = nlight.m_stokes_vector[3] * f;
		
		return nlight;
	}

	template<class T>
	const tPolarizedLight<D, Radiance>& operator*=(const T& f)
	{
		m_stokes_vector[0] *= f;

		if (m_polarized) {
			m_stokes_vector[1] *= f;
			m_stokes_vector[2] *= f;
			m_stokes_vector[3] *= f;
		}

		return *this;
	}

	template<class T>
	tPolarizedLight<D, Radiance> operator/(const T& f) const
	{
		T inv_f = T(1.) / f;
		
		if (!m_polarized)
			return tPolarizedLight<D, Radiance>(m_stokes_vector[0] * inv_f);

		tPolarizedLight<D, Radiance> nlight(*this);
		
		nlight.m_stokes_vector[0] *= inv_f;
		nlight.m_stokes_vector[1] *= inv_f;
		nlight.m_stokes_vector[2] *= inv_f;
		nlight.m_stokes_vector[3] *= inv_f;

		return nlight;
	}

	template<class T>
	const tPolarizedLight<D, Radiance>& operator/=(const T& f)
	{
		T inv_f = T(1.) / f;
		m_stokes_vector[0] *= inv_f;

		if (m_polarized) {
			m_stokes_vector[1] *= inv_f;
			m_stokes_vector[2] *= inv_f;
			m_stokes_vector[3] *= inv_f;
		}

		return *this;
	}

	// ------------------------------------------------------------------------------------
	// Operations with Polarized light
	tPolarizedLight<D, Radiance> operator+(const tPolarizedLight<D, Radiance>& r) const
	{
		if (!r.m_polarized) {
			tPolarizedLight<D, Radiance> rr(*this);
			rr.m_stokes_vector[0] += r.m_stokes_vector[0];

			return rr;
		} else if (!m_polarized) {
			tPolarizedLight<D, Radiance> rr(r);
			rr.m_stokes_vector[0] += m_stokes_vector[0];

			return rr;
		} else {
			tPolarizedLight<D, Radiance> rr = r.aligned(m_frame);

			for (size_t i = 0; i < 4; ++i)
				rr.m_stokes_vector[i] += m_stokes_vector[i];

			return rr;
		}
	}
	
	const tPolarizedLight<D, Radiance>& operator+=(const tPolarizedLight<D, Radiance>& r)
	{
		if (!r.m_polarized) {
			m_stokes_vector[0] += r.m_stokes_vector[0];
		} else if (!m_polarized) {
			m_frame = r.m_frame;
			m_polarized = true;

			for (size_t i = 0; i < 4; ++i)
				m_stokes_vector[i] += r.m_stokes_vector[i];
		} else {
			tPolarizedLight<D, Radiance> rr = r.aligned(m_frame);

			for (unsigned int i = 0; i < 4; ++i)
				m_stokes_vector[i] += rr.m_stokes_vector[i];
		}

		return *this;
	}

	tPolarizedLight<D, Radiance> operator-(const tPolarizedLight<D, Radiance>& r) const
	{
		if (!r.m_polarized) {
			tPolarizedLight<D, Radiance> rr(*this);
			rr.m_stokes_vector[0] -= r.m_stokes_vector[0];

			return rr;
		} else if (!m_polarized) {
			tPolarizedLight<D, Radiance> rr(r);
			rr.m_stokes_vector[0] -= m_stokes_vector[0];

			return rr;
		} else {
			tPolarizedLight<D, Radiance> rr = r.aligned(m_frame);

			for (size_t i = 0; i < 4; ++i)
				rr.m_stokes_vector[i] -= m_stokes_vector[i];

			return rr;
		}
	}

	const tPolarizedLight<D, Radiance>& operator-=(const tPolarizedLight<D, Radiance>& r)
	{
		if (!r.m_polarized) {
			m_stokes_vector[0] -= r.m_stokes_vector[0];
		} else if (!m_polarized) {
			m_frame = r.m_frame;
			m_polarized = true;

			for (size_t i = 0; i < 4; ++i)
				m_stokes_vector[i] -= r.m_stokes_vector[i];
		} else {
			tPolarizedLight<D, Radiance> rr = r.aligned(m_frame);

			for (size_t i = 0; i < 4; ++i)
				m_stokes_vector[i] -= rr.m_stokes_vector[i];
		}

		return *this;
	}

	tPolarizedLight<D, Radiance> align(const PolarizationFrame<D>& f)
	{
		if (!m_polarized) {
			m_frame = f;
			return *this;
		}

		Real cos_alpha, sin_alpha;
		f.get_cos_sin(m_frame, cos_alpha, sin_alpha);

		Real cos_2alpha = 2. * cos_alpha * cos_alpha - 1.;
		Real sin_2alpha = 2. * cos_alpha * sin_alpha;

		m_stokes_vector[1] = cos_2alpha * m_stokes_vector[1] + sin_2alpha * m_stokes_vector[2];
		m_stokes_vector[2] = cos_2alpha * m_stokes_vector[2] - sin_2alpha * m_stokes_vector[1];
		m_frame = f;

		return *this;
	}

	tPolarizedLight<D, Radiance> aligned(const PolarizationFrame<D>& f) const
	{
		if (!m_polarized) {
			tPolarizedLight<D, Radiance> v(*this);
			v.m_frame = f;

			return v;
		} else {
			Real cos_alpha, sin_alpha;
			f.get_cos_sin(m_frame, cos_alpha, sin_alpha);

			Real cos_2alpha = 2. * cos_alpha * cos_alpha - 1.;
			Real sin_2alpha = 2. * cos_alpha * sin_alpha;

			tPolarizedLight<D, Radiance> v(
				m_stokes_vector[0],
				cos_2alpha * m_stokes_vector[1] + sin_2alpha * m_stokes_vector[2],
				cos_2alpha * m_stokes_vector[2] - sin_2alpha * m_stokes_vector[1],
				m_stokes_vector[3],
				f
			);

			return v;
		}
	}

	inline const PolarizationFrame<D>& get_frame() const
	{
		return m_frame;
	}

	inline void set_frame(const PolarizationFrame<D>& frame)
	{
		m_frame = frame;
	}

	void print(FILE* stream) const
	{
		Vector3 I = m_stokes_vector[0].to_rgb();
		Vector3 Q = m_stokes_vector[1].to_rgb();
		Vector3 U = m_stokes_vector[2].to_rgb();
		Vector3 V = m_stokes_vector[3].to_rgb();

		fprintf(stream, "- I: [%f, %f, %f]\n", I[0], I[1], I[2]);
		fprintf(stream, "- Q: [%f, %f, %f]\n", Q[0], Q[1], Q[2]);
		fprintf(stream, "- U: [%f, %f, %f]\n", U[0], U[1], U[2]);
		fprintf(stream, "- V: [%f, %f, %f]\n", V[0], V[1], V[2]);

		m_frame.print(stream);
	}

	Vector3 to_rgb() const
	{
		return m_stokes_vector[0].to_rgb();
	}
};
// tPolarizedLight

template<unsigned D, class Radiance>
tPolarizedLight<D, Radiance> exp(const tPolarizedLight<D, Radiance> &v)
{
	return v.exp();
}

template<unsigned D, class Radiance>
tPolarizedLight<D, Radiance> operator*(const Real &s, const tPolarizedLight<D, Radiance> &f)
{
	return f * s;
}

template<unsigned D>
using PolarizedLight = tPolarizedLight<D, Spectrum>;

/*************************************************************************************/
/* Orientation functions */
template<class Radiance, unsigned D>
inline void align_to_frame(const PolarizationFrame<D>& f, tPolarizedLight<D, Radiance>& light)
{
	PolarizationFrame<3> frame(f[0], light.get_frame()[1]);
	light = light.aligned(frame);
}

template<class Radiance, unsigned D>
inline void set_direction(const VectorN<D>& dir, tPolarizedLight<D, Radiance>& light)
{
	PolarizationFrame<D> new_frame = transform_frame(light.get_frame(), dir);
	light.set_frame(new_frame);
}

#endif //_POLARIZED_SPECTRUM_H_
