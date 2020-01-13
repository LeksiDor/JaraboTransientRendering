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

#ifndef _VECTORN_H_
#define _VECTORN_H_

#include "bunnykiller.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <initializer_list>

template<unsigned N, class T>
class tVectorN
{
protected:
	std::array<T, N> m_data;
public:
	// Constructors
	tVectorN(const T f = 0.0)
	{
		m_data.fill(f);
	}

	/** The following constructors are very useful in the case
		of Vecotor<2> and Vector<3>, which are actually the most
		common types of vectors in rendering					*/
	tVectorN(const T x, const T y)
	{
		if (N == 2) {
			m_data[0]=x;
			m_data[1]=y;
		} else {
			m_data[0] = x;

			if (N > 1) {
				m_data[1] = y;

				for(size_t i = 2; i < N; ++i) {
					m_data[i] = 0.;
				}
			}
		}
	}

	tVectorN(const T x, const T y, const T z)
	{
		if (N == 3) {
			m_data[0]=x;
			m_data[1]=y;
			m_data[2]=z;
		} else {
			m_data[0] = x;
			if (N > 1) {
				m_data[1] = y;
				if (N > 2) {
					m_data[2] = z;

					for (size_t i=2; i<N; ++i)
						m_data[i] = 0.;
				}
			}
		}
	}

	constexpr tVectorN(const tVectorN<N, T> &v) :
		m_data(v.m_data)
	{}

	// Arithmetic operators
	tVectorN<N, T> operator-() const
	{
		tVectorN<N, T> nspectrum;
		for (size_t i = 0; i < N; ++i) {
			nspectrum[i] = -m_data[i];
		}
		return nspectrum;
	}

	tVectorN<N, T> operator+(const tVectorN<N, T>& v) const
	{
		tVectorN<N, T> nspectrum;
		for (size_t i = 0; i < N; ++i) {
			nspectrum[i] = m_data[i] + v.m_data[i];
		}
		return nspectrum;
	}

	tVectorN<N, T> operator-(const tVectorN<N, T>& v) const
	{
		tVectorN<N, T> nspectrum;
		for (size_t i = 0; i < N; ++i) {
			nspectrum[i] = m_data[i] - v.m_data[i];
		}
		return nspectrum;
	}

	tVectorN<N, T> operator*(const T f) const
	{
		tVectorN<N, T> nspectrum;
		for (size_t i = 0; i < N; ++i) {
			nspectrum[i] = m_data[i]*f;
		}
		return nspectrum;
	}

	tVectorN<N, T> operator*(const tVectorN<N, T>& v) const
	{
		tVectorN<N, T> nspectrum;
		for (size_t i = 0; i < N; ++i) {
			nspectrum[i] = m_data[i]*v.m_data[i];
		}
		return nspectrum;
	}

	tVectorN<N, T> operator/(const T f) const
	{
		tVectorN<N, T> nspectrum;
		T inv_f = 1./f;
		for (size_t i = 0; i < N; ++i) {
			nspectrum[i] = m_data[i]*inv_f;
		}
		return nspectrum;
	}

	tVectorN<N, T> operator/(const tVectorN<N, T>& v) const
	{
		tVectorN<N, T> nspectrum;
		for (size_t i = 0; i < N; ++i) {
			nspectrum[i] = m_data[i] / v.m_data[i];
		}
		return nspectrum;
	}

	const tVectorN<N, T> &operator+=(const tVectorN<N, T>& v)
	{
		for (size_t i = 0; i < N; ++i) {
			m_data[i] += v.m_data[i];
		}
		return *this;
	}

	const tVectorN<N, T> &operator-=(const tVectorN<N, T>& v)
	{
		for (size_t i = 0; i < N; ++i) {
			m_data[i] -= v.m_data[i];
		}
		return *this;
	}

	const tVectorN<N, T> &operator*=(const tVectorN<N, T>& v)
	{
		for (size_t i = 0; i < N; ++i) {
			m_data[i] *= v.m_data[i];
		}
		return *this;
	}

	const tVectorN<N, T> &operator/=(const tVectorN<N, T>& v)
	{
		for (size_t i = 0; i < N; ++i) {
			m_data[i] /= v.m_data[i];
		}
		return *this;
	}

	const tVectorN<N, T> &operator*=(const T f)
	{
		for (size_t i = 0; i < N; ++i) {
			m_data[i] *= f;
		}
		return *this;
	}

	const tVectorN<N, T> &operator/=(const T f)
	{
		T inv_f = 1./f;
		for (size_t i = 0; i < N; ++i) {
			m_data[i] *= inv_f;
		}
		return *this;
	}

	// Assignation operators
	const tVectorN<N, T> &operator=(const tVectorN<N, T>& v)
	{
		m_data = v.m_data;
		return *this;
	}

	const tVectorN<N, T> &operator=(const T f)
	{
		m_data.fill(f);

		return *this;
	}

	// Comparison operators
	bool operator==(const tVectorN<N, T>& v) const
	{
		for (size_t i = 0; i < N; ++i) {
			if (m_data[i] != v[i]) {
				return false;
			}
		}
		return true;
	}

   	bool operator!= (const tVectorN<N, T> & v) const
	{
		for (size_t i = 0; i < N; ++i) {
			if (m_data[i] == v[i]) {
				return false;
			}
		}
		return true;
	}

	// Access functions. Define _SAFE_CHECK_ at 'Utils/globals.h'
	// to check if the dimension is in the range [0,2]
	inline T operator[](const int dimension) const
	{
		#ifdef _SAFE_CHECK_		
			if (dimension > N-1 || dimension < 0) throw("Out-of-range");
		#endif	
		return m_data[dimension];	
	}

	inline T &operator[](const int dimension)
	{
		#ifdef _SAFE_CHECK_		
			if (dimension > N-1 || dimension < 0) throw("Out-of-range");
		#endif	
		return m_data[dimension];	
	}

	/* Numeric Functions */
	inline T total()const
	{
		T value = 0.;
		for (size_t i = 0; i < N; ++i) {
			value += m_data[i];
		}
		return value;
	}

	inline T avg()const
	{
		return total() / T(N);
	}

	inline bool isZero(void) const
	{
		for (size_t i = 0; i < N; ++i) {
			if (m_data[i] != 0.)
				return false;
		}
		return true;
	}

	inline bool has_nans() const
	{
		for (size_t i = 0; i<N; ++i) {
			if (std::isnan(m_data[i]) != 0)
				return true;
		}

		return false;
	}

	/* Geometric functions */
	inline T length() const
	{
		return std::sqrt(length2());
	}

	inline T length2() const
	{
		T value = 0.;
		for (size_t i = 0; i < N; ++i) {
			value += (m_data[i] * m_data[i]);
		}
		return value;
	}

	unsigned int dominant() const
	{
		Real dominant_value = std::abs(m_data[0]);
		int dominant_idx = 0;
		
		for (size_t i = 0; i < N; ++i) {
			if (std::abs(m_data[i]) > dominant_value) {
				dominant_value = std::abs(m_data[i]);
				dominant_idx = i;
			}
		}
		return dominant_idx;
	}

	inline void normalize()
	{
		(*this) /= length();
	}

	inline tVectorN<N, T> normalized(void) const
	{
		return (tVectorN<N, T>(*this) /= length());
	}

	inline T min() const
	{
		T min_v = m_data[0];
		for (size_t i = 1; i < N; ++i) {
			min_v = std::min(min_v, m_data[i]);
		}
		return min_v;
	}

	inline tVectorN<N, T> min(const tVectorN<N, T> &v) const
	{
		tVectorN<N, T> min_v;
		for (size_t i = 0; i < N; ++i) {
			min_v[i] = std::min(m_data[i], v.m_data[i]);
		}
		return min_v;
	}
	
	inline T max() const
	{
		T max_v = m_data[0];
		for (size_t i = 1; i < N; ++i) {
			max_v = std::max(max_v, m_data[i]);
		}
		return max_v;
	}

	inline tVectorN<N, T> max(const tVectorN<N, T> &v) const
	{
		tVectorN<N, T> max_v;
		for (size_t i = 0; i < N; ++i) {
			max_v[i] = std::max(m_data[i], v.m_data[i]);
		}
		return max_v;
	}
	inline T dot(const tVectorN<N, T> &v) const
	{
		T value = 0.;
		for (size_t i = 0; i < N; ++i) {
			value += m_data[i] * v.m_data[i];
		}
		return value;
	}

	inline T dot_clamped(const tVectorN<N, T> &v) const
	{
		return std::max(dot(v), Real(0.));
	}

	inline T dot_abs(const tVectorN<N, T> &v) const
	{
		return std::abs(dot(v));
	}

	inline tVectorN<N, T> cross(const tVectorN<N, T> &v) const;

	//Function that returns the reflected direction, depending on normal
	tVectorN<N, T> reflect(const tVectorN<N, T> &normal)const
	{
		return  ((*this) - 2.*dot(normal)*normal);
	}

	//Function that returns the reflected direction, depending on normal and refractive_index
	tVectorN<N, T> refract(const tVectorN<N, T> &normal, T eta) const
	{
		T cos_i = dot(normal);
		T cos_o = sqrt(std::max(T(0.), T(1.) - eta*eta*(T(1.) - cos_i*cos_i)));

		return ((eta*(*this) + (eta*std::abs(cos_i) - cos_o)*(cos_i<0.f?1.f:-1.f)*normal));
	}

	//Functions that change (direction) vector to other coordinate system
	tVectorN<N, T> transform_matrix_to(const tVectorN<N, T>& from, const tVectorN<N, T>& to) const;

	tVectorN<N, T> transform_matrix_to_z(const tVectorN<N, T>& from, const tVectorN<N, T>& from_up) const;
}; // tVectorN

// /** Template specializations */
// Cross product
template<unsigned N, class T>
tVectorN<N, T> tVectorN<N, T>::cross(const tVectorN<N, T> &v) const
{
	static_assert((N == 3) || (N == 2), "No cross product operation defined for vectors of this dimension.");

	if (N == 2) {
		return tVectorN<N, T>(m_data[0] * v.m_data[1] - m_data[1] * v.m_data[0], 0.);
	}

	if (N == 3) {
		return tVectorN<N, T>((m_data[1] * v.m_data[2] - m_data[2] * v.m_data[1]),
							  (m_data[2] * v.m_data[0] - m_data[0] * v.m_data[2]),
							  (m_data[0] * v.m_data[1] - m_data[1] * v.m_data[0]));
	}

	return tVectorN<N, T>();
}

// transform_matrix_to
template<unsigned N, class T>
tVectorN<N, T> tVectorN<N, T>::transform_matrix_to(const tVectorN<N, T>& from, const tVectorN<N, T>& to) const
{
	static_assert((N == 3) || (N == 2), "No \"transform_matrix_to\" operation defined for vectors of this dimension");

	if (N == 2) {
		T mtx[2][2];
	
		T cos_theta = from.dot(to);
		T sin_theta = from.cross(to)[0];//sqrt(1.-cos_theta*cos_theta);
			
		tVectorN<N, T> sol;
		
		mtx[0][0] = cos_theta;  // cos    -sin
		mtx[0][1] = -sin_theta; 
		mtx[1][0] = sin_theta;  
		mtx[1][1] = cos_theta;  // sin     cos
		
		for (size_t i=0; i<2; i++)
			sol[i] = mtx[i][0]*m_data[0]+mtx[i][1]*m_data[1];

		return sol;
	}

	if (N == 3) {
		T mtx[3][3];
		const T epsilon = 1.e-4;
		tVectorN<N, T> v =from.cross(to);
		T e = from.dot(to);
		T f = (e<0.0)?-e:e;
		if (f > 1.0 - epsilon) {
			tVectorN<N, T> u, x;
			T c1, c2, c3;
			size_t i, j;

			x[0] = (from[0] > 0.0)? from[0] : -from[0];
			x[1] = (from[1] > 0.0)? from[1] : -from[1];
			x[2] = (from[2] > 0.0)? from[2] : -from[2];
		
			if (x[0] < x[1]) {
				if (x[0] < x[2]) {
					x[0] = 1.0; x[1] = x[2] = 0.0;
				} else {
					x[2] = 1.0; x[0] = x[1] = 0.0;
				}
			} else {
				if (x[1] < x[2]) {
					x[1] = 1.0; x[0] = x[2] = 0.0;
				} else {
					x[2] = 1.0; x[0] = x[1] = 0.0;
				}
			}

			u[0] = x[0] - from[0]; u[1] = x[1] - from[1]; u[2] = x[2] - from[2];
			v[0] = x[0] - to[0];   v[1] = x[1] - to[1];   v[2] = x[2] - to[2];

			c1 = 2.0 / u.dot(u);
			c2 = 2.0 / v.dot(v);
			c3 = c1 * c2  * u.dot(v);

			for (i = 0; i < 3; i++) {
				for (j = 0; j < 3; j++)
					mtx[i][j] = -c1*u[i]*u[j] - c2*v[i]*v[j] + c3*v[i]*u[j];
				mtx[i][i] += 1.0;
			}
		} else { /* the most common case, unless "from"="to", or "from"=-"to" */
			/* ...otherwise use this hand optimized version (9 mults less) */
			T h, hvx, hvz, hvxy, hvxz, hvyz;
			/* h = (1.0 - e)/DOT(v, v); old code */
			h = 1.0/(1.0 + e);      /* optimization by Gottfried Chen */
			hvx = h * v[0];
			hvz = h * v[2];
			hvxy = hvx * v[1];
			hvxz = hvx * v[2];
			hvyz = hvz * v[1];
			mtx[0][0] = e + hvx * v[0];
			mtx[0][1] = hvxy - v[2];
			mtx[0][2] = hvxz + v[1];

			mtx[1][0] = hvxy + v[2];
			mtx[1][1] = e + h * v[1] * v[1];
			mtx[1][2] = hvyz - v[0];

			mtx[2][0] = hvxz - v[1];
			mtx[2][1] = hvyz + v[0];
			mtx[2][2] = e + hvz * v[2];
		}

		tVectorN<N, T> sol;
		for (size_t i=0; i<3; i++)
			sol[i] = mtx[i][0]*m_data[0] + mtx[i][1]*m_data[1] + mtx[i][2]*m_data[2];
		sol.normalize();

		return sol;
	}

	return tVectorN<N, T>();
}

// transform_matrix_to_z
template<unsigned N, class T>
tVectorN<N, T> tVectorN<N, T>::transform_matrix_to_z(const tVectorN<N, T>& from, const tVectorN<N, T>& from_up)const
{
	static_assert((N == 3), "No \"transform_matrix_to_z\" operation defined for vectors of this dimension");

	if (N == 3) {
		tVectorN<N, T> f_up = from_up;
		T dotup = from.dot(f_up);
		if (abs(dotup) > .999)
			f_up = tVectorN<N, T>(0.,0.,dotup);

		tVectorN<N, T> right = f_up.cross(from).normalized();
		
		tVectorN<N, T> s1(this->dot(right), this->dot(f_up), this->dot(from));
		s1.normalized();

		return s1;
	}

	return tVectorN<N, T>();
}

/** Non-member operations involving Vectors	*/
template<unsigned N, class T>
inline tVectorN<N, T> normalize(const tVectorN<N, T> &v)
{
	return v.normalized();
}

template<unsigned N, class T>
inline tVectorN<N, T> minimum(const tVectorN<N, T> &v1, const tVectorN<N, T> &v2)
{
	return v1.min(v2);
}

template<unsigned N, class T>
inline tVectorN<N, T> maximum(const tVectorN<N, T> &v1, const tVectorN<N, T> &v2)
{
	return v1.max(v2);
}

template<unsigned N, class T>
inline T dot(const tVectorN<N, T> &v1, const tVectorN<N, T> &v2)
{
	return v1.dot(v2);
}

template<unsigned N, class T>
inline T dot_clamped(const tVectorN<N, T> &v1, const tVectorN<N, T> &v2)
{
	return v1.dot_clamped(v2);
}

template<unsigned N, class T>
inline T dot_abs(const tVectorN<N, T> &v1, const tVectorN<N, T> &v2)
{
	return v1.dot_abs(v2);
}

template<unsigned N, class T>
inline tVectorN<N, T> cross(const tVectorN<N, T> &v1, const tVectorN<N, T> &v2)
{
	return v1.cross(v2);
}

template<unsigned N, class T, class F>
inline tVectorN<N, T> operator+(const F f, const tVectorN<N, T> &v)
{
	return v + f;
}

template<unsigned N, class T, class F>
inline tVectorN<N, T> operator-(const F f, const tVectorN<N, T> &v)
{
	return tVectorN<N, T>(f) - v;
}

template<unsigned N, class T, class F>
inline tVectorN<N, T> operator*(const F f, const tVectorN<N, T> &v)
{
	return v*f;
}

template<unsigned N, class T, class F>
inline tVectorN<N, T> operator/(const F f, const tVectorN<N, T> &v)
{
	return tVectorN<N, T>(f) / v;
}

/* Define Real as default type for easy of use */
template<unsigned N>
using VectorN = tVectorN<N, Real>;

#endif //_VECTORN_H_
