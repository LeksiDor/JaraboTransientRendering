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

#ifndef _MATRIX_N_H_
#define _MATRIX_N_H_

#include <stdio.h>

#include "bunnykiller.h"
#include "LinearAlgebra/VectorN.h"

template<unsigned N>
class MatrixN
{
public:
	Real m_data[N][N];
public:
	MatrixN()
	{	
		for(size_t i = 0; i < N; i++)
			for(size_t j = 0; j < N; j++)
				m_data[i][j] = 0.0f;
	}

	MatrixN(const MatrixN<N> &m)
	{	
		for(size_t i = 0; i < N; i++)
			for(size_t j = 0; j < N; j++)
				m_data[i][j] = m.m_data[i][j];
	}

	void setColumn(const int index, const VectorN<N>& v)
	{
		#ifdef _SAFE_CHECK_		
			if (index > N-1 || index < 0) throw("Out-of-range");
		#endif	

		for(size_t i = 0; i < N; i++)
			m_data[i][index] = v[i];
	}
	void setColumn(const int index, const VectorN<N-1>& v)
	{
		#ifdef _SAFE_CHECK_		
			if (index > N-1 || index < 0) throw("Out-of-range");
		#endif	

		for(size_t i = 0; i < N-1; i++)
			m_data[i][index] = v[i];

		m_data[N-1][index] = (index == (N-1)) ? 1.0f : 0.0f;
	}

	
	/** Operations */
	const MatrixN<N> &invert();

	MatrixN<N> inverse()const
	{
		MatrixN<N> m(*this);
		m.invert();
		return m;
	}

	const MatrixN<N> &transpose()
	{
		for( int i=1; i<N; ++i )
		for( int j=0; j<i; ++j )
			m_data[i][j] = m_data[j][i];
		
		return (*this);
	}

	MatrixN<N> transposed()const
	{
		MatrixN<N> m(*this);
		m.transpose();
		return m;
	}

	MatrixN<N> operator*(const Real f) const
	{
		MatrixN<N> nm;
		for( int i=0; i<N; ++i )
		for( int j=0; j<N; ++j )
			nm.m_data[i][j] = m_data[i][j]*f;

		return nm;
	}

	VectorN<N> operator*(const VectorN<N>& v) const
	{
		VectorN<N> nvector(0.);
		for( int i=0; i<N; ++i )
		for( int j=0; j<N; ++j )
			nvector[i] += m_data[i][j]*v[j];

		return nvector;
	}

	VectorN<N-1> operator*(const VectorN<N-1>& v) const
	{
		VectorN<N> nvector(0.);
		Real w = 0.;
		for( int i=0; i<N-1; ++i )
		{
			for( int j=0; j<N-1; ++j )
				nvector[i] += m_data[i][j]*v[j];
			
			nvector[i] += m_data[i][N-1];
		}		
		
		for( int j=0; j<N-1; ++j )
			w += m_data[N-1][j];


		return nvector/w;
	}

	MatrixN<N> operator*(const MatrixN<N> &m) const
	{
		MatrixN<N> nm;
		for( int i=0; i<N; ++i )
		for( int j=0; j<N; ++j )
		for( int k=0; k<N; ++k )
			nm.m_data[i][j] += m_data[i][k]*m.m_data[k][j];

		return nm;
	}

	const MatrixN<N> &operator*=(const Real f)
	{
		for( int i=0; i<N; ++i )
		for( int j=0; j<N; ++j )
			m_data[i][j] *= f;

		return *this;
	}

	const MatrixN<N> &operator*=(const MatrixN<N> &m)
	{ 
		MatrixN<N> nm(*this);
		for (size_t i=0; i<N; ++i ) {
			for (size_t j=0; j<N; ++j ) {
				m_data[i][j] = 0.;
				for (size_t k=0; k<N; ++k) {
					m_data[i][j] += nm.m_data[i][k]*m.m_data[k][j];
				}
			}
		}

		return nm;
	}
	
	/** Predefined matrices */
	static MatrixN<N> identity()
	{
		MatrixN<N> m;
		for (size_t i=0; i<N; ++i)
			m.m_data[i][i] = 1;

		return m;
	}

	static MatrixN<N> rotation(Real angle, const VectorN<N-1>& v = VectorN<N-1>(1.));

	static MatrixN<N> scale(const VectorN<N-1>& v)
	{
		MatrixN<N> m=MatrixN<N>::identity();

		for( int i=0; i<N-1; ++i )
			m.m_data[i][i] = v[i];

		return m;
	}

	static MatrixN<N> translate(const VectorN<N-1>& v)
	{
		MatrixN<N> m = MatrixN<N>::identity();

		for( int i=0; i<N-1; ++i )
			m.m_data[i][N-1] = v[i];

		return m;
	}
}; //MatrixN

/**	Specialization for the 'invert' function.
	Note that there's only implementation for N=3 and N=4,
	which are the most common matrices in ray-tracing.	*/
template<unsigned N>
const MatrixN<N> &MatrixN<N>::invert()
{
	throw std::runtime_error("No function \"invert\" for this dimension");
}
template<>
const MatrixN<4> &MatrixN<4>::invert()
{
		Real determinant =	(m_data[0][3] * m_data[1][2] * m_data[2][1] * m_data[3][0] - m_data[0][2] * m_data[1][3] * m_data[2][1] * m_data[3][0] - m_data[0][3] * m_data[1][1] * m_data[2][2] * m_data[3][0] + m_data[0][1] * m_data[1][3] * m_data[2][2] * m_data[3][0]+
							 m_data[0][2] * m_data[1][1] * m_data[2][3] * m_data[3][0] - m_data[0][1] * m_data[1][2] * m_data[2][3] * m_data[3][0] - m_data[0][3] * m_data[1][2] * m_data[2][0] * m_data[3][1] + m_data[0][2] * m_data[1][3] * m_data[2][0] * m_data[3][1]+
							 m_data[0][3] * m_data[1][0] * m_data[2][2] * m_data[3][1] - m_data[0][0] * m_data[1][3] * m_data[2][2] * m_data[3][1] - m_data[0][2] * m_data[1][0] * m_data[2][3] * m_data[3][1] + m_data[0][0] * m_data[1][2] * m_data[2][3] * m_data[3][1]+
							 m_data[0][3] * m_data[1][1] * m_data[2][0] * m_data[3][2] - m_data[0][1] * m_data[1][3] * m_data[2][0] * m_data[3][2] - m_data[0][3] * m_data[1][0] * m_data[2][1] * m_data[3][2] + m_data[0][0] * m_data[1][3] * m_data[2][1] * m_data[3][2]+
							 m_data[0][1] * m_data[1][0] * m_data[2][3] * m_data[3][2] - m_data[0][0] * m_data[1][1] * m_data[2][3] * m_data[3][2] - m_data[0][2] * m_data[1][1] * m_data[2][0] * m_data[3][3] + m_data[0][1] * m_data[1][2] * m_data[2][0] * m_data[3][3]+
							 m_data[0][2] * m_data[1][0] * m_data[2][1] * m_data[3][3] - m_data[0][0] * m_data[1][2] * m_data[2][1] * m_data[3][3] - m_data[0][1] * m_data[1][0] * m_data[2][2] * m_data[3][3] + m_data[0][0] * m_data[1][1] * m_data[2][2] * m_data[3][3]);
		
		Real rcpDeterminant = 1.0 / determinant;
		if (determinant == 0.0 || determinant == -0.0) {
			throw std::runtime_error("Matrix non-invertible");
		}
		
		Real inv[4][4];
		inv[0][0] = m_data[1][2]*m_data[2][3]*m_data[3][1] - m_data[1][3]*m_data[2][2]*m_data[3][1] + m_data[1][3]*m_data[2][1]*m_data[3][2] - m_data[1][1]*m_data[2][3]*m_data[3][2] - m_data[1][2]*m_data[2][1]*m_data[3][3] + m_data[1][1]*m_data[2][2]*m_data[3][3];
		inv[0][1] = m_data[0][3]*m_data[2][2]*m_data[3][1] - m_data[0][2]*m_data[2][3]*m_data[3][1] - m_data[0][3]*m_data[2][1]*m_data[3][2] + m_data[0][1]*m_data[2][3]*m_data[3][2] + m_data[0][2]*m_data[2][1]*m_data[3][3] - m_data[0][1]*m_data[2][2]*m_data[3][3];
		inv[0][2] = m_data[0][2]*m_data[1][3]*m_data[3][1] - m_data[0][3]*m_data[1][2]*m_data[3][1] + m_data[0][3]*m_data[1][1]*m_data[3][2] - m_data[0][1]*m_data[1][3]*m_data[3][2] - m_data[0][2]*m_data[1][1]*m_data[3][3] + m_data[0][1]*m_data[1][2]*m_data[3][3];
		inv[0][3] = m_data[0][3]*m_data[1][2]*m_data[2][1] - m_data[0][2]*m_data[1][3]*m_data[2][1] - m_data[0][3]*m_data[1][1]*m_data[2][2] + m_data[0][1]*m_data[1][3]*m_data[2][2] + m_data[0][2]*m_data[1][1]*m_data[2][3] - m_data[0][1]*m_data[1][2]*m_data[2][3];
		
		inv[1][0] = m_data[1][3]*m_data[2][2]*m_data[3][0] - m_data[1][2]*m_data[2][3]*m_data[3][0] - m_data[1][3]*m_data[2][0]*m_data[3][2] + m_data[1][0]*m_data[2][3]*m_data[3][2] + m_data[1][2]*m_data[2][0]*m_data[3][3] - m_data[1][0]*m_data[2][2]*m_data[3][3];
		inv[1][1] = m_data[0][2]*m_data[2][3]*m_data[3][0] - m_data[0][3]*m_data[2][2]*m_data[3][0] + m_data[0][3]*m_data[2][0]*m_data[3][2] - m_data[0][0]*m_data[2][3]*m_data[3][2] - m_data[0][2]*m_data[2][0]*m_data[3][3] + m_data[0][0]*m_data[2][2]*m_data[3][3];
		inv[1][2] = m_data[0][3]*m_data[1][2]*m_data[3][0] - m_data[0][2]*m_data[1][3]*m_data[3][0] - m_data[0][3]*m_data[1][0]*m_data[3][2] + m_data[0][0]*m_data[1][3]*m_data[3][2] + m_data[0][2]*m_data[1][0]*m_data[3][3] - m_data[0][0]*m_data[1][2]*m_data[3][3];
		inv[1][3] = m_data[0][2]*m_data[1][3]*m_data[2][0] - m_data[0][3]*m_data[1][2]*m_data[2][0] + m_data[0][3]*m_data[1][0]*m_data[2][2] - m_data[0][0]*m_data[1][3]*m_data[2][2] - m_data[0][2]*m_data[1][0]*m_data[2][3] + m_data[0][0]*m_data[1][2]*m_data[2][3];
		
		inv[2][0] = m_data[1][1]*m_data[2][3]*m_data[3][0] - m_data[1][3]*m_data[2][1]*m_data[3][0] + m_data[1][3]*m_data[2][0]*m_data[3][1] - m_data[1][0]*m_data[2][3]*m_data[3][1] - m_data[1][1]*m_data[2][0]*m_data[3][3] + m_data[1][0]*m_data[2][1]*m_data[3][3];
		inv[2][1] = m_data[0][3]*m_data[2][1]*m_data[3][0] - m_data[0][1]*m_data[2][3]*m_data[3][0] - m_data[0][3]*m_data[2][0]*m_data[3][1] + m_data[0][0]*m_data[2][3]*m_data[3][1] + m_data[0][1]*m_data[2][0]*m_data[3][3] - m_data[0][0]*m_data[2][1]*m_data[3][3];
		inv[2][2] = m_data[0][1]*m_data[1][3]*m_data[3][0] - m_data[0][3]*m_data[1][1]*m_data[3][0] + m_data[0][3]*m_data[1][0]*m_data[3][1] - m_data[0][0]*m_data[1][3]*m_data[3][1] - m_data[0][1]*m_data[1][0]*m_data[3][3] + m_data[0][0]*m_data[1][1]*m_data[3][3];
		inv[2][3] = m_data[0][3]*m_data[1][1]*m_data[2][0] - m_data[0][1]*m_data[1][3]*m_data[2][0] - m_data[0][3]*m_data[1][0]*m_data[2][1] + m_data[0][0]*m_data[1][3]*m_data[2][1] + m_data[0][1]*m_data[1][0]*m_data[2][3] - m_data[0][0]*m_data[1][1]*m_data[2][3];
		
		inv[3][0] = m_data[1][2]*m_data[2][1]*m_data[3][0] - m_data[1][1]*m_data[2][2]*m_data[3][0] - m_data[1][2]*m_data[2][0]*m_data[3][1] + m_data[1][0]*m_data[2][2]*m_data[3][1] + m_data[1][1]*m_data[2][0]*m_data[3][2] - m_data[1][0]*m_data[2][1]*m_data[3][2];
		inv[3][1] = m_data[0][1]*m_data[2][2]*m_data[3][0] - m_data[0][2]*m_data[2][1]*m_data[3][0] + m_data[0][2]*m_data[2][0]*m_data[3][1] - m_data[0][0]*m_data[2][2]*m_data[3][1] - m_data[0][1]*m_data[2][0]*m_data[3][2] + m_data[0][0]*m_data[2][1]*m_data[3][2];
		inv[3][2] = m_data[0][2]*m_data[1][1]*m_data[3][0] - m_data[0][1]*m_data[1][2]*m_data[3][0] - m_data[0][2]*m_data[1][0]*m_data[3][1] + m_data[0][0]*m_data[1][2]*m_data[3][1] + m_data[0][1]*m_data[1][0]*m_data[3][2] - m_data[0][0]*m_data[1][1]*m_data[3][2];
		inv[3][3] = m_data[0][1]*m_data[1][2]*m_data[2][0] - m_data[0][2]*m_data[1][1]*m_data[2][0] + m_data[0][2]*m_data[1][0]*m_data[2][1] - m_data[0][0]*m_data[1][2]*m_data[2][1] - m_data[0][1]*m_data[1][0]*m_data[2][2] + m_data[0][0]*m_data[1][1]*m_data[2][2];
		
		for(int i = 0; i < 4; i++)
		for(int j = 0; j < 4; j++)
			m_data[i][j] = (inv[i][j] * rcpDeterminant);
		
		return (*this);
}

template<>
const MatrixN<3> &MatrixN<3>::invert()
{
	MatrixN<3> m(*this);
	Real det = m_data[0][0]*m_data[1][1]*m_data[2][2] +
			   m_data[0][1]*m_data[1][2]*m_data[2][0] +
			   m_data[0][2]*m_data[1][0]*m_data[2][1] -
			   m_data[0][0]*m_data[1][2]*m_data[2][1] -
			   m_data[0][1]*m_data[1][0]*m_data[2][2] -
			   m_data[0][2]*m_data[1][1]*m_data[2][0];

	m_data[0][0] = (m.m_data[1][1]*m.m_data[2][2] - m.m_data[1][2]*m.m_data[2][1])/det;
	m_data[0][1] = (m.m_data[0][2]*m.m_data[2][1] - m.m_data[0][1]*m.m_data[2][2])/det;
	m_data[0][2] = (m.m_data[0][1]*m.m_data[1][2] - m.m_data[0][2]*m.m_data[1][1])/det;
	m_data[1][0] = (m.m_data[1][2]*m.m_data[2][0] - m.m_data[1][0]*m.m_data[2][2])/det;
	m_data[1][1] = (m.m_data[0][0]*m.m_data[2][2] - m.m_data[0][2]*m.m_data[2][0])/det;
	m_data[1][2] = (m.m_data[0][2]*m.m_data[1][0] - m.m_data[0][0]*m.m_data[1][2])/det;
	m_data[2][0] = (m.m_data[1][0]*m.m_data[2][1] - m.m_data[1][1]*m.m_data[2][0])/det;
	m_data[2][1] = (m.m_data[0][1]*m.m_data[2][0] - m.m_data[0][0]*m.m_data[2][1])/det;
	m_data[2][2] = (m.m_data[0][0]*m.m_data[1][1] - m.m_data[0][1]*m.m_data[1][0])/det;

	return (*this);
}

/**	Specialization for the 'rotation' function.
	Note that there's only implementation for N=3 and N=4,
	which are the most common matrices in ray-tracing.	*/
template<unsigned N>
MatrixN<N> rotation(Real angle, const VectorN<N-1>& v)
{
	throw std::runtime_error("No \"rotatiom\" matrix implemented for this dimension");
}

template<>
MatrixN<4> rotation(Real angle, const VectorN<3>& v)
{
	VectorN<3> vn = v.normalized();
	MatrixN<4> m = MatrixN<4>::identity();

	float s = std::sin(angle);
	float c = std::cos(angle);

	m.m_data[0][0] = vn[0] * vn[0] + (1.f - vn[0] * vn[0]) * c;
	m.m_data[0][1] = vn[0] * vn[1] * (1.f - c) - vn[2] * s;
	m.m_data[0][2] = vn[0] * vn[2] * (1.f - c) + vn[1] * s;
	m.m_data[1][0] = vn[0] * vn[1] * (1.f - c) + vn[2] * s;
	m.m_data[1][1] = vn[1] * vn[1] + (1.f - vn[1] * vn[1]) * c;
	m.m_data[1][2] = vn[1] * vn[2] * (1.f - c) - vn[0] * s;
	m.m_data[2][0] = vn[0] * vn[2] * (1.f - c) - vn[1] * s;
	m.m_data[2][1] = vn[1] * vn[2] * (1.f - c) + vn[0] * s;
	m.m_data[2][2] = vn[2] * vn[2] + (1.f - vn[2] * vn[2]) * c;

	return m;
}

/** Here we are rotating in the 3D-axis defined by the plane z, which
	actually does not exist in 2D domain. Thus, vector v is being ignored. */
template<>
MatrixN<3> rotation(Real angle, const VectorN<2>& v)
{
	MatrixN<3> m = MatrixN<3>::identity();

	float s = std::sin(angle);
	float c = std::cos(angle);

	m.m_data[0][0] = c;
	m.m_data[0][1] = -s;
	m.m_data[1][0] = s;
	m.m_data[1][1] = c;

	return m;
}

/**	Typedef for Matrices 4x4 and 3x3	*/
typedef MatrixN<4> Matrix4x4;
typedef MatrixN<3> Matrix3x3;

#endif //_MATRIX_N_H_
