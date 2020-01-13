/*
 * Copyright (C) 2017, Victor Arellano (http://giga.cps.unizar.es/~varella/)
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

#ifndef _FLUORESCENT_ATTENUATION_H_
#define _FLUORESCENT_ATTENUATION_H_

#include "bunnykiller.h"

#include <cstdio>
#include <cstring>

#include "LinearAlgebra/VectorN.h"
#include "RayTracing/RayTraceDirection.h"
#include "Color/PolarizedSpectrum.h"
#include "Color/Spectrum.h"
#include "Color/FluorescentMatrix.h"
#include "Color/PolarizationFrame.h"

#ifndef _RAYTRACE_DIRECTION_H_
enum TraceDirection;
#endif

template<unsigned D, class Radiance>
class tFluorescentAttenuation
{
public:
	static constexpr unsigned spectral_samples = Radiance::spectral_samples;
	static constexpr unsigned components = Radiance::components;
protected:
	tFluorescentMatrix<Radiance> m_matrix;
public:
	tFluorescentAttenuation() {}

	tFluorescentAttenuation(tFluorescentMatrix<Radiance> mat) :
		m_matrix(mat)
	{}

	tFluorescentAttenuation(const Radiance &intensity) :
		m_matrix(intensity)
	{}

	tFluorescentAttenuation(const Real &intensity) :
		m_matrix(intensity)
	{}

	tFluorescentAttenuation(const Real inMatrix[]) :
		m_matrix(inMatrix)
	{}

	tFluorescentAttenuation(const tFluorescentAttenuation<D, Radiance> &f) :
		m_matrix(f.m_matrix)
	{}

	// Set operators
	void set(const Spectrum &intensity)
	{
		m_matrix.set(intensity);
	}

	void set(const Real inMatrix[])
	{
		m_matrix.set(inMatrix);
	}

	// Access operators and functions
	inline const Spectrum& operator()(unsigned int i, unsigned int j) const
	{
		return m_matrix(i, j);
	}

	inline bool isFluorescent() const
	{
		return m_matrix.isFluorescent();
	}

	// ------------------------------------------------------------------------------------
	// Operations with Spectrum or Scalars
	tFluorescentAttenuation<D, Radiance> operator*(const Real &f) const
	{
		tFluorescentAttenuation<D, Radiance> natt(*this);
		natt.m_matrix *= f;

		return natt;
	}

	const tFluorescentAttenuation<D, Radiance> &operator*=(const Real &f)
	{
		m_matrix *= f;

		return *this;
	}

	tFluorescentAttenuation<D, Radiance> operator/(const Real &f) const
	{
		tFluorescentAttenuation<D, Radiance> natt(*this);
		natt.m_matrix /= f;

		return natt;
	}

	const tFluorescentAttenuation<D, Radiance> &operator/=(const Real &f)
	{
		m_matrix /= f;

		return *this;
	}

	// ------------------------------------------------------------------------------------
	// Operations with spectrum
	Radiance operator*(const Radiance &f)const
	{
		return m_matrix * f;
	}

	tFluorescentAttenuation<D, Radiance> &operator*=(const Radiance &f)
	{
		m_matrix * f;

		return (*this);
	}

	// ------------------------------------------------------------------------------------
	// Operations with itself: 
	tFluorescentAttenuation<D, Radiance> operator*(const tFluorescentAttenuation<D, Radiance> &f) const
	{
		tFluorescentMatrix<Radiance> mat = m_matrix * f.m_matrix;

		return tFluorescentAttenuation<D, Radiance>(mat);
	}

	tFluorescentAttenuation<D, Radiance> operator/(const tFluorescentAttenuation<D, Radiance> &f) const
	{
		throw("Invalid operation"); //YES
	}

	// ------------------------------------------------------------------------------------
	// Assign operators
	tFluorescentAttenuation<D, Radiance> &operator=(const Radiance &intensity)
	{
		m_matrix = intensity;

		return (*this);
	}

	const tFluorescentAttenuation<D, Radiance> &operator=(const Real f)
	{
		m_matrix = f;

		return *this;
	}

	tFluorescentAttenuation<D, Radiance> &operator=(const tFluorescentAttenuation<D, Radiance> &f)
	{
		m_matrix = f.m_matrix;

		return *this;
	}

	void print(FILE* stream) const
	{
		fprintf(stream, "Fluorescence Added: %s\n", m_matrix.isFluorescent() ? "true" : "false");
		m_matrix.print(stream);
	}
}; // tFluorescentAttenuation

template<unsigned D>
using FluorescentAttenuation = tFluorescentAttenuation<D, Spectrum>;

#endif //_FLUORESCENT_ATTENUATION_H_
