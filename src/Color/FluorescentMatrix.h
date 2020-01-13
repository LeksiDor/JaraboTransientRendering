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

#ifndef _FLUORESCENT_MATRIX_H_
#define _FLUORESCENT_MATRIX_H_

#include "bunnykiller.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdio>

#include "LinearAlgebra/VectorN.h"
#include "RayTracing/RayTraceDirection.h"
#include "Color/PolarizedSpectrum.h"
#include "Color/Spectrum.h"
#include "Color/PolarizationFrame.h"

template<class Radiance>
class tFluorescentMatrix
{
public:
	static constexpr unsigned spectral_samples = Radiance::spectral_samples;
	static constexpr unsigned components = Radiance::components;
protected:
	/* Flurescence matrix stored as an array */
	std::array<Real, spectral_samples*spectral_samples> m_fluorescence;
	bool isFluorescenceAdded = false;
protected:
	void setup(const Real inMatrix[])
	{
		std::copy(&inMatrix[0], &inMatrix[spectral_samples*spectral_samples], m_fluorescence.begin());
		isFluorescenceAdded = true;
	}

	void setup()
	{
		m_fluorescence.fill(0.);
		isFluorescenceAdded = false;
	}
public:
	constexpr tFluorescentMatrix() :
		m_fluorescence(), isFluorescenceAdded(false)
	{}

	tFluorescentMatrix(const Radiance &intensity)
	{
		setup();

		for (size_t i = 0; i < spectral_samples; i++) {
			m_fluorescence[i + spectral_samples*i] = intensity[i];
		}
	}

	tFluorescentMatrix(const Real intensity)
	{
		setup();

		for (size_t i = 0; i < spectral_samples; i++) {
			m_fluorescence[i + spectral_samples*i] = intensity;
		}
	}

	tFluorescentMatrix(const Real inMatrix[])
	{
		setup(inMatrix);
	}

	tFluorescentMatrix(const tFluorescentMatrix<Radiance> &f)
	{
		set(f);
	}

	// ------------------------------------------------------------------------------------
	// Access operators and functions
	tFluorescentMatrix<Radiance> &operator=(const Radiance &intensity)
	{
		if (isFluorescenceAdded)
			setup();

		for (size_t i = 0; i < spectral_samples; i++) {
			m_fluorescence[i + spectral_samples*i] = intensity[i];
		}

		return (*this);
	}

	const tFluorescentMatrix<Radiance> &operator=(const Real f)
	{
		if (isFluorescenceAdded)
			setup();

		for (size_t i = 0; i < spectral_samples; i++) {
			m_fluorescence[i + spectral_samples*i] = f;
		}

		return *this;
	}

	const Real* getFluorescenceMatrix()
	{
		return m_fluorescence;
	}

	void set(const tFluorescentMatrix<Radiance>& f)
	{
		m_fluorescence = f.m_fluorescence;
		isFluorescenceAdded = f.isFluorescenceAdded;
	}

	void set(const Real n)
	{
		m_fluorescence.fill(n);
		isFluorescenceAdded = n != 0.;
	}

	void set(const Spectrum & intensity)
	{
		setup();

		for (size_t i = 0; i < spectral_samples; i++) {
			m_fluorescence[i + spectral_samples*i] = intensity[i];
		}
	}

	void set(const Real inMatrix[])
	{
		setup(inMatrix);
	}

	inline const Spectrum& operator()(unsigned int i, unsigned int j) const
	{
		return m_fluorescence[i*spectral_samples + j];
	}

	inline bool isFluorescent() const
	{
		return isFluorescenceAdded;
	}

	// ------------------------------------------------------------------------------------
	// Operations with Real
	tFluorescentMatrix<Radiance> operator*(const Real &f) const
	{
		tFluorescentMatrix<Radiance> natt(*this);
		for (Real& elem: natt.m_fluorescence) {
			elem *= f;
		}

		return natt;
	}

	const tFluorescentMatrix<Radiance> &operator*=(const Real &f)
	{
		for (Real& elem: m_fluorescence) {
			elem *= f;
		}

		return *this;
	}

	tFluorescentMatrix<Radiance> operator/(const Real &f) const
	{
		Real inv_f = 1. / f;

		tFluorescentMatrix<Radiance> natt(*this);
		for (Real& elem: natt.m_fluorescence) {
			elem *= inv_f;
		}

		return natt;
	}

	const tFluorescentMatrix<Radiance>& operator/=(const Real &f)
	{
		Real inv_f = 1. / f;
		for (Real& elem: m_fluorescence) {
			elem *= inv_f;
		}

		return *this;
	}

	// ------------------------------------------------------------------------------------
	// Operations with spectrum
	Radiance operator*(const Radiance &f)const
	{
		Radiance toRet(0.);

		if (!isFluorescenceAdded) {
			for (size_t i = 0; i < spectral_samples; i++) {
				toRet[i] = f[i] * m_fluorescence[i + spectral_samples*i];
			}

		} else {
			for (size_t i = 0; i < spectral_samples; i++) {
				for (size_t w = 0; w < spectral_samples; w++) {
					toRet[i] += m_fluorescence[i*spectral_samples + w] * f[w];
				}
			}
		}

		return toRet;
	}

	tFluorescentMatrix<Radiance> &operator*=(const Radiance &f)
	{
		if (!isFluorescenceAdded) {
			for (size_t i = 0; i < spectral_samples; i++) {
				m_fluorescence[i + spectral_samples*i] *= f[i];
			}
		} else {
			for (size_t i = 0; i < spectral_samples; ++i) {
				for (size_t w = 0; w < spectral_samples; ++w) {
					m_fluorescence[i*spectral_samples + w] *= f[w];
				}
			}
		}

		return (*this);
	}

	// ------------------------------------------------------------------------------------
	// Operations with itself:
	tFluorescentMatrix<Radiance> operator*(const tFluorescentMatrix<Radiance> &f) const
	{
		if (!isFluorescenceAdded && !f.isFluorescenceAdded) {
			Radiance toRet(0.);
			for (size_t i = 0; i < spectral_samples; i++) {
				toRet[i] = f.m_fluorescence[spectral_samples*i + i] * m_fluorescence[spectral_samples*i + i];
			}

			return tFluorescentMatrix<Radiance>(toRet);
		}

		tFluorescentMatrix<Radiance> res(0.);
		for (size_t i = 0; i < spectral_samples; i++) {
			for (size_t j = 0; j < spectral_samples; j++) {
				for (size_t k = 0; k < spectral_samples; k++) {
					res.m_fluorescence[i*spectral_samples + j] +=
							m_fluorescence[i*spectral_samples + k] * f.m_fluorescence[k*spectral_samples + j];
				}
			}
		}
		res.isFluorescenceAdded = true;

		return res;
	}

	tFluorescentMatrix<Radiance> operator*=(const tFluorescentMatrix<Radiance> &f)
	{
		if (!isFluorescenceAdded && !f.isFluorescenceAdded) {
			for (size_t i = 0; i < spectral_samples; i++) {
				m_fluorescence[spectral_samples*i + i] *= f.m_fluorescence[spectral_samples*i + i];
			}

			return (*this);
		}

		tFluorescentMatrix<Radiance> res(0.);
		for (size_t i = 0; i < spectral_samples; i++) {
			for (size_t j = 0; j < spectral_samples; j++) {
				for (size_t k = 0; k < spectral_samples; k++) {
					res.m_fluorescence[i*spectral_samples + j] +=
							m_fluorescence[i*spectral_samples + k] * f.m_fluorescence[k*spectral_samples + j];
				}
			}
		}
		res.isFluorescenceAdded = true;
		set(res);

		return (*this);
	}

	tFluorescentMatrix<Radiance> &operator+=(const tFluorescentMatrix<Radiance> &f)
	{
		isFluorescenceAdded = isFluorescenceAdded || f.isFluorescenceAdded;
		
		for (size_t i = 0; i < spectral_samples*spectral_samples; i++) {
			m_fluorescence[i] += f.m_fluorescence[i];
		}

		return (*this);
	}

	tFluorescentMatrix<Radiance> operator-(const tFluorescentMatrix<Radiance> &f) const
	{
		bool fluorAdded = false;

		tFluorescentMatrix<Radiance> toRet(0.);
		for (size_t i = 0; i < spectral_samples*spectral_samples; i++) {
			toRet.m_fluorescence[i] = m_fluorescence[i] - f.m_fluorescence[i];

			if (toRet.m_fluorescence[i] != 0 && i % 4 == 0) {
				fluorAdded = true;
			}
		}
		toRet.isFluorescenceAdded = fluorAdded;

		return toRet;
	}

	tFluorescentMatrix<Radiance> operator+(const tFluorescentMatrix<Radiance> &f) const
	{
		bool fluorAdded = false;

		tFluorescentMatrix<Radiance> toRet(0.);
		for (size_t i = 0; i < spectral_samples*spectral_samples; i++) {
			toRet.m_fluorescence[i] = m_fluorescence[i] + f.m_fluorescence[i];

			if (toRet.m_fluorescence[i] != 0 && i % 4 == 0) {
				fluorAdded = true;
			}
		}
		toRet.isFluorescenceAdded = fluorAdded;

		return toRet;
	}

	tFluorescentMatrix<Radiance> operator/(const tFluorescentMatrix<Radiance> &f) const
	{
		static_assert(!std::is_same<decltype(f), decltype(*this)>::value,
				"Division of fluorescent matrices is not a valid operation");
		return tFluorescentMatrix<Radiance>();
	}

	Real avg() const
	{
		Real res = 0.;
		for (size_t i = 0; i < spectral_samples*spectral_samples; i++) {
			res += m_fluorescence[i];
		}
		return res / spectral_samples*spectral_samples;
	}

	void print(FILE* stream) const
	{
		for (size_t x = 0; x < spectral_samples; x++) {
			fprintf(stream, "[");
			for (size_t y = 0; y < spectral_samples; y++) {
				fprintf(stream, "%f ", m_fluorescence[x*spectral_samples + y]);
			}
			fprintf(stream, "]\n");
		}
	}
}; // tFluorescentMatrix

template<class Radiance>
Radiance operator*(const Radiance &s, const tFluorescentMatrix<Radiance> &f)
{
	return f*s;
}

template<class Radiance>
tFluorescentMatrix<Radiance> operator*(const Real &s, const tFluorescentMatrix<Radiance> &f)
{
	return f*s;
}

using FluorescentMatrix = tFluorescentMatrix<Spectrum>;

#endif //_POLARIZED_ATTENUATION_H_
