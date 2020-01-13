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

#ifndef _SPECTRUM_H_
#define _SPECTRUM_H_

#include "bunnykiller.h"

#include <array>
#include <cmath>
#include <cfloat>
#include <cstddef>
#include <cstdio>

#include "Color/Color.h"
#include "RayTracing/RayTraceDirection.h"
#include "Utils/Utils.h"

template<unsigned N, int WL_MIN=380, int WL_MAX=780>
class tSpectrum
{
public:
	static constexpr unsigned spectral_samples = N;
	static constexpr unsigned components = 1;
protected:
	std::array<Real, N> data;
public:
	tSpectrum(const Real p = 0.)
	{
		for (size_t i = 0; i<N; ++i)
			data[i] = p;
	}

	tSpectrum(const tSpectrum<N, WL_MIN, WL_MAX> &s)
	{
		for (size_t i = 0; i<N; ++i)
			data[i] = s.data[i];
	}

	/**	This function transforms general RGB into spectral
		color. It bases on:
		Smits, B. "An RGB-to-spectrum conversion for 
		reflectances". In Journal of Graphics Tools Vol.4(4).
		http://www.cs.utah.edu/~bes/papers/color/		*/	
	//	[NEED TO BE DONE]								
	tSpectrum(const Real r, const Real g, const Real b)
	{
	}

	tSpectrum(const std::vector<Real> &lambda, const Real *values)
	{
		extract_values(lambda, values);
	}
	
	virtual ~tSpectrum() {}

	// Access operations
	inline Real operator[](const unsigned dimension) const
	{
		#ifdef _SAFE_CHECK_		
			if (dimension > N-1 || dimension < 0) throw("Out-of-range");
		#endif	
		return data[dimension];	
	}

	inline Real& operator[](const unsigned dimension)
	{
		#ifdef _SAFE_CHECK_		
			if (dimension > N-1 || dimension < 0) throw("Out-of-range");
		#endif	
		return data[dimension];	
	}

	inline Real get_max() const
	{
		Real value = -1.;
		for (size_t i = 0; i<N; ++i) {
			if (value < data[i])
				value = data[i];
		}
		return value;
	}

	inline Real get_min() const
	{
		Real value = std::numeric_limits<Real>::infinity();
		for (size_t i = 0; i<N; ++i) {
			if (value > data[i])
				value = data[i];
		}
		return value;
	}

	virtual void get_lambda(Real* lambda) const
	{
		Real delta_lambda = static_cast<Real>(WL_MAX - WL_MIN);
		
		if (N == 1)
			delta_lambda = 0.;
		else
			delta_lambda /= static_cast<Real>(N-1);

		Real current_lambda = WL_MIN;
		for (size_t i = 0; i<N; ++i, current_lambda += delta_lambda) {
			lambda[i] = current_lambda;
		}
	}

	virtual void extract_values(const std::vector<Real>& lambda, const Real* values)
	{
		Real delta_lambda = static_cast<Real>(WL_MAX - WL_MIN);
		
		if (N == 1)
			delta_lambda = 0.;
		else
			delta_lambda /= static_cast<Real>(N-1);

		Real current_lambda = WL_MIN;
		for (size_t i = 0; i<N; ++i, current_lambda += delta_lambda) {
			data[i] = Utils::interpolate(lambda, values, current_lambda);
		}

	}
	
	// Comparison 
	inline bool operator==(const tSpectrum<N,WL_MIN,WL_MAX> &v)const
	{
		for (size_t i = 0; i<N; ++i) {
			if (data[i] != v.data[i])
				return false;
		}

		return true;
	}

	inline bool is_negative() const
	{
		for (size_t i = 0; i<N; ++i) {
			if (data[i] < 0.)
				return true;
		}

		return false;
	}

	inline tSpectrum<N,WL_MIN,WL_MAX> max(const tSpectrum<N,WL_MIN,WL_MAX> &v)const
	{
		tSpectrum<N,WL_MIN,WL_MAX> nspectrum;
		for (size_t i = 0; i<N; ++i)
			nspectrum[i] = std::max(data[i], v.data[i]);

		return nspectrum;
	}
	
	inline tSpectrum<N,WL_MIN,WL_MAX> min(const tSpectrum<N,WL_MIN,WL_MAX> &v)const
	{
		tSpectrum<N,WL_MIN,WL_MAX> nspectrum;
		for (size_t i = 0; i<N; ++i)
			nspectrum[i] = std::min(data[i], v.data[i]);

		return nspectrum;
	}

	// Operations
	inline Real sum() const
	{
		Real value = 0.;
		for (size_t i = 0; i<N; ++i)
			value += data[i];

		return value;
	}

	inline Real avg() const
	{
		Real value = 0.;
		for (size_t i = 0; i<N; ++i)
			value += data[i];

		return value / static_cast<Real>(N);
	}

	inline Real maxChannel() const
	{
		Real max = 0.;
		for (size_t i = 0; i<N; ++i) {
			if (max < data[i])
				max = data[i];
		}

		return max;
	}

	inline bool is_zero() const
	{
		for (size_t i = 0; i<N; ++i) {
			if (data[i] != 0.)
				return false;
		}


		return true;
	}

	inline bool is_valid() const
	{
		for (size_t i = 0; i<N; ++i) {
			if (std::isnan(data[i]) != 0)
				return false;
		}

		return true;
	}

	inline bool is_finite() const
	{
		for (size_t i = 0; i<N; ++i) {
			if (std::isfinite(data[i]) == 0)
				return false;
		}

		return true;
	}

	inline Real length() const
	{
		return std::sqrt(length2());
	}

	inline Real length2() const
	{
		Real value = 0.;
		for (size_t i=0; i<N; ++i)
			value += data[i]*data[i];
		
		return value;
	}

	tSpectrum<N,WL_MIN,WL_MAX> exp() const
	{
		tSpectrum<N,WL_MIN,WL_MAX> nspectrum;
		for (size_t i=0; i<N; ++i)
			nspectrum[i] = std::exp(data[i]);

		return nspectrum;
	}

	tSpectrum<N,WL_MIN,WL_MAX> operator-()const
	{
		tSpectrum<N,WL_MIN,WL_MAX> nspectrum;
		for (size_t i=0; i<N; ++i)
			nspectrum[i] = -data[i];

		return nspectrum;
	}

	tSpectrum<N,WL_MIN,WL_MAX> operator+(const tSpectrum<N,WL_MIN,WL_MAX>& v) const
	{
		tSpectrum<N,WL_MIN,WL_MAX> nspectrum;
		for (size_t i=0; i<N; ++i)
			nspectrum[i] = data[i]+v[i];

		return nspectrum;
	}

	tSpectrum<N,WL_MIN,WL_MAX> operator-(const tSpectrum<N,WL_MIN,WL_MAX>& v) const
	{
		tSpectrum<N,WL_MIN,WL_MAX> nspectrum;
		for (size_t i=0; i<N; ++i)
			nspectrum[i] = data[i]-v[i];

		return nspectrum;
	}

	tSpectrum<N,WL_MIN,WL_MAX> operator*(const Real f) const
	{
		tSpectrum<N,WL_MIN,WL_MAX> nspectrum;
		for (size_t i=0; i<N; ++i)
			nspectrum[i] = data[i]*f;

		return nspectrum;
	}

	tSpectrum<N,WL_MIN,WL_MAX> operator*(const tSpectrum<N,WL_MIN,WL_MAX>& v) const
	{
		tSpectrum<N,WL_MIN,WL_MAX> nspectrum;
		for (size_t i=0; i<N; ++i)
			nspectrum[i] = data[i]*v[i];

		return nspectrum;
	}

	tSpectrum<N,WL_MIN,WL_MAX> operator/(const Real f) const
	{
		tSpectrum<N,WL_MIN,WL_MAX> nspectrum;
		Real inv_f = 1./f;
		for (size_t i=0; i<N; ++i)
			nspectrum[i] = data[i]*inv_f;

		return nspectrum;
	}

	tSpectrum<N,WL_MIN,WL_MAX> operator/(const tSpectrum<N,WL_MIN,WL_MAX>& v) const
	{
		tSpectrum<N,WL_MIN,WL_MAX> nspectrum;
		for (size_t i=0; i<N; ++i)
			nspectrum[i] = data[i]/v[i];

		return nspectrum;
	}

	const tSpectrum<N,WL_MIN,WL_MAX> &operator+=(const tSpectrum<N,WL_MIN,WL_MAX>& v)
	{
		for (size_t i=0; i<N; ++i)
			data[i] += v[i];

		return *this;
	}

	const tSpectrum<N,WL_MIN,WL_MAX> &operator-=(const tSpectrum<N,WL_MIN,WL_MAX>& v)
	{
		for (size_t i=0; i<N; ++i)
			data[i]-=v[i];

		return *this;
	}

	const tSpectrum<N,WL_MIN,WL_MAX> &operator*=(const tSpectrum<N,WL_MIN,WL_MAX>& v)
	{
		for (size_t i=0; i<N; ++i)
			data[i] *= v[i];

		return *this;
	}

	const tSpectrum<N,WL_MIN,WL_MAX> &operator/=(const tSpectrum<N,WL_MIN,WL_MAX>& v)
	{
		for (size_t i=0; i<N; ++i)
			data[i]/=v[i];

		return *this;
	}

	/*const tSpectrum<N, WL_MIN, WL_MAX> &operator*=(const FluorescenceMatrix<N>& f)
	{
		if (!f.isFluorescenceAdded) {
			for (int i = 0; i < N; i++) data[i] *= f.fluorescenceMatrix[i + N*i];
		}
		else {
			Spectrum fluorescent(0);
			for (unsigned int i = 0; i < N; ++i)
			{
				for (unsigned int w = 0; w< N; ++w)
				{
					fluorescent.data[i] += f.fluorescenceMatrix[i*SPECTRUM_COMPONENTS + w] * data[w];
				}
			}
			for (int i = 0; i < N; i++) data[i] = fluorescent[i];
		}
		return *this;
	}
	const tSpectrum<N, WL_MIN, WL_MAX> &operator/=(const FluorescenceMatrix<N>& v)
	{
		throw("Invalid operation");
	}*/

	const tSpectrum<N,WL_MIN,WL_MAX> &operator*=(const Real f)
	{
		for (size_t i=0; i<N; ++i)
			data[i]*=f;

		return *this;
	}

	const tSpectrum<N,WL_MIN,WL_MAX> &operator/=(const Real f)
	{
		Real inv_f = 1./f;
		for (size_t i=0; i<N; ++i)
			data[i]*=inv_f;

		return *this;
	}

	//Operators with fluorescent light
	
	// Assignation operators
	const tSpectrum<N,WL_MIN,WL_MAX> &operator=(const tSpectrum<N,WL_MIN,WL_MAX>& v)
	{
		for (size_t i=0; i<N; ++i)
			data[i]=v[i];

		return *this;
	}

	const tSpectrum<N,WL_MIN,WL_MAX> &operator=(const Real f)
	{
		for (size_t i=0; i<N; ++i)
			data[i]=f;

		return *this;
	}

	// Transformation to RGB functions
	Vector3 to_rgb() const;

	void print(FILE* stream) const;
}; //tSpectrum

template<unsigned N, int WL_MIN, int WL_MAX>
inline tSpectrum<N,WL_MIN,WL_MAX> operator+(const Real f, const tSpectrum<N,WL_MIN,WL_MAX> &v)
{
	return v + f;
}

template<unsigned N, int WL_MIN, int WL_MAX>
inline tSpectrum<N,WL_MIN,WL_MAX> operator-(const Real f, const tSpectrum<N,WL_MIN,WL_MAX> &v)
{
	return tSpectrum<N,WL_MIN,WL_MAX>(f) - v;
}

template<unsigned N, int WL_MIN, int WL_MAX>
inline tSpectrum<N,WL_MIN,WL_MAX> operator*(const Real f, const tSpectrum<N,WL_MIN,WL_MAX> &v)
{
	return v*f;
}

template<unsigned N, int WL_MIN, int WL_MAX>
inline tSpectrum<N,WL_MIN,WL_MAX> operator/(const Real f, const tSpectrum<N,WL_MIN,WL_MAX> &v)
{
	return tSpectrum<N,WL_MIN,WL_MAX>(f)/v;
}

template<unsigned N, int WL_MIN, int WL_MAX>
inline tSpectrum<N,WL_MIN,WL_MAX> exp(const tSpectrum<N,WL_MIN,WL_MAX> &v)
{
	return v.exp();
}

template<unsigned N, int WL_MIN, int WL_MAX>
inline tSpectrum<N,WL_MIN,WL_MAX> max(const tSpectrum<N,WL_MIN,WL_MAX> &v1, const tSpectrum<N,WL_MIN,WL_MAX> &v2)
{
	return v1.max(v2);
}

template<unsigned N, int WL_MIN, int WL_MAX>
inline tSpectrum<N,WL_MIN,WL_MAX> min(const tSpectrum<N,WL_MIN,WL_MAX> &v1, const tSpectrum<N,WL_MIN,WL_MAX> &v2)
{
	return v1.min(v2);
}

template<unsigned N, int WL_MIN, int WL_MAX>
void tSpectrum<N, WL_MIN, WL_MAX>::print(FILE* stream)const
{}

//------------------------------------------------------------------------
// To RGB
template<>
Vector3 tSpectrum<81,380,780>::to_rgb() const
{
	Vector3 xyz;
	for (size_t i=0; i<3; ++i)
		xyz += Vector3(Color::WL2XYZ[i][0]*data[i],Color::WL2XYZ[i][1]*data[i],Color::WL2XYZ[i][2]*data[i]); 

	return Color::xyz2rgb(xyz);
}

// For now, need to be done... Oh well...
template<unsigned N, int WL_MIN,int WL_MAX>
Vector3 tSpectrum<N,WL_MIN,WL_MAX>::to_rgb() const
{
	//Vector3 xyz;

	//int step_size = N / (WL_MAX-WL_MIN);
	//int init = WL_MIN + 
	//Real delta = static_cast<Real>(DELTA_SAMPLES_SPECTRUM)/static_cast<Real>(step_size);

	//
	//Real delta_wl = static_cast<Real>(N)/static_cast<Real>(WL_MAX-WL_MIN);
	//Real current = static_cast<Real>(WL_MIN);
	//

	//for( int i=0; i<N; ++i, current += delta_wl )
	//	
	//return Vector3(data[0],data[1],data[2]);
	return Vector3(0.);
}

class Spectrum3: public tSpectrum<3, 380, 780>
{ 
public:
	Spectrum3(const Real p = 0.) : tSpectrum<3,380,780>(0.)
	{
		data[0] =  data[1] = data[2] = p;
	}

	Spectrum3(const tSpectrum<3,380,780> &s)
	{
		data[0] = s[0];
		data[1] = s[1];
		data[2] = s[2];
	}

	// Need to be done, base on: http://www.cs.utah.edu/~bes/papers/color/
	Spectrum3(const Real r, const Real g, const Real b)
	{
		data[0] = r;
		data[1] = g;
		data[2] = b;
	}
	
	Spectrum3(const std::vector<Real>& lambda, const Real* values)
	{
		extract_values(lambda, values);
	}

	virtual ~Spectrum3()
	{}

	Vector3 to_rgb() const
	{
		return Vector3(data[0], data[1], data[2]);
	}

	virtual void get_lambda(Real *lambda) const
	{
		lambda[0] = 650;
		lambda[1] =	520;
		lambda[2] =	465;
	}

	virtual void extract_values(const std::vector<Real>& lambda, const Real* values)
	{
		data[0] = Utils::interpolate(lambda, values, (Real)650.);
		data[1] = Utils::interpolate(lambda, values, (Real)520.);
		data[2] = Utils::interpolate(lambda, values, (Real)465.);
	}

	void print(FILE* stream) const
	{
		fprintf(stream, "- I: [%f, %f, %f]\n", data[0], data[1], data[2]);
	}
}; // Spectrum3

class Spectrum1: public tSpectrum<1, 380, 380>
{ 
public:
	Spectrum1(const Real p=0.) : tSpectrum<1,380,380>(p)
	{}

	Spectrum1(const tSpectrum<1, 380, 380> &s) : tSpectrum<1,380,380>(s[0])
	{}

	/** Transforms RGB to luminances using the values given
		by Stokes et al. 1996. "A Standard Default Color Space 
		for the Internet - sRGB".
		http://www.w3.org/Graphics/Color/sRGB			*/
	Spectrum1(const Real r, const Real g, const Real b)
	{
		data[0] = 0.2126*r + 0.7152*g + 0.0722*b;
	}
	
	virtual ~Spectrum1()
	{}

	Vector3 to_rgb() const
	{
		return Vector3(data[0]);
	}

	void print(FILE* stream) const
	{
		fprintf(stream, "- I: [%f, %f, %f]\n", data[0], data[0], data[0]);
	}
}; // Spectrum1

#if SPECTRUM_CHANNELS == 1
using Spectrum = Spectrum1;
#elif SPECTRUM_CHANNELS == 3
using Spectrum = Spectrum3;
#else
using Spectrum = tSpectrum<SPECTRUM_CHANNELS>;
#endif

/*************************************************************************************/
/* Orientation functions (dummies for orientation-agnostic spectrum classes) */
template<class Radiance, unsigned D = DIM>
inline void set_attenuation_tracing_direction(const TraceDirection dir, Radiance &R)
{

}

template<class Radiance, unsigned D = DIM>
inline void set_direction(const VectorN<D>&, Radiance&)
{

}

template<class Frame, class Radiance, unsigned D = DIM>
inline void align_to_frame(const Frame&, Radiance&)
{

}

template<class Radiance, class RadianceAtt, unsigned D = DIM>
inline void set_coordinate_system(const Radiance&, RadianceAtt&)
{

}

#endif //_SPECTRUM_H_
