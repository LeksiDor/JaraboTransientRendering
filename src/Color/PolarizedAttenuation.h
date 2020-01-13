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

#ifndef _POLARIZED_ATTENUATION_H_
#define _POLARIZED_ATTENUATION_H_

#include "bunnykiller.h"

#include <cmath>
#include <array>
#include <cstring>
#include <cstdio>
#include <type_traits>

#include "LinearAlgebra/VectorN.h"
#include "RayTracing/RayTraceDirection.h"
#include "Color/PolarizedSpectrum.h"
#include "Color/Spectrum.h"
#include "Color/PolarizationFrame.h"

/**
 * Implemented following Mueller calculus according to [1]
 *
 * [1] A. Wilkie, A. Weidlich 2012
 * 	   Polarised light in computer graphics
 *     SIGGRAPH Asia 2012 Courses
 */
template<unsigned D, class RadianceAtt>
class tPolarizedAttenuation
{
public:
	static constexpr unsigned spectral_samples = RadianceAtt::spectral_samples;
	static constexpr unsigned components = 4;
public:
	enum AttenuationType
	{
		UNKNOWN = 0,
		DEPOLARISER = 1,
		POLARISER = 2,
		ATTENUATION = 4,
		NON_POLARISER = DEPOLARISER | ATTENUATION
	};
private:
	using MuellerMatrixRow = std::array<RadianceAtt, 4>;
	using MuellerMatrix = std::array<MuellerMatrixRow, 4>;

	/* Mueller matrix */
	MuellerMatrix m_mueller;
	/* Polarization reference frames */
	PolarizationFrame<D> m_input_frame, m_output_frame;
	/* Polarization state */
	AttenuationType m_type;

	inline void setup(
		const RadianceAtt &s00, const RadianceAtt &s01, const RadianceAtt &s02, const RadianceAtt &s03,
		const RadianceAtt &s10, const RadianceAtt &s11, const RadianceAtt &s12, const RadianceAtt &s13,
		const RadianceAtt &s20, const RadianceAtt &s21, const RadianceAtt &s22, const RadianceAtt &s23,
		const RadianceAtt &s30, const RadianceAtt &s31, const RadianceAtt &s32, const RadianceAtt &s33,
		const PolarizationFrame<D> &input_frame, const PolarizationFrame<D> &output_frame,
		const AttenuationType type)
	{
		m_mueller[0][0] = s00; m_mueller[0][1] = s01; m_mueller[0][2] = s02; m_mueller[0][3] = s03;
		m_mueller[1][0] = s10; m_mueller[1][1] = s11; m_mueller[1][2] = s12; m_mueller[1][3] = s13;
		m_mueller[2][0] = s20; m_mueller[2][1] = s21; m_mueller[2][2] = s22; m_mueller[2][3] = s23;
		m_mueller[3][0] = s30; m_mueller[3][1] = s31; m_mueller[3][2] = s32; m_mueller[3][3] = s33;

		m_input_frame = input_frame;
		m_output_frame = output_frame;

		m_type = type;
	}

	inline void fill_matrix(const RadianceAtt& f)
	{
		m_mueller.fill(MuellerMatrixRow({f, f, f, f}));
	}

	tPolarizedAttenuation<D, RadianceAtt> align(const PolarizationFrame<D>& frame)
	{
		if (!m_type || (m_type & ATTENUATION)) {
			m_output_frame = frame;
			return *this;
		}
		
		Real cos_alpha, sin_alpha;
		frame.get_cos_sin(m_output_frame, cos_alpha, sin_alpha);

		Real C = 2. * cos_alpha * cos_alpha - 1.;
		Real S = 2. * cos_alpha * sin_alpha;

		const MuellerMatrix& M = m_mueller;
		setup(
			M[0][0],               C*M[0][1] - S*M[0][2],                                 C*M[0][2] + S*M[0][1],                                 M[0][3],
			C*M[1][0] - S*M[2][0], C*(C*M[1][1] - S*M[2][1]) - S*(C*M[1][2] - S*M[2][2]), C*(C*M[1][2] - S*M[2][2]) + S*(C*M[1][1] - S*M[2][1]), C*M[1][3] - S*M[2][3],
			C*M[2][0] + S*M[1][0], C*(C*M[2][1] + S*M[1][1]) - S*(C*M[2][2] + S*M[1][2]), C*(C*M[2][2] + S*M[1][2]) + S*(C*M[2][1] + S*M[1][1]), C*M[2][3] + S*M[1][3],
			M[3][0],               C*M[3][1] - S*M[3][2],                                 C*M[3][2] + S*M[3][1],                                 M[3][3],
			m_input_frame, frame, POLARISER
		);

		return *this;
	}
public:
	tPolarizedAttenuation()
	{
		fill_matrix(RadianceAtt(0.));
		m_type = UNKNOWN;
	}

	template<class S00>
	tPolarizedAttenuation(const S00& s00)
	{
		fill_matrix(RadianceAtt(0.));
		m_mueller[0][0] = m_mueller[1][1] = m_mueller[2][2] = m_mueller[3][3] = RadianceAtt(s00);

		m_type = ATTENUATION;
	}

	template<class S00>
	tPolarizedAttenuation(const S00& s00,
			const PolarizationFrame<D>& input_frame,
			const PolarizationFrame<D>& output_frame)
	{
		fill_matrix(RadianceAtt(0.));
		m_mueller[0][0] = RadianceAtt(s00);

		m_input_frame = input_frame;
		m_output_frame = output_frame;

		m_type = DEPOLARISER;
	}

	template <class S00, class S01, class S02, class S03,
	          class S10, class S11, class S12, class S13,
			  class S20, class S21, class S22, class S23,
			  class S30, class S31, class S32, class S33>
	tPolarizedAttenuation(
		const S00 &s00, const S01 &s01, const S02 &s02, const S03 &s03,
		const S10 &s10, const S11 &s11, const S12 &s12, const S13 &s13,
		const S20 &s20, const S21 &s21, const S22 &s22, const S23 &s23,
		const S30 &s30, const S31 &s31, const S32 &s32, const S33 &s33,
		const PolarizationFrame<D>& input_frame, const PolarizationFrame<D>& output_frame)
	{
		using R = RadianceAtt;

		m_mueller[0][0] = R(s00); m_mueller[0][1] = R(s01); m_mueller[0][2] = R(s02); m_mueller[0][3] = R(s03);
		m_mueller[1][0] = R(s10); m_mueller[1][1] = R(s11); m_mueller[1][2] = R(s12); m_mueller[1][3] = R(s13);
		m_mueller[2][0] = R(s20); m_mueller[2][1] = R(s21); m_mueller[2][2] = R(s22); m_mueller[2][3] = R(s23);
		m_mueller[3][0] = R(s30); m_mueller[3][1] = R(s31); m_mueller[3][2] = R(s32); m_mueller[3][3] = R(s33);

		m_input_frame = input_frame;
		m_output_frame = output_frame;

		m_type = POLARISER;
	}

	tPolarizedAttenuation(const tPolarizedAttenuation<D, RadianceAtt> &f)
	{
		m_mueller = f.m_mueller;
		m_input_frame = f.m_input_frame;
		m_output_frame = f.m_output_frame;
		m_type = f.m_type;
	}
	
	/*************************************************************************************/
	/* Assign operators */
	template<class S00>
	tPolarizedAttenuation<D, RadianceAtt> operator=(const S00 &f)
	{
		if (m_type) {
			fill_matrix(RadianceAtt(0.));
		}
		m_mueller[0][0] = m_mueller[1][1] = m_mueller[2][2] = m_mueller[3][3] = RadianceAtt(f);

		m_type = ATTENUATION;

		return *this;
	}

	const tPolarizedAttenuation<D, RadianceAtt>& operator=(const tPolarizedAttenuation<D, RadianceAtt> &f)
	{
		m_mueller = f.m_mueller;
		m_input_frame = f.m_input_frame;
		m_output_frame = f.m_output_frame;
		m_type = f.m_type;

		return *this;
	}

	/*************************************************************************************/
	/* Set functions */
	template<class S00>
	void set(const S00& s00)
	{
		if (m_type) {
			fill_matrix(RadianceAtt(0.));
		}
		m_mueller[0][0] = m_mueller[1][1] = m_mueller[2][2] = m_mueller[3][3] = RadianceAtt(s00);
		
		m_type = ATTENUATION;
	}

	template<class S00>
	void set(const S00 &s00,
			 const PolarizationFrame<D>& input_frame,
			 const PolarizationFrame<D>& output_frame)
	{
		if (m_type) {
			fill_matrix(RadianceAtt(0.));
		}
		m_mueller[0][0] = RadianceAtt(s00);

		m_input_frame = input_frame;
		m_output_frame = output_frame;

		m_type = DEPOLARISER;
	}

	template <class S00, class S01, class S02, class S03,
			  class S10, class S11, class S12, class S13,
			  class S20, class S21, class S22, class S23,
			  class S30, class S31, class S32, class S33>
	void set(const S00 &s00, const S01 &s01, const S02 &s02, const S03 &s03,
			 const S10 &s10, const S11 &s11, const S12 &s12, const S13 &s13,
			 const S20 &s20, const S21 &s21, const S22 &s22, const S23 &s23,
			 const S30 &s30, const S31 &s31, const S32 &s32, const S33 &s33,
		     const PolarizationFrame<D> &input_frame, const PolarizationFrame<D> &output_frame)
	{
		using R = RadianceAtt;

		setup(R(s00), R(s01), R(s02), R(s03),
			  R(s10), R(s11), R(s12), R(s13),
			  R(s20), R(s21), R(s22), R(s23),
			  R(s30), R(s31), R(s32), R(s33),
			  input_frame, output_frame, POLARISER);
	}

	/*************************************************************************************/
	/* Access operators and functions */
	inline const RadianceAtt& operator()(unsigned i, unsigned j) const
	{
		return m_mueller[i][j];
	}

	inline const PolarizationFrame<D> &get_input_frame() const
	{
		return m_input_frame;
	}

	inline void set_input_frame(const PolarizationFrame<D> &frame)
	{
		m_input_frame = frame;
	}

	inline const PolarizationFrame<D> &get_output_frame() const
	{
		return m_output_frame;
	}

	inline void set_output_frame(const PolarizationFrame<D> &frame)
	{
		m_output_frame = frame;
	}

	inline bool polarized() const
	{
		return m_type & POLARISER;
	}

	/*************************************************************************************/
	/* Operations with Scalars or Spectrum */
	template<class T>
	tPolarizedAttenuation<D, RadianceAtt> operator*(const T& f) const
	{
		tPolarizedAttenuation<D, RadianceAtt> natt(*this);

		for (size_t i = 0; i < 4; ++i) {
			MuellerMatrixRow& natt_row = natt.m_mueller[i];
			for (size_t j = 0; j < 4; ++j) {
				natt_row[j] *= f;
			}
		}

		return natt;
	}

	template<class T>
	const tPolarizedAttenuation<D, RadianceAtt> &operator*=(const T& f)
	{
		for (size_t i = 0; i < 4; ++i) {
			MuellerMatrixRow& m_row = m_mueller[i];
			for (size_t j = 0; j < 4; ++j) {
				m_row[j] *= f;
			}
		}

		return *this;
	}

	template<class T>
	tPolarizedAttenuation<D, RadianceAtt> operator/(const T& f) const
	{
		tPolarizedAttenuation<D, RadianceAtt> natt(*this);
		T inv_f = 1. / f;

		for (size_t i = 0; i < 4; ++i) {
			MuellerMatrixRow& natt_row = natt.m_mueller[i];
			for (size_t j = 0; j < 4; ++j) {
				natt_row[j] *= inv_f;
			}
		}

		return natt;
	}

	template<class T>
	const tPolarizedAttenuation<D, RadianceAtt> &operator/=(const T& f)
	{
		T inv_f = 1. / f;

		for (size_t i = 0; i < 4; ++i) {
			MuellerMatrixRow& m_row = m_mueller[i];
			for (size_t j = 0; j < 4; ++j) {
				m_row[j] *= inv_f;
			}
		}

		return *this;
	}

	/*************************************************************************************/
	/* Operations with Polarized Light */
	template<class Radiance>
	tPolarizedLight<D, Radiance> operator*(const tPolarizedLight<D, Radiance> &f)const
	{
		if (m_type & DEPOLARISER) {
			return tPolarizedLight<D, Radiance>(f[0]*m_mueller[0][0]); /* Careful! */
		}

		if (m_type & ATTENUATION) {
			return f*m_mueller[0][0]; /* Careful! */
		}
		
		if (!f.polarized()) {
			return tPolarizedLight<D, Radiance>(f[0]*m_mueller[0][0], f[0]*m_mueller[1][0],
					                            f[0]*m_mueller[2][0], f[0]*m_mueller[3][0],
												m_output_frame);
		}
		
		tPolarizedLight<D, Radiance> rlight = f.aligned(m_input_frame);
		Radiance nlight[4];

		for (size_t i = 0; i < 4; ++i) {
			const MuellerMatrixRow& m_row = m_mueller[i];

			nlight[i] = rlight[0] * m_row[0] + rlight[1] * m_row[1]
			          + rlight[2] * m_row[2] + rlight[3] * m_row[3];
		}

		return tPolarizedLight<D, Radiance>(nlight[0], nlight[1], nlight[2], nlight[3],
				                            m_output_frame);
	}

	/*************************************************************************************/
	/* Operations with itself */
	/* Note that following Mueller calculus not all operations are defined */
	tPolarizedAttenuation<D, RadianceAtt> operator*(
			const tPolarizedAttenuation<D, RadianceAtt> &f) const
	{
		if (m_type & ATTENUATION)
			return f*m_mueller[0][0];
		
		if (f.m_type & ATTENUATION)
			return (*this)*f.m_mueller[0][0];

		/* If multiplying with a non-polarized attenuation... */
		if ((m_type & DEPOLARISER) && (f.m_type & DEPOLARISER)) {
			return tPolarizedAttenuation<D, RadianceAtt>(f.m_mueller[0][0] * m_mueller[0][0],
					                                     f.m_input_frame, m_output_frame);
		}
		
		tPolarizedAttenuation<D, RadianceAtt> natt;
		natt.m_input_frame = f.m_input_frame;
		natt.m_output_frame = m_output_frame;
		natt.m_type = m_type;

		if (f.m_type & DEPOLARISER) {
			const RadianceAtt& f00 = f.m_mueller[0][0];

			for (size_t i = 0; i < 4; ++i) {
				natt.m_mueller[i][0] = m_mueller[i][0] * f00;
			}

			return natt;
		}

		tPolarizedAttenuation<D, RadianceAtt> nrot(f);
		nrot.align(m_input_frame); 

		/* If multiplying with a polarized attenuation... */
		if (m_type & DEPOLARISER) {
			const RadianceAtt& m00 = m_mueller[0][0];

			for (size_t i = 0; i < 4; ++i) {
				natt.m_mueller[0][i] = nrot.m_mueller[0][i] * m00;
			}

			return natt;
		}

		/* Multiply Mueller matrices */
		for (size_t i = 0; i < 4; ++i) {
			const MuellerMatrixRow& m_row = m_mueller[i];
			MuellerMatrixRow& natt_row = natt.m_mueller[i];

			for (size_t j = 0; j < 4; ++j) {
				RadianceAtt res = RadianceAtt(0.);

				for (size_t k = 0; k < 4; ++k) {
					res += m_row[k] * nrot.m_mueller[k][j];
				}

				natt_row[j] = res;
			}
		}

		return natt;
	}
	
	tPolarizedAttenuation<D, RadianceAtt> operator/(const tPolarizedAttenuation<D, RadianceAtt> &f) const
	{
		/* Forbid division by itself */
		static_assert(!std::is_same<decltype(f), decltype(*this)>::value,
				"Division of Mueller matrices is not a valid operation");
		return tPolarizedAttenuation<D, RadianceAtt>();
	}

	/*************************************************************************************/
	/* Reorientation */
	void frames_to_light_direction(const TraceDirection dir)
	{
		if (dir & TraceDirection::FROM_EYE) {
			PolarizationFrame<D> aux_frame(m_output_frame);

			m_output_frame = PolarizationFrame<D>(m_input_frame.m_frame[0], -m_input_frame.m_frame[1]);
			m_input_frame = PolarizationFrame<D>(aux_frame.m_frame[0], -aux_frame.m_frame[1]);
		}
	}

	/*************************************************************************************/
	/* Debug functions */
	void print(FILE* stream) const
	{
		print_type(stream);
		m_input_frame.print(stream);
		m_output_frame.print(stream);
		fprintf(stream, "[%0.10f \t %0.10f \t %0.10f \t %0.10f]\n",
			m_mueller[0][0].avg(), m_mueller[0][1].avg(), m_mueller[0][2].avg(), m_mueller[0][3].avg());
		fprintf(stream, "[%0.10f \t %0.10f \t %0.10f \t %0.10f]\n",
			m_mueller[1][0].avg(), m_mueller[1][1].avg(), m_mueller[1][2].avg(), m_mueller[1][3].avg());
		fprintf(stream, "[%0.10f \t %0.10f \t %0.10f \t %0.10f]\n",
			m_mueller[2][0].avg(), m_mueller[2][1].avg(), m_mueller[2][2].avg(), m_mueller[2][3].avg());
		fprintf(stream, "[%0.10f \t %0.10f \t %0.10f \t %0.10f]\n",
			m_mueller[3][0].avg(), m_mueller[3][1].avg(), m_mueller[3][2].avg(), m_mueller[3][3].avg());
	}

	void printFull(FILE* stream) const
	{
		print_type(stream);
		m_input_frame.print(stream);
		m_output_frame.print(stream);

		for (size_t i = 0; i < 4; ++i) {
			for (size_t j = 0; j < 4; ++j) {
				m_mueller[i][j].print(stream);
			}
		}
	}

	void print_type(FILE* stream) const
	{
		switch (m_type) {
			case DEPOLARISER:
				fprintf(stream, "Interaction - Depolariser\n");
				break;
			case POLARISER:
				fprintf(stream, "Interaction - Polariser\n");
				break;
			case ATTENUATION:
				fprintf(stream, "Interaction - Attenuation\n");
				break;
			default:
				break;
		}
	}
}; /* tPolarizedAttenuation */

template<unsigned D>
using PolarizedAttenuation = tPolarizedAttenuation<D, Spectrum>;

/*************************************************************************************/
/* Orientation functions */
template<class RadianceAtt, unsigned D>
inline void set_attenuation_tracing_direction(const TraceDirection dir, tPolarizedAttenuation<D, RadianceAtt> &R)
{
	R.frames_to_light_direction(dir);
}

template<class Radiance, class RadianceAtt, unsigned D = DIM>
inline void set_coordinate_system(const tPolarizedLight<D, Radiance> &R0, tPolarizedAttenuation<D, RadianceAtt> &R1)
{
	R1.set_input_frame(R0.get_frame());
}

template<class Radiance, class RadianceAtt, unsigned D = DIM>
inline void set_coordinate_system(const tPolarizedAttenuation<D, Radiance> &R0, tPolarizedAttenuation<D, RadianceAtt> &R1)
{
	R1.set_input_frame(R0.get_output_frame());
}

#endif /* _POLARIZED_ATTENUATION_H_ */
