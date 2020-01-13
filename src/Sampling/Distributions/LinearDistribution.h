/*
 * LinearDistribution.h
 *
 *  Created on: Jul 6, 2017
 *      Author: ibon
 */

#ifndef _LINEAR_DISTRIBUTION_H_
#define _LINEAR_DISTRIBUTION_H_

#include "bunnykiller.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <memory>
#include <utility>
#include <vector>

#include "Sampling/Distributions/TabulatedDistribution.h"

/** Stores a picewise linear distribution with fixed intervals. */
class LinearDistribution : public TabulatedDistribution {
protected:
	typedef TabulatedDistribution TDTR;
public:
	LinearDistribution() :
		TabulatedDistribution()
	{}

	/**
	 * Creates a normalized linear distribution for sampling
	 * 		f  : values of the function
	 * 		N  : number of values
	 * 		x0 : start of the function domain
	 * 		xN : end of the function domain
	 *		C  : Normalization factor (usually 1.0)
	 */
	LinearDistribution(const Real* f, size_t N, Real x0 = 0., Real xN = 1.)
	{
		TDTR::m_N = N;
		TDTR::m_x0 = x0;

		const Real xT = (xN - x0);
		TDTR::m_delta = xT / Real(N-1);
		TDTR::m_inv_delta = Real(N-1) / xT;

		TDTR::m_data = new DistributionEntry[N];
		TDTR::m_data[0] = DistributionEntry(
			f[0], Real(0.)
		);

		for (size_t i = 1; i < N; i++) {
			TDTR::m_data[i] = DistributionEntry(
				f[i], m_data[i-1].CDF + (f[i] + f[i-1])*(TDTR::m_delta)*Real(0.5)
			);
		}

		Real norm= 1./TDTR::m_data[N-1].CDF;

		for (size_t i = 0; i < N; i++) {
			TDTR::m_data[i] = DistributionEntry(
				m_data[i].pdf*norm, m_data[i].CDF*norm
			);
		}

		// Sanity check
		TDTR::m_data[N-1].CDF = 1.;
	}

	virtual ~LinearDistribution() {}

	LinearDistribution& operator=(LinearDistribution&& cd)
	{
		TabulatedDistribution::operator=(cd);
		return *this;
	}

	LinearDistribution(LinearDistribution&& cd)
	{
		*this=cd;
	};

	LinearDistribution& operator=(LinearDistribution& cd)
	{
		TabulatedDistribution::operator=(cd);
		return *this;
	}

	LinearDistribution(LinearDistribution& cd) :
		TabulatedDistribution()
	{
		*this=cd;
	};
public:
	virtual inline Real pdf(Real x) const {
		x = x - TDTR::m_x0;

		const size_t idx = size_t(std::floor(x * TDTR::m_inv_delta));

		const Real f0 = TDTR::m_data[idx].pdf;
		const Real f1 = TDTR::m_data[idx+1].pdf;

		const Real x0 = Real(idx)*TDTR::m_delta;

		// First order interpolation
		const Real t = (x - x0)*TDTR::m_inv_delta;

		return f0*(Real(1.) - t) + f1*t;
	}

	virtual inline Real CDF(Real x) const {
		x = x - TDTR::m_x0;

		const size_t idx = size_t(std::floor(x * TDTR::m_inv_delta));

		const Real f0 = TDTR::m_data[idx].pdf;
		const Real f1 = TDTR::m_data[idx+1].pdf;

		const Real F0 = TDTR::m_data[idx].CDF;

		const Real x0 = Real(idx)*TDTR::m_delta;

		// Second order interpolation
		const Real t = (x - x0);

		return F0 + (f1 - f0)*t*t*TDTR::m_inv_delta*0.5 + f0*t;
	}

	virtual inline Real sample(Real xi, Real& pdf) const {
		const DistributionEntry* entry = std::lower_bound(&TDTR::m_data[0], &TDTR::m_data[m_N-1], xi,
			[](const DistributionEntry& e, const Real& r) { return e.CDF < r; });

		const size_t idx = size_t(std::max(ptrdiff_t(0), entry - &TDTR::m_data[0] - 1));

		const Real f0 = TDTR::m_data[idx].pdf;
		const Real f1 = TDTR::m_data[idx+1].pdf;

		const Real F0 = TDTR::m_data[idx].CDF;

		const Real x0 = TDTR::m_x0 + Real(idx)*TDTR::m_delta;

		// Second order CDF, solve quadratic system
		const Real a = (f1 - f0)*TDTR::m_inv_delta*Real(0.5);
		const Real b = f0;
		const Real c = F0 - xi;

		Real x = (-b + std::sqrt(b*b - 4.*a*c))/(2.*a);

		if (x < Real(0.)) {
			x = (-b - std::sqrt(b*b - 4.*a*c))/(2.*a);
		}

		const Real t = x*TDTR::m_inv_delta;

		pdf = f0*(1. - t) + f1*t;
		return x0 + x;
	}
};

#endif // _LINEAR_DISTRIBUTION_H_
