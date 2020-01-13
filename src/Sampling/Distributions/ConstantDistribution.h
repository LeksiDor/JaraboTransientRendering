/*
 * ConstantDistribution.h
 *
 *  Created on: Jul 6, 2017
 *      Author: ibon
 */

#ifndef _CONSTANTDISTRIBUTION_H_
#define _CONSTANTDISTRIBUTION_H_

#include "bunnykiller.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <memory>
#include <utility>
#include <vector>

#include "Sampling/Distributions/TabulatedDistribution.h"

/** Stores a picewise constant distribution with fixed intervals. */
class ConstantDistribution : public TabulatedDistribution {
protected:
	typedef TabulatedDistribution TDTR;
public:
	ConstantDistribution() :
		TabulatedDistribution()
	{}

	/**
	 * Creates a normalized constant distribution for sampling
	 * 		f  : values of the function
	 * 		N  : number of values
	 * 		x0 : start of the function domain
	 * 		xN : end of the function domain
	 *		C  : Normalization factor (usually 1.0)
	 */
	ConstantDistribution(const Real* f, size_t N, Real x0 = 0., Real xN = 1.)
	{
		TDTR::m_N = N;
		TDTR::m_x0 = x0;

		const Real xT = (xN - x0);
		TDTR::m_delta = xT / Real(N);
		TDTR::m_inv_delta = Real(N) / xT;

		TDTR::m_data = new DistributionEntry[N+1];

		TDTR::m_data[0] = DistributionEntry(
			f[0], Real(0.)
		);

		/* Sum rectangles */
		for (size_t i = 1; i < N; i++) {
			TDTR::m_data[i] = DistributionEntry(
				f[i], m_data[i-1].CDF + f[i-1] * TDTR::m_delta
			);
		}

		TDTR::m_data[N] = DistributionEntry(
			f[N-1], m_data[N-1].CDF + f[N-1] * TDTR::m_delta
		);

		Real norm = Real(1.)/TDTR::m_data[N].CDF;

		for (size_t i = 0; i < N+1; i++) {
			TDTR::m_data[i] = DistributionEntry(
				m_data[i].pdf*norm, m_data[i].CDF*norm
			);
		}

		/* Sanity check */
		TDTR::m_data[N].CDF = Real(1.);
	}

	virtual ~ConstantDistribution() {}

	ConstantDistribution& operator=(ConstantDistribution&& cd)
	{
		TabulatedDistribution::operator=(cd);
		return *this;
	}

	ConstantDistribution(ConstantDistribution&& cd)
	{
		*this=cd;
	}

	ConstantDistribution& operator=(ConstantDistribution& cd)
	{
		TabulatedDistribution::operator=(cd);
		return *this;
	}

	ConstantDistribution(ConstantDistribution& cd) :
		TabulatedDistribution()
	{
		*this=cd;
	}
public:
	virtual inline Real pdf(Real x) const {
		x = x - TDTR::m_x0;

		const size_t idx = size_t(floor(x * TDTR::m_inv_delta));

		return TDTR::m_data[idx].pdf;
	}

	virtual inline Real CDF(Real x) const {
		x = x - TDTR::m_x0;

		const size_t idx = size_t(floor(x * TDTR::m_inv_delta));

		const Real F0 = TDTR::m_data[idx].CDF;
		const Real F1 = TDTR::m_data[idx+1].CDF;

		const Real x0 = Real(idx)*m_delta;

		/* First order interpolation */
		const Real t = (x - x0)*m_inv_delta;

		return F0*(Real(1.) - t) + F1*t;
	}

	virtual inline Real sample(Real xi, Real& pdf) const {
		const DistributionEntry* entry = std::lower_bound(&TDTR::m_data[0], &TDTR::m_data[m_N+1], xi,
			[](const DistributionEntry& e, const Real& r) { return e.CDF < r; });

		size_t idx = size_t(std::max(ptrdiff_t(0), entry - &TDTR::m_data[0] - 1));

		Real f0 = TDTR::m_data[idx].pdf;

		/* Avoid zero probabilities */
		while (f0 <= std::numeric_limits<Real>::epsilon()) {
			f0 = TDTR::m_data[++idx].pdf;
		}

		const Real F0 = TDTR::m_data[idx].CDF;
		const Real F1 = TDTR::m_data[idx+1].CDF;

		const Real x0 = TDTR::m_x0 + Real(idx)*TDTR::m_delta;

		/* First order CDF, solve linear system */
		const Real t = (xi - F0)/(F1 - F0);

		const Real x = t*TDTR::m_delta;

		pdf = f0;
		return x0 + x;
	}
}; /* ConstantDistribution */

#endif /* _CONSTANTDISTRIBUTION_H_ */
