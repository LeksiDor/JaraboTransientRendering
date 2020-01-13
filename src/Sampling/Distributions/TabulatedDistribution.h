/*
 * TabulatedDistribution.h
 *
 *  Created on: Jul 10, 2017
 *      Author: ibon
 */

#ifndef _TABULATED_DISTRIBUTION_H_
#define _TABULATED_DISTRIBUTION_H_

#include "bunnykiller.h"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <memory>
#include <utility>
#include <vector>

/** Stores a tabulated distribution with fixed intervals. */
class TabulatedDistribution {
protected:
	struct DistributionEntry {
		Real pdf;
		Real CDF;

		DistributionEntry() :
			pdf(0.), CDF(0.)
		{}

		DistributionEntry(Real _pdf, Real _CDF) :
			pdf(_pdf), CDF(_CDF)
		{}
	};
protected:
	DistributionEntry* m_data;
	size_t m_N; // Number of steps
	Real m_x0; // Interval start
	Real m_delta, m_inv_delta; // Interval step
public:
	TabulatedDistribution() :
		m_data(nullptr), m_N(0), m_x0(0.), m_delta(0.), m_inv_delta(0.)
	{}

	virtual ~TabulatedDistribution()
	{
		if (m_data)
			delete[] m_data;
	}

	TabulatedDistribution& operator=(TabulatedDistribution&& ld)
	{
		m_data = std::move(ld.m_data);
		ld.m_data = nullptr;

		m_N = ld.m_N;
		m_delta = ld.m_delta;
		m_x0 = ld.m_x0;
		m_inv_delta = ld.m_inv_delta;

		return *this;
	}

	TabulatedDistribution(TabulatedDistribution&& ld)
	{
		*this=ld;
	}

	TabulatedDistribution& operator=(TabulatedDistribution& ld)
	{
		m_data = std::move(ld.m_data);
		ld.m_data = nullptr;

		m_N = ld.m_N;
		m_delta = ld.m_delta;
		m_x0 = ld.m_x0;
		m_inv_delta = ld.m_inv_delta;

		return *this;
	}

	TabulatedDistribution(TabulatedDistribution& ld)
	{
		*this=ld;
	}
public:
	virtual inline Real pdf(Real x) const = 0;

	virtual inline Real CDF(Real x) const = 0;

	virtual inline Real sample(Real xi, Real& pdf) const = 0;
};

#endif /* _TABULATED_DISTRIBUTION_H_ */
