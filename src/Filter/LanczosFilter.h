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

#ifndef _LANCZOS_FILTER_H_
#define _LANCZOS_FILTER_H_

#include "bunnykiller.h"

#include <cmath>

#include "Filter.h"

class LanczosFilter : public Filter
{
protected:
	Real m_tau;

	inline Real Sinc1D(Real x) const
	{
		x = std::abs(x);
		if (x < 1e-5)
			return 1.;
		if (x > 1.)
			return 0.;

		x *= M_PI;
		// Get sinc
		Real sinc = std::sin(x * m_tau) / (x * m_tau);
		// Get lanczos windowing
		Real lanczos = std::sin(x) / x;
		return sinc * lanczos;
	}
public:
	LanczosFilter(Real tau) :
			Filter(),
			m_tau(tau)
	{
	}

	LanczosFilter(int size, Real tau) :
			Filter(size),
			m_tau(tau)
	{
	}

	~LanczosFilter()
	{
	}

	virtual Real evaluate(Real x) const
	{
		return Sinc1D(x * inv_size);
	}

	virtual Real evaluate(Real x, Real y) const
	{
		return Sinc1D(x * inv_size) * Sinc1D(y * inv_size);
	}

	virtual Real evaluate(Real x, Real y, Real z) const
	{
		return Sinc1D(x * inv_size) * Sinc1D(y * inv_size) * Sinc1D(z * inv_size);
	}

	virtual Real evaluate(Real x, Real y, Real z, Real w) const
	{
		return Sinc1D(x * inv_size) * Sinc1D(y * inv_size) * Sinc1D(z * inv_size)
				* Sinc1D(w * inv_size);
	}

	virtual Real evaluate(std::vector<Real> &v) const
	{
		Real sol(1.);
		for (size_t i = 0; i < v.size(); i++)
			sol *= Sinc1D(v[i] * inv_size);

		return sol;
	}
	
	virtual Filter *get_subfilter(const unsigned int d) const
	{
		return new LanczosFilter(size, m_tau);
	}
};
// LanczosFilter
#endif // _LANCZOS_FILTER_H_
