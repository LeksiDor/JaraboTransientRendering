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

#ifndef _GAUSSIAN_FILTER_H_
#define _GAUSSIAN_FILTER_H_

#include "bunnykiller.h"

#include <cmath>

#include "Filter.h"

class GaussianFilter : public Filter
{
	std::vector<Real> _1_sigma_2;

	inline Real evaluate_d(Real x, unsigned int d) const
	{
		if (d < _1_sigma_2.size())
			return std::exp(-x * x * _1_sigma_2[d]);
		else
			return std::exp(-x * x * _1_sigma_2[0]);
	}
public:
	GaussianFilter(Real sigma = 1.) :
			Filter()
	{
		_1_sigma_2.push_back(1. / (sigma * sigma));
	}

	GaussianFilter(Real sigma_x, Real sigma_y) :
			Filter()
	{
		_1_sigma_2.push_back(1. / (sigma_x * sigma_x));
		_1_sigma_2.push_back(1. / (sigma_y * sigma_y));
	}

	GaussianFilter(Real sigma_x, Real sigma_y, Real sigma_z) :
			Filter()
	{
		_1_sigma_2.push_back(1. / (sigma_x * sigma_x));
		_1_sigma_2.push_back(1. / (sigma_y * sigma_y));
		_1_sigma_2.push_back(1. / (sigma_z * sigma_z));
	}

	GaussianFilter(int _size, Real sigma = 1.) :
			Filter(_size)
	{
		_1_sigma_2.push_back(1. / (sigma * sigma));
	}

	GaussianFilter(int _size, Real sigma_x, Real sigma_y) :
			Filter(_size)
	{
		_1_sigma_2.push_back(1. / (sigma_x * sigma_x));
		_1_sigma_2.push_back(1. / (sigma_y * sigma_y));
	}

	GaussianFilter(int _size, Real sigma_x, Real sigma_y, Real sigma_z) :
			Filter(_size)
	{
		_1_sigma_2.push_back(1. / (sigma_x * sigma_x));
		_1_sigma_2.push_back(1. / (sigma_y * sigma_y));
		_1_sigma_2.push_back(1. / (sigma_z * sigma_z));
	}

	virtual ~GaussianFilter()
	{
	}

	Real get_sigma(unsigned int d = 0) const
	{
		return 1. / std::sqrt(_1_sigma_2[(d < _1_sigma_2.size()) ? d : 0]);
	}

	Real evaluate(Real x) const
	{
		return evaluate_d(x, 0);
	}

	Real evaluate(Real x, Real y) const
	{
		return evaluate_d(x, 0) * evaluate_d(y, 1);
	}

	Real evaluate(Real x, Real y, Real z) const
	{
		return evaluate_d(x, 0) * evaluate_d(y, 1) * evaluate_d(z, 2);
	}

	Real evaluate(Real x, Real y, Real z, Real w) const
	{
		return evaluate_d(x, 0) * evaluate_d(y, 1) * evaluate_d(z, 2) * evaluate_d(w, 3);
	}

	Real evaluate(std::vector<Real> &v) const
	{
		Real sol(1.);
		for (unsigned int i = 0; i < v.size(); ++i)
			sol *= evaluate_d(v[i], i);
		return sol;
	}

	Filter* get_subfilter(const unsigned int d) const
	{
		return new GaussianFilter(size, get_sigma(d));
	}
};
//GaussianFilter
#endif //_GAUSSIAN_FILTER_H_
