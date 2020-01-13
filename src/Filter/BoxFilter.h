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

#ifndef _BOX_FILTER_H_
#define _BOX_FILTER_H_

#include "bunnykiller.h"

#include "Filter/Filter.h"

class BoxFilter: public Filter
{
public:
	BoxFilter() :
		Filter()
	{}
	BoxFilter(Real size) :
		Filter(size)
	{}

	virtual ~BoxFilter()
	{}

	Real evaluate(Real x) const
	{
		return 1.;
	}

	Real evaluate(Real x, Real y) const
	{
		return 1.;
	}

	Real evaluate(Real x, Real y, Real z) const
	{
		return 1.;
	}

	Real evaluate(Real x, Real y, Real z, Real w) const
	{
		return 1.;
	}

	Real evaluate(std::vector<Real> &v) const
	{
		return 1.;
	}
	
	Filter *get_subfilter(const unsigned int d)const
	{
		return new BoxFilter(size);
	}
}; // BoxFilter

#endif // _BOX_FILTER_H
