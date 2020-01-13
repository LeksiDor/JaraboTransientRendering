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

#ifndef _FILTER_H_
#define _FILTER_H_

#include "bunnykiller.h"

#include <vector>

#include "LinearAlgebra/Vector2.h"
#include "LinearAlgebra/Vector3.h"

#define _DEFAULT_FILTER_SIZE 3

class Filter
{
protected:
	int size;
	Real inv_size;
public:
	Filter(int _size = _DEFAULT_FILTER_SIZE) :
		size(_size), inv_size(1. / Real(_size))
	{}

	virtual ~Filter()
	{}

	int get_size() const
	{
		return size;
	}

	Real get_inv_size() const
	{
		return inv_size;
	}

	virtual Real evaluate(Real x) const = 0;
	virtual Real evaluate(Real x, Real y) const = 0;
	virtual Real evaluate(Real x, Real y, Real z) const = 0;
	virtual Real evaluate(Real x, Real y, Real z, Real w) const = 0;
	virtual Real evaluate(std::vector<Real> &v) const = 0;

	virtual Real evaluate(const Vector2 &v) const
	{
		return evaluate(v[0], v[1]);
	}

	virtual Real evaluate(const Vector3 &v) const
	{
		return evaluate(v[0], v[1], v[2]);
	}
	
	virtual Filter *get_subfilter(const unsigned int d) const = 0;
}; // Filter
#endif // _FILTER_H
