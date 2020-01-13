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

#ifndef _UBER_FILTER_H_
#define _UBER_FILTER_H_

#include "Filter.h"
#include <vector>

class UberFilter: public Filter
{
	std::vector<const Filter*> filters;

	Real evaluate_d( Real x, unsigned int d )const
	{
		if (d < filters.size())
			return filters[d]->evaluate(x);
		else
			return filters[0]->evaluate(x);
	}
public:
	UberFilter( const Filter* f)
		:Filter() 
	{	filters.push_back(f);		}	
	UberFilter( const Filter* fx, const Filter* fy )
		:Filter() 
	{	filters.push_back(fx);
		filters.push_back(fy);	}
	UberFilter( const Filter* fx, const Filter* fy, const Filter* fz )
		:Filter()
	{	filters.push_back(fx);
		filters.push_back(fy);
		filters.push_back(fz);	}
	UberFilter(int _size, Filter* f)
		:Filter(_size)
	{	filters.push_back(f);		}

	UberFilter(int _size, Filter* fx, const Filter* fy )
		:Filter(_size) 
	{	filters.push_back(fx);
		filters.push_back(fy);	}

	UberFilter(int _size, Filter* fx, const Filter* fy, const Filter* fz )
		:Filter(_size)
	{	filters.push_back(fx);
		filters.push_back(fy);
		filters.push_back(fz);	}

	~UberFilter(){}

	Real evaluate(Real x)const { return filters[0]->evaluate(x); }
	Real evaluate(Real x, Real y)const { return evaluate(x)*evaluate_d(y,1); }
	Real evaluate(Real x, Real y, Real z)const { return evaluate(x)*evaluate_d(y,1)*evaluate_d(z,2); }
	Real evaluate(Real x, Real y, Real z, Real w)const {return evaluate(x)*evaluate_d(y,1)*evaluate_d(z,2)*evaluate_d(w,3);};
	Real evaluate(std::vector<Real> &v)const{Real sol(1.); for(unsigned int i=0; i<v.size(); ++i) sol*=evaluate_d(v[i],i); return sol;}

	Filter *get_subfilter(const unsigned int d)const
	{
		if (d < filters.size())
			return filters[d]->get_subfilter(d);
		else
			return filters[0]->get_subfilter(d);
	}
}; //UberFilter
#endif //_UBER_FILTER_H
