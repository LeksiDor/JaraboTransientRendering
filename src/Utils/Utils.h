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

#ifndef _UTILS_H_
#define _UTILS_H_

#include "bunnykiller.h"

#include <vector>
#include <algorithm>

namespace Utils
{
	template<class K, class V>
	V interpolate(const std::vector<K> keys, const V* values, const K &x)
	{
		if (x < keys[0])
			return values[0];

		if (x >= (*keys.end()))
			return values[keys.size() - 1];

		typename std::vector<K>::const_iterator low = std::lower_bound(keys.begin(), keys.end(), x);
		unsigned int idx = low - keys.begin();
		Real weight = (x - keys[idx]) / static_cast<Real>(keys[idx + 1] - keys[idx]);

		return (1. - weight) * values[idx] + weight * values[idx + 1];
	}

	template<class T0, class T1, class T2>
	inline T0 clamp(const T0 value, const T1 min, const T2 max)
	{
		return std::max(T0(min), std::min(value, T0(max)));
	}
}
;
//Utils

#endif //_UTILS_H_
