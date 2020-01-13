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

#ifndef _DENSITY_ESTIMATION_KERNEL_H_
#define _DENSITY_ESTIMATION_KERNEL_H_

#include "bunnykiller.h"

#include <cmath>

namespace DensityEsimationKernels
{
	//=============================================================
	template<unsigned D>
	class Perlin
	{
	public:
		Real operator()(const Real t, const Real bandwidth) const
		{
			Real tn = t / bandwidth;
			Real k_t = 1. - 6.*pow(tn, 5) + 15.*pow(tn, 4) - 10.*pow(tn, 3);
			return c(bandwidth)*k_t;
		}
		Real c(const Real bandwidth) const;
	};
	//-------------------------------------------------------------
	template<>
	Real Perlin<3>::c(const Real bandwidth) const
	{
		return 1. / (M_PI);
	}
	//-------------------------------------------------------------
	template<>
	Real Perlin<2>::c(const Real bandwidth) const
	{
		return 7. / (2.*M_PI*bandwidth*bandwidth);
	}
	//-------------------------------------------------------------
	template<>
	Real Perlin<1>::c(const Real bandwidth) const
	{
		return 1. / bandwidth;
	}
	//=============================================================
	template<unsigned D>
	class Box
	{
	public:
		Real operator()(const Real t, const Real bandwidth) const
		{
			return c(bandwidth);
		}
		Real c(const Real bandwidth) const;
	};
	//-------------------------------------------------------------
	template<>
	Real Box<3>::c(const Real bandwidth)const
	{
		return 3. / (4.*M_PI*bandwidth*bandwidth*bandwidth);
	}
	//-------------------------------------------------------------
	template<>
	Real Box<2>::c(const Real bandwidth)const
	{
		return 1. / (M_PI*bandwidth*bandwidth);
	}
	//-------------------------------------------------------------
	template<>
	Real Box<1>::c(const Real bandwidth)const
	{
		return 1. / (2.*bandwidth);
	}
}; // DensityEsimationKernels

#endif // _DENSITY_ESTIMATION_KERNEL_H_